#!/usr/bin/env python3
"""motif_caller — anchor-free, iterative R-M methylation motif finder.

Drop-in for STREME inside MicrobeMod. Reads centered methylation windows
(pos.fa) and a random control (neg.fa); writes streme.xml.

Usage:
    motif_caller.py . -p pos.fa --n neg.fa -o out_dir
"""
from __future__ import annotations
import math, os, sys, time
import numpy as np

# ── Algorithm constants ──────────────────────────────────────────────────────
SEQ_LEN, MIN_W, MAX_W, MAX_MOTIFS = 26, 4, 15, 25
IC_THRESH = 0.7
FREQ_CUTOFF, FREQ_CUTOFF_TRI = 0.8, 0.95
MIN_SEED_COUNT, MIN_MOTIF_SITES, MIN_MOTIF_SPECIFIC = 5, 10, 4
MAX_LOG10_EVAL, N_TESTS, PSEUDO = -1.0, float(MAX_MOTIFS), 0.5
MAX_NEG_MATCH_RATE = 0.2
MIN_ENRICH_RATIO, ENRICH_NEG_FLOOR = 3.0, 0.001
MIN_SPACER, MAX_SPACER = 1, 12
MIN_FAR_HALF, MAX_BIPARTITE_TOTAL = 2, 17
BIP_FAR_IC_THRESH, MIN_FAR_IC_SUM = 0.5, 2.5
RESCUE_ACTIVE_FLOOR, RESCUE_MIN_RESULTS = 80, 6
RESCUE_ACTIVE_FRAC, RESCUE_SEED_SCORE_FLOOR = 0.03, 0.0
SEED_ACTIVE_N_THRESH = 0.265
METH_CENTER, CENTER_TOL = SEQ_LEN // 2, 2

# Encoded bases: A=0 C=1 G=2 T=3 N/other=4
_ENC = np.full(256, 4, dtype=np.uint8)
for c, v in [('A',0),('a',0),('C',1),('c',1),('G',2),('g',2),
             ('T',3),('t',3),('U',3),('u',3)]:
    _ENC[ord(c)] = v
COMP = np.array([3, 2, 1, 0, 4], dtype=np.uint8)

# IUPAC bit-encoding (A=1, C=2, G=4, T=8); union of bits = degenerate code
IUPAC_BITS = {'A':1,'C':2,'G':4,'T':8,'M':3,'R':5,'W':9,'S':6,'Y':10,'K':12,
              'V':7,'H':11,'D':13,'B':14,'N':15}
BITS_IUPAC = {v: k for k, v in IUPAC_BITS.items()}
COMP_IUPAC = dict(zip("ACGTRYWSKMBDHVN", "TGCAYRWSMKVHDBN"))

def iupac_bits(c): return IUPAC_BITS.get(c, 0)
def revcomp_iupac(s): return "".join(COMP_IUPAC.get(c, 'N') for c in reversed(s))
def comp_iupac(c): return COMP_IUPAC.get(c, 'N')

# Match table[bits, encoded_base] → bool
_MATCH_TABLE = np.zeros((16, 5), dtype=bool)
for b in IUPAC_BITS.values():
    if b & 1: _MATCH_TABLE[b, 0] = True
    if b & 2: _MATCH_TABLE[b, 1] = True
    if b & 4: _MATCH_TABLE[b, 2] = True
    if b & 8: _MATCH_TABLE[b, 3] = True

# ── FASTA loader: returns [n_seqs, SEQ_LEN] uint8 array padded with 4 ───────
def load_fasta(path):
    seqs, cur = [], []
    with open(path, "rb") as f:
        for ln in f:
            ln = ln.rstrip(b"\r\n")
            if not ln: continue
            if ln[:1] == b">":
                if cur: seqs.append(b"".join(cur)); cur = []
            else:
                cur.append(ln)
        if cur: seqs.append(b"".join(cur))
    if not seqs:
        return np.zeros((0, SEQ_LEN), dtype=np.uint8)
    out = np.full((len(seqs), SEQ_LEN), 4, dtype=np.uint8)
    for i, s in enumerate(seqs):
        L = min(len(s), SEQ_LEN)
        if L:
            out[i, :L] = _ENC[np.frombuffer(s[:L], dtype=np.uint8)]
    return out

# ── Core numerical primitives ───────────────────────────────────────────────
def freq_matrix(seqs, start, end):
    """[width, 4] Laplace-smoothed frequencies; rows with N in window excluded."""
    w = end - start
    if seqs.shape[0] == 0 or w <= 0:
        return np.full((max(0, w), 4), 0.25)
    win = seqs[:, start:end]
    valid = ~(win == 4).any(axis=1)
    if not valid.any():
        return np.full((w, 4), 0.25)
    sub = win[valid]
    cnt = np.stack([(sub == j).sum(axis=0) for j in range(4)], axis=1)
    s = cnt.sum(axis=1, keepdims=True).astype(float)
    return (cnt + PSEUDO) / (s + 4 * PSEUDO)

def ic(pos, bg):
    """Information content (bits) of a position vs background."""
    mask = (pos > 1e-12) & (bg > 1e-12)
    return float(np.sum(pos[mask] * np.log(pos[mask] / bg[mask])) / math.log(2))

# ── Fisher's exact (one-tailed, log10 p-value) ──────────────────────────────
def _log_factorial(n): return math.lgamma(n + 1)

def _log_binom(n, k):
    if k < 0 or k > n: return -math.inf
    if k == 0 or k == n: return 0.0
    return _log_factorial(n) - _log_factorial(k) - _log_factorial(n - k)

def _log_hyper(a, b, c, d):
    n, K, sn = a + b + c + d, a + c, a + b
    return _log_binom(K, a) + _log_binom(n - K, sn - a) - _log_binom(n, sn)

def fisher_log10_pvalue(a, b, c, d):
    """One-tailed Fisher's exact: P(X >= a | hypergeometric). Returns log10."""
    if a == 0: return 0.0
    upper = min(a + b, a + c)
    log_obs = _log_hyper(a, b, c, d)
    terms, cur, ka = [log_obs], log_obs, a
    while True:
        ka += 1
        if ka > upper: break
        num = (a + c - ka + 1) * (a + b - ka + 1)
        den = ka * (d + ka - a)
        if den <= 0: break
        cur += math.log(num / den)
        terms.append(cur)
        if cur < log_obs - 40 * math.log(10): break
    mx = max(terms)
    return (mx + math.log(sum(math.exp(t - mx) for t in terms))) / math.log(10)

# ── IUPAC consensus ─────────────────────────────────────────────────────────
_TWO_LETTER = {(0,1):'M',(0,2):'R',(0,3):'W',(1,2):'S',(1,3):'Y',(2,3):'K'}
_THREE_LETTER = {0:'B', 1:'D', 2:'H', 3:'V'}  # missing-base index → B/D/H/V

def iupac_char(freq):
    """Map [pA, pC, pG, pT] → IUPAC char.
    >=80% one base → A/C/G/T; >=80% two bases → 2-letter; >=95% three → 3-letter."""
    a, c, g, t = float(freq[0]), float(freq[1]), float(freq[2]), float(freq[3])
    if a >= FREQ_CUTOFF: return 'A'
    if c >= FREQ_CUTOFF: return 'C'
    if g >= FREQ_CUTOFF: return 'G'
    if t >= FREQ_CUTOFF: return 'T'
    ranked = sorted([(a, 0), (c, 1), (g, 2), (t, 3)], reverse=True)
    if ranked[0][0] + ranked[1][0] >= FREQ_CUTOFF:
        idx = tuple(sorted([ranked[0][1], ranked[1][1]]))
        return _TWO_LETTER[idx]
    if ranked[0][0] + ranked[1][0] + ranked[2][0] >= FREQ_CUTOFF_TRI:
        return _THREE_LETTER[ranked[3][1]]
    return 'N'

# Background-corrected IUPAC consensus.  Each candidate base only counts toward
# the consensus if its observed frequency is enriched at least BG_ENR_THRESH-
# fold over its frequency in the reference genome.  This suppresses spurious
# IUPAC over-calls in genomes with skewed composition, where the absolute-
# frequency gate in iupac_char treats positions that just track local genome
# composition as informative.  Falls back to iupac_char when bg is near
# uniform (e.g. E. coli, GC ~50%) since the per-base gate then only tightens
# 2-letter calls without adding bg-correction value.
#
# 3-letter codes (B/D/H/V) are NOT emitted by this function: max per-base
# enrichment over a 3-base group at uniform bg is 4/3≈1.33, below the 1.5x
# gate.  Asymmetric depletion tests at flank positions over-fire on local-bg
# skew, generating spurious BGGCA-style boundary extensions.
BG_ENR_THRESH = 1.5
BG_UNIFORM_TOL = 0.05
def iupac_char_bg(freq, bg):
    bg_skew = max(abs(float(b) - 0.25) for b in bg)
    if bg_skew < BG_UNIFORM_TOL:
        return iupac_char(freq)
    a, c, g, t = float(freq[0]), float(freq[1]), float(freq[2]), float(freq[3])
    bga, bgc, bgg, bgt = max(float(bg[0]), 1e-6), max(float(bg[1]), 1e-6), max(float(bg[2]), 1e-6), max(float(bg[3]), 1e-6)
    sa = a if a >= BG_ENR_THRESH * bga else 0.0
    sc = c if c >= BG_ENR_THRESH * bgc else 0.0
    sg = g if g >= BG_ENR_THRESH * bgg else 0.0
    st = t if t >= BG_ENR_THRESH * bgt else 0.0
    if sa >= FREQ_CUTOFF: return 'A'
    if sc >= FREQ_CUTOFF: return 'C'
    if sg >= FREQ_CUTOFF: return 'G'
    if st >= FREQ_CUTOFF: return 'T'
    ranked = sorted([(sa, 0), (sc, 1), (sg, 2), (st, 3)], reverse=True)
    if ranked[0][0] + ranked[1][0] >= FREQ_CUTOFF:
        idx = tuple(sorted([ranked[0][1], ranked[1][1]]))
        return _TWO_LETTER[idx]
    return 'N'

def compute_genome_bg(fasta_path):
    """Compute background base frequencies [pA, pC, pG, pT] from a reference
    genome FASTA.  Returns [0.25, 0.25, 0.25, 0.25] if no parsable bases."""
    cnt = [0, 0, 0, 0]
    code = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'): continue
            for ch in line.strip().upper():
                i = code.get(ch)
                if i is not None: cnt[i] += 1
    total = sum(cnt)
    if total == 0: return [0.25, 0.25, 0.25, 0.25]
    return [cnt[i] / total for i in range(4)]

def compute_bg_from_seqs(seqs):
    """Background frequencies from an encoded seq array (0=A,1=C,2=G,3=T,4=N)."""
    import numpy as _np
    if seqs is None or len(seqs) == 0:
        return [0.25, 0.25, 0.25, 0.25]
    flat = _np.asarray(seqs).ravel()
    valid = flat[flat < 4]
    if valid.size == 0:
        return [0.25, 0.25, 0.25, 0.25]
    counts = _np.bincount(valid, minlength=4)
    total = float(counts.sum())
    return [float(counts[i] / total) for i in range(4)]

def consensus_seed_aware(filtered, active_freq_full, start, end, seed_start, seed_end):
    """IUPAC consensus from filtered set; positions inside the seed where the
    pre-filter active distribution looks uniform get N (the seed locked them
    by accident, not biology)."""
    f_filt = freq_matrix(filtered, start, end)
    out_freq, out_iupac = [], []
    for p in range(start, end):
        if seed_start <= p < seed_end and float(np.max(active_freq_full[p])) < SEED_ACTIVE_N_THRESH:
            f = np.array([0.25] * 4)
        else:
            f = f_filt[p - start]
        out_freq.append(f)
        out_iupac.append(iupac_char(f))
    return np.array(out_freq), "".join(out_iupac)

# ── IUPAC string utilities ──────────────────────────────────────────────────
def iupac_specificity(c):
    if c in "ACGT":  return 1.0
    if c in "MRWSYK": return 0.5
    if c in "BDHV":  return 0.20751874963942  # -log4(3/4)
    return 0.0

def iupac_levenshtein(a, b):
    """Edit distance; substitution costs 0 if IUPAC bit-sets overlap, else 1."""
    n, m = len(a), len(b)
    prev = list(range(m + 1))
    for i in range(1, n + 1):
        cur = [i] + [0] * m
        ai = iupac_bits(a[i - 1])
        for j in range(1, m + 1):
            sub = 0 if (ai & iupac_bits(b[j - 1])) else 1
            cur[j] = min(cur[j - 1] + 1, prev[j] + 1, prev[j - 1] + sub)
        prev = cur
    return prev[m]

def _is_specific(b): return b and bin(b).count('1') <= 2  # ≤2-base IUPAC code

def specific_overlap(a, b):
    """Count positions where both motifs are ≤2-letter AND bits overlap."""
    n = 0
    for x, y in zip(a, b):
        bx, by = iupac_bits(x), iupac_bits(y)
        if _is_specific(bx) and _is_specific(by) and (bx & by): n += 1
    return n

def n_specific_positions(a):
    return sum(1 for c in a if _is_specific(iupac_bits(c)))

def n_identical_specific(a, b):
    """Count positions where both are ≤2-letter AND have IDENTICAL bits."""
    return sum(1 for x, y in zip(a, b)
               if _is_specific(iupac_bits(x)) and iupac_bits(x) == iupac_bits(y))

def n_specific_bit_disjoint(a, b):
    """Count positions where both are specific but bits disjoint (= true mismatch)."""
    n = 0
    for x, y in zip(a, b):
        bx, by = iupac_bits(x), iupac_bits(y)
        if _is_specific(bx) and _is_specific(by) and (bx & by) == 0: n += 1
    return n

def merge_iupac(a, b):
    """Position-wise union of IUPAC bit-sets."""
    return "".join(x if x == y else BITS_IUPAC.get(iupac_bits(x) | iupac_bits(y), 'N')
                   for x, y in zip(a, b))

def motifs_alignable(a, b, max_h):
    """True if a and b are essentially the same motif (or its RC).
    Gates: lengths within 1; IUPAC-Lev <= max_h in fwd or rc; specific-base
    overlap covers max(4, ceil(0.7 * min_specific)); same-length lev>=2 must
    have <=1 disjoint specific positions."""
    if abs(len(a) - len(b)) > 1: return False
    b_rc = revcomp_iupac(b)
    lev_f, lev_r = iupac_levenshtein(a, b), iupac_levenshtein(a, b_rc)
    if lev_f > max_h and lev_r > max_h: return False
    spec_a = n_specific_positions(a)
    same_len = len(a) == len(b)
    for cand, lev in [(b, lev_f), (b_rc, lev_r)]:
        if lev > max_h: continue
        floor = max(4, math.ceil(min(spec_a, n_specific_positions(cand)) * 0.7))
        if specific_overlap(a, cand) < floor: continue
        if same_len and lev >= 2 and n_specific_bit_disjoint(a, cand) > 1: continue
        return True
    return False

# ── Sequence-matching primitives ────────────────────────────────────────────
def iupac_to_bits_array(s):
    return np.array([iupac_bits(c) for c in s], dtype=np.uint8)

def matches_centered_mask(seqs, pat_bits, pat_rc_bits):
    """Boolean mask: which seqs have pat or pat_rc anywhere with the matched
    span overlapping METH_CENTER ± CENTER_TOL."""
    m = len(pat_bits)
    n = seqs.shape[0]
    if m > seqs.shape[1] or n == 0:
        return np.zeros(n, dtype=bool)
    lo, hi = max(0, METH_CENTER - CENTER_TOL), METH_CENTER + CENTER_TOL
    found = np.zeros(n, dtype=bool)
    for s in range(seqs.shape[1] - m + 1):
        if s > hi: break
        if s + m <= lo: continue
        win = seqs[:, s:s + m]
        for pat in (pat_bits, pat_rc_bits):
            ok = np.ones(n, dtype=bool)
            for j in range(m):
                ok &= _MATCH_TABLE[pat[j], win[:, j]]
                if not ok.any(): break
            found |= ok
    return found

def matches_at_mask(seqs, pat_bits, start):
    """Boolean mask: which seqs match pat at exactly position `start`."""
    m = len(pat_bits)
    if start + m > seqs.shape[1]: return np.zeros(seqs.shape[0], dtype=bool)
    win = seqs[:, start:start + m]
    out = np.ones(seqs.shape[0], dtype=bool)
    for j in range(m):
        out &= _MATCH_TABLE[pat_bits[j], win[:, j]]
        if not out.any(): break
    return out

# ── Canonical k-mer counting ────────────────────────────────────────────────
def count_canonical_kmers(seqs, start, k):
    """{canonical_int: count} for k-mers at window [start, start+k).
    Canonical = min(kmer, rc_kmer) so palindromes count once."""
    end = start + k
    if seqs.shape[1] < end or seqs.shape[0] == 0: return {}
    win = seqs[:, start:end]
    valid = ~(win == 4).any(axis=1)
    if not valid.any(): return {}
    sub = win[valid].astype(np.int64)
    shifts = np.arange(k, dtype=np.int64) * 2
    kms = (sub << shifts).sum(axis=1)
    rcs = np.zeros_like(kms)
    for i in range(k):
        rcs |= (3 - ((kms >> (i * 2)) & 3)) << ((k - 1 - i) * 2)
    canon = np.minimum(kms, rcs)
    unique, counts = np.unique(canon, return_counts=True)
    return dict(zip(unique.tolist(), counts.tolist()))

def decode_kmer(km, k):
    return np.array([(km >> (i * 2)) & 3 for i in range(k)], dtype=np.uint8)

# ── Motif data class ────────────────────────────────────────────────────────
class Motif:
    __slots__ = ("idx", "iupac", "width", "evalue", "pvalue", "total_sites", "freq")
    def __init__(self, idx, iupac, width, evalue, pvalue, total_sites, freq):
        self.idx, self.iupac, self.width = idx, iupac, width
        self.evalue, self.pvalue, self.total_sites = evalue, pvalue, total_sites
        self.freq = freq

# ── Dedup ───────────────────────────────────────────────────────────────────
def dedup_results(motifs, max_h):
    """Drop or merge near-duplicate motifs.
    Two paths:
      1) Same length, IUPAC-Lev=1: merge to bit-union (with specific-overlap
         + identical-specific-overlap floors + ≥0.5 support ratio).
      2) motifs_alignable (lev<=max_h fwd or rc): drop the later motif.
         For RC-aligned same-length pairs differing at one position, broaden
         the kept motif's IUPAC at that position toward the bit-union."""
    kept = []
    for m in motifs:
        handled = False
        for i, k in enumerate(kept):
            ac, kc = m.iupac, k.iupac
            if len(ac) == len(kc) and iupac_levenshtein(ac, kc) == 1:
                min_spec = min(n_specific_positions(ac), n_specific_positions(kc))
                so_floor = max(4, math.ceil(min_spec * 0.70))
                id_floor = max(4, math.ceil(min_spec * 0.60))
                lo, hi = sorted([m.total_sites, k.total_sites])
                if (specific_overlap(ac, kc) >= so_floor
                        and n_identical_specific(ac, kc) >= id_floor
                        and hi > 0 and lo / hi >= 0.5):
                    kept[i].iupac = merge_iupac(ac, kc)
                    kept[i].total_sites = hi
                    handled = True
                    break
            if motifs_alignable(ac, kc, max_h):
                # RC-merge: broaden kept[i] at the single differing position
                if len(ac) == len(kc):
                    krc = revcomp_iupac(kc)
                    if iupac_levenshtein(ac, krc) == 1:
                        diffs = [p for p in range(len(ac))
                                 if (iupac_bits(ac[p]) & iupac_bits(krc[p])) == 0]
                        if len(diffs) == 1:
                            p = diffs[0]
                            kept_pos = len(kc) - 1 - p
                            cur_bits = iupac_bits(kc[kept_pos])
                            var_bits = iupac_bits(comp_iupac(ac[p]))
                            if (_is_specific(cur_bits) and (var_bits & ~cur_bits)):
                                lst = list(kc)
                                lst[kept_pos] = BITS_IUPAC.get(cur_bits | var_bits, 'N')
                                kept[i].iupac = "".join(lst)
                handled = True
                break
        if not handled:
            kept.append(m)
    return kept

# ── Main motif finder ───────────────────────────────────────────────────────
def find_motifs(pos_seqs, neg_seqs, bg=None):
    """Iterative: find best seed → filter → extend → check for bipartite far-half →
    consensus → garbage gates → Fisher → mask matched → repeat."""
    active = pos_seqs.copy()
    results = []
    pos_total = pos_seqs.shape[0]
    neg_freq_full = (np.full((SEQ_LEN, 4), 0.25) if neg_seqs.shape[0] == 0
                     else freq_matrix(neg_seqs, 0, SEQ_LEN))

    seed_blacklist, seed_canon_blacklist = set(), set()
    attempts = 0
    while attempts < MAX_MOTIFS * 12 and len(results) < MAX_MOTIFS:
        attempts += 1
        if active.shape[0] < MIN_MOTIF_SITES: break
        active_freq_full = freq_matrix(active, 0, SEQ_LEN)
        pos_n, neg_n = float(active.shape[0]), float(neg_seqs.shape[0])

        # Rescue mode: deeper iteration when residual signal remains
        active_frac = active.shape[0] / pos_total if pos_total > 0 else 0.0
        rescue_active = (active.shape[0] > RESCUE_ACTIVE_FLOOR
                         and (len(results) < RESCUE_MIN_RESULTS
                              or active_frac >= RESCUE_ACTIVE_FRAC))

        # ── Step 1: seed scan ──
        best_score, best_seed = -math.inf, None
        for k in range(3 if rescue_active else MIN_W, MAX_W + 1):
            for start in range(SEQ_LEN - k + 1):
                pos_win = count_canonical_kmers(active, start, k)
                if not pos_win or all(c < MIN_SEED_COUNT for c in pos_win.values()):
                    continue
                neg_all = count_canonical_kmers(neg_seqs, start, k) if neg_n else {}
                for can, pc in pos_win.items():
                    if (pc < MIN_SEED_COUNT
                            or (k, start, can) in seed_blacklist
                            or (k, can) in seed_canon_blacklist):
                        continue
                    nc = neg_all.get(can, 0)
                    log_odds = math.log2((pc / pos_n) / ((nc + PSEUDO) / (neg_n + PSEUDO)))
                    score = pc * log_odds + k * 0.01
                    if score > best_score:
                        best_score, best_seed = score, (k, start, can)

        seed_floor = RESCUE_SEED_SCORE_FLOOR if rescue_active else 1.0
        if best_seed is None or best_score <= seed_floor: break
        seed_k, seed_start, seed_can = best_seed
        seed_end = seed_start + seed_k

        # ── Step 2: filter to seqs containing the seed at the seed window ──
        seed_fwd = decode_kmer(seed_can, seed_k)
        seed_rc = COMP[seed_fwd[::-1]]
        seed_palindromic = np.array_equal(seed_fwd, seed_rc)
        if seed_end > active.shape[1]: break
        win = active[:, seed_start:seed_end]
        match_fwd = (win == seed_fwd[None, :]).all(axis=1)
        match_rc  = (win == seed_rc[None, :]).all(axis=1)
        filt_fwd, filt_rc = active[match_fwd], active[match_rc]
        filt_union = active[match_fwd | match_rc]

        # Strand-collision split: for non-palindromic seeds with both strands
        # populating the seed window, partitioning to the dominant strand can
        # recover strand-asymmetric flank specificity (e.g., CCCAGG fwd vs
        # CCTGGG rc → without partition we get CCWGG; with partition CCCAGG).
        if (seed_palindromic
                or filt_fwd.shape[0] < MIN_MOTIF_SITES
                or filt_rc.shape[0] < MIN_MOTIF_SITES):
            filtered = filt_union
        else:
            min_part = min(filt_fwd.shape[0], filt_rc.shape[0])
            if min_part / (filt_fwd.shape[0] + filt_rc.shape[0]) < 0.25:
                filtered = filt_union
            else:
                dom = filt_fwd if filt_fwd.shape[0] >= filt_rc.shape[0] else filt_rc
                dom_freq = freq_matrix(dom, 0, SEQ_LEN)
                uni_freq = freq_matrix(filt_union, 0, SEQ_LEN)
                keep = False
                # (a) flanking position has IC crossing threshold only in partition
                for p in (seed_start - 1, seed_end):
                    if 0 <= p < SEQ_LEN:
                        dic, uic = ic(dom_freq[p], neg_freq_full[p]), ic(uni_freq[p], neg_freq_full[p])
                        if dic >= IC_THRESH and uic < IC_THRESH:
                            keep = True; break
                # (b) seed-internal position becomes specific in partition
                if not keep:
                    for p in range(seed_start, seed_end):
                        if iupac_char(dom_freq[p]) in "ACGT" and iupac_char(uni_freq[p]) not in "ACGT":
                            keep = True; break
                filtered = dom if keep else filt_union

        if filtered.shape[0] < MIN_MOTIF_SITES:
            if rescue_active:
                seed_blacklist.add((seed_k, seed_start, seed_can))
                continue
            break

        # ── Step 3: extend seed bounds via per-position IC ──
        filt_freq_full = freq_matrix(filtered, 0, SEQ_LEN)
        pos_ic_full = np.array([ic(filt_freq_full[i], neg_freq_full[i]) for i in range(SEQ_LEN)])
        ms, me = seed_start, seed_end
        # IC walk with one-position look-ahead: skip a single mid-IC base
        # (>=0.35 but <0.7) if the next base is a strong anchor (>=1.0)
        IC_LA_FLOOR, IC_LA_ANCHOR = 0.35, 1.0
        while ms > 0:
            cur = pos_ic_full[ms - 1]
            if cur >= IC_THRESH:
                ms -= 1
            elif ms >= 2 and cur >= IC_LA_FLOOR and pos_ic_full[ms - 2] >= IC_LA_ANCHOR:
                ms -= 2
            else: break
        while me < SEQ_LEN:
            cur = pos_ic_full[me]
            if cur >= IC_THRESH:
                me += 1
            elif me + 1 < SEQ_LEN and cur >= IC_LA_FLOOR and pos_ic_full[me + 1] >= IC_LA_ANCHOR:
                me += 2
            else: break
        # Trim to MAX_W keeping the higher-IC end
        while me - ms > MAX_W:
            if pos_ic_full[ms] <= pos_ic_full[me - 1]: ms += 1
            else: me -= 1

        # ── Step 4: bipartite far-half search (Type I R-M motifs) ──
        close_len = me - ms
        bipartite, best_far_ic = None, 0.0
        # Scan right of close
        for sp in range(MIN_SPACER, MAX_SPACER + 1):
            far_start = me + sp
            if far_start >= SEQ_LEN: break
            max_far = max(0, MAX_BIPARTITE_TOTAL - close_len - sp)
            far_end = far_start
            while far_end < SEQ_LEN and far_end - far_start < max_far:
                if pos_ic_full[far_end] < BIP_FAR_IC_THRESH: break
                far_end += 1
            far_len = far_end - far_start
            if far_len >= MIN_FAR_HALF:
                ic_sum = float(pos_ic_full[far_start:far_end].sum())
                if ic_sum >= MIN_FAR_IC_SUM and ic_sum > best_far_ic:
                    best_far_ic, bipartite = ic_sum, (True, sp, far_start, far_end)
        # Scan left of close
        for sp in range(MIN_SPACER, MAX_SPACER + 1):
            if ms < sp + MIN_FAR_HALF: continue
            far_end_pos = ms - sp
            max_far = max(0, MAX_BIPARTITE_TOTAL - close_len - sp)
            far_start = far_end_pos
            while far_start > 0 and far_end_pos - far_start < max_far:
                if pos_ic_full[far_start - 1] < BIP_FAR_IC_THRESH: break
                far_start -= 1
            far_len = far_end_pos - far_start
            if far_len >= MIN_FAR_HALF:
                ic_sum = float(pos_ic_full[far_start:far_end_pos].sum())
                if ic_sum >= MIN_FAR_IC_SUM and ic_sum > best_far_ic:
                    best_far_ic, bipartite = ic_sum, (False, sp, far_start, far_end_pos)

        # ── Step 5: build IUPAC consensus ──
        if bipartite is None:
            motif_freq, iupac_str = consensus_seed_aware(
                filtered, active_freq_full, ms, me, seed_start, seed_end)
        else:
            right, sp, far_start, far_end = bipartite
            close_freq, close_iupac = consensus_seed_aware(
                filtered, active_freq_full, ms, me, seed_start, seed_end)
            far_freq = freq_matrix(filtered, far_start, far_end)
            far_iupac = "".join(iupac_char(f) for f in far_freq)
            spacer = np.full((sp, 4), 0.25)
            spacer_str = "N" * sp
            if right:
                motif_freq = np.concatenate([close_freq, spacer, far_freq])
                iupac_str = close_iupac + spacer_str + far_iupac
            else:
                motif_freq = np.concatenate([far_freq, spacer, close_freq])
                iupac_str = far_iupac + spacer_str + close_iupac

        # ── Step 6: garbage gates + Fisher ──
        if sum(1 for c in iupac_str if c in "ACGT") < MIN_MOTIF_SPECIFIC:
            seed_blacklist.add((seed_k, seed_start, seed_can))
            seed_canon_blacklist.add((seed_k, seed_can))
            continue

        pat_bits = iupac_to_bits_array(iupac_str)
        pat_rc_bits = iupac_to_bits_array(revcomp_iupac(iupac_str))
        pos_match_mask = matches_centered_mask(active, pat_bits, pat_rc_bits)
        pos_match = int(pos_match_mask.sum())
        if pos_match < MIN_MOTIF_SITES:
            if rescue_active:
                seed_blacklist.add((seed_k, seed_start, seed_can)); continue
            break
        neg_match = int(matches_centered_mask(neg_seqs, pat_bits, pat_rc_bits).sum())

        neg_rate = neg_match / neg_n if neg_n else 0.0
        cap = 0.15 if len(iupac_str) <= 5 else MAX_NEG_MATCH_RATE
        if neg_rate > cap:
            seed_blacklist.add((seed_k, seed_start, seed_can))
            seed_canon_blacklist.add((seed_k, seed_can))
            continue

        eff_neg = max(neg_rate, ENRICH_NEG_FLOOR)
        if (pos_match / active.shape[0]) / eff_neg < MIN_ENRICH_RATIO:
            seed_blacklist.add((seed_k, seed_start, seed_can))
            seed_canon_blacklist.add((seed_k, seed_can))
            continue

        # Raw-rate gate for short motifs (catches genome-frequency artifacts
        # that pass active-rate enrichment because active is depleted post-mask)
        stripped = iupac_str.strip("N")
        if len(stripped) <= 5:
            raw = (pos_match / pos_total) / eff_neg
            min_raw = 3.0 if len(stripped) <= 4 else 2.0
            if raw < min_raw:
                seed_blacklist.add((seed_k, seed_start, seed_can))
                seed_canon_blacklist.add((seed_k, seed_can))
                continue

        log10_p = fisher_log10_pvalue(pos_match, active.shape[0] - pos_match,
                                      neg_match, neg_seqs.shape[0] - neg_match)
        log10_e = log10_p + math.log10(N_TESTS)
        if log10_e > MAX_LOG10_EVAL:
            if rescue_active:
                seed_blacklist.add((seed_k, seed_start, seed_can)); continue
            break

        evalue = 0.0 if log10_e < -300 else 10 ** log10_e
        pvalue = 0.0 if log10_p < -300 else 10 ** log10_p
        idx = len(results) + 1
        print(f"  Motif {idx}: {iupac_str} (evalue={evalue:.2e}, "
              f"pos_match={pos_match}/{active.shape[0]}, neg_match={neg_match}/{neg_n:.0f})",
              file=sys.stderr)
        results.append(Motif(idx, iupac_str, motif_freq.shape[0],
                             evalue, pvalue, pos_match, motif_freq))

        # ── Step 7: mask matched seqs; reset blacklists ──
        active = active[~pos_match_mask]
        seed_blacklist.clear()
        seed_canon_blacklist.clear()

    # Fallback bg: if caller didn't supply, derive from neg_seqs.
    if bg is None:
        bg = compute_bg_from_seqs(neg_seqs)

    # ── Post-loop: dedup + iterated refinement → palindromize → final dedup ──
    refined = refine_motifs(dedup_results(results, 2), pos_seqs, bg)
    for _ in range(2):
        prev = [m.iupac for m in refined]
        refined = refine_motifs(refined, pos_seqs, bg)
        if all(m.iupac == p for m, p in zip(refined, prev)): break
    palindromized = palindromize(refined)
    return dedup_results(palindromized, 2)

def palindromize(motifs):
    """Merge each motif with its RC. Accept the merge if specificity loss
    <= 1.0 bit (catches palindromic R-M motifs reported as one strand,
    e.g. CCTGG → CCWGG; rejects non-palindromic motifs like GAGNNNNNGGG
    where the merge collapses to junk)."""
    out = []
    for m in motifs:
        rc_str = revcomp_iupac(m.iupac)
        if rc_str == m.iupac:
            out.append(m); continue
        merged = merge_iupac(m.iupac, rc_str)
        old_eff = sum(iupac_specificity(c) for c in m.iupac)
        new_eff = sum(iupac_specificity(c) for c in merged)
        if old_eff - new_eff > 1.0 + 1e-6:
            out.append(m); continue
        # Update freq at positions that changed; keep original freq elsewhere.
        new_freq = []
        for i, c in enumerate(merged):
            if i < len(m.iupac) and c == m.iupac[i]:
                new_freq.append(m.freq[i])
            else:
                bits = iupac_bits(c)
                n = bin(bits).count('1') if bits else 0
                if n > 0:
                    f = np.array([1.0 / n if (bits >> j) & 1 else 0.0
                                  for j in range(4)])
                else:
                    f = np.array([0.25] * 4)
                new_freq.append(f)
        out.append(Motif(m.idx, merged, len(merged), m.evalue, m.pvalue,
                          m.total_sites, np.array(new_freq)))
    return out

# ── Refinement: re-derive consensus from ALL pos_seqs that match centered ──
def refine_motifs(motifs, pos_seqs, bg):
    """For each motif, re-derive consensus from every centered match (not the
    seed-filtered subset). Optionally extend by up to 2 flank positions per
    side, trim leading/trailing N, trim borderline 2-letter flanks.

    `bg` is the genome-wide background base frequency [pA, pC, pG, pT],
    used to suppress IUPAC over-calls at positions that just track local
    genome composition."""
    out = []
    lo, hi = max(0, METH_CENTER - CENTER_TOL), METH_CENTER + CENTER_TOL
    for m in motifs:
        m_len = len(m.iupac)
        if not (0 < m_len <= SEQ_LEN):
            out.append(m); continue
        pat_bits = iupac_to_bits_array(m.iupac)
        pat_rc_bits = iupac_to_bits_array(revcomp_iupac(m.iupac))

        cnt = np.zeros((m_len, 4), dtype=np.int64)
        left_flanks = np.zeros((2, 4), dtype=np.int64)
        right_flanks = np.zeros((2, 4), dtype=np.int64)
        left_totals, right_totals = np.zeros(2, dtype=np.int64), np.zeros(2, dtype=np.int64)

        # Vectorized: for each start s in [lo-m_len+1, hi+1), check fwd/rc
        # match against all seqs. For each seq, take the LOWEST-offset hit
        # (fwd preferred over rc at the same offset), record its (s, is_rc).
        n = pos_seqs.shape[0]
        offset = np.full(n, -1, dtype=np.int64)
        is_rc_arr = np.zeros(n, dtype=bool)
        starts = []
        for s in range(pos_seqs.shape[1] - m_len + 1):
            if s > hi: break
            if s + m_len <= lo: continue
            starts.append(s)
        for s in starts:
            unset = (offset == -1)
            if not unset.any(): break
            sub = pos_seqs[unset]
            mfwd = matches_at_mask(sub, pat_bits, s)
            mrc  = matches_at_mask(sub, pat_rc_bits, s)
            # Reject windows containing N (4) — already handled by MATCH_TABLE
            # which never matches base 4, but double-check span has no N
            no_n = ~(sub[:, s:s+m_len] == 4).any(axis=1)
            mfwd &= no_n
            mrc  &= no_n
            new_fwd = unset.copy()
            new_fwd[unset] = mfwd
            new_rc = unset.copy()
            new_rc[unset] = mrc & ~mfwd  # fwd preferred at same offset
            offset[new_fwd | new_rc] = s
            is_rc_arr[new_rc] = True

        matched_mask = (offset != -1)
        total = int(matched_mask.sum())

        # Now accumulate cnt + flanks. Group matched seqs by (s, is_rc).
        slen = pos_seqs.shape[1]
        for s in starts:
            for rc in (False, True):
                idx = matched_mask & (offset == s) & (is_rc_arr == rc)
                if not idx.any(): continue
                sub = pos_seqs[idx]
                window = sub[:, s:s + m_len]  # [n_sub, m_len]
                if rc:
                    # Complement and reverse → motif orientation
                    rwin = COMP[window][:, ::-1]
                    for i in range(m_len):
                        col = rwin[:, i]
                        for j in range(4):
                            cnt[i, j] += int((col == j).sum())
                    # Left flank in motif orientation = pos s+m_len, s+m_len+1 (complemented)
                    for k in range(2):
                        p = s + m_len + k
                        if p < slen:
                            col = COMP[sub[:, p]]
                            valid = col < 4
                            for j in range(4):
                                left_flanks[k, j] += int(((col == j) & valid).sum())
                            left_totals[k] += int(valid.sum())
                    # Right flank = pos s-1, s-2 (complemented)
                    for k in range(2):
                        if s >= k + 1:
                            col = COMP[sub[:, s - 1 - k]]
                            valid = col < 4
                            for j in range(4):
                                right_flanks[k, j] += int(((col == j) & valid).sum())
                            right_totals[k] += int(valid.sum())
                else:
                    for i in range(m_len):
                        col = window[:, i]
                        for j in range(4):
                            cnt[i, j] += int((col == j).sum())
                    for k in range(2):
                        if s >= k + 1:
                            col = sub[:, s - 1 - k]
                            valid = col < 4
                            for j in range(4):
                                left_flanks[k, j] += int(((col == j) & valid).sum())
                            left_totals[k] += int(valid.sum())
                    for k in range(2):
                        p = s + m_len + k
                        if p < slen:
                            col = sub[:, p]
                            valid = col < 4
                            for j in range(4):
                                right_flanks[k, j] += int(((col == j) & valid).sum())
                            right_totals[k] += int(valid.sum())

        if total < MIN_MOTIF_SITES:
            out.append(m); continue

        def freq(c):
            s = float(c.sum())
            return (c + PSEUDO) / (s + 4 * PSEUDO)
        new_freq = [freq(cnt[i]) for i in range(m_len)]
        new_iupac = "".join(iupac_char_bg(f, bg) for f in new_freq)

        # Flank extension (greedy nearest-out, capped at SEQ_LEN)
        for k in range(2):
            if len(new_iupac) >= SEQ_LEN or left_totals[k] < MIN_MOTIF_SITES: break
            f = freq(left_flanks[k])
            c = iupac_char_bg(f, bg)
            if c == 'N': break
            new_freq.insert(0, f); new_iupac = c + new_iupac
        for k in range(2):
            if len(new_iupac) >= SEQ_LEN or right_totals[k] < MIN_MOTIF_SITES: break
            f = freq(right_flanks[k])
            c = iupac_char_bg(f, bg)
            if c == 'N': break
            new_freq.append(f); new_iupac += c

        # Trim flanking N and borderline 2-letter codes (top2_sum < 0.85)
        while new_iupac.startswith('N') and len(new_iupac) > MIN_W:
            new_iupac, new_freq = new_iupac[1:], new_freq[1:]
        while new_iupac.endswith('N') and len(new_iupac) > MIN_W:
            new_iupac, new_freq = new_iupac[:-1], new_freq[:-1]

        TOP2 = 0.85
        def top2_sum(f): s = sorted(f, reverse=True); return s[0] + s[1]
        did_trim = False
        while len(new_iupac) > MIN_W and new_iupac[0] in "MRWSYK" and top2_sum(new_freq[0]) < TOP2:
            new_iupac, new_freq = new_iupac[1:], new_freq[1:]; did_trim = True
        while len(new_iupac) > MIN_W and new_iupac[-1] in "MRWSYK" and top2_sum(new_freq[-1]) < TOP2:
            new_iupac, new_freq = new_iupac[:-1], new_freq[:-1]; did_trim = True

        # Accept refined motif only if strictly better, same-eff-shorter, or borderline-trim
        old_eff = sum(iupac_specificity(c) for c in m.iupac)
        new_eff = sum(iupac_specificity(c) for c in new_iupac)
        if (new_eff > old_eff + 0.001
                or (abs(new_eff - old_eff) < 0.001 and len(new_iupac) < len(m.iupac))
                or did_trim):
            out.append(Motif(m.idx, new_iupac, len(new_iupac),
                             m.evalue, m.pvalue, m.total_sites,
                             np.array(new_freq) if new_freq else m.freq))
        else:
            out.append(m)
    return out

# ── XML output ──────────────────────────────────────────────────────────────
def _format_sci(v):
    if v == 0.0: return "0.0"
    if math.isnan(v) or math.isinf(v): return "1.0"
    e = math.floor(math.log10(abs(v)))
    return f"{v / 10**e:.1f}e{e:+d}"

# IUPAC -> set of allowed bases (A,C,G,T index order)
_IUPAC_ALLOWED = {
    "A": (0,), "C": (1,), "G": (2,), "T": (3,),
    "W": (0, 3), "S": (1, 2), "M": (0, 1), "K": (2, 3), "R": (0, 2), "Y": (1, 3),
    "B": (1, 2, 3), "D": (0, 2, 3), "H": (0, 1, 3), "V": (0, 1, 2),
    "N": (0, 1, 2, 3),
}

def _iupac_pwm_row(ch):
    """Sharp 4-vector PWM for one IUPAC char, [pA,pC,pG,pT].

    Puts ~0.97 of the mass on the allowed base(s) so that MicrobeMod's
    assign_motifs() round-trips the call faithfully: a 1-letter code keeps
    max(freq) > 0.8, a 2-letter code keeps top-two-sum > 0.8.  Without this,
    assign_motifs would re-apply a raw-frequency gate and N-mask exactly the
    background-corrected positions this caller intentionally called."""
    allowed = _IUPAC_ALLOWED.get(ch, (0, 1, 2, 3))
    row = [0.01, 0.01, 0.01, 0.01]
    if len(allowed) == 4:                       # N: leave uniform
        return [0.25, 0.25, 0.25, 0.25]
    share = 0.97 / len(allowed)
    rest = 0.03 / (4 - len(allowed))
    for b in range(4):
        row[b] = share if b in allowed else rest
    return row

def write_xml(path, motifs, pos_n, neg_n):
    with open(path, "w") as f:
        f.write('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n')
        f.write('<STREME version="5.5.3" release="motif_caller-py-1.0">\n')
        f.write('  <model>\n')
        f.write('    <command_line>motif_caller</command_line>\n')
        f.write(f'    <train_positives count="{pos_n}" positions="{pos_n*SEQ_LEN}" '
                f'maxlen="{SEQ_LEN}" file="pos.fasta"/>\n')
        f.write(f'    <train_negatives count="{neg_n}" positions="{neg_n*SEQ_LEN}" '
                f'from="file" file="control.fasta"/>\n')
        f.write('  </model>\n  <motifs>\n')
        for m in motifs:
            log_e = math.log10(m.evalue) if m.evalue > 0 else -300.0
            f.write(f'    <motif id="{m.idx}-{m.iupac}" alt="STREME-{m.idx}" '
                    f'width="{m.width}" test_pvalue="{_format_sci(m.pvalue)}" '
                    f'test_evalue="{_format_sci(m.evalue)}" '
                    f'test_log_evalue="{log_e:.4f}" total_sites="{m.total_sites}">\n')
            for ch in m.iupac:
                fr = _iupac_pwm_row(ch)
                f.write(f'      <pos A="{fr[0]:.6f}" C="{fr[1]:.6f}" '
                        f'G="{fr[2]:.6f}" T="{fr[3]:.6f}"/>\n')
            f.write('    </motif>\n')
        f.write('  </motifs>\n</STREME>\n')

# ── TSV output (stand-alone subcommand) ────────────────────────────────────
def write_tsv(path, motifs, mod_type=None):
    """Compact tab-separated motif report. Columns:
        motif  width  total_sites  evalue  pvalue  methylation_type
    """
    with open(path, "w") as f:
        cols = ["motif", "width", "total_sites", "evalue", "pvalue"]
        if mod_type is not None: cols.append("methylation_type")
        f.write("\t".join(cols) + "\n")
        for m in motifs:
            row = [m.iupac, str(m.width), str(m.total_sites),
                   _format_sci(m.evalue), _format_sci(m.pvalue)]
            if mod_type is not None: row.append(mod_type)
            f.write("\t".join(row) + "\n")

# ── In-process API for MicrobeMod's call_methylation pipeline ──────────────
def run_from_fastas(pos_path, neg_path, out_dir, output_type="xml",
                    genome_path=None):
    """Read pos/neg FASTAs, find motifs, write streme.xml (or motifs.tsv).

    Returns the output directory if motifs were called, else None.
    Used by call_methylation as a drop-in for run_streme().

    If `genome_path` is given, the genome's base composition is used as the
    background for IUPAC consensus calls (suppresses spurious 2-letter IUPAC
    codes at positions that just track local genome bias).  Otherwise the
    background is derived from neg_seqs.
    """
    os.makedirs(out_dir, exist_ok=True)
    pos_seqs = load_fasta(pos_path)
    neg_seqs = load_fasta(neg_path)
    if pos_seqs.shape[0] < MIN_MOTIF_SITES:
        return None
    bg = compute_genome_bg(genome_path) if genome_path else None
    motifs = find_motifs(pos_seqs, neg_seqs, bg=bg)
    if output_type == "xml":
        write_xml(os.path.join(out_dir, "streme.xml"), motifs,
                  pos_seqs.shape[0], neg_seqs.shape[0])
    else:
        write_tsv(os.path.join(out_dir, "motifs.tsv"), motifs)
    return out_dir

# ── Stand-alone subcommand: bedmethyl + reference → motifs ────────────────
def _build_pos_neg_from_bed(bed_path, fasta_path, mod_code, percent_cutoff=0.66,
                            min_coverage=10, n_control=100000, seed=13):
    """Walk a modkit bedmethyl, extract SEQ_LEN windows centered on every
    methylated position passing the cutoffs, plus n_control random windows
    from the reference for negatives. Returns (pos_seqs, neg_seqs)."""
    import random
    from Bio import SeqIO
    contigs = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}
    half = SEQ_LEN // 2
    pos_strs = []
    with open(bed_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 18:                      # newer modkit space-delimits
                parts = parts[:-1] + parts[-1].split()  # the trailing stat columns
            if len(parts) < 18: continue
            try:
                contig, start = parts[0], int(parts[1])
                code = parts[3]
                cov = int(parts[9]); frac = float(parts[10]) / 100.0
                n_mod = int(parts[11])
            except (ValueError, IndexError):
                continue
            if code != mod_code: continue
            if cov < min_coverage: continue
            if frac < percent_cutoff: continue
            seq = contigs.get(contig)
            if seq is None: continue
            lo, hi = start - half, start + (SEQ_LEN - half)
            if lo < 0 or hi > len(seq): continue
            win = seq[lo:hi].upper()
            if "N" in win: continue
            pos_strs.append(win)
    rng = random.Random(seed)
    neg_strs = []
    keys = list(contigs.keys())
    weights = [len(contigs[k]) for k in keys]
    total_len = sum(weights)
    if total_len > 2 * SEQ_LEN:
        tries = 0
        while len(neg_strs) < n_control and tries < n_control * 5:
            tries += 1
            r = rng.randrange(total_len)
            acc = 0
            for k, w in zip(keys, weights):
                acc += w
                if r < acc: contig = k; break
            seq = contigs[contig]
            p = rng.randrange(half, len(seq) - SEQ_LEN + half)
            win = seq[p - half:p + (SEQ_LEN - half)].upper()
            if "N" not in win:
                neg_strs.append(win)

    def to_array(strs):
        if not strs: return np.zeros((0, SEQ_LEN), dtype=np.uint8)
        arr = np.full((len(strs), SEQ_LEN), 4, dtype=np.uint8)
        for i, s in enumerate(strs):
            L = min(len(s), SEQ_LEN)
            if L:
                arr[i, :L] = _ENC[np.frombuffer(s[:L].encode("ascii"), dtype=np.uint8)]
        return arr
    return to_array(pos_strs), to_array(neg_strs)

# Map modkit codes (file column 4) to user-facing methylation labels
_MOD_CODES = {"a": "6mA", "m": "5mC", "21839": "4mC", "h": "5hmC"}

def subcommand_main(bed_path, fasta_path, threads=1, output_type="tsv",
                    output_prefix=None, methylation_types="6mA,5mC,4mC",
                    percent_cutoff=0.66, min_coverage=10):
    """Stand-alone `MicrobeMod microbe_motif` entry point.

    Iterates over methylation types, builds pos/neg sequence sets directly
    from the bedmethyl + reference, runs find_motifs, and writes one TSV
    (or XML) per methylation type.
    """
    if output_prefix is None:
        output_prefix = os.path.splitext(os.path.basename(bed_path))[0]

    type_list = [t.strip() for t in methylation_types.split(",") if t.strip()]
    name_to_code = {"6mA": "a", "5mC": "m", "4mC": "21839", "5hmC": "h"}

    rows = []  # for combined TSV across mod types
    any_emitted = False
    for mod_label in type_list:
        code = name_to_code.get(mod_label)
        if code is None:
            print(f"Skipping unknown methylation type: {mod_label}", file=sys.stderr)
            continue
        print(f"[microbe_motif] building windows for {mod_label} ({code})…",
              file=sys.stderr)
        pos_seqs, neg_seqs = _build_pos_neg_from_bed(
            bed_path, fasta_path, code, percent_cutoff, min_coverage)
        n_pos = int(pos_seqs.shape[0])
        n_neg = int(neg_seqs.shape[0])
        print(f"  positives: {n_pos}  negatives: {n_neg}", file=sys.stderr)
        if n_pos < MIN_MOTIF_SITES:
            print(f"  too few positives (<{MIN_MOTIF_SITES}); skipping {mod_label}",
                  file=sys.stderr)
            continue
        bg = compute_genome_bg(fasta_path)
        motifs = find_motifs(pos_seqs, neg_seqs, bg=bg)
        any_emitted = any_emitted or bool(motifs)
        if output_type == "xml":
            xml_dir = f"{output_prefix}_{mod_label}_microbe_motif"
            os.makedirs(xml_dir, exist_ok=True)
            write_xml(os.path.join(xml_dir, "streme.xml"), motifs, n_pos, n_neg)
            print(f"  wrote {xml_dir}/streme.xml ({len(motifs)} motifs)",
                  file=sys.stderr)
        else:
            for m in motifs:
                rows.append((mod_label, m.iupac, m.width, m.total_sites,
                             m.evalue, m.pvalue))

    if output_type == "tsv":
        out_path = f"{output_prefix}_motifs.tsv"
        with open(out_path, "w") as f:
            f.write("methylation_type\tmotif\twidth\ttotal_sites\tevalue\tpvalue\n")
            for mod_label, iup, w, n, e, p in rows:
                f.write(f"{mod_label}\t{iup}\t{w}\t{n}\t"
                        f"{_format_sci(e)}\t{_format_sci(p)}\n")
        print(f"[microbe_motif] wrote {out_path} ({len(rows)} motifs across "
              f"{len(type_list)} methylation types)", file=sys.stderr)

# ── CLI ─────────────────────────────────────────────────────────────────────
def main():
    args = sys.argv[1:]
    pos_file = neg_file = out_dir = genome_file = None
    i = 0
    while i < len(args):
        a = args[i]
        if a == "-p" and i + 1 < len(args): pos_file = args[i+1]; i += 2
        elif a in ("-n", "--n") and i + 1 < len(args): neg_file = args[i+1]; i += 2
        elif a in ("-g", "--genome") and i + 1 < len(args): genome_file = args[i+1]; i += 2
        elif a == "-o" and i + 1 < len(args): out_dir = args[i+1]; i += 2
        elif a.startswith("--") and i + 1 < len(args) and not args[i+1].startswith('-'): i += 2
        else: i += 1

    if not (pos_file and neg_file and out_dir):
        sys.exit("Usage: motif_caller.py . -p pos.fa --n neg.fa [-g genome.fa] "
                 "-o out_dir")
    os.makedirs(out_dir, exist_ok=True)

    print("microbe_motif: loading sequences…", file=sys.stderr)
    t0 = time.time()
    pos_seqs, neg_seqs = load_fasta(pos_file), load_fasta(neg_file)
    print(f"  Positive: {pos_seqs.shape[0]}  Negative: {neg_seqs.shape[0]}", file=sys.stderr)

    if pos_seqs.shape[0] < MIN_MOTIF_SITES:
        write_xml(os.path.join(out_dir, "streme.xml"), [],
                  pos_seqs.shape[0], neg_seqs.shape[0])
        return

    bg = compute_genome_bg(genome_file) if genome_file else None
    if bg is not None:
        print(f"  bg (from genome): A={bg[0]:.3f} C={bg[1]:.3f} G={bg[2]:.3f} T={bg[3]:.3f}",
              file=sys.stderr)

    print("microbe_motif: finding motifs…", file=sys.stderr)
    motifs = find_motifs(pos_seqs, neg_seqs, bg=bg)
    print(f"microbe_motif: found {len(motifs)} motif(s)  ({time.time()-t0:.1f}s)",
          file=sys.stderr)
    write_xml(os.path.join(out_dir, "streme.xml"), motifs,
              pos_seqs.shape[0], neg_seqs.shape[0])

if __name__ == "__main__":
    main()
