"""Regression tests for the built-in (default) motif caller, microbe_motif.

These cover the pieces that the legacy STREME path never exercised: the IUPAC
consensus logic, the Fisher exact p-value, k-mer counting, the STREME-format
XML round-trip, and an end-to-end synthetic recovery of a known motif.

The motif caller is fully deterministic, so end-to-end output can be asserted
directly.
"""
import math
import os
import random
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from MicrobeMod import microbe_motif as mm


# ── IUPAC string utilities ──────────────────────────────────────────────────
def test_revcomp_iupac_palindromes_and_degenerate():
    assert mm.revcomp_iupac("GATC") == "GATC"          # palindrome
    assert mm.revcomp_iupac("CCWGG") == "CCWGG"         # W is self-complement
    assert mm.revcomp_iupac("GAATTC") == "GAATTC"       # EcoRI palindrome
    assert mm.revcomp_iupac("AACGTT") == "AACGTT"
    assert mm.revcomp_iupac("GGTGA") == "TCACC"         # non-palindrome
    # R (A/G) complements to Y (C/T)
    assert mm.revcomp_iupac("R") == "Y"


def test_iupac_char_thresholds():
    # >=80% one base -> that base
    assert mm.iupac_char([0.9, 0.04, 0.03, 0.03]) == "A"
    assert mm.iupac_char([0.03, 0.03, 0.04, 0.9]) == "T"
    # two bases summing >=80% -> 2-letter code (A+G -> R)
    assert mm.iupac_char([0.45, 0.05, 0.45, 0.05]) == "R"
    # C+T -> Y
    assert mm.iupac_char([0.05, 0.45, 0.05, 0.45]) == "Y"
    # uniform -> N
    assert mm.iupac_char([0.25, 0.25, 0.25, 0.25]) == "N"


def test_merge_iupac_bit_union():
    assert mm.merge_iupac("A", "G") == "R"      # 1 | 4 = 5 = R
    assert mm.merge_iupac("C", "T") == "Y"      # 2 | 8 = 10 = Y
    assert mm.merge_iupac("GATC", "GATC") == "GATC"
    # A vs C -> M, identical positions preserved
    assert mm.merge_iupac("AAAA", "CAAA") == "MAAA"


def test_iupac_levenshtein_overlap_is_free():
    # identical -> 0
    assert mm.iupac_levenshtein("GATC", "GATC") == 0
    # R overlaps A (bit sets share A) -> substitution is free
    assert mm.iupac_levenshtein("AATC", "RATC") == 0
    # disjoint specific bases -> cost 1
    assert mm.iupac_levenshtein("AATC", "CATC") == 1
    # one indel
    assert mm.iupac_levenshtein("GATC", "GAT") == 1


# ── Fisher's exact (one-tailed, log10) ──────────────────────────────────────
def test_fisher_log10_known_values():
    # a=0 is defined to return 0.0 (p = 1)
    assert mm.fisher_log10_pvalue(0, 5, 5, 5) == 0.0

    # 2x2 = [[1,0],[0,1]]: P(X>=1) for hypergeom(n=2,K=1,draws=1) = 0.5
    assert mm.fisher_log10_pvalue(1, 0, 0, 1) == pytest.approx(math.log10(0.5), abs=1e-9)

    # 2x2 = [[2,0],[0,2]]: P(X>=2) = C(2,2)C(2,0)/C(4,2) = 1/6
    assert mm.fisher_log10_pvalue(2, 0, 0, 2) == pytest.approx(
        math.log10(1.0 / 6.0), abs=1e-9
    )


def test_fisher_log10_monotonic_and_bounded():
    # p-value is a probability -> log10 <= 0 always
    assert mm.fisher_log10_pvalue(20, 5, 5, 20) <= 0.0
    # a stronger enrichment is more significant (more negative log10)
    weak = mm.fisher_log10_pvalue(12, 8, 8, 12)
    strong = mm.fisher_log10_pvalue(20, 0, 0, 20)
    assert strong < weak


# ── FASTA + k-mer primitives ────────────────────────────────────────────────
def test_load_fasta_encoding(tmp_path):
    p = tmp_path / "x.fa"
    p.write_text(">1\nACGTN\n>2\nTTTT\n")
    arr = mm.load_fasta(str(p))
    assert arr.shape == (2, mm.SEQ_LEN)
    # A=0 C=1 G=2 T=3 N=4
    assert list(arr[0, :5]) == [0, 1, 2, 3, 4]
    assert list(arr[1, :4]) == [3, 3, 3, 3]
    # unfilled positions are padded with 4
    assert arr[1, 4] == 4


def test_count_canonical_kmers_palindrome_counts_once():
    # GATC at positions 0..3 in two identical seqs; canonical == its own rc
    seqs = np.array([[2, 0, 3, 1] + [4] * (mm.SEQ_LEN - 4)] * 2, dtype=np.uint8)
    counts = mm.count_canonical_kmers(seqs, 0, 4)
    assert sum(counts.values()) == 2
    # decoding the single canonical key back yields GATC (its own rc)
    (canon,) = counts.keys()
    assert list(mm.decode_kmer(canon, 4)) == [2, 0, 3, 1]  # GATC


# ── XML round-trip ──────────────────────────────────────────────────────────
def test_write_xml_roundtrip(tmp_path):
    freq = np.array([mm._iupac_pwm_row(c) for c in "GATC"])
    motif = mm.Motif(1, "GATC", 4, evalue=1e-10, pvalue=1e-12, total_sites=42, freq=freq)
    out = tmp_path / "streme.xml"
    mm.write_xml(str(out), [motif], pos_n=100, neg_n=1000)

    root = ET.parse(str(out)).getroot()
    assert root.tag == "STREME"
    motif_nodes = root.find("motifs").findall("motif")
    assert len(motif_nodes) == 1
    node = motif_nodes[0]
    assert node.get("id") == "1-GATC"
    assert node.get("total_sites") == "42"
    pos_rows = node.findall("pos")
    assert len(pos_rows) == 4
    # each PWM row sums to ~1.0 and is sharply peaked on the called base
    for row, base in zip(pos_rows, "GATC"):
        vals = {k: float(v) for k, v in row.attrib.items()}
        assert sum(vals.values()) == pytest.approx(1.0, abs=1e-6)
        assert max(vals, key=vals.get) == base


# ── End-to-end synthetic recovery ───────────────────────────────────────────
def _embed_motif_fastas(tmp_path, motif, n_pos=200, n_neg=3000, seed=0):
    """Write pos/neg FASTAs: positives carry `motif` centered on METH_CENTER,
    negatives are uniform random. Returns (pos_path, neg_path)."""
    rng = random.Random(seed)
    bases = "ACGT"
    half = len(motif) // 2
    lo = mm.METH_CENTER - half  # place motif so its middle sits on the center

    def rand_seq():
        return [rng.choice(bases) for _ in range(mm.SEQ_LEN)]

    pos_path = tmp_path / "pos.fasta"
    with open(pos_path, "w") as f:
        for i in range(n_pos):
            s = rand_seq()
            s[lo : lo + len(motif)] = list(motif)
            f.write(">{}\n{}\n".format(i, "".join(s)))

    neg_path = tmp_path / "neg.fasta"
    with open(neg_path, "w") as f:
        for i in range(n_neg):
            f.write(">{}\n{}\n".format(i, "".join(rand_seq())))

    return str(pos_path), str(neg_path)


def test_find_motifs_recovers_palindrome(tmp_path):
    pos_path, neg_path = _embed_motif_fastas(tmp_path, "GATC", seed=1)
    pos = mm.load_fasta(pos_path)
    neg = mm.load_fasta(neg_path)
    motifs = mm.find_motifs(pos, neg)
    called = [m.iupac for m in motifs]
    assert any("GATC" in c for c in called), called


def test_find_motifs_recovers_longer_palindrome(tmp_path):
    # GANTC (HinfI) — palindromic with a degenerate center
    pos_path, neg_path = _embed_motif_fastas(tmp_path, "GAATTC", seed=2)
    pos = mm.load_fasta(pos_path)
    neg = mm.load_fasta(neg_path)
    motifs = mm.find_motifs(pos, neg)
    called = [m.iupac for m in motifs]
    assert any("GAATTC" in c for c in called), called


def test_find_motifs_no_signal_returns_nothing(tmp_path):
    # positives and negatives drawn from the same random distribution -> no motif
    rng = random.Random(3)
    bases = "ACGT"

    def write(path, n):
        with open(path, "w") as f:
            for i in range(n):
                seq = "".join(rng.choice(bases) for _ in range(mm.SEQ_LEN))
                f.write(">{}\n{}\n".format(i, seq))

    pos_path = tmp_path / "pos.fasta"
    neg_path = tmp_path / "neg.fasta"
    write(pos_path, 300)
    write(neg_path, 3000)
    motifs = mm.find_motifs(
        mm.load_fasta(str(pos_path)), mm.load_fasta(str(neg_path))
    )
    assert motifs == []


# ── Methylated-center constraint (issue #52) ────────────────────────────────
def test_constrain_center_narrows_to_modifiable():
    a_bits = mm._allowed_center_bits("a")      # A|T (6mA)
    c_bits = mm._allowed_center_bits("m")      # C|G (5mC)
    assert mm._constrain_center("GMCGKC", 1, a_bits) == "GACGKC"   # M -> A
    assert mm._constrain_center("GMATTC", 1, c_bits) == "GCATTC"   # M -> C
    assert mm._constrain_center("GACGGC", 1, a_bits) == "GACGGC"   # already A
    assert mm._constrain_center("GMCGKC", 1, None) == "GMCGKC"     # disabled

    # every modkit code / label maps to the right allowed + modifiable base, so
    # a typo in an alias can't silently disable #52 for a modification type.
    W, S = mm.IUPAC_BITS["A"] | mm.IUPAC_BITS["T"], mm.IUPAC_BITS["C"] | mm.IUPAC_BITS["G"]
    assert mm._allowed_center_bits("a") == mm._allowed_center_bits("6ma") == W
    assert (mm._allowed_center_bits("m") == mm._allowed_center_bits("5mc")
            == mm._allowed_center_bits("21839") == mm._allowed_center_bits("4mc")
            == mm._allowed_center_bits("h") == mm._allowed_center_bits("5hmc") == S)
    assert mm._allowed_center_bits("zzz") is None and mm._allowed_center_bits(None) is None
    assert mm._modifiable_base_bits("a") == mm._modifiable_base_bits("6ma") == mm.IUPAC_BITS["A"]
    assert (mm._modifiable_base_bits("m") == mm._modifiable_base_bits("21839")
            == mm._modifiable_base_bits("h") == mm.IUPAC_BITS["C"])
    assert mm._modifiable_base_bits("zzz") is None and mm._modifiable_base_bits(None) is None


def test_palindromize_rejects_widening_methyl_center():
    """#52: merging GACGGC with its RC GCCGTC gives GMCGKC (A/C at the 6mA).
    Without a methylation type palindromize still merges; with 6mA it rejects
    the merge and keeps the clean motif."""
    freq = np.array([mm._iupac_pwm_row(c) for c in "GACGGC"])
    m = mm.Motif(1, "GACGGC", 6, 1e-10, 1e-12, 100, freq, meth_index=1)
    (loose,) = mm.palindromize([m], allowed_center_bits=None)
    assert loose.iupac == "GMCGKC"
    m2 = mm.Motif(1, "GACGGC", 6, 1e-10, 1e-12, 100, freq, meth_index=1)
    (kept,) = mm.palindromize([m2], allowed_center_bits=mm._allowed_center_bits("a"))
    assert kept.iupac == "GACGGC"


def test_find_motifs_keeps_methyl_center_modifiable(tmp_path):
    """#52 end-to-end: without a methylation type the caller widens the methyl
    center (GACGGC -> GMCGKC); passing 6mA keeps that column an A."""
    pos_path, neg_path = _embed_motif_fastas(tmp_path, "GACGGC", seed=1)
    pos, neg = mm.load_fasta(pos_path), mm.load_fasta(neg_path)
    loose = [m.iupac for m in mm.find_motifs(pos, neg)]
    typed = [m.iupac for m in mm.find_motifs(
        pos, neg, allowed_center_bits=mm._allowed_center_bits("a"),
        modifiable_bits=mm._modifiable_base_bits("a"))]
    assert "GACGGC" in typed, typed          # methyl center kept as the A
    assert "GACGGC" not in loose, loose       # untyped widens it (the bug)


def test_dedup_results_constrains_methyl_center():
    """#52: the final dedup_results merge must not re-widen the methylated
    center past the modifiable base (mirrors the palindromize guard — dedup is
    the last transform before emit)."""
    def M(s, mi):
        return mm.Motif(1, s, len(s), 1e-10, 1e-12, 100,
                        np.array([mm._iupac_pwm_row(c) for c in s]), meth_index=mi)
    # merging GGGAGGG + GGGCGGG bit-unions the center (index 3) A -> M=[A/C]
    assert mm.dedup_results([M("GGGAGGG", 3), M("GGGCGGG", 3)], 2)[0].iupac == "GGGMGGG"
    kept = mm.dedup_results([M("GGGAGGG", 3), M("GGGCGGG", 3)], 2,
                            mm._allowed_center_bits("a"))
    assert kept[0].iupac == "GGGAGGG"
    assert kept[0].iupac[kept[0].meth_index] == "A"


def _embed_windows(fragments, n_pos=300, n_neg=4000, seed=0):
    """Build encoded pos/neg window arrays. `fragments` = list of
    (frag_or_list, offset); a list of fragments alternates across windows (to
    make a column genuinely mixed). The methylated base is expected at
    METH_CENTER by the caller, so place fragments accordingly."""
    rng = random.Random(seed)
    bases = "ACGT"

    def arr(strs):
        a = np.full((len(strs), mm.SEQ_LEN), 4, dtype=np.uint8)
        for i, x in enumerate(strs):
            a[i, : len(x)] = mm._ENC[np.frombuffer(x.encode(), np.uint8)]
        return a

    pos = []
    for k in range(n_pos):
        s = [rng.choice(bases) for _ in range(mm.SEQ_LEN)]
        for frag, off in fragments:
            f = frag[k % len(frag)] if isinstance(frag, list) else frag
            s[off : off + len(f)] = list(f)
        pos.append("".join(s))
    neg = ["".join(rng.choice(bases) for _ in range(mm.SEQ_LEN)) for _ in range(n_neg)]
    return arr(pos), arr(neg)


def _call_typed(pos, neg):
    return mm.find_motifs(
        pos, neg, bg=mm.compute_bg_from_seqs(neg),
        allowed_center_bits=mm._allowed_center_bits("a"),
        modifiable_bits=mm._modifiable_base_bits("a"))


def test_find_motifs_bipartite_meth_center_tracked():
    """#52: meth_index is tracked correctly through the BIPARTITE consensus
    (close half + N spacer + far half) — the every other find_motifs test uses
    short contiguous motifs, so the far/close/spacer win_pos math is otherwise
    unexercised. The constraint must land on the methyl-A in the close half."""
    pos, neg = _embed_windows([("TGACC", 11), ("GGTT", 21)], seed=0)  # methyl-A at 13
    bip = [m for m in _call_typed(pos, neg) if "N" in m.iupac.strip("N")]
    assert bip, [m.iupac for m in _call_typed(pos, neg)]
    m = bip[0]
    assert m.iupac == "TGACCNNNNNGGTT"
    assert m.meth_index == 2 and m.iupac[m.meth_index] == "A"


def test_find_motifs_narrows_mixed_center():
    """#52: a genuinely mixed methyl center (A in half the windows, C in the
    other -> consensus M) is narrowed to the modifiable base for 6mA; untyped
    keeps the degenerate M (confirms the column is really mixed)."""
    pos, neg = _embed_windows([(["GACGGC", "GCCGGC"], 12)], seed=0)  # index1 -> pos 13
    typed = [m.iupac for m in _call_typed(pos, neg)]
    untyped = [m.iupac for m in mm.find_motifs(pos, neg, bg=mm.compute_bg_from_seqs(neg))]
    assert "GACGGC" in typed and not any("M" in u for u in typed), typed
    assert "GACGGC" not in untyped and any("M" in u for u in untyped), untyped


# ── Orient emitted motifs to the modified base (issue #51 comment) ──────────
def test_orient_to_modifiable_flips_rc_emitted():
    # 6mA: modifiable base A. A motif emitted in RC orientation, with the methyl
    # position on the complement (T), is flipped so it reads as the A.
    A = mm.IUPAC_BITS["A"]
    freq = np.array([mm._iupac_pwm_row(c) for c in "ACCTGA"])
    m = mm.Motif(1, "ACCTGA", 6, 1e-10, 1e-12, 100, freq, meth_index=3)  # T @3
    (out,) = mm._orient_to_modifiable([m], A)
    assert out.iupac == "TCAGGT"
    assert out.meth_index == 2
    assert out.iupac[out.meth_index] == "A"
    # the freq array is RC'd too (reverse + A<->T, C<->G swap) so the emitted
    # PWM stays consistent with the flipped IUPAC (write_xml derives from IUPAC,
    # but keep freq correct for any consumer / round-trip)
    L = len(m.freq)
    expected = np.array([m.freq[L - 1 - i][[3, 2, 1, 0]] for i in range(L)])
    assert np.allclose(out.freq, expected)
    assert "".join("ACGT"[r.argmax()] for r in out.freq) == "TCAGGT"


def test_orient_to_modifiable_leaves_modifiable_center():
    A = mm.IUPAC_BITS["A"]
    freq = np.array([mm._iupac_pwm_row(c) for c in "GACGGC"])
    m = mm.Motif(1, "GACGGC", 6, 1e-10, 1e-12, 100, freq, meth_index=1)  # A @1
    (out,) = mm._orient_to_modifiable([m], A)
    assert out.iupac == "GACGGC"       # already shows the modified base
    assert out.meth_index == 1


def test_orient_to_modifiable_skips_degenerate_and_untracked():
    A = mm.IUPAC_BITS["A"]
    # W center (methylated on both strands) is ambiguous -> leave as-is
    mW = mm.Motif(1, "GWTGC", 5, 1e-10, 1e-12, 100,
                  np.array([mm._iupac_pwm_row(c) for c in "GWTGC"]), meth_index=1)
    # no tracked methyl center -> leave as-is
    mN = mm.Motif(2, "ACCTGA", 6, 1e-10, 1e-12, 100,
                  np.array([mm._iupac_pwm_row(c) for c in "ACCTGA"]), meth_index=None)
    outW, outN = mm._orient_to_modifiable([mW, mN], A)
    assert outW.iupac == "GWTGC"
    assert outN.iupac == "ACCTGA"


def test_find_motifs_orients_center_end_to_end():
    """#52 orientation wiring: find_motifs re-orients an emitted motif whose
    methyl center would read as the complement (T) so it reads as the modified
    base (A). This embed reads T at the center, so the caller emits ACCTGA
    (center T) natively and _orient_to_modifiable flips it to TCAGGT (center A);
    without that wiring the emitted motif keeps the T center."""
    pos, neg = _embed_windows([("ACCTGA", 10)], seed=0)  # center reads T at METH_CENTER
    iupacs = [m.iupac for m in _call_typed(pos, neg)]
    assert "TCAGGT" in iupacs and "ACCTGA" not in iupacs, iupacs
    m = next(m for m in _call_typed(pos, neg) if m.iupac == "TCAGGT")
    assert m.iupac[m.meth_index] == "A"


def test_run_from_fastas_writes_streme_xml(tmp_path):
    pos_path, neg_path = _embed_motif_fastas(tmp_path, "GATC", seed=4)
    out_dir = tmp_path / "out"
    result = mm.run_from_fastas(
        pos_path, neg_path, str(out_dir), output_type="xml"
    )
    assert result == str(out_dir)
    xml_path = os.path.join(str(out_dir), "streme.xml")
    assert os.path.isfile(xml_path)
    root = ET.parse(xml_path).getroot()
    ids = [m.get("id") for m in root.find("motifs").findall("motif")]
    assert any("GATC" in i for i in ids), ids


def test_run_from_fastas_applies_methylation_type(tmp_path):
    """#52 plumbing through the production entry point: run_from_fastas with a
    methylation type constrains the emitted center (GACGGC, not GMCGKC), and
    without it the degenerate center persists."""
    pos_path, neg_path = _embed_motif_fastas(tmp_path, "GACGGC", seed=5)

    def emitted(mod_type):
        out_dir = tmp_path / ("out_" + (mod_type or "none"))
        mm.run_from_fastas(pos_path, neg_path, str(out_dir),
                           output_type="xml", mod_type=mod_type)
        root = ET.parse(os.path.join(str(out_dir), "streme.xml")).getroot()
        return [m.get("id") for m in root.find("motifs").findall("motif")]

    typed, untyped = emitted("a"), emitted(None)
    assert any("GACGGC" in i for i in typed), typed
    assert not any("GMCGKC" in i for i in typed), typed
    assert any("GMCGKC" in i for i in untyped), untyped   # confirms the case is live


# ── refine_motifs meth_index shift, in isolation (issue #52) ────────────────
def test_refine_motifs_tracks_meth_center_through_left_prepend():
    """#52: refine_motifs carries the methyl-center index through a flank edit.

    Every other #52 test drives the center constraint through the full
    find_motifs → palindromize → dedup stack, all of which re-apply the center
    constraint — so a wrong +1/-1 shift *inside* refine_motifs is masked. Here we
    call refine_motifs directly: the input core "TGATM" has its methyl center at
    index 4 (the mixed A/C column); refine re-derives it, prepends a consensus
    "C" left flank (shifting the center to index 5) and extends "GG" on the
    right. The 6mA constraint must then land on the *shifted* column: typed
    narrows M->A at index 5 (CTGATAGG), untyped keeps M (CTGATMGG). If the +1
    left-prepend shift were dropped, the constraint would fall on index 4 (a
    fixed 'T', a no-op) and the M center would survive typed."""
    ENC = {"A": 0, "C": 1, "G": 2, "T": 3}
    n = 300
    pos = np.full((n, mm.SEQ_LEN), ENC["G"], dtype=np.uint8)   # 'G' filler
    for i in range(n):
        pos[i, 7] = i % 4                                       # uniform -> N (blocks 2nd prepend)
        pos[i, 8] = ENC["C"]                                    # consensus C -> left prepend (+1)
        pos[i, 9], pos[i, 10], pos[i, 11], pos[i, 12] = ENC["T"], ENC["G"], ENC["A"], ENC["T"]
        pos[i, 13] = ENC["A"] if i % 2 == 0 else ENC["C"]       # mixed -> M at METH_CENTER
    bg = np.array([0.25, 0.25, 0.25, 0.25])
    core = mm.Motif(1, "TGATM", 5, 1e-10, 1e-12, 100,
                    np.array([mm._iupac_pwm_row(c) for c in "TGATM"]), meth_index=4)

    untyped = mm.refine_motifs([core], pos, bg, allowed_center_bits=None)[0]
    assert untyped.iupac == "CTGATMGG" and untyped.meth_index == 5
    assert untyped.iupac[untyped.meth_index] == "M"            # center un-narrowed

    typed = mm.refine_motifs([core], pos, bg,
                             allowed_center_bits=mm._allowed_center_bits("a"))[0]
    assert typed.iupac == "CTGATAGG" and typed.meth_index == 5
    assert typed.iupac[typed.meth_index] == "A"                # narrowed on the shifted column


# ── subcommand + CLI plumbing (issues #52 / #50) ────────────────────────────
def _write_bed_and_fasta(tmp_path, motif="GACGGC", meth_off=1, code="a", n=40):
    """Synthetic reference + 18-column modkit bedmethyl with `n` methylated
    sites (one per embedded motif occurrence). Returns (bed_path, fasta_path)."""
    gap, start0 = 25, 50
    length = start0 + n * gap + 50
    genome = list("G" * length)
    positions = []
    for i in range(n):
        ms = start0 + i * gap
        genome[ms:ms + len(motif)] = list(motif)
        positions.append(ms + meth_off)                         # methylated base
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">contig1\n" + "".join(genome) + "\n")
    bed = tmp_path / "calls.bed"
    with open(bed, "w") as f:
        for p in positions:
            cols = ["contig1", str(p), str(p + 1), code, "1000", "+",
                    str(p), str(p + 1), "255,0,0", "20", "90.0", "18",
                    "2", "0", "0", "0", "0", "0"]          # 18 cols; cov=20, frac=90%
            f.write("\t".join(cols) + "\n")
    return str(bed), str(fasta)


def test_subcommand_main_threads_methylation_constraint(tmp_path, monkeypatch):
    """#52 plumbing through the stand-alone subcommand: subcommand_main must
    pass the methylation type's center constraint (allowed_center_bits +
    modifiable_bits) into find_motifs. Capture the call rather than the emit so
    the assertion pins the plumbing itself; a dropped kwarg would capture None."""
    bed, fasta = _write_bed_and_fasta(tmp_path)
    captured = {}

    def fake_find_motifs(pos_seqs, neg_seqs, bg=None,
                         allowed_center_bits=None, modifiable_bits=None):
        captured["n_pos"] = int(pos_seqs.shape[0])
        captured["allowed"] = allowed_center_bits
        captured["modifiable"] = modifiable_bits
        return []

    monkeypatch.setattr(mm, "find_motifs", fake_find_motifs)
    monkeypatch.chdir(tmp_path)
    mm.subcommand_main(bed, fasta, output_type="tsv",
                       output_prefix="probe", methylation_types="6mA")

    assert captured.get("n_pos", 0) >= mm.MIN_MOTIF_SITES     # find_motifs actually reached
    assert captured["allowed"] == mm._allowed_center_bits("a")
    assert captured["modifiable"] == mm._modifiable_base_bits("a")
    assert (tmp_path / "probe_motifs.tsv").exists()


def test_main_cli_parses_args_and_writes_xml(tmp_path):
    """#50 plumbing: the stand-alone main() CLI parses -p/--n/-o, runs the
    caller, and writes streme.xml with the recovered motif."""
    pos_path, neg_path = _embed_motif_fastas(tmp_path, "GATC", seed=1)
    out_dir = tmp_path / "cli_out"
    old_argv = sys.argv
    sys.argv = ["motif_caller.py", ".", "-p", pos_path, "--n", neg_path,
                "-o", str(out_dir)]
    try:
        mm.main()
    finally:
        sys.argv = old_argv

    xml_path = out_dir / "streme.xml"
    assert xml_path.exists()
    root = ET.parse(str(xml_path)).getroot()
    ids = [m.get("id") for m in root.find("motifs").findall("motif")]
    assert any("GATC" in i for i in ids), ids
