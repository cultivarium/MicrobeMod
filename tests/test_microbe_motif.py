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
