"""Regression tests for the built-in (default) motif caller, microbe_motif.

These cover the pieces that the legacy STREME path never exercised: the IUPAC
consensus logic, the Fisher exact p-value, k-mer counting, the STREME-format
XML round-trip, and an end-to-end synthetic recovery of a known motif.

The end-to-end tests pin opt_max_seconds=0 so the wall-clock optimizer
is disabled and the output is fully deterministic for assertion.
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
    motifs = mm.find_motifs(pos, neg, opt_max_seconds=0)
    called = [m.iupac for m in motifs]
    assert any("GATC" in c for c in called), called


def test_find_motifs_recovers_longer_palindrome(tmp_path):
    # GANTC (HinfI) — palindromic with a degenerate center
    pos_path, neg_path = _embed_motif_fastas(tmp_path, "GAATTC", seed=2)
    pos = mm.load_fasta(pos_path)
    neg = mm.load_fasta(neg_path)
    motifs = mm.find_motifs(pos, neg, opt_max_seconds=0)
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
        mm.load_fasta(str(pos_path)), mm.load_fasta(str(neg_path)), opt_max_seconds=0
    )
    assert motifs == []


def test_run_from_fastas_writes_streme_xml(tmp_path):
    pos_path, neg_path = _embed_motif_fastas(tmp_path, "GATC", seed=4)
    out_dir = tmp_path / "out"
    result = mm.run_from_fastas(
        pos_path, neg_path, str(out_dir), output_type="xml", opt_max_seconds=0
    )
    assert result == str(out_dir)
    xml_path = os.path.join(str(out_dir), "streme.xml")
    assert os.path.isfile(xml_path)
    root = ET.parse(xml_path).getroot()
    ids = [m.get("id") for m in root.find("motifs").findall("motif")]
    assert any("GATC" in i for i in ids), ids


# ── Optimizer (wall-clock budget; warns when time-truncated) ────────────────
def test_optimize_motifs_respects_zero_seconds(tmp_path):
    # max_seconds<=0 disables the optimizer (used by the end-to-end tests). The
    # input motif set is returned untouched and nothing is tried.
    pos_path, neg_path = _embed_motif_fastas(tmp_path, "GATC", seed=6)
    pos = mm.load_fasta(pos_path)
    neg = mm.load_fasta(neg_path)
    freq = np.array([mm._iupac_pwm_row(c) for c in "GANC"])
    start = [mm.Motif(1, "GANC", 4, 1e-5, 1e-6, 100, freq)]
    out, info = mm.optimize_motifs(start, pos, neg, max_seconds=0)
    assert out is start
    assert info["perturbations_tried"] == 0
    assert info["converged"] is True


def test_optimize_motifs_warns_when_time_truncated(tmp_path, capsys, monkeypatch):
    # When the hill-climb exhausts its time budget before converging, the run is
    # machine-/load-dependent: it must flag converged=False and warn on stderr,
    # never produce a non-reproducible result silently. A fake monotonic clock
    # makes "time runs out mid-search" deterministic, independent of real speed.
    pos_path, neg_path = _embed_motif_fastas(tmp_path, "GATC", seed=7)
    pos = mm.load_fasta(pos_path)
    neg = mm.load_fasta(neg_path)
    # A deliberately too-broad seed gives the climb real work to do.
    freq = np.array([mm._iupac_pwm_row(c) for c in "GANC"])
    start = [mm.Motif(1, "GANC", 4, 1e-5, 1e-6, 100, freq)]

    # Clock reads 0.0 first (t0), then advances 1s per call. With max_seconds=5
    # the search enters, tries a few perturbations, then crosses the budget.
    ticks = iter([0.0] + [float(n) for n in range(1, 100000)])
    monkeypatch.setattr(mm.time, "time", lambda: next(ticks))

    _out, info = mm.optimize_motifs(start, pos, neg, max_seconds=5)
    assert info["perturbations_tried"] > 0
    assert info["converged"] is False
    assert "WARNING" in capsys.readouterr().err
