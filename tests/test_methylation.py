import pytest
import numpy as np
import pandas as pd
import warnings
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from MicrobeMod import microbemod
from MicrobeMod import microbe_motif as mm
from MicrobeMod.microbemod import read_modkit
from MicrobeMod.microbemod import write_to_fasta
from MicrobeMod.microbemod import assign_motifs
from MicrobeMod.microbemod import make_motif_table

from pandas.testing import assert_frame_equal


def test_read_modkit():
    min_coverage = 10

    test_fasta = "./tests/test_data/EcoliCVM05_GCF_000005845.2_ASM584v2_genomic.fna"
    for record in SeqIO.parse(test_fasta, "fasta"):
        microbemod.REF[record.id] = record.seq
    r = read_modkit("./tests/test_data/test.bed", min_coverage).round(2)
    print(r)
    read_modkit_result = pd.read_csv("./tests/test_data/read_modkit_result.csv")
    assert r.equals(read_modkit_result)


def test_write_to_fasta():
    percent_cutoff_streme = 0.9
    min_coverage = 10
    methylation_type = "a"
    output_prefix = "./tests/test_data/test"

    modkit_table_tmp = pd.read_csv("./tests/test_data/read_modkit_result.csv")

    fasta_out = write_to_fasta(
        modkit_table_tmp,
        output_prefix,
        methylation_type,
        percent_cutoff_streme,
        min_coverage,
    )

    # Read in the FASTA that we just wrote
    ref_seqs = [
        str(record.seq)
        for record in SeqIO.parse("./tests/test_data/test_fasta.fasta", "fasta")
    ]
    new_seqs = [str(record.seq) for record in SeqIO.parse(fasta_out, "fasta")]

    assert new_seqs == ref_seqs


def test_assign_motifs():
    modkit_table = pd.read_csv("./tests/test_data/read_modkit_result.csv")
    streme_output = "./tests/test_data/test_streme_output/"

    test_fasta = "./tests/test_data/EcoliCVM05_GCF_000005845.2_ASM584v2_genomic.fna"
    for record in SeqIO.parse(test_fasta, "fasta"):
        microbemod.REF[record.id] = record.seq

    result = assign_motifs(modkit_table, streme_output)

    expected = pd.Series(index=["GATC", "CCWGG"], data=[131, 94]).to_dict()

    assert result.motif.value_counts().to_dict() == expected


def test_assign_motifs_both_strands_nonpalindromic(tmp_path):
    """Issue #51: the partner-strand methyl site of a NON-palindromic motif is
    assigned to that motif rather than dropped into "No Motif Assigned"; a
    palindrome is assigned on both strands either way; and make_motif_table
    coverage stays <= 1.0 (the both-strand keying does not inflate it)."""
    # Only-C filler so GATGC / GATC / their RCs each occur exactly once.
    genome = "C" * 50 + "GATGC" + "C" * 50 + "GATC" + "C" * 50
    microbemod.REF.clear()
    microbemod.REF["c"] = Seq(genome)

    # GATGC: methyl-A at +1 on '+', partner methyl-A at +2 on '-'.
    # GATC:  methyl-A at +1 on '+', partner at +2 on '-'. Plus a control site.
    rows = [("c:51", "+"), ("c:52", "-"),      # non-palindromic, both strands
            ("c:106", "+"), ("c:107", "-"),    # palindrome, both strands
            ("c:10", "+")]                     # control (no motif here)
    df = pd.DataFrame({
        "SNP_Position": [r[0] for r in rows],
        "Strand": [r[1] for r in rows],
        "Modification": ["a"] * len(rows),
        "Percent_modified": [0.95] * len(rows),
    })

    # streme.xml with exactly the two planted motifs (forward orientation)
    motifs = [mm.Motif(1, "GATGC", 5, 1e-20, 1e-22, 100,
                       np.array([mm._iupac_pwm_row(c) for c in "GATGC"])),
              mm.Motif(2, "GATC", 4, 1e-20, 1e-22, 100,
                       np.array([mm._iupac_pwm_row(c) for c in "GATC"]))]
    mm.write_xml(str(tmp_path / "streme.xml"), motifs, 100, 1000)

    result = assign_motifs(df, str(tmp_path))
    assigned = dict(zip(result.SNP_Position + result.Strand, result.motif))

    assert assigned["c:51+"] == "GATGC"
    assert assigned["c:52-"] == "GATGC"     # partner strand — NaN before the fix
    assert assigned["c:106+"] == "GATC"
    assert assigned["c:107-"] == "GATC"
    assert pd.isna(assigned["c:10+"])       # control stays unassigned

    # coverage is recomputed independently of the motif column and stays <= 1.0
    tbl = make_motif_table(result)
    cov = pd.to_numeric(
        tbl[tbl.Motif != "No Motif Assigned"].Methylation_coverage,
        errors="coerce")
    assert (cov <= 1.0).all()
