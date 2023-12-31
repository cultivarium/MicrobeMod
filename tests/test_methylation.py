import pytest
import pandas as pd
import warnings
import subprocess
from Bio import SeqIO
from MicrobeMod import microbemod
from MicrobeMod.microbemod import read_modkit
from MicrobeMod.microbemod import write_to_fasta
from MicrobeMod.microbemod import assign_motifs

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
