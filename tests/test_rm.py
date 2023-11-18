import pytest
import os
import pandas as pd
import warnings
import subprocess
from Bio import SeqIO
from MicrobeMod.restriction_modification import parse_hmmer
from MicrobeMod.restriction_modification import read_blast
from MicrobeMod.restriction_modification import create_gene_table

from pandas.testing import assert_frame_equal


def test_parse_hmmer():
    gene_hits, gene_locations, evalues = parse_hmmer(
        "./tests/test_data/EcoliCVM05_GCF_000005845.resolved.hits"
    )
    expected = {
        "NC_000913.3_1130": ["Type_IV_05-RM_Type_IV__Type_IV_REases"],
        "NC_000913.3_1932": ["Type_II_MTases_FAM_2"],
        "NC_000913.3_2152": ["Type_I_REases_FAM_2.einsi_trimmed"],
        "NC_000913.3_3194": ["Type_II_MTases_FAM_1"],
        "NC_000913.3_3316": ["Type_II_MTases_FAM_4"],
        "NC_000913.3_4262": ["FAM_2-RM_Type_IV__Type_IV_REases"],
        "NC_000913.3_4263": ["FAM_1-RM_Type_IV__Type_IV_REases"],
        "NC_000913.3_4265": ["Type_I_S_52"],
        "NC_000913.3_4266": ["Type_I_MTases_FAM_2"],
        "NC_000913.3_4267": ["Type_I_REases_FAM_2.einsi_trimmed"],
        "NC_000913.3_4268": [
            "FAM_0-RM_Type_IV__Type_IV_REases",
            "FAM_0-RM_Type_IV__Type_IV_REases",
        ],
        "NC_000913.3_464": ["Type_II_REase06"],
    }
    expected_locations = {
        "NC_000913.3_1130": ("NC_000913.3", 1130),
        "NC_000913.3_1932": ("NC_000913.3", 1932),
        "NC_000913.3_2152": ("NC_000913.3", 2152),
        "NC_000913.3_3194": ("NC_000913.3", 3194),
        "NC_000913.3_3316": ("NC_000913.3", 3316),
        "NC_000913.3_4262": ("NC_000913.3", 4262),
        "NC_000913.3_4263": ("NC_000913.3", 4263),
        "NC_000913.3_4265": ("NC_000913.3", 4265),
        "NC_000913.3_4266": ("NC_000913.3", 4266),
        "NC_000913.3_4267": ("NC_000913.3", 4267),
        "NC_000913.3_4268": ("NC_000913.3", 4268),
        "NC_000913.3_464": ("NC_000913.3", 464),
    }
    expected_evalues = {
        "NC_000913.3_1130": {"Type_IV_05-RM_Type_IV__Type_IV_REases": 6.8e-27},
        "NC_000913.3_1932": {"Type_II_MTases_FAM_2": 1.2e-118},
        "NC_000913.3_2152": {"Type_I_REases_FAM_2.einsi_trimmed": 5.1e-31},
        "NC_000913.3_3194": {"Type_II_MTases_FAM_1": 5.1e-93},
        "NC_000913.3_3316": {"Type_II_MTases_FAM_4": 1.4e-87},
        "NC_000913.3_4262": {"FAM_2-RM_Type_IV__Type_IV_REases": 8.6e-137},
        "NC_000913.3_4263": {"FAM_1-RM_Type_IV__Type_IV_REases": 4e-149},
        "NC_000913.3_4265": {"Type_I_S_52": 1.8e-81},
        "NC_000913.3_4266": {"Type_I_MTases_FAM_2": 4.7e-252},
        "NC_000913.3_4267": {"Type_I_REases_FAM_2.einsi_trimmed": 0.0},
        "NC_000913.3_4268": {"FAM_0-RM_Type_IV__Type_IV_REases": 2.1e-76},
        "NC_000913.3_464": {"Type_II_REase06": 7.9e-10},
    }
    assert (
        (gene_hits == expected)
        and (gene_locations == expected_locations)
        and (evalues == expected_evalues)
    )


def test_read_blast():
    blast_hits = read_blast("./tests/test_data/EcoliCVM05_GCF_000005845.blast")
    expected = {
        "NC_000913.3_1130": ("SsoSE61ORF22640P", 100.0, "", "YCGR"),
        "NC_000913.3_1932": ("M.SflLIN6DcmP", 100.0, "m5C", "CCWGG"),
        "NC_000913.3_3194": ("M.Eco4792LORF2734P", 100.0, "", "ATGCAT"),
        "NC_000913.3_3316": ("M.UbaC1152DamP", 100.0, "m6A", "GATC"),
        "NC_000913.3_4262": ("SenHNK130McrBCP_(SenHNK130McrCP)", 100.0, "", ""),
        "NC_000913.3_4263": ("Eco1655dMcrBCP_(Eco1655dMcrBP)", 100.0, "", ""),
        "NC_000913.3_4265": ("S.SenHNK130ORF17125P", 100.0, "", "AACNNNNNNGTGC"),
        "NC_000913.3_4266": ("M.SenHNK130ORF17125P", 100.0, "m6A", "AACNNNNNNGTGC"),
        "NC_000913.3_4267": ("Msa17ORFC2P", 100.0, "", "AACNNNNNNGTGC"),
        "NC_000913.3_4268": ("EcoZK126MrrP", 100.0, "", ""),
    }

    assert blast_hits == expected


def test_create_gene_table():
    metadata_file = (
        "/"
        + os.path.dirname(__file__)[:-6]
        + "/MicrobeMod/db/restriction_metadata.csv"
    )
    metadata = pd.read_csv(metadata_file)
    system_types = {}
    for index, row in metadata.iterrows():
        system_types[row["Name"]] = (row["Enzyme_type"], row["System"])

    gene_hits, gene_locations, evalues = parse_hmmer(
        "./tests/test_data/EcoliCVM05_GCF_000005845.resolved.hits"
    )
    blast_hits = read_blast("./tests/test_data/EcoliCVM05_GCF_000005845.blast")

    gene_table = create_gene_table(
        gene_hits, gene_locations, system_types, evalues, blast_hits
    )

    expected_json = '{"Operon":{"0":"RM Operon #1","6":"RM Operon #2","10":"RM Operon #2","8":"RM Operon #2","9":"RM Operon #2","7":"RM Operon #2","5":"RM Operon #2","11":"Singleton #1","1":"Singleton #2","2":"Singleton #3","3":"Singleton #4","4":"Singleton #5"},"Gene":{"0":"NC_000913.3_1130","6":"NC_000913.3_4263","10":"NC_000913.3_4268","8":"NC_000913.3_4266","9":"NC_000913.3_4267","7":"NC_000913.3_4265","5":"NC_000913.3_4262","11":"NC_000913.3_464","1":"NC_000913.3_1932","2":"NC_000913.3_2152","3":"NC_000913.3_3194","4":"NC_000913.3_3316"},"System Type":{"0":"RM_Type_IV","6":"RM_Type_IV","10":"RM_Type_IV","8":"RM_Type_I","9":"RM_Type_I","7":"RM_Type_I","5":"RM_Type_IV","11":"RM_Type_II","1":"RM_Type_II","2":"RM_Type_I","3":"RM_Type_II","4":"RM_Type_II"},"Gene type":{"0":"RE","6":"RE","10":"RE","8":"MT","9":"RE","7":"SP","5":"RE","11":"RE","1":"MT","2":"RE","3":"MT","4":"MT"},"HMM":{"0":"Type_IV_05-RM_Type_IV__Type_IV_REases","6":"FAM_1-RM_Type_IV__Type_IV_REases","10":"FAM_0-RM_Type_IV__Type_IV_REases","8":"Type_I_MTases_FAM_2","9":"Type_I_REases_FAM_2.einsi_trimmed","7":"Type_I_S_52","5":"FAM_2-RM_Type_IV__Type_IV_REases","11":"Type_II_REase06","1":"Type_II_MTases_FAM_2","2":"Type_I_REases_FAM_2.einsi_trimmed","3":"Type_II_MTases_FAM_1","4":"Type_II_MTases_FAM_4"},"Evalue":{"0":6.8e-27,"6":4e-149,"10":2.1e-76,"8":4.7e-252,"9":0.0,"7":1.8e-81,"5":8.6e-137,"11":0.0000000008,"1":1.2e-118,"2":5.1e-31,"3":5.1e-93,"4":1.4e-87},"REBASE homolog":{"0":"SsoSE61ORF22640P","6":"Eco1655dMcrBCP_(Eco1655dMcrBP)","10":"EcoZK126MrrP","8":"M.SenHNK130ORF17125P","9":"Msa17ORFC2P","7":"S.SenHNK130ORF17125P","5":"SenHNK130McrBCP_(SenHNK130McrCP)","11":"","1":"M.SflLIN6DcmP","2":"","3":"M.Eco4792LORF2734P","4":"M.UbaC1152DamP"},"Homolog identity(%)":{"0":100.0,"6":100.0,"10":100.0,"8":100.0,"9":100.0,"7":100.0,"5":100.0,"11":"","1":100.0,"2":"","3":100.0,"4":100.0},"Homolog methylation":{"0":"","6":"","10":"","8":"m6A","9":"","7":"","5":"","11":"","1":"m5C","2":"","3":"","4":"m6A"},"Homolog motif":{"0":"YCGR","6":"","10":"","8":"AACNNNNNNGTGC","9":"AACNNNNNNGTGC","7":"AACNNNNNNGTGC","5":"","11":"","1":"CCWGG","2":"","3":"ATGCAT","4":"GATC"}}'

    assert gene_table.to_json() == expected_json
