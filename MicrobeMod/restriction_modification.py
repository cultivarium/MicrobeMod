import os
import sys
import logging
import subprocess
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)-8s %(message)s",
    datefmt="%y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)


def run_prodigal(output_prefix, input_fasta):
    """Runs prodigal on a genome.
    Args:
            output_prefix: prefix of output file.
            input_fasta: path to FASTA file.

    Returns:
            The name of the prodigal FAA file generated.
            gene_locations: a dictionary of genes to (contig, gene_num)
    """

    cmd = "prodigal -i {} -a {}.faa"
    cmd = cmd.format(input_fasta, output_prefix)

    if os.path.isfile(input_fasta):
        logging.info("Running prodigal: %s", cmd)

        subprocess.run(
            cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True, check=True
        )

        gene_locations = {}  # values are (contig, gene_num)
        for record in SeqIO.parse(output_prefix + ".faa", "fasta"):
            gene_locations[record.id] = (
                "_".join(record.id.split("_")[:-1]),
                int(record.id.split("_")[-1]),
            )
        return output_prefix + ".faa", gene_locations

    else:
        logging.error("ERROR: FASTA file doesn't exist")
        raise SystemExit("Cannot find FASTA file- quitting....")


def convert_genbank(output_prefix, input_genbank):
    """Converts a Genbank file to resemble a .FAA file from prodigal.
    Args:
            output_prefix: prefix of output file.
            input_genbank: path to Genbank file.

    Returns:
            The name of the FAA file converted.
            gene_locations: a dictionary of genes to (contig, gene_num)
    """

    if os.path.isfile(input_genbank):
        logging.info("Reading genbank file: %s", input_genbank)

        f = open(output_prefix + ".faa", "w+",  encoding="UTF-8")
        gene_locations = {}  # values are (contig, gene_num)
        for record in SeqIO.parse(input_genbank, "genbank"):
            i = 1
            for feature in record.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    f.write(
                        ">{} # {} # {} # {} # {}_{}\n".format(
                            feature.qualifiers["locus_tag"][0],
                            feature.location.start,
                            feature.location.end,
                            feature.location.strand,
                            record.id,
                            i
                        )
                    )
                    f.write(feature.qualifiers["translation"][0] + "\n")
                    gene_locations[feature.qualifiers["locus_tag"][0]] = (record.id, i)
                    i += 1

        return output_prefix + ".faa", gene_locations

    else:
        logging.error("ERROR: FASTA file doesn't exist")
        raise SystemExit("Cannot find FASTA file- quitting....")


def run_hmmer(output_prefix, prodigal_fasta, threads):
    """Runs HMMER to identify RM genes on a set of protein predictions.
    Args:
            output_prefix: prefix of output file.
            prodigal_fasta: path to prodigal FAA file.
            threads: number of threads to pass to HMMER

    Returns:
            The name of the HMMER file created.
    """

    cmd = "hmmsearch --cut_ga --cpu {} --domtblout {}.hits {} {}"

    rm_hmm_file = os.path.dirname(__file__) + "/db/HMMs/RM_HMMs.hmm"

    cmd = cmd.format(threads, output_prefix, rm_hmm_file, prodigal_fasta)

    if os.path.isfile(prodigal_fasta):
        logging.info("Running HMMER: %s", cmd)
        subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT
        )
        return "{}.hits".format(output_prefix)
    else:
        logging.error("ERROR: Prodigal output file not found")
        raise SystemExit("Cannot find FAA file- quitting....")


def extract_genes(hits, prodigal_fasta, output_prefix):
    """Creates a FASTA file of genes involved in RM.
    Args:
            hits: a dictionary of genes to their corresponding HMM hits.
            prodigal_fasta: path to the FAA file created by prodigal.
            output_prefix: prefix of output file.

    Returns:
            Name of RM protein sequence file created.
    """

    rm_gene_file = output_prefix + ".rm.genes.faa"
    f = open(rm_gene_file, "w+",  encoding="UTF-8")

    for record in SeqIO.parse(prodigal_fasta, "fasta"):
        if record.id in hits:
            f.write(">" + record.description + "\n")
            f.write(str(record.seq) + "\n")
    f.close()

    return rm_gene_file


def resolve_hits(hmmer_output, output_prefix):
    """Runths cath-resolve-hits on the output of HMMER. This resolves
    the best non-overlapping set of HMMs for each gene.

    Args:
            hmmer_output: the path of the file output by HMMER
            output_prefix: prefix of output file.

    Returns:
            Name of cath resolve hits file created.
    """

    cmd = "cath-resolve-hits --input-format hmmer_domtblout {} --hits-text-to-file {}.resolved.hits"

    cmd = cmd.format(hmmer_output, output_prefix)
    logging.info("Running cath: %s", cmd)
    subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    return "{}.resolved.hits".format(output_prefix)


def parse_hmmer(resolved_hits):
    """Reads the output of cath resolve hits and converts it into dictionaries
    linking each gene to information about that gene.

    Args:
            resolved_hits: the path of the file output by cath-resolve-hits

    Returns:
            hits: a defaultdictionary where keys are genes and values are lists of their HMM names.
            evalues: a defaultdictionary where keys are genes, values are dictionaries with HMMs as keys and values as evalues.
    """

    offtarget_file = os.path.dirname(__file__) + "/db/HMMs/off_target.txt"
    offtarget = []
    f = open(offtarget_file,  encoding="UTF-8")
    for line in f.readlines():
        offtarget.append(line.strip())
    f.close()

    f = open(resolved_hits,  encoding="UTF-8")
    hits = defaultdict(list)
    evalues = defaultdict(dict)
    for line in f.readlines():
        if not line.startswith("#"):
            gene = line.split()[0]
            hmm = line.split()[1]
            evalue = float(line.split()[-1])

            if hmm not in offtarget:
                hits[gene].append(hmm)
                evalues[gene][hmm] = evalue
    f.close()
    return hits, evalues


def create_gene_table(
    hits, gene_locations, system_types, evalues, blast_hits, gene_window=10
):
    """Creates a final table of RM genes organized by their types and their operon status.
    Args:
            hits: a defaultdictionary where keys are genes and values are lists of their HMM names.
            gene_locations: a dictionary where keys are genes and values are tuples of (contig, gene number).
            system_types: Dictionary of HMMs to RM system types
            evalues: a defaultdictionary where keys are genes, values are dictionaries with HMMs as keys and values as evalues.
            blast_hits: dictionary of genes-> blast hits, from read_blast()
            gene_window: the window size to use to call operons.
    Returns:
            gene_table: A final pandas table of RM genes organized by type and operon.

    """

    gene_table = []

    for hit in hits:
        ## Get BLAST info
        best_blast = blast_pid = meth_type = motif = ""
        if hit in blast_hits:
            best_blast = blast_hits[hit][0]
            blast_pid = blast_hits[hit][1]
            meth_type = blast_hits[hit][2]
            motif = blast_hits[hit][3]

        gene_types = set()  # So that each HMM MT/RE hit only counts once per gene

        for hmm in hits[hit]:
            ## Get the system type of this HMM
            if system_types[hmm][0] not in gene_types:
                gene_types.add(system_types[hmm][0])
                gene_table.append(
                    {
                        "Gene": hit,
                        "Contig": gene_locations[hit][0],
                        "Gene Position": gene_locations[hit][1],
                        "System Type": system_types[hmm][1],
                        "Gene type": system_types[hmm][0],
                        "HMM": hmm,
                        "Evalue": evalues[hit][hmm],
                        "REBASE homolog": best_blast,
                        "Homolog identity(%)": blast_pid,
                        "Homolog methylation": meth_type,
                        "Homolog motif": motif,
                    }
                )

    gene_table = pd.DataFrame(gene_table)
    gene_table = gene_table.sort_values(["Contig", "Gene Position"], ascending=True)

    ## Assign operons
    gene2operon = {}
    operon2gene = defaultdict(list)
    operon_number = 0
    pos = 0
    contig = ""

    for _, row in gene_table.iterrows():
        ## New contig = new operon
        if row["Contig"] != contig:
            operon_number += 1
            contig = row["Contig"]
            pos = row["Gene Position"]
            gene2operon[row["Gene"]] = operon_number
            operon2gene[operon_number].append(
                row["System Type"] + ":" + row["Gene type"]
            )

        # Same contig, gene within window = same operon
        # Update position (extend operons)
        elif row["Contig"] == contig and row["Gene Position"] <= pos + gene_window:
            pos = row["Gene Position"]
            gene2operon[row["Gene"]] = operon_number
            operon2gene[operon_number].append(
                row["System Type"] + ":" + row["Gene type"]
            )

        # Same contig, gene outside of window = new operon
        else:
            operon_number += 1
            contig = row["Contig"]
            pos = row["Gene Position"]
            gene2operon[row["Gene"]] = operon_number
            operon2gene[operon_number].append(
                row["System Type"] + ":" + row["Gene type"]
            )

    ## Only operons with 1 MT, 1 RE are kept - others are singletons
    ## Keep Type IVs

    operon_number = 0
    singleton_number = 0
    operon_labels = {}

    for _, row in gene_table.iterrows():
        found_MT = False
        found_RE = False
        found_RE_IV = False
        operon = gene2operon[row["Gene"]]

        if operon not in operon_labels:
            for genetype in operon2gene[operon]:
                if "MT" in genetype:
                    found_MT = True
                if "RE" in genetype:
                    found_RE = True
                if "RM_Type_IV" in genetype or "IIG" in genetype:
                    found_RE_IV = True

            if found_RE_IV or (found_MT and found_RE):
                operon_number += 1
                operon_labels[operon] = "RM Operon #" + str(operon_number)
            else:
                singleton_number += 1
                operon_labels[operon] = "Singleton #" + str(singleton_number)

    gene_table.insert(
        0, "Operon", gene_table.Gene.map(gene2operon).map(operon_labels)
    )
    gene_table = gene_table.sort_values(["Operon", "REBASE homolog", "Gene Position"])
    del gene_table["Gene Position"]
    del gene_table["Contig"]

    return gene_table


def run_rebase_blast(rm_gene_file, output_prefix, threads):
    """Run a BLAST of RM gene files against REBASE
    Args:
            rm_gene_file: a defaultdictionary where keys are genes and values are lists of their HMM names.
            output_prefix: a dictionary where keys are genes and values are tuples of (contig, gene number).
    Returns:
            file path of BLAST output.
    """

    cmd = "blastp -query {} -db {} -outfmt 6 -evalue 1e-5 -num_threads {} > {}.blast"

    blast_db_file = (
        os.path.dirname(__file__) + "/db/rebase_blast/all_rebase_proteins.faa"
    )

    cmd = cmd.format(rm_gene_file, blast_db_file, threads, output_prefix)

    logging.info("Running BLASTP against REBASE: %s", cmd)
    subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)

    return "{}.blast".format(output_prefix)


def read_blast(blast_file):
    """Read the output of BLAST
    Args:
            blast_file: path to a file output by BLAST against REBASE database.
    Returns:
            blast_hits: dictionary of BLAST hits
    """

    f = open(blast_file,  encoding="UTF-8")
    blast_hits = {}
    prots = set()

    for line in f.readlines():
        hit = line.split()[0]
        pid = float(line.split()[2])

        if hit not in prots:  # top hit only
            prots.add(hit)
            if pid > 50:
                best_hit = line.split("\t")[1].split(":")[1].split("-")[0]
                meth_type = line.split("\t")[1].split("-")[-2]
                motif = line.split("\t")[1].split("-")[-1]

                blast_hits[hit] = (best_hit, pid, meth_type, motif)

    return blast_hits


def main(fasta, genbank, output_prefix, threads):
    metadata_file = os.path.dirname(__file__) + "/db/restriction_metadata.csv"
    metadata = pd.read_csv(metadata_file)
    system_types = {}
    for _, row in metadata.iterrows():
        system_types[row["Name"]] = (row["Enzyme_type"], row["System"])

    if genbank:
        logging.info("Using Genbank file: %s", genbank)
        if (
            genbank.split(".")[-1] != "gbk"
            and genbank.split(".")[-1] != "gb"
            and genbank.split(".")[-1] != "gbff"
        ):
            logging.info(
                "WARNING: is your -g genbank file really a genbank file? It does not end in .gb, .gbk, or .gbff"
            )
        prodigal_fasta, gene_locations = convert_genbank(output_prefix, genbank)
    elif fasta:
        if (
            fasta.split(".")[-1] != "fasta"
            and fasta.split(".")[-1] != "fna"
            and fasta.split(".")[-1] != "fa"
        ):
            logging.info(
                "WARNING: is your -f fasta file really a genomic FASTA file? It does not end in .fa, .fna, or .fasta."
            )
        logging.info("Calling prodigal on FASTA file: %s", fasta)
        prodigal_fasta, gene_locations = run_prodigal(output_prefix, fasta)
    else:
        raise SystemExit(
            "Error: You must supply either a WGS FASTA file with --fasta or a Genbank file with --genbank."
        )

    hmmer_output = run_hmmer(output_prefix, prodigal_fasta, threads)

    resolved_hits = resolve_hits(hmmer_output, output_prefix)
    gene_hits, evalues = parse_hmmer(resolved_hits)

    rm_gene_file = extract_genes(gene_hits, prodigal_fasta, output_prefix)

    blast_file = run_rebase_blast(rm_gene_file, output_prefix, threads)
    blast_hits = read_blast(blast_file)

    gene_table = create_gene_table(
        gene_hits, gene_locations, system_types, evalues, blast_hits
    )
    gene_table.to_csv(output_prefix + ".rm.genes.tsv", sep="\t", index=False)
