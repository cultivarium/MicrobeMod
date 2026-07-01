import logging
import os
import random
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)-8s %(message)s",
    datefmt="%y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)

METHYLATION_TYPES = {"6mA": "a", "5mC": "m", "4mC": "21839", "5hmC": "h"}
METHYLATION_TYPES_REV = {"a": "6mA", "m": "5mC", "21839": "4mC", "h": "5hmC"}
REF = {}
# 13 -> get_seq yields 26 bp windows with the methylated base at index 13,
# matching the python caller's SEQ_LEN=26 / METH_CENTER=13 convention (and the
# motif_caller benchmark). STREME is unaffected by the extra 2 bp of flank.
WINDOW_SIZE = 13
# Fixed seed for the genomic background control sampling so motif-calling
# results are reproducible run-to-run. Matches the standalone microbe_motif
# bed path, which also seeds with 13.
BACKGROUND_SEED = 13
MIN_EVALUE = 0.1
MOTIF_FREQ_CUTOFF = (
    0.8  # If a nucleotide isn't 80% of sites in a motif, converted to an N
)


def get_seq(position):
    """Returns a sequence of a fixed size at a particular position in the reference genome.
    Args:
        position: A string in the form of contig:position specifically a genomic position.
    Returns:
        A string of the sequence in the REF global variable of global variable WINDOW_SIZE length.
    """

    contig = position.split(":")[0]
    pos = int(position.split(":")[1])
    return str(REF[contig][pos - WINDOW_SIZE : pos + WINDOW_SIZE])


def get_ref_pos(position):
    """Returns the base pair at a particular position in the reference genome.
    Args:
        position: A string in the form of contig:position specifically a genomic position
    Returns:
        A string of the base pair in the REF global variable.
    """

    contig = position.split(":")[0]
    pos = int(position.split(":")[1])
    return str(REF[contig][pos])


def run_modkit(
    prefix, bamfile, fasta_file, methylation_confidence_threshold=0.6, threads=14
):
    """Runs Modkit on a given BAM to identify methylated sites; created methylation BED file.
    Args:
        prefix: string prefix of the output
        bamfile: the input mapped BAM to run Modkit on
        fasta_file: The FASTA file for the reference genome associated with the mapped BAM.
        methylation_confidence_threshold: The confidence threshold to pass to Modkit. Essentially a probability of methylation.
        threads: The number of threads to run Modkit with.
    Returns:
        low_modkit_file: path to the output file of modkit
    """

    low_modkit_file = prefix + "_low.bed"

    ## Run first iteration
    cmd = "modkit pileup -t {threads} {bam} {output} -r {fasta} --suppress-progress --filter-threshold {threshold}"
    cmd = cmd.format(
        threads=threads,
        bam=bamfile,
        output=low_modkit_file,
        fasta=fasta_file,
        threshold=methylation_confidence_threshold,
    )
    logging.info("Running modkit: %s", cmd)
    subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode()

    return low_modkit_file


def _read_bedmethyl(path):
    """Read a modkit bedMethyl into the canonical 18-column frame, tolerant of
    the trailing-delimiter change across modkit versions.

    The 9 standard BED columns are always tab-delimited. The 9 trailing
    methylation-stat columns are tab-delimited in some modkit versions
    (18 tab fields total) but space-delimited in newer ones (10 tab fields,
    with the stats packed into the last field). Reading with sep="\\t" and,
    when that yields fewer than 18 columns, splitting the packed final field on
    whitespace recovers the same 18 columns either way. The RGB colour field
    (e.g. "255,0,0") is comma-joined, so it is never split. The re-split stat
    columns are coerced back to numeric so downstream arithmetic stays numeric
    (otherwise "0" + "30" would string-concatenate).
    """
    d = pd.read_csv(path, sep="\t", header=None)
    if d.shape[1] >= 18:
        return d                                   # already all-tab (fast path)
    head = d.iloc[:, :-1].reset_index(drop=True)
    tail = d.iloc[:, -1].astype(str).str.split(expand=True)
    tail = tail.apply(pd.to_numeric, errors="coerce").reset_index(drop=True)
    out = pd.concat([head, tail], axis=1)
    out.columns = range(out.shape[1])
    return out


def read_modkit(low_modkit_output, min_coverage=10):
    """Read the output of Modkit
    Args:
        low_modkit_output: path to a the output file of modkit (returned by run_modkit())
        min_coverage: the minimum coverage required to consider a site for methylation
    Returns:
        final_table: A pandas dataframe containing information about modified sites in the genome.
    """

    logging.info("Reading Modkit table...")
    d = _read_bedmethyl(low_modkit_output)

    ## Rename columns
    d.columns = [
        "Contig",
        "Position",
        "drop-1",
        "Modification",
        "drop0",
        "Strand",
        "drop1",
        "drop2",
        "drop3",
        "drop4",
        "drop5",
        "Modified_bases",
        "Unmodified_bases",
        "Other_mod_base",
        "drop6",
        "Modification_below_threshold",
        "Other_bases",
        "drop7",
    ]
    d = d[d.columns[~d.columns.str.contains("drop")]]

    d["Total_coverage"] = (
        d["Modified_bases"]
        + d["Unmodified_bases"]
        + d["Other_mod_base"]
        + d["Modification_below_threshold"]
        + d["Other_bases"]
    )
    d["Percent_modified"] = d["Modified_bases"] / d["Total_coverage"]
    d["SNP_Position"] = d["Contig"].astype(str) + ":" + d["Position"].astype(str)

    d["SNP_Position_Strand"] = d["SNP_Position"] + "-" + d["Strand"]

    ## Merge tables
    final_table = d

    logging.info("Retrieving sequences")
    final_table["RefBase"] = final_table.SNP_Position.map(get_ref_pos)
    final_table["Sequence"] = final_table.SNP_Position.map(get_seq)

    ## Summary: mean coverage
    logging.info(
        "Total number of sites with any possible calls (including low quality): %s bases",
        final_table.shape[0],
    )

    logging.info(
        "Median coverage (stranded): %sx", final_table["Total_coverage"].median()
    )

    ## Sites > min coverage
    min_coverage_table = final_table[final_table.Total_coverage >= min_coverage]
    logging.info(
        "Total number of sites with greater than min coverage (%sx): %s bases",
        min_coverage,
        min_coverage_table.shape[0],
    )

    for meth_type in d["Modification"].unique():
        m_table = min_coverage_table[min_coverage_table.Modification == meth_type]
        logging.info(
            "Potential sites (>33%%) for methylation type %s: %s bases",
            meth_type,
            m_table[m_table.Percent_modified >= 0.33].shape[0],
        )
        logging.info(
            "High quality sites (>66%%) for methylation type %s: %s bases",
            meth_type,
            m_table[m_table.Percent_modified >= 0.66].shape[0],
        )
        logging.info(
            "Very high quality sites (>90%%) for methylation type %s: %s bases",
            meth_type,
            m_table[m_table.Percent_modified >= 0.90].shape[0],
        )

    return final_table


def write_to_fasta(
    modkit_table, prefix, mod_type, percent_cutoff, min_coverage, subsample=1,
    motif_caller="python",
):
    """Write kmers of length WINDOW_SIZE around methylated positions in reference genome REF.
    Args:
        modkit_table: pandas dataframe containing information about each methylated site.
        prefix: prefix of output files
        mod_type: either "a" or "m", referring to m5C or 6mA
        percent_cutoff: The cutoff of % of reads that has to be methylated to be written out.
        min_coverage: Minimum coverage of a site required to write that site to the output.
        subsample: optional parameter, the portion of sites to randomly subsample down.
        motif_caller: the active motif caller ("python" or "streme"), for the log line only.
    Returns:
        fn_name: file name of the FASTA file created.
    """

    logging.info("Modkit table size: %s", modkit_table.shape[0])
    modkit_table = modkit_table[modkit_table.Percent_modified >= percent_cutoff]
    modkit_table = modkit_table[modkit_table.Total_coverage >= min_coverage]
    modkit_table = modkit_table[modkit_table.Modification == mod_type]

    if subsample != 1:
        modkit_table = modkit_table.sample(int(modkit_table.shape[0] * subsample))

    logging.info(
        "%s sites for motif calling (caller: %s) with cutoff %s and min coverage %s",
        modkit_table.shape[0],
        motif_caller,
        percent_cutoff,
        min_coverage,
    )

    if modkit_table.shape[0] >= 10:
        fn_name = prefix + "_" + mod_type + "_pos.fasta"
        f = open(fn_name, "w+", encoding="utf-8")
        i = 0
        for _, row in modkit_table.iterrows():
            i += 1
            f.write(">" + str(i) + "\n")
            f.write(str(row["Sequence"]) + "\n")

        f.close()

        background_file = fn_name.replace("_pos.fasta", "_control.fasta")
        logging.info("Creating background control distribution: %s", background_file)
        reference_genome = ""
        for _, contig in REF.items():
            reference_genome += str(contig)
        sites = []
        rng = random.Random(BACKGROUND_SEED)
        for i in range(0, 100000):
            random_pos = rng.randint(200, len(reference_genome) - 200)
            sites.append(
                reference_genome[random_pos - WINDOW_SIZE : random_pos + WINDOW_SIZE]
            )
        f = open(background_file, "w+", encoding="utf-8")
        i = 0
        for s in sites:
            f.write(">" + str(i) + "\n")
            f.write(s + "\n")
        f.close()

        return fn_name
    else:
        logging.info(
            "Note: Fewer than 10 methylated sites identified, no motif finding was run."
        )
        return None


def run_streme(kmer_file, streme_path="streme"):
    """Run STREME on a set of files
    Args:
        kmer_file: The FASTA file of methylated kmers, returned by write_to_fasta
        streme_path: Path to the STREME program.
    Returns:
        output_dir: path to output directory of STREME
    """

    # STREME is optional since #49 (python is the default caller). Fail with a
    # clear message instead of a bare exit-127 CalledProcessError when the user
    # asks for --motif_caller streme without the MEME suite installed.
    if shutil.which(streme_path) is None:
        raise SystemExit(
            "Error: STREME executable '%s' not found on PATH. Install the MEME "
            "suite (which provides STREME), or use --motif_caller python (the "
            "default)." % streme_path
        )

    ## Run STREME
    output_dir = kmer_file.split(".fasta")[0] + "_streme"
    cmd = "{streme} . --n {background_input} -p {input_file} -o {output}"

    subprocess.check_output(
        "rm -rf " + output_dir, shell=True, stderr=subprocess.STDOUT
    )

    cmd = cmd.format(
        streme=streme_path,
        output=output_dir,
        input_file=kmer_file,
        background_input=kmer_file.replace("_pos.fasta", "_control.fasta"),
    )
    logging.info("Running STREME: %s", cmd)

    subprocess.check_output(
        cmd, shell=True, stderr=subprocess.DEVNULL, universal_newlines=True
    )
    return output_dir


def assign_motifs(modkit_table, streme_output):
    """Read the output of STREME and assign methylated sites to any identified significant motifs.
    Args:
        modkit_table: The output of read_modkit()
        streme_output: The directory output of STREME
    Returns:
        modkit_table2: A modified version of modkit_table with the motif column added for motif assignments.
    """

    tree = ET.parse(streme_output + "/streme.xml")
    root = tree.getroot()

    motifs = []
    motif_new_to_original = {}  # store original motif names
    motif_sites = {}

    for motif in root[1]:
        # STREME >=5.4.0 reports "test_evalue"; older versions (e.g. 5.3.0)
        # only report "test_pvalue". Fall back so both formats work.
        evalue = motif.get("test_evalue")
        if evalue is None:
            evalue = motif.get("test_pvalue")
        evalue = float(evalue)
        if evalue < MIN_EVALUE:
            ### Evalue cutoff
            motif_old = motif.get("id").split("-")[1]
            motif_new = ""

            ### Convert to N's
            i = 0
            for pos in motif:
                freqs = [float(x) for x in pos.attrib.values()]
                if motif_old[i] in ["A", "C", "G", "T"]:  # If it's a single base
                    if max(freqs) <= MOTIF_FREQ_CUTOFF:  # If it's not >80%, make it N
                        motif_new += "N"
                    else:
                        motif_new += motif_old[i]
                elif motif_old[i] in ["W", "S", "M", "K", "R", "Y"]:  # If it's a double
                    top_two = sum(sorted(freqs)[-2:])
                    if top_two <= MOTIF_FREQ_CUTOFF:  # If it's not >80%, make it N
                        motif_new += "N"
                    else:
                        motif_new += motif_old[i]
                elif motif_old[i] in ["B", "D", "H", "V"]:  # All triplies become N
                    motif_new += "N"
                else:
                    motif_new += motif_old[i]
                i += 1

            ### trim N's
            motif_new = motif_new.strip("N")
            # STREME >=5.4.0 reports "total_sites"; older versions (e.g. 5.3.0)
            # don't, so derive it from the train/test positive site counts.
            total_sites = motif.get("total_sites")
            if total_sites is None:
                total_sites = str(
                    int(motif.get("train_pos_count", 0))
                    + int(motif.get("test_pos_count", 0))
                )
            logging.info(
                "Found and trimmed motif: "
                + motif.get("id")
                + "\t"
                + motif_new
                + "\t"
                + total_sites
            )
            motifs.append(motif_new)
            motif_new_to_original[motif_new] = motif.get("id")

            ## Find motif occurrences in reference
            motif_len = len(motif_new)

            # A motif occurrence physically covers both strands, so mark BOTH
            # the + and - keys for every hit. Keying only the search strand
            # dropped the partner-strand methylation of non-palindromic motifs
            # into "No Motif Assigned" (issue #51). make_motif_table recomputes
            # coverage independently from a strand-agnostic methylated set
            # (methylated_sites is keyed by SNP_Position, no strand), so this
            # does not change Genome_sites / Methylated_sites / coverage — only
            # the per-site motif column and the No-Motif bin.
            for r, contig in REF.items():
                for site in nt_search(str(contig), motif_new)[1:]:
                    for i in range(site, site + motif_len):
                        motif_sites[r + ":" + str(i) + "+"] = motif_new
                        motif_sites[r + ":" + str(i) + "-"] = motif_new

                for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
                    1:
                ]:
                    for i in range(site, site + motif_len):
                        motif_sites[r + ":" + str(i) + "+"] = motif_new
                        motif_sites[r + ":" + str(i) + "-"] = motif_new

    modkit_table2 = modkit_table.copy()
    modkit_table2["site_strand"] = modkit_table2.SNP_Position + modkit_table2.Strand

    modkit_table2["motif"] = modkit_table2["site_strand"].map(motif_sites)
    modkit_table2["motif_raw"] = modkit_table2["motif"].map(motif_new_to_original)

    return modkit_table2

def make_motif_table(modkit_table):
    """Creates the motif table summarizing the data for each motif.
    Args:
        modkit_table: Pandas dataframe of information for each methylated site.
    Returns:
        motif_table: Pandas dataframe of information for each motif.
    """

    motif_table = []

    for meth in modkit_table.Modification.unique():
        modtable_temp = modkit_table[modkit_table.Modification == meth]

        motif_counts = modtable_temp.motif.value_counts()
        methylated_sites = set(modtable_temp["SNP_Position"])

        ### Get motif methylation statistics
        for motif in motif_counts.index:
            motif_data = modtable_temp[modtable_temp.motif == motif]

            ### motif raw version
            motif_raw = motif_data.motif_raw.iloc[0]

            motif_methylated_positions = []

            ### Get Average_Percent_methylation_per_site
            average = motif_data.Percent_modified.mean()

            logging.info("Counting genome sites for: %s", motif)
            logging.info("Motif has %s methylated instances", motif_data.shape[0])
            ### Get genome counts
            genome_counts = 0
            methylated_counts = 0
            for contig, contig_seq in REF.items():
                for motif_site in nt_search(str(contig_seq), Seq(motif))[1:]:
                    genome_counts += 1
                    methylated = False
                    pos_in_motif = 0

                    for i in range(motif_site, motif_site + len(motif)):
                        seq_pos = contig + ":" + str(i)
                        if seq_pos in methylated_sites:
                            methylated = True
                            motif_methylated_positions.append(pos_in_motif + 1)
                        pos_in_motif += 1

                    if methylated:
                        methylated_counts += 1
                for motif_site in nt_search(
                    str(contig_seq), Seq(motif).reverse_complement()
                )[1:]:
                    genome_counts += 1
                    methylated = False
                    pos_in_motif = len(motif)

                    for i in range(motif_site, motif_site + len(motif)):
                        pos_in_motif -= 1
                        seq_pos = contig + ":" + str(i)
                        if seq_pos in methylated_sites:
                            methylated = True
                            motif_methylated_positions.append(pos_in_motif + 1)

                    if methylated:
                        methylated_counts += 1
            logging.info("Methylated positions in motif %s:", motif)
            m_pos_counts = pd.Series(motif_methylated_positions).value_counts()
            m_pos_counts = m_pos_counts / m_pos_counts.sum() * 100
            logging.info(m_pos_counts)

            ## Get top two positions and their frequencies
            if len(m_pos_counts.index) >= 1:
                methylated_position_1 = m_pos_counts.index[0]
                methylated_position_1_percent = m_pos_counts.values[0]
            else:
                methylated_position_1 = "NA"
                methylated_position_1_percent = 0

            if len(m_pos_counts.index) >= 2:
                methylated_position_2 = m_pos_counts.index[1]
                methylated_position_2_percent = m_pos_counts.values[1]
            else:
                methylated_position_2 = "NA"
                methylated_position_2_percent = 0

            motif_table.append(
                {
                    "Motif": motif,
                    "Motif_raw": motif_raw,
                    "Methylation_type": meth,
                    "Genome_sites": genome_counts,
                    "Methylated_sites": methylated_counts,
                    "Methylation_coverage": round(methylated_counts / genome_counts, 3),
                    "Average_Percent_Methylation_per_site": average,
                    "Methylated_position_1": methylated_position_1,
                    "Methylated_position_1_percent": round(
                        methylated_position_1_percent, 3
                    ),
                    "Methylated_position_2": methylated_position_2,
                    "Methylated_position_2_percent": round(
                        methylated_position_2_percent, 3
                    ),
                }
            )

    ## Get number of methylated sites without a motif
    for meth in modkit_table.Modification.unique():
        no_motif = modkit_table[modkit_table.Modification == meth]
        no_motif = no_motif[no_motif.motif.isna()]
        average = no_motif.Percent_modified.mean()

        motif_table.append(
            {
                "Motif": "No Motif Assigned",
                "Motif_raw": "NA",
                "Methylation_type": meth,
                "Genome_sites": "NA",
                "Methylated_sites": no_motif.shape[0],
                "Methylation_coverage": "NA",
                "Average_Percent_Methylation_per_site": average,
                "Methylated_position_1": "NA",
                "Methylated_position_1_percent": "NA",
                "Methylated_position_2": "NA",
                "Methylated_position_2_percent": "NA",
            }
        )

    motif_table = pd.DataFrame(motif_table)

    return motif_table


def main(
    bam_file,
    fasta_file,
    methylation_types,
    output_prefix,
    threads,
    streme_path,
    min_coverage=10,
    percent_cutoff=0.66,
    percent_cutoff_streme=0.9,
    methylation_confidence_threshold=0.6,
    motif_caller="python",
):
    if (
        percent_cutoff > 1
        or percent_cutoff_streme > 1
        or methylation_confidence_threshold > 1
    ):
        raise SystemExit(
            "Error: Methylation fractions and scores should be a decimal between 0 and 1."
        )

    if not os.path.isfile(bam_file + ".bai"):
        raise SystemExit(
            "Error: Bam file index not found (.bai file). Please run samtools index on your Bam."
        )
    if not output_prefix:
        output_prefix = bam_file.split("/")[-1].split(".bam")[0]

    for record in SeqIO.parse(fasta_file, "fasta"):
        REF[record.id] = record.seq

    ### Step 1: Run Modkit on the BAM
    low_modkit_output = run_modkit(
        output_prefix,
        bam_file,
        fasta_file,
        threads=threads,
        methylation_confidence_threshold=methylation_confidence_threshold,
    )

    ### Step 2: Process modkit table
    modkit_table = read_modkit(low_modkit_output, min_coverage)

    # Step 3: For each methylation, create FASTA and call motifs
    logging.info("Motif caller: %s", motif_caller)
    reference_fasta = fasta_file  # preserve genome path (loop var below is the pos FASTA)
    if motif_caller == "python":
        from MicrobeMod import microbe_motif
    i = 0
    final_table = None
    final_motif_table = None
    for methylation_type in methylation_types.split(","):
        logging.info("Running for Methylation type: %s", methylation_type)
        methylation = METHYLATION_TYPES[methylation_type]

        modkit_table_tmp = modkit_table[modkit_table.Modification == methylation]

        # Write methylated-site windows (and a genomic background) to FASTA
        pos_fasta = write_to_fasta(
            modkit_table_tmp,
            output_prefix,
            methylation,
            percent_cutoff_streme,
            min_coverage,
            motif_caller=motif_caller,
        )
        # Call motifs with the selected engine. Both write a STREME-format
        # streme.xml that assign_motifs() parses, so downstream is identical.
        if pos_fasta:
            if motif_caller == "python":
                out_dir = pos_fasta.split(".fasta")[0] + "_microbe_motif"
                streme_out = microbe_motif.run_from_fastas(
                    pos_fasta,
                    pos_fasta.replace("_pos.fasta", "_control.fasta"),
                    out_dir,
                    output_type="xml",
                    genome_path=reference_fasta,
                    mod_type=methylation,
                )
            else:
                streme_out = run_streme(pos_fasta, streme_path)
        else:
            streme_out = None

        # If there is STREME output, then you have motifs!
        if streme_out:
            modkit_table_tmp = assign_motifs(
                modkit_table_tmp, streme_out
            )  # Assign motifs to sites
            modkit_table_tmp = modkit_table_tmp[
                modkit_table_tmp.Percent_modified >= percent_cutoff
            ]  # Just get sites above percent cutoff
            modkit_table_tmp = modkit_table_tmp[
                modkit_table_tmp.Total_coverage >= min_coverage
            ]  # Just get sites above minimum coverage
            motif_table = make_motif_table(modkit_table_tmp)

            if i == 0:  # first loop
                final_table = modkit_table_tmp
                final_motif_table = motif_table
                i += 1
            else:
                final_table = pd.concat([final_table, modkit_table_tmp])
                final_motif_table = pd.concat([final_motif_table, motif_table])

    if final_table is not None and final_motif_table is not None:
        final_table.Modification = final_table.Modification.map(METHYLATION_TYPES_REV)
        final_motif_table.Methylation_type = final_motif_table.Methylation_type.map(METHYLATION_TYPES_REV)
        final_table.reset_index().round(2).to_csv(
            output_prefix + "_methylated_sites.tsv", sep="\t"
        )
        final_motif_table.reset_index().round(2).to_csv(
            output_prefix + "_motifs.tsv", sep="\t"
        )
        logging.info("Complete!")
        logging.info(
            "Saving methylated site table to: %s",
            output_prefix + "_methylated_sites.tsv",
        )
        logging.info("Saving motif output to: %s", output_prefix + "_motifs.tsv")
    else:
        logging.info("Note: No methylated sites identified, so there is no output.")
