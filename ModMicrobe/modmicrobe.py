import pandas as pd
from Bio import SeqIO
import tqdm
import numpy as np
import random
import argparse
from collections import defaultdict
import subprocess
from Bio.SeqUtils import nt_search
from Bio.Seq import Seq
import xml.etree.ElementTree as ET

METHYLATION_TYPES = {"6mA": "a", "5mC": "m"}

REF = {}
WINDOW_SIZE = 12

def get_seq(position):
    contig = position.split(":")[0]
    pos = int(position.split(":")[1])
    return(str(REF[contig][pos-WINDOW_SIZE:pos+WINDOW_SIZE]))

def get_ref_pos(position):
    contig = position.split(":")[0]
    pos = int(position.split(":")[1])
    return(str(REF[contig][pos]))


def run_modkit(prefix, bamfile, fasta_file, methylation_confidence_threshold = 0.6, threads=14):
    
    low_modkit_file  = prefix + "_low.bed"
    high_modkit_file = prefix + "_high.bed"
    
    ## Run first iteration
    print("Running modkit")
    cmd = "modkit pileup -t {threads} {bam} {output} -r {fasta} --only-tabs --filter-threshold {threshold}"
    cmd = cmd.format(threads=threads,
                     bam = bamfile,
                     output = low_modkit_file,
                     fasta = fasta_file,
                     threshold = methylation_confidence_threshold)
    
    output = subprocess.check_output(
        cmd, shell=True, stderr=subprocess.STDOUT).decode()
    
    return low_modkit_file

def read_modkit(low_modkit_output, min_coverage = 20):
    
    print("Reading low cutoff Modkit table...")
    d = pd.read_csv(low_modkit_output, sep="\t", header=None)
    
    ## Rename columns
    d.columns = ['Contig','Position','drop-1','Modification','drop0','Strand',
     'drop1','drop2','drop3','drop4','drop5','Modified_bases','Unmodified_bases','Other_mod_base',
     'drop6','Modification_below_threshold','Other_bases','drop7']
    d = d[d.columns[~d.columns.str.contains("drop")]]
    
    d['Total_coverage'] = d['Modified_bases'] + d['Unmodified_bases'] + d['Other_mod_base'] + d['Modification_below_threshold'] + d['Other_bases']
    d['Percent_modified'] = d['Modified_bases'] / d['Total_coverage']
    d['SNP_Position'] = d['Contig'].astype(str) + ":" + d['Position'].astype(str)
        
    d['SNP_Position_Strand'] = d['SNP_Position'] + '-' + d['Strand']

    ## Merge tables     
    final_table = d

    print("Retrieving sequences")
    final_table['RefBase'] = final_table.SNP_Position.map(get_ref_pos)
    final_table['Sequence'] = final_table.SNP_Position.map(get_seq)

    ## Summary: mean coverage
    print("Total number of sites with any possible calls (including low quality): {} bases".format(final_table.shape[0]))
    print("Median coverage: {}x".format(final_table['Total_coverage'].median()))

    ## Sites > min coverage
    min_coverage_table = final_table[final_table.Total_coverage >= min_coverage]
    print("Total number of sites with greater than min coverage ({}x): {} bases".format(min_coverage, min_coverage_table.shape[0]))
    
    for meth_type in d['Modification'].unique():
        m_table = min_coverage_table[min_coverage_table.Modification == meth_type]
        print("Potential sites (>33%) for methylation type {}: {} bases".format(meth_type, m_table[m_table.Percent_modified >= 0.33].shape[0]))
        print("High quality sites (>66%) for methylation type {}: {} bases".format(meth_type, m_table[m_table.Percent_modified >= 0.66].shape[0]))
        print("Very high quality sites (>90%) for methylation type {}: {} bases".format(meth_type, m_table[m_table.Percent_modified >= 0.90].shape[0]))

    return final_table

def write_to_fasta(modkit_table, prefix, mod_type, percent_cutoff, min_coverage, subsample = 1):
    
    print("Modkit table size: {}".format(modkit_table.shape[0]))
    modkit_table = modkit_table[modkit_table.Percent_modified >= percent_cutoff]
    modkit_table = modkit_table[modkit_table.Total_coverage >= min_coverage]
    modkit_table = modkit_table[modkit_table.Modification == mod_type]


    modkit_table = modkit_table.sample(int(modkit_table.shape[0]*subsample))

    print("{} sites sent to STREME with cutoff {} and min coverage {}".format(modkit_table.shape[0],percent_cutoff,min_coverage))
    
    fn_name = prefix + "_" + mod_type + "_pos.fasta"
    f  = open(fn_name, 'w+')
    i = 0 
    for index, row in modkit_table.iterrows():
        i += 1
        f.write(">" + str(i) + "\n")
        f.write(str(row['Sequence']) + "\n")
        
    f.close()

    background_file = fn_name.replace("_pos.fasta", "_control.fasta")
    print("Creating background control distribution: {}".format(background_file))
    reference_genome = ''
    for r in REF:
        reference_genome += str(REF[r])
    sites = []
    for i in range(0, 100000):
        random_pos = random.randint(200, len(reference_genome)-200)
        sites.append(reference_genome[random_pos-WINDOW_SIZE:random_pos+WINDOW_SIZE])
    f =  open(background_file, "w+")
    i = 0
    for s in sites:
        f.write(">" + str(i) + "\n")
        f.write(s + "\n")
    f.close()

    return fn_name

def run_streme(kmer_file, prefix, streme_path = 'streme'):
    ## Run MEME
    output_dir = kmer_file.split(".fasta")[0] + "_streme"
    print("Running STREME")
    cmd = "{streme} --minw 4 --n {background_input} -p {input_file} -o {output}"

    output = subprocess.check_output("rm -rf " + output_dir, shell=True, stderr=subprocess.STDOUT).decode()

    cmd = cmd.format(streme = streme_path,
                     output = output_dir,
                     input_file = kmer_file, 
                     background_input = kmer_file.replace("_pos.fasta","_control.fasta"))

    print(cmd)
    
    try:
        output = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.DEVNULL, universal_newlines=True)
        return output_dir
    except:
        return False
    

def assign_motifs(modkit_table, streme_output):
    
    tree = ET.parse(streme_output + "/streme.xml")
    root = tree.getroot()

    motifs = []
    motif_new_to_original = {} # store original motif names
    motif_sites = {}
    motif_total = defaultdict(int)
    for motif in root[1]:
        evalue = float(motif.get("test_evalue"))
        if evalue < 0.1:
                
            ### Evalue cutoff
            motif_old = motif.get("id").split("-")[1]
            motif_new = ''

            ### Convert to N's
            i = 0
            for pos in motif:
                freqs = [float(x) for x in pos.attrib.values()]
                if motif_old[i] in ['A','C','G','T']: # If it's a single base
                    if max(freqs) <= 0.8: # If it's not >80%, make it N
                        motif_new += 'N'
                    else:
                        motif_new += motif_old[i]
                elif motif_old[i] in ['W','S','M','K','R','Y']: # If it's a double
                    top_two = sum(sorted(freqs)[-2:])
                    if top_two <= 0.8: #If it's not >80%, make it N
                        motif_new += 'N'
                    else:
                        motif_new += motif_old[i]
                elif motif_old[i] in ['B','D','H','V']: # All triplies become N
                    motif_new += 'N'
                else:
                    motif_new += motif_old[i]
                i += 1

            ### trim N's
            motif_new = motif_new.strip("N")
            print("Found and trimmed motif: " + motif.get("id") + "\t" + motif_new + "\t" + motif.get("npassing"))
            motifs.append(motif_new)
            motif_new_to_original[motif_new] = motif.get("id")

            ## Find motif occurrences in reference
            motif_len = len(motif_new)

            for contig in REF:
                for site in nt_search(str(REF[contig]), motif_new)[1:]:
                    for i in range(site,site+motif_len):
                        motif_sites[contig + ":" + str(i) + "+"] = motif_new

                for site in nt_search(str(REF[contig]), Seq(motif_new).reverse_complement())[1:]:
                    for i in range(site,site+motif_len):
                        motif_sites[contig + ":" + str(i) + "-"] = motif_new


    modkit_table2 = modkit_table.copy()
    modkit_table2['site_strand'] = modkit_table2.SNP_Position + modkit_table2.Strand
    
    modkit_table2['motif'] = modkit_table2['site_strand'].map(motif_sites)
    modkit_table2['motif_raw'] = modkit_table2['motif'].map(motif_new_to_original)


    return modkit_table2
    
def make_motif_table(modkit_table):
    
    motif_table = []

    for meth in modkit_table.Modification.unique():
        modtable_temp = modkit_table[modkit_table.Modification == meth]

        motif_counts = modtable_temp.motif.value_counts()
        methylated_sites = set(modtable_temp['SNP_Position'])

        ### Get motif methylation statistics
        for motif in motif_counts.index:
            motif_data = modtable_temp[modtable_temp.motif == motif]
            
            ### motif raw version
            motif_raw = motif_data.motif_raw.iloc[0]

            motif_methylated_positions = []

            ### Get Average_Percent_methylation_per_site
            average = motif_data.Percent_modified.mean()

            print("Counting genome sites for: " + motif)
            print("Motif has {} methylated instances".format(motif_data.shape[0]))
            ### Get genome counts
            genome_counts = 0
            methylated_counts = 0
            for contig in REF:
                for motif_site in nt_search(str(REF[contig]), motif)[1:]:
                    genome_counts += 1
                    methylated = False
                    pos_in_motif = 0

                    for i in range(motif_site, motif_site + len(motif)):
                        seq_pos = contig + ":" + str(i)
                        if seq_pos in methylated_sites:
                            methylated = True
                            motif_methylated_positions.append(pos_in_motif+1)
                        pos_in_motif += 1

                    if methylated:
                        methylated_counts += 1
                for motif_site in nt_search(str(Seq(REF[contig]).reverse_complement()), motif)[1:]:
                    genome_counts += 1
                    methylated = False
                    pos_in_motif = 0

                    for i in range(motif_site, motif_site + len(motif)):
                        seq_pos = contig + ":" + str(i)
                        if seq_pos in methylated_sites:
                            methylated = True
                            motif_methylated_positions.append(pos_in_motif+1)
                        pos_in_motif += 1

                    if methylated:
                        methylated_counts += 1
            print("Methylated positions in motif {}:".format(motif))
            m_pos_counts = pd.Series(motif_methylated_positions).value_counts()
            m_pos_counts = m_pos_counts / m_pos_counts.sum() * 100
            print(m_pos_counts)

            ## Get top two positions and their frequencies
            motif_table.append({"Motif":motif,"Motif_raw": motif_raw, "Methylation_type":meth,"Genome_sites":genome_counts,"Methylated_sites":methylated_counts,
                                "Methylation_coverage":methylated_counts/genome_counts,"Average_Percent_Methylation_per_site":average,
                                "Methylated_position_1":m_pos_counts.index[0],'Methylated_position_1_percent': m_pos_counts.values[0],
                                "Methylated_position_2":m_pos_counts.index[1],'Methylated_position_2_percent': m_pos_counts.values[1],})
    
    ## Get number of methylated sites without a motif
    for meth in modkit_table.Modification.unique():
        no_motif = modkit_table[modkit_table.Modification == meth]
        no_motif = no_motif[no_motif.motif.isna()]
        average = no_motif.Percent_modified.mean()

        motif_table.append({"Motif":"No Motif Assigned", "Motif_raw": "NA", "Methylation_type":meth,"Genome_sites":'NA',"Methylated_sites":no_motif.shape[0],
                            "Methylation_coverage":'NA',"Average_Percent_Methylation_per_site":average,
                                "Methylated_position_1":'NA','Methylated_position_1_percent': 'NA',
                                "Methylated_position_2":'NA','Methylated_position_2_percent': 'NA',})

    motif_table = pd.DataFrame(motif_table)
    
    return motif_table

def main(bam_file, fasta_file, methylation_types, output_prefix, threads, streme_path, 
    min_coverage = 20, percent_cutoff = 0.66, percent_cutoff_streme = 0.9, methylation_confidence_threshold = 0.6):

    if not output_prefix:
        output_prefix = bam_file.split("/")[-1].split(".bam")[0]

    for record in SeqIO.parse(fasta_file, 'fasta'):
        REF[record.id] = record.seq

    ### Step 1: Run Modkit on the BAM
    low_modkit_output  = run_modkit(output_prefix, 
                                    bam_file, 
                                    fasta_file,
                                    threads = threads, 
                                    methylation_confidence_threshold = methylation_confidence_threshold)

    ### Step 2: Process modkit table
    modkit_table = read_modkit(low_modkit_output, min_coverage)

    # Step 3: For each methylation, create FASTA and run STREME
    i = 0
    for methylation_type in methylation_types.split(","):
        print("Running for Methylation type: {}".format(methylation_type))
        methylation = METHYLATION_TYPES[methylation_type]
        
        modkit_table_tmp = modkit_table[modkit_table.Modification == methylation]

        # Write sites to FASTA files
        fasta_file = write_to_fasta(modkit_table_tmp, output_prefix, methylation, percent_cutoff_streme, min_coverage)
        # Run STREME
        streme_out = run_streme(fasta_file, output_prefix, streme_path)

        # If there is STREME output, then you have motifs!
        if streme_out:
            modkit_table_tmp = assign_motifs(modkit_table_tmp, streme_out) # Assign motifs to sites
            modkit_table_tmp = modkit_table_tmp[modkit_table_tmp.Percent_modified >= percent_cutoff] # Just get sites above percent cutoff
            modkit_table_tmp = modkit_table_tmp[modkit_table_tmp.Total_coverage >= min_coverage] # Just get sites above minimum coverage
            motif_table = make_motif_table(modkit_table_tmp)

            if i == 0: # first loop
                final_table = modkit_table_tmp
                final_motif_table = motif_table
                i += 1
            else:
                final_table = pd.concat([final_table, modkit_table_tmp])
                final_motif_table = pd.concat([final_motif_table, motif_table])

    final_table.reset_index().round(2).to_csv(output_prefix + "_methylated_sites.tsv", sep="\t")            
    final_motif_table.reset_index().round(2).to_csv(output_prefix + "_motifs.tsv", sep="\t")

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Pipeline for plasmid sequence verification.')

    parser.add_argument("-b", "--bam_file", action="store", default=None, required=True, \
        help='BAM file of nanopore reads mapped to reference genome with the MM and ML tags preserved.')

    parser.add_argument("-r", "--reference_fasta", action="store", default=None, required=True, \
        help='Reference genome FASTA file.')
    parser.add_argument("-m", "--methylation_types", action="store", default="6mA,5mC", required=False, \
        help='Methylation types to profile')
    parser.add_argument("-o", "--output_prefix", action="store", default=None, required=False, \
        help='Output prefix. Default is based on the BAM filename.')
    parser.add_argument("-s", "--streme_path", action="store", default='streme', required=False, \
        help='Path to streme executable.')

    parser.add_argument("--min_coverage", action="store", default=20, required=False, type=int, \
        help='Minimum coverage required to call a site as methylated.')
    parser.add_argument("--methylation_confidence_threshold", action="store", default=0.66, required=False, type=float, \
        help='The minimum confidence score to call a base on a read as methylated. Passed to modkit. Default: 0.66')    
    parser.add_argument("--percent_methylation_cutoff", action="store", default=0.66, required=False, type=float, \
        help='The fraction of methylated reads mapping to a site to count that site as methylated. Default: 0.66')
    parser.add_argument("--percent_cutoff_streme", action="store", default=0.9, required=False, type=float, \
        help='The fraction of methylated reads mapping to a site to pass that site to motif calling. Default: 0.9')


    parser.add_argument("-t", "--threads", action="store", default=12, \
        help='Number of threads to use. Only the first step (modkit) is multithreaded.')
    args = parser.parse_args()

    if args.percent_methylation_cutoff > 1 or args.percent_cutoff_streme > 1 or args.methylation_confidence_threshold > 1:
        raise SystemExit('Error: Methylation fractions and scores should be a decimal between 0 and 1.')

    main(args.bam_file, args.reference_fasta, args.methylation_types, args.output_prefix, args.threads, args.streme_path, 
        args.min_coverage, args.percent_methylation_cutoff, args.percent_cutoff_streme, args.methylation_confidence_threshold)