import pandas as pd
from Bio import SeqIO
import tqdm
import argparse
import subprocess
from pymemesuite.common import MotifFile
from Bio.SeqUtils import nt_search
from Bio.Seq import Seq

METHYLATION_TYPES = {"6mA": "a", "5mC": "m"}
REF = {}
WINDOW_SIZE = 11

REF = {}
WINDOW_SIZE = 11

def get_seq(position):
    contig = position.split(":")[0]
    pos = int(position.split(":")[1])
    return(str(REF[contig][pos-WINDOW_SIZE:pos+WINDOW_SIZE]))

def get_ref_pos(position):
    contig = position.split(":")[0]
    pos = int(position.split(":")[1])
    return(str(REF[contig][pos]))

def run_modkit(prefix, bamfile, fasta_file, low_threshold = 0.6, high_threshold = 0.90, threads=12):
    
    low_modkit_file  = prefix + "_low.bed"
    high_modkit_file = prefix + "_high.bed"
    
    ## Run first iteration
    print("Running modkit round 1")
    cmd = "modkit pileup -t {threads} {bam} {output} -r {fasta} --only-tabs"
    cmd = cmd.format(threads=threads,
                     bam = bamfile,
                     output = low_modkit_file,
                     fasta = fasta_file)
    
    output = subprocess.check_output(
        cmd, shell=True, stderr=subprocess.STDOUT).decode()
    
    ## Run second iteration
    print("Running modkit round 2")
    cmd = "modkit pileup -t {threads} {bam} {output} -r {fasta} --only-tabs --filter-threshold 0.95"
    cmd = cmd.format(threads=threads,
                     bam = bamfile,
                     output = high_modkit_file,
                     fasta = fasta_file)
    output = subprocess.check_output(
        cmd, shell=True, stderr=subprocess.STDOUT).decode()
    
    
    return low_modkit_file, high_modkit_file

def filter_modkit(low_modkit_output, high_modkit_output):
    
    print("Reading low cutoff table")
    d = pd.read_csv(low_modkit_output, sep="\t", header=None)
    
    ## Rename columns
    d.columns = ['Contig','Position','drop-1','Modification','drop0','Strand',
     'drop1','drop2','drop3','drop4','drop5','Modified_bases','Unmodified_bases','Other_mod_base',
     'drop6','Modification_below_threshold','Other_bases','drop7']
    d = d[d.columns[~d.columns.str.contains("drop")]]
    
    d['Total_coverage'] = d['Modified_bases'] + d['Unmodified_bases'] + d['Other_mod_base'] + d['Modification_below_threshold'] + d['Other_bases']
    d['Percent_modified'] = d['Modified_bases'] / d['Total_coverage']
    d['SNP_Position'] = d['Contig'] + ":" + d['Position'].astype(str)
    
    ## Get 2nd table
    print("Reading high cutoff table")
    d2 = pd.read_csv(high_modkit_output, sep="\t", header=None)
    
    ## Rename columns
    d2.columns = ['Contig','Position','drop-1','Modification','drop0','Strand',
     'drop1','drop2','drop3','drop4','drop5','Modified_bases','Unmodified_bases','Other_mod_base',
     'drop6','Modification_below_threshold','Other_bases','drop7']
    d2 = d2[d2.columns[~d2.columns.str.contains("drop")]]
    
    d2['Total_coverage'] = d2['Modified_bases'] + d2['Unmodified_bases'] + d2['Other_mod_base'] + d2['Modification_below_threshold'] + d2['Other_bases']
    d2['Percent_modified'] = d2['Modified_bases'] / d2['Total_coverage']
    d2['SNP_Position'] = d2['Contig'] + ":" + d2['Position'].astype(str)
    
    d['SNP_Position_Strand'] = d['SNP_Position'] + '-' + d['Strand']
    d2['SNP_Position_Strand'] = d2['SNP_Position'] + '-' + d2['Strand']

    ## Merge tables 
    print("Merging")
    d2 = d2[['SNP_Position_Strand','Percent_modified']]
    d2.columns = ['SNP_Position_Strand','Percent_modified_90']    
    final_table = pd.merge(d,d2, on='SNP_Position_Strand')
    del final_table['SNP_Position_Strand']
    
    print("Retrieving sequences")
    final_table['RefBase'] = final_table.SNP_Position.map(get_ref_pos)
    final_table['Sequence'] = final_table.SNP_Position.map(get_seq)

    return final_table

def write_to_fasta(modkit_table, prefix, mod_type, percent_cutoff):
    
    modkit_table = modkit_table[modkit_table.Percent_modified_90 >= percent_cutoff]
    
    modkit_table = modkit_table[modkit_table.Modification == mod_type]
    
    fn_name = prefix + "_" + mod_type + "_pos.fasta"
    f  = open(fn_name, 'w+')
    i = 0 
    for index, row in modkit_table.iterrows():
        i += 1
        f.write(">" + str(i) + "\n")
        f.write(str(row['Sequence']) + "\n")
        
    f.close()
    
    return fn_name

def run_meme(kmer_file, prefix, meme_path = 'meme', threads = 12):
    ## Run MEME
    output_dir = kmer_file.split(".fasta")[0] + "_meme"
    print("Running MEME")
    cmd = "{meme} -dna -mod zoops -nmotifs 10 -minw 3 -maxw 15 -maxsize 0 -p {threads} -o {output} {input_file}"
    cmd = cmd.format(meme = meme_path,
                     threads=threads,
                     output = output_dir,
                     input_file = kmer_file)
    
    output = subprocess.check_output(
        cmd, shell=True, stderr=subprocess.STDOUT).decode()
    
    return output_dir

def assign_motifs(modkit_table, meme_output):
    motif_sites = {}
    
    fn_output = meme_output + "/meme.txt"
    with MotifFile(fn_output) as motif_file:
        while True:
            motif = motif_file.read()
            
            if not motif or motif.evalue >= 1:
                break

            print(motif.consensus)
            motif = motif.consensus.strip("NBVHD")
            print(motif)

            ### count all occcurrences in original sequence
            motif_len = len(motif)
            for contig in REF:
                print(contig)
                for site in nt_search(str(REF[contig]), motif)[1:]:
                    for i in range(site,site+motif_len):
                        motif_sites[contig + ":" + str(i) + "+"] = motif

                for site in nt_search(str(REF[contig]), Seq(motif).reverse_complement())[1:]:
                    for i in range(site,site+motif_len):
                        motif_sites[contig + ":" + str(i) + "-"] = motif

    print(len(motif_sites))
    modkit_table['site_strand'] = modkit_table.SNP_Position + modkit_table.Strand
    modkit_table['motif'] = modkit_table['site_strand'].map(motif_sites)
    return modkit_table

def assign_motifs_round2(modkit_table, meme_output):
    motif_sites = {}
    
    fn_output = meme_output + "/meme.txt"
    with MotifFile(fn_output) as motif_file:
        while True:
            motif = motif_file.read()
            
            if not motif or motif.evalue >= 1:
                break

            print(motif.consensus)
            motif = motif.consensus.strip("NBVHD")
            print(motif)

            ### count all occcurrences in original sequence
            motif_len = len(motif)
            for contig in REF:
                print(contig)
                for site in nt_search(str(REF[contig]), motif)[1:]:
                    for i in range(site,site+motif_len):
                        motif_sites[contig + ":" + str(i) + "+"] = motif

                for site in nt_search(str(REF[contig]), Seq(motif).reverse_complement())[1:]:
                    for i in range(site,site+motif_len):
                        motif_sites[contig + ":" + str(i) + "-"] = motif
    
    print(len(motif_sites))
    
    motif = modkit_table.site_strand.map(motif_sites)
    modkit_table['motif'] = modkit_table['motif'].fillna(motif)
    
    return modkit_table
    
def make_motif_table(modkit_table):
    
    motif_table = []
    
    for meth in modkit_table.Modification.unique():
        modtable_temp = modkit_table[modkit_table.Modification == meth]
        motif_counts = modtable_temp.query("Percent_modified > 0.66").motif.value_counts()

        ### Get motif methylation statistics
        for motif in motif_counts.index:
            motif_data = modtable_temp[modtable_temp.motif == motif]
            
            ### Get Average_Percent_methylation_per_site
            average = motif_data.Percent_modified.mean()

            ### Get those >66%
            print(motif_data.shape)
            methylated_sites = set(motif_data.query("Percent_modified_90 > 0.66")['SNP_Position'])
            print(motif)

            ### Get genome counts
            genome_counts = 0
            methylated_counts = 0
            for contig in REF:
                for motif_site in nt_search(str(REF[contig]), motif)[1:]:
                    genome_counts += 1
                    methylated = False
                    for i in range(motif_site, motif_site + len(motif)):
                        seq_pos = contig + ":" + str(i)
                        if seq_pos in methylated_sites:
                            methylated = True
                            break
                    if methylated:
                        methylated_counts += 1
                for motif_site in nt_search(str(REF[contig]), Seq(motif).reverse_complement())[1:]:
                    genome_counts += 1
                    methylated = False
                    for i in range(motif_site, motif_site + len(motif)):
                        seq_pos = contig + ":" + str(i)
                        if seq_pos in methylated_sites:
                            methylated = True
                            break
                    if methylated:
                        methylated_counts += 1
                    
            motif_table.append({"Motif":motif,"Methylation_type":meth,"Genome_sites":genome_counts,"Methylated_sites":methylated_counts,
                                "Methylation_coverage":methylated_counts/genome_counts,
                                "Average_Percent_Methylation_per_site":average})
    
    motif_table = pd.DataFrame(motif_table)
    
    return motif_table

def main(bam_file, fasta_file, methylation_types, output_prefix, threads):
    
    if not output_prefix:
        output_prefix = bam_file.split("/")[-1].split(".bam")[0]

    for record in SeqIO.parse(fasta_file, 'fasta'):
        REF[record.id] = record.seq


    ### Step 1: Run Modkit on the BAM
    print("Running Modkit")
    low_modkit_output, high_modkit_output = run_modkit(output_prefix, 
                                                   bam_file, 
                                                   fasta_file)
    ### Step 2: Process modkit table

    print("Processing Modkit")
    modkit_table = filter_modkit(low_modkit_output, high_modkit_output)

    ### Step 3: For each methylation, create FASTA and run MEME
    for methylation_type in methylation_types.split(","):
        print("Running for Methylation type: {}".format(methylation_type))
        methylation = METHYLATION_TYPES[methylation_type]
        
        modkit_table_tmp = modkit_table[modkit_table.Modification == methylation]

        print("{} sites".format(modkit_table_tmp.shape[0]))

        ## MEME ROUND 1
        fasta_file = write_to_fasta(modkit_table_tmp, output_prefix, methylation, 0.66)
        meme1 = run_meme(fasta_file, output_prefix + "_1", '/home/alex/miniconda3/envs/meme/bin/meme')
        modkit_table_tmp = assign_motifs(modkit_table_tmp, meme1)

        ## MEME ROUND 2
        fasta_file = write_to_fasta(modkit_table_tmp[modkit_table_tmp.motif.isna()], output_prefix + "_2", methylation, 0.66)
        meme2 = run_meme(fasta_file, output_prefix + "_2", '/home/alex/miniconda3/envs/meme/bin/meme')
    
        modkit_table_tmp = assign_motifs_round2(modkit_table_tmp, meme2)

        modkit_table_tmp.to_csv(output_prefix + "_methylated_sites_" + methylation + ".tsv", sep="\t")

        motif_table = make_motif_table(modkit_table_tmp)
        make_motif_table.to_csv(output_prefix + "_motif_" + methylation + ".tsv", sep="\t")

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

    parser.add_argument("-t", "--threads", action="store", default=12, \
        help='Number of threads to use.')
    args = parser.parse_args()
    
    main(args.bam_file, args.reference_fasta, args.methylation_types, args.output_prefix, args.threads)