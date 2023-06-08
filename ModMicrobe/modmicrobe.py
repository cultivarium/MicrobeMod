import argparse
import numpy as np
import pandas as pd
import sys
import tqdm
from Bio import SeqIO
from modbampy import ModBam
from collections import defaultdict
from multiprocessing import Pool

METHYLATION_TYPES = {"6mA": "a", "5mC": "m"}

def load_fasta(fasta, bam_file, methylation_type, chunk_size = 100000):
    '''
    Loads a FASTA file in chunks
    '''
    refs = {}
    for record in SeqIO.parse(fasta, 'fasta'):
        refs[record.id] = record.seq
        
    segments = []
    i = 0
    for ref in refs:
        while True:
            if i+chunk_size <= len(refs[ref]):
                segments.append((ref,bam_file,i,i+chunk_size-1,methylation_type))
                i += chunk_size
            else:
                segments.append((ref,bam_file,i,len(refs[ref]),methylation_type))
                break 

    return refs, segments

def create_mod_table(args):
    '''
    Modification codes: 5mC, 6mA.
    '''

    (ref, bam_file, start, stop, mod) = args
    print(str(start) + ",", end='')
    
    with ModBam(bam_file) as bam:
        positions, counts = bam.pileup(
            ref, start, stop,
            low_threshold=0.33, high_threshold=0.66, mod_base=METHYLATION_TYPES[mod])        

    table = []
    i = 0
    rev = {0:"A",1:"C",2:"G",3:"T"}
    fwd = {4:"A",5:"C",6:"G",7:"T"}
    
    if mod == '6mA':
        fwd_i = 4 # A
        rev_i = 0 # A
    elif mod == '5mC':
        fwd_i = 5 # C
        rev_i = 1 # C
        
    fwd_m = 11
    rev_m = 10
    
    results = []
    for pos in positions:
        mod_forward = counts[i][fwd_m]
        mod_rev = counts[i][rev_m]
        if mod_forward > 0: #Forward strand
            total_coverage = np.sum(counts[i][[4,5,6,7,fwd_m]])
            modified_coverage = counts[i][fwd_m]
            base_coverage = counts[i][fwd_i] + modified_coverage
            
            data = {"Contig": ref, "Position": pos, "Strand": '+', "Base": fwd[fwd_i], "Modification": mod, 
                    "Modified_coverage": modified_coverage, "Base_coverage": base_coverage, "Total_coverage": total_coverage}
            results.append(data)
            
        if mod_rev > 0: # Reverse strand
            total_coverage = np.sum(counts[i][[0,1,2,3,rev_m]])
            modified_coverage = counts[i][rev_m]
            base_coverage = counts[i][rev_i] + modified_coverage
            
            data = {"Contig": ref, "Position": pos, "Strand": '-', "Base": rev[rev_i], "Modification": mod, 
                    "Modified_coverage": modified_coverage, "Base_coverage": base_coverage, "Total_coverage": total_coverage}
            results.append(data)

        i += 1
        
    return(pd.DataFrame(results))

def main(bam_file, fasta_file, methylation_types, threads):
    
    ## Load sequence and segments
    pool = Pool(threads)
    
    for methylation_type in methylation_types.split(","):
        print(methylation_type)
        print("Loading FASTA")
        refs, segments = load_fasta(fasta_file, bam_file, methylation_type)
        print("Creating table")
        with Pool(threads) as pool:
            results = list(tqdm.tqdm(pool.imap(create_mod_table,segments), total=len(segments)))
        print(results)

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Pipeline for plasmid sequence verification.')

    parser.add_argument("-b", "--bam_file", action="store", default=None, required=True, \
        help='BAM file of nanopore reads mapped to reference genome with the MM and ML tags preserved.')
    parser.add_argument("-r", "--reference_fasta", action="store", default=None, required=True, \
        help='Reference genome FASTA file.')
    parser.add_argument("-m", "--methylation_types", action="store", default="6mA,5mC", required=False, \
        help='Methylation types to profile')

    parser.add_argument("-t", "--threads", action="store", default=12, \
        help='Number of threads to use.')
    args = parser.parse_args()
    
    main(args.bam_file, args.reference_fasta, args.methylation_types, args.threads)