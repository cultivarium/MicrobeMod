import os
import sys
import glob
from collections import defaultdict
import subprocess
import pandas as pd
from Bio import SeqIO
import argparse
import subprocess
from Bio.SeqUtils import nt_search
from Bio.Seq import Seq
from pathlib import Path

### Run HMMs
def run_hmmer(output_prefix, input_fasta, threads):
	cmd = 'hmmsearch --cpu {} --domtblout {}.hits RM_HMMs_From_DefenseFinder.hmm {}'

	cmd = cmd.format(threads, output_prefix, input_fasta)

	if os.path.isfile(input_fasta):
		print("Running HMMER:")
		print(cmd)
		subprocess.run(cmd, stdout=subprocess.DEVNULL,
					stderr=subprocess.STDOUT, shell=True)
		return "{}.hits".format(output_prefix)
	else:
		print("ERROR: FASTA file doesn't exist")

def extract_genes(hits, input_fasta, output_prefix):
	rm_gene_file = output_prefix + "_RM_genes.faa"
	f = open(rm_gene_file, 'w+')

	for record in SeqIO.parse(input_fasta, 'fasta'):
		if record.id in hits:
			f.write(">" + record.description + "\n")
			f.write(str(record.seq) + "\n")
	f.close()

	return rm_gene_file

def resolve_hits(hmmer_output, output_prefix):
	cmd = "cath-resolve-hits --input-format hmmer_domtblout {} --hits-text-to-file {}.resolved.hits"

	cmd = cmd.format(hmmer_output, output_prefix)
	print("Running cath:")
	print(cmd)
	subprocess.run(cmd, stdout=subprocess.DEVNULL,
				stderr=subprocess.STDOUT, shell=True)
	return "{}.resolved.hits".format(output_prefix)


def parse_hmmer(resolved_hits):
	f = open(resolved_hits)
	hits = defaultdict(list)
	locations = {}
	for line in f.readlines():
		if not line.startswith("#"):
			gene = line.split()[0]
			hmm = line.split()[1]
			hits[gene].append(hmm)
			contig = "_".join(gene.split("_")[:-1])
			gene_num = int(gene.split("_")[-1])
			locations[gene] = (contig, gene_num)
	f.close()
	return hits, locations

def find_operons(hits, gene_locations, system_types, gene_window = 10):
	links = defaultdict(list)

	for hit in hits:
		### Get types for this gene
		# hit_types = []
		# for hmm in hits[hit]:
		# 	hit_types.append()

		### Find any close genes
		for hit2 in hits:
			if hit != hit2:
				if gene_locations[hit][0] == gene_locations[hit2][0]: # same contig
					if abs(gene_locations[hit][1] - gene_locations[hit2][1]) <= gene_window: # Within 10 genes
						## Found match: get the types of this gene
						links[hit].append(hit2)

	print(links)		

	## Now define operons
	ME_operons = []

	## For each RE
	for hit in hits:
		gene_types = [system_types[hmm][0] for hmm in hits[hit]]
		
		## Is it both RE and ME? if so, done
		if 'RE' in gene_types and 'MT' in gene_types:
			ME_operons.append({"RE": hit, "MT": hit})
			continue

		## Is it a singleton? if so, done
		if len(links[hit]) == 0:
			if 'IIG' in gene_types:
				ME_operons.append({"RE":hit,"MT":hit})
			if 'RE' in gene_types:
				ME_operons.append({"RE":hit, "MT":""})
			if 'MT' in gene_types:
				ME_operons.append({"RE":"","MT":hit})
			continue

		# Does it have any links? 

		if len(links[hit]) > 0:

			# If it is an ME, does it have any RE links?
			if 'MT' in gene_types: 
				REs = ''
				for hit2 in links[hit]:
					gene_types2 = [system_types[hmm][0] for hmm in hits[hit2]]
					if 'RE' in gene_types2:
						if REs == '':
							REs = hit2
						else:
							REs = REs + "," + hit2
				ME_operons.append({"RE":REs,"MT":hit})

			# If it is an RE, and just has RE links, then it is a singleton
			elif 'RE' in gene_types:
				REs = hit
				for hit2 in links[hit]:
					gene_types2 = [system_types[hmm][0] for hmm in hits[hit2]]
					if 'RE' in gene_types2:
							REs = REs + "," + hit2
				ME_operons.append({"RE":REs,"MT":""})

	## Create final operon table
	operon_table = pd.DataFrame(ME_operons)
	return operon_table

def assign_systems_and_filter(operon_table, system_types):
	##

	pass

def main(fasta, output_prefix, threads):

	p = Path(__file__).with_name('restriction_metadata.csv')
	metadata_file = p.absolute()
	metadata = pd.read_csv(metadata_file)
	system_types = {}
	for index, row in metadata.iterrows():
		system_types[row['Name']] = (row['Enzyme_type'], row['System'])

	
	hmmer_output = run_hmmer(output_prefix, fasta, threads)
	resolved_hits = resolve_hits(hmmer_output, output_prefix)
	gene_hits, gene_locations = parse_hmmer(resolved_hits)

	rm_gene_file = extract_genes(gene_hits, fasta, output_prefix)
	operon_table = find_operons(gene_hits, gene_locations, system_types)
	operon_table.to_csv(output_prefix + ".all.tsv", sep="\t", index=False)
	operon_table[(operon_table.RE != '') & (operon_table.MT != '') ].to_csv(output_prefix + ".operon.tsv", sep="\t", index=False)

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Identify Restriction Modification operons in a genome.')

    parser.add_argument("-f", "--fasta", action="store", default=None, required=True, \
        help='Prodigal FAA amino acid file (prodigal -a) for a genome.')
    parser.add_argument("-o", "--output_prefix", action="store", default=None, required=False, \
        help='Output prefix.')

    parser.add_argument("-t", "--threads", action="store", default=12, \
        help='Number of threads to use.')
    args = parser.parse_args()

    if not args.output_prefix:
    	args.output_prefix = args.fasta.split("/")[-1].split(".")[0]
    
    main(args.fasta, args.output_prefix, args.threads)