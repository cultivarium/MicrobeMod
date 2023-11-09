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

def run_prodigal(output_prefix, input_fasta):
	''' Runs prodigal on a genome.
	Args: 
		output_prefix: prefix of output file.
		input_fasta: path to FASTA file.

	Returns:
		The name of the prodigal FAA file generated. 
	'''

	cmd = 'prodigal -i {} -a {}.faa'
	cmd = cmd.format(input_fasta, output_prefix)
	
	if os.path.isfile(input_fasta):
		print("Running prodigal: {}".format(cmd))

		subprocess.run(cmd, stdout=subprocess.DEVNULL,
				stderr=subprocess.STDOUT, shell=True)
	
		return output_prefix + ".faa"

	else:
		print("ERROR: FASTA file doesn't exist")

def run_hmmer(output_prefix, prodigal_fasta, threads):
	''' Runs HMMER to identify RM genes on a set of protein predictions.
	Args: 
		output_prefix: prefix of output file.
		prodigal_fasta: path to prodigal FAA file.
		threads: number of threads to pass to HMMER

	Returns:
		The name of the HMMER file created. 
	'''

	cmd = 'hmmsearch --cut_ga --cpu {} --domtblout {}.hits {} {}'

	rm_hmm_file = os.path.dirname(__file__)   + "/db/HMMs/RM_HMMs.hmm"

	cmd = cmd.format(threads, output_prefix, rm_hmm_file, prodigal_fasta)

	if os.path.isfile(prodigal_fasta):
		print("Running HMMER: " + cmd)
		subprocess.run(cmd, stdout=subprocess.DEVNULL,
					stderr=subprocess.STDOUT, shell=True)
		return "{}.hits".format(output_prefix)
	else:
		print("ERROR: Prodigal output file not found")

def extract_genes(hits, prodigal_fasta, output_prefix):
	'''Creates a FASTA file of genes involved in RM.
	Args:
		hits: a dictionary of genes to their corresponding HMM hits.
		prodigal_fasta: path to the FAA file created by prodigal.
		output_prefix: prefix of output file.

	Returns:
		Name of RM protein sequence file created.
	'''

	rm_gene_file = output_prefix + ".rm.genes.faa"
	f = open(rm_gene_file, 'w+')

	for record in SeqIO.parse(prodigal_fasta, 'fasta'):
		if record.id in hits:
			f.write(">" + record.description + "\n")
			f.write(str(record.seq) + "\n")
	f.close()

	return rm_gene_file

def resolve_hits(hmmer_output, output_prefix):
	'''Runths cath-resolve-hits on the output of HMMER. This resolves
	the best non-overlapping set of HMMs for each gene.

	Args:
		hmmer_output: the path of the file output by HMMER
		output_prefix: prefix of output file.

	Returns:
		Name of cath resolve hits file created.
	'''

	cmd = "cath-resolve-hits --input-format hmmer_domtblout {} --hits-text-to-file {}.resolved.hits"

	cmd = cmd.format(hmmer_output, output_prefix)
	print("Running cath: " + cmd)
	subprocess.run(cmd, stdout=subprocess.DEVNULL,
				stderr=subprocess.STDOUT, shell=True)
	return "{}.resolved.hits".format(output_prefix)


def parse_hmmer(resolved_hits):
	'''Reads the output of cath resolve hits and converts it into dictionaries
	linking each gene to information about that gene.

	Args:
		resolved_hits: the path of the file output by cath-resolve-hits

	Returns:
		hits: a defaultdictionary where keys are genes and values are lists of their HMM names.
		locations: a dictionary where keys are genes and values are tuples of (contig, gene number).
		evalues: a defaultdictionary where keys are genes, values are dictionaries with HMMs as keys and values as evalues.
	'''	

	offtarget_file = os.path.dirname(__file__)   + "/db/HMMs/off_target.txt"
	offtarget = []
	f = open(offtarget_file)
	for line in f.readlines():
		offtarget.append(line.strip())
	f.close()

	f = open(resolved_hits)
	hits = defaultdict(list)
	locations = {}
	evalues = defaultdict(dict)
	for line in f.readlines():
		if not line.startswith("#"):
			gene = line.split()[0]
			hmm = line.split()[1]
			contig = "_".join(gene.split("_")[:-1])
			gene_num = int(gene.split("_")[-1])
			evalue = float(line.split()[-1])

			if hmm not in offtarget:
				hits[gene].append(hmm)
				evalues[gene][hmm] = evalue
				locations[gene] = (contig, gene_num)
	f.close()
	return hits, locations, evalues

def create_gene_table(hits, gene_locations, system_types, evalues, blast_hits, gene_window = 10):
	'''Creates a final table of RM genes organized by their types and their operon status.
	Args:
		hits: a defaultdictionary where keys are genes and values are lists of their HMM names.
		gene_locations: a dictionary where keys are genes and values are tuples of (contig, gene number).
		system_types: Dictionary of HMMs to RM system types
		evalues: a defaultdictionary where keys are genes, values are dictionaries with HMMs as keys and values as evalues.
		blast_hits: dictionary of genes-> blast hits, from read_blast()
		gene_window: the window size to use to call operons.
	Returns:
		gene_table: A final pandas table of RM genes organized by type and operon.
		
	'''	

	gene_table = []

	for hit in hits:

		## Get BLAST info
		best_blast=blast_pid=meth_type=motif = ''
		if hit in blast_hits:
			best_blast = blast_hits[hit][0]
			blast_pid = blast_hits[hit][1]
			meth_type = blast_hits[hit][2]
			motif = blast_hits[hit][3]

		gene_types = set() # So that each HMM MT/RE hit only counts once per gene

		for hmm in hits[hit]:
			## Get the system type of this HMM
			if system_types[hmm][0] not in gene_types:
				gene_types.add(system_types[hmm][0])
				gene_table.append({
					"Gene": hit,
					"Contig": "_".join(hit.split("_")[:-1]),
					"Gene Position": int(hit.split("_")[-1]),
					"System Type": system_types[hmm][1],
					"Gene type": system_types[hmm][0],
					"HMM": hmm,
					"Evalue": evalues[hit][hmm],
					"REBASE homolog": best_blast,
					"Homolog identity(%)": blast_pid,
					"Homolog methylation": meth_type,
					"Homolog motif": motif
					})
			
	gene_table = pd.DataFrame(gene_table)
	gene_table = gene_table.sort_values(["Contig","Gene Position"], ascending=True)

	## Assign operons
	gene2operon = {}
	operon2gene = defaultdict(list)
	operon_number = 0
	pos = 0
	contig = ''

	for index, row in gene_table.iterrows():	
		## New contig = new operon
		if row['Contig'] != contig:
			operon_number += 1
			contig = row['Contig']
			pos = row['Gene Position']
			gene2operon[row['Gene']] = operon_number
			operon2gene[operon_number].append(row['System Type'] + ":" +row['Gene type'])

		# Same contig, gene within window = same operon
		# Update position (extend operons)
		elif row['Contig'] == contig and row['Gene Position'] <= pos + gene_window:
			pos = row['Gene Position']
			gene2operon[row['Gene']] = operon_number
			operon2gene[operon_number].append(row['System Type'] + ":" +row['Gene type'])			

		# Same contig, gene outside of window = new operon
		else:
			operon_number += 1
			contig = row['Contig']
			pos = row['Gene Position']
			gene2operon[row['Gene']] = operon_number
			operon2gene[operon_number].append(row['System Type'] + ":" +row['Gene type'])			


	## Only operons with 1 MT, 1 RE are kept - others are singletons
	## Keep Type IVs

	operon_number = 0
	singleton_number = 0
	operon_labels = {}

	for index, row in gene_table.iterrows():
		found_MT = False
		found_RE = False
		found_RE_IV = False 
		operon = gene2operon[row['Gene']]

		if operon not in operon_labels:
			for genetype in operon2gene[operon]:
				if 'MT' in genetype:
					found_MT = True
				if 'RE' in genetype:
					found_RE = True
				if 'RM_Type_IV' in genetype or 'IIG' in genetype:
					found_RE_IV = True

			if found_RE_IV or (found_MT and found_RE):
				operon_number += 1
				operon_labels[operon] = "RM Operon #" + str(operon_number)
			else:
				singleton_number += 1
				operon_labels[operon] = "Singleton #" + str(singleton_number)

	gene_table.insert(0,'Operon', gene_table['Gene'].map(gene2operon).map(operon_labels))
	gene_table = gene_table.sort_values(["Operon","REBASE homolog", "Gene Position"])
	del gene_table['Gene Position']
	del gene_table['Contig']

	return gene_table

def run_rebase_blast(rm_gene_file, output_prefix, threads):
	'''Run a BLAST of RM gene files against REBASE
	Args:
		rm_gene_file: a defaultdictionary where keys are genes and values are lists of their HMM names.
		output_prefix: a dictionary where keys are genes and values are tuples of (contig, gene number).
	Returns:
		file path of BLAST output.
	'''	

	cmd = 'blastp -query {} -db {} -outfmt 6 -evalue 1e-5 -num_threads {} > {}.blast'

	blast_db_file = os.path.dirname(__file__)  + "/db/rebase_blast/all_rebase_proteins.faa"

	cmd = cmd.format(rm_gene_file, blast_db_file, threads, output_prefix)

	print("Running BLASTP against REBASE: " + cmd)
	subprocess.run(cmd, stdout=subprocess.DEVNULL,
				stderr=subprocess.STDOUT, shell=True)

	return "{}.blast".format(output_prefix)

def read_blast(blast_file):
	''' Read the output of BLAST
	Args:
		blast_file: path to a file output by BLAST against REBASE database.
	Returns:
		blast_hits: dictionary of BLAST hits
	'''
	
	f = open(blast_file)
	blast_hits = {}
	prots = set()

	for line in f.readlines():
		hit = line.split()[0]
		pid = float(line.split()[2])

		if hit not in prots: # top hit only
			prots.add(hit)
			if pid > 50:
				best_hit = line.split("\t")[1].split(":")[1].split("-")[0]
				meth_type = line.split("\t")[1].split('-')[-2]
				motif = line.split("\t")[1].split("-")[-1]

				blast_hits[hit] = (best_hit, pid, meth_type, motif)

	return blast_hits


def main(fasta, output_prefix, threads):

	metadata_file = os.path.dirname(__file__) + "/db/restriction_metadata.csv"
	metadata = pd.read_csv(metadata_file)
	system_types = {}
	for index, row in metadata.iterrows():
		system_types[row['Name']] = (row['Enzyme_type'], row['System'])

	prodigal_fasta = run_prodigal(output_prefix, fasta)
	hmmer_output = run_hmmer(output_prefix, prodigal_fasta, threads)

	resolved_hits = resolve_hits(hmmer_output, output_prefix)
	gene_hits, gene_locations, evalues = parse_hmmer(resolved_hits)

	rm_gene_file = extract_genes(gene_hits, prodigal_fasta, output_prefix)

	blast_file = run_rebase_blast(rm_gene_file, output_prefix, threads)
	blast_hits = read_blast(blast_file)

	gene_table = create_gene_table(gene_hits, gene_locations, system_types, evalues, blast_hits)
	gene_table.to_csv(output_prefix + "_RM_genes.tsv", sep="\t", index=False)