import sys
import glob
import pysam
import random
import numpy as np
import operator
from collections import defaultdict
from Bio import SeqIO

def quick_coverage(fasta_file, bam_file):

	iterations = 100
	reference_fasta = sys.argv[1]
	bam_file = sys.argv[2]

	seqs = {}
	for record in SeqIO.parse(reference_fasta, 'fasta'):
		seqs[record.id] = record.seq


	samfile = pysam.AlignmentFile(bam_file, "rb")

	print("Contig\tContig_length\tMedian_coverage\tMean_coverage\tMedian_read_ANI")
	for seq in seqs:
		if len(seqs[seq]) >= 500:
			covs = []
			avg_nm = []
			for i in range(0, iterations):
				random_pos = random.randint(1+100,len(seqs[seq])-100)
				cov = 0
				for pileupcolumn in samfile.pileup(seq, start=random_pos, end=random_pos+1, truncate=True): ## Contig name
					for pileupread in pileupcolumn.pileups:
						cov += 1
						read_len = len(pileupread.alignment.get_reference_positions())
						read_pid = 1 - (int(pileupread.alignment.get_tag("NM")) / read_len)
						avg_nm.append(read_pid)

				covs.append(cov)
			print(seq + "\t" + str(len(seqs[seq])) + "\t" + str(np.median(covs)) + "\t" + str(np.mean(covs)) + "\t" + str(round(np.median(avg_nm),3)*100))

if __name__ == '__main__':
	main()