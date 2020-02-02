#!/usr/bin/env python

######################################################################
#
# SAM Coverage Filter 0.1
# Helper script to take merged SAM files from Meiji and filter them by
# the coverage results calcualated by pileup.sh in the bbmap suite.
# Default percent coverage cutoff = 1.0
#
# David Levy-Booth, 2020
#
######################################################################

import numpy as np
import pandas as pd
import sys


# Arguments ##########################################################

cut_off = snakemake.params[0]
cov_path = snakemake.input[0]
sam_path = snakemake.input[1]
out_path = snakemake.output[0]

print("\nSAM Coverage 0.1\n")

# Definitions #######################################################
def main():
	filter_coverage(cov_path, sam_path, out_path, cut_off)

def filter_coverage(cov_path, sam_path, out_path, cut_off):
	#open and filter COV file by cut off
	cov = read_cov(cov_path, cut_off)
	n = 0
	kill_list = []
	kill_reads = []
	#open and filter processed SAM file by COV
	with open(out_path, 'w') as newfile:
		with open(sam_path,"r") as fi:
			for ln in fi:
				if not ln.startswith("@"):
					lnrow = ln.split()
					read = lnrow[0]
					genome = lnrow[2]
					n+=1
					if n % 1000000 == 0:
						print("Processed " + str(n) + " alignments. " + str(len(kill_list)) + " genomes < " + str(cut_off) + "% coverage.")
					if genome in cov:
						ID = identity(lnrow)
						newfile.write(read+"\t"+genome+"\t"+str(ID)+"\n")
					else:
						if genome not in kill_list:
							kill_list.append(genome)

def read_cov(cov_path, cut_off):
	cov = pd.read_csv(cov_path, sep="\t", header = 0)
	cov = cov[cov.Covered_percent > cut_off]
	cov_list = cov['#ID'].tolist()
	print("Retaining " + str(len(cov_list)) + " genomes > " + str(cut_off) + "% coverage.")
	return(cov_list)

def identity(lnrow):
	mismatch = lnrow[5].count('1X')
	nlen = len(lnrow[9])
	ID = 100.0-(float(mismatch)/float(nlen)*100)
	ID = round(ID,1)
	return(ID)


#Run Script ###########################################################
if __name__ == '__main__':
	main()
