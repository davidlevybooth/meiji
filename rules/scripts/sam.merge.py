#!/usr/bin/env python

######################################################################
#
# Merge SAMs 0.1
# Helper script to take SAM files from result from alignment of reads
# to multiple indices. Built for use in Meiji 0.1
#
# David Levy-Booth, 2020
#
######################################################################

import os

# Arguments ##########################################################

dirpath = snakemake.input[0]
outpath = snakemake.output[0]
mapping = True
print("\nMerge SAM 0.1\n")
if mapping:
	print("Ony writing mapped reads to:\n" + outpath)
else:
	print("Writing mapped and unmapped reads to:\n" + outpath)

# Definitions ########################################################

def main():
	sam_merge(dirpath, outpath)

def sam_merge(dirpath, outpath):
	sam_list = sam_collect(dirpath)
	headers = header_collect(sam_list)
	genomes = genomes_collect(headers)

	with open(outpath, "w") as fo:
		for header in headers:
			fo.write(header)
		for sam in sam_list:
			print("Merging " + sam)
			with open(os.path.join(dirpath,sam),"r") as fi:
				for ln in fi:
					lnrow = ln.split() #Add error checking by length of list
					if not ln.startswith("@"):
						if mapping:
							if "*" not in lnrow[4]:
								fo.write(ln)
						else:
							fo.write(ln)

def genomes_collect(headers):
	genomes = [header.split('\t')[1] for header in headers]
	genomes = [genome.strip('SN:') for genome in genomes]
	return(genomes)

def sam_collect(dirpath):
	sam_list = os.listdir(dirpath)
	sam_list = [f for f in sam_list if f[-4:] == ".sam"]
	sam_list.sort()
	if len(sam_list) > 0:
		print("Merging " + str(len(sam_list)) + " SAM files")
	if len(sam_list) == 0:
		print("No SAM files found in directory.")
	return(sam_list)


def header_collect(sam_list):
	header = []
	PG = True
	for sam in sam_list:
		print("Collecting " + sam + " header")
		with open(os.path.join(dirpath,sam),"r") as fi:
			for ln in fi:
				if not ln.startswith("@"):
					break
				if ln.startswith("@PG") and PG:
					header.append(ln)
					PG = False
				elif ln.startswith("@"):
					if ln not in header:
						header.append(ln)
	return(header)


# Run Script ########################################################
if __name__ == '__main__':
	main()
