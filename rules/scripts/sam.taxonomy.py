#!/usr/bin/env python

######################################################################
#
# Map Taxonomy 0.1
# Takes in a two files:
#   1. input file with three columns: Read, Genome, Identity
#   2. taxonomy table: Genome, taxonomy string
# Adds taxonomy string at fourth column to input file and
# writes to output file
#
# David Levy-Booth, 2020
#
######################################################################

import pandas as pd

# Arguments ##########################################################

taxonomy_path = snakemake.input[1]
input_path =  snakemake.input[0]
output_path = snakemake.output[0]

print("\nMap Taxonomy 0.1\n")

# Definitions ########################################################
def main():
	output_taxonomy = add_taxonomy(taxonomy_path, input_path)
	print("Complete! Writing output file: ")
	print(output_path)
	output_taxonomy.to_csv(output_path, sep = "\t", index=False)

def add_taxonomy(taxonomy_path, input_path):
	print("Parsing Taxonomy: " + taxonomy_path)
	taxonomy = pd.read_csv(taxonomy_path, sep="\t", header = None, names=['Genome', 'Species'])
	print("Reading Parsed SAM file: ")
	print(input_path)
	output_file = pd.read_csv(input_path, sep="\t", header = None, names=['Read', 'Genome', 'Identity'])
	print("Adding taxonomy to SAM file")
	try:
		output_taxonomy = output_file.merge(taxonomy, on='Genome', how='left')
	except MergeError:
		print("Could not merge Taxonomy and SAM files")
		print("Ensure SAM is formatted properly: ")
		print("input file with three columns: Read, Genome, Identity")
	return(output_taxonomy)


# Run Script ########################################################
if __name__ == '__main__':
	main()
