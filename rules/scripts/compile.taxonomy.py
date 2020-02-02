#!/usr/bin/env python

######################################################################
#
# Compile taxonomy from processed SAM file 0.1
#
# Lots of pandas stuff in here. Basically, we calculate fractional
# counts for reads that map to multiple genomes. Then find the sum
# of reads for each species in each sample. Finally pivot to create
# the species x sample table. Built for use in Meiji
#
# David Levy-Booth, 2020
#
######################################################################

import numpy as np
import pandas as pd

# arguments ##########################################################

taxonomy_path = snakemake.input[0]
taxonomy_count_path = snakemake.output[0]

print("\nCompile Taxonomy Table 0.1\n")


# Functions ###########################################################
def main():
	#Read and format taxonomy
	read_taxonomy = pd.read_csv(taxonomy_path, sep="\t")
	read_taxonomy = read_taxonomy.dropna(subset=['Species'])
	#Create taxonomy count table
	taxonomy_table = pd.DataFrame(read_taxonomy.groupby('Read').Species.value_counts(normalize = True)) #.drop_duplicates(keep=False)
	taxonomy_table.columns = ['Count']
	taxonomy_table.reset_index(inplace = True)
	taxonomy_table["Sample"] = taxonomy_table["Read"].str.split("_", n = 1, expand = True)[0]
	taxonomy_table = taxonomy_table.groupby(['Sample', 'Species']).sum()[['Count']]
	taxonomy_table.reset_index(inplace = True)
	#Create species x sample table
	taxonomy_table = taxonomy_table.pivot_table(index=['Species'], columns='Sample', values='Count', aggfunc=np.sum)
	taxonomy_table = taxonomy_table.fillna(0).astype(int)
	taxonomy_table = taxonomy_table.reset_index()
	#Write Taxonomy Count Table
	taxonomy_table.to_csv(taxonomy_count_path, sep = "\t", index=False)

#Run Script #####################################################################
if __name__ == '__main__':
	main()
