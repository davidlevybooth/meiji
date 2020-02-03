#!/usr/bin/env python

######################################################################
#
# Community Detection 0.1
# Find ANI-based community centroids using python-louvain
# https://python-louvain.readthedocs.io/en/latest/
# https://github.com/taynaud/python-louvain
# Depends on http://networkx.lanl.gov/
#
# Cluster reads to cluster centroids to avoid "genome splitting"
# Clustering performed by default at ANI > 94 (See link below)
# https://microscope.readthedocs.io/en/stable/content/compgenomics/genoclust.html#clustering-genomes
#
# Cluster centroids calculated with total counts for each genome.
# If you want to use percent coverage to calculate cluster centroids, you have to provide
# the coverage file (produced by pileup.sh) with the -v flag
# (They're almost always the same though)
#
# Build for use with Meiji 0.1
#
# David Levy-Booth, 2020
#
######################################################################

import pandas as pd
import numpy as np
import community
import networkx as nx
import matplotlib.pyplot as plt

# Arguments ##########################################################

cut_off = snakemake.params[0]
filter_cutoff = snakemake.params[1]
write_plot = snakemake.params[2]

input_path = snakemake.input[0]
cov_path =  snakemake.input[1]
ani_path = snakemake.input[2]

output_path = snakemake.output[0]


print("\nCommunity Detection 0.1")
print("Collapsing genomes to ANI centriods at " + str(cut_off) + " ANI.\n")


# Definitions #######################################################
def main():
	ANI = pd.read_csv(ani_path, sep="\t", header = None, names=['Genome1', 'Genome2', "ANI", "X1", "X2", "Taxonomy1", "Taxonomy2"])
	COUNTS = pd.read_csv(input_path, sep="\t")


	#Lots of pandas and buckle up because none of this is documented yet.
	ANI_COUNTS = ANI[ANI['Taxonomy1'].isin(COUNTS["Species"]) & ANI['Taxonomy2'].isin(COUNTS["Species"])]
	ANI_COUNTS = ANI_COUNTS[ANI_COUNTS["ANI"] > cut_off]

	COUNTS['total_counts'] = COUNTS.sum(axis=1)
	TOTAL_COUNTS = COUNTS[['Species', 'total_counts']]

	#### COV Location
	if cov_path:
			taxonomy_string = "Taxonomy1"
	# 	TOTAL_COUNTS = cov_filter(cov_path, taxonomy_string, ANI, TAXONOMY)
	# else:
	# 	taxonomy_string = "Species"
	COV = pd.read_csv(cov_path, sep="\t")
	TAXONOMY = ANI[["Genome1", taxonomy_string]]
	COV = pd.merge(TAXONOMY, COV, left_on='Genome1', right_on='#ID').drop_duplicates()
	COV = COV[[taxonomy_string, "Covered_percent"]]
	TOTAL_COUNTS = pd.merge(TOTAL_COUNTS, COV, how="outer", left_on='Species', right_on=taxonomy_string)

	# networkx graph loading using pandas dataframe format
	G = nx.from_pandas_edgelist(ANI_COUNTS, "Taxonomy1", "Taxonomy2", ['ANI'])
	partition, rev_partition = partitioning(G)

	#first compute the best partition
	# partition = community.best_partition(G)
	# rev_partition = {}
	# for key, value in partition.items():
	# 	rev_partition.setdefault(value, set()).add(key)

	#Create dataframe of clusters
	#COUNTS = collapse_centroids(COUNTS, TOTAL_COUNTS, rev_partition, taxonomy_string)

	clusters = pd.DataFrame(columns=['Species', 'total_counts', 'Covered_percent'])
	for key in rev_partition:
		cluster = TOTAL_COUNTS[TOTAL_COUNTS["Species"].isin(rev_partition[key])]
		genome_max_counts = cluster.loc[cluster['total_counts'].idxmax()][taxonomy_string]
		SUM_COUNTS = COUNTS[COUNTS["Species"].isin(rev_partition[key])].sum(axis=0)
		SUM_COUNTS["Species"] = genome_max_counts

		COUNTS = COUNTS[~COUNTS["Species"].isin(rev_partition[key])] # ~ not in
		COUNTS = COUNTS.append(SUM_COUNTS, ignore_index=True)

	COUNTS.sort_values(by=['Species'])
	COUNTS = COUNTS[COUNTS.total_counts > filter_cutoff]
	COUNTS = COUNTS.drop(['total_counts'], axis=1)
	COUNTS.to_csv(output_path, sep = "\t", index=False)

	#Call plotting of network in matplotlib if flagged
	if write_plot:
		plotting(partition, output_path, G)


#drawing
def plotting(partition, output_path, G):
	#show_plot = False #Debugging (and it looks cool)
	size = float(len(set(partition.values())))
	pos = nx.spring_layout(G)
	count = 0.
	for com in set(partition.values()) :
		count = count + 1.
		list_nodes = [nodes for nodes in partition.keys()
									if partition[nodes] == com]
		nx.draw_networkx_nodes(G, pos, list_nodes, node_size = 20,
									node_color = str(count / size))

	nx.draw_networkx_edges(G, pos, alpha=0.5)
	plt.savefig(output_path+".pdf")
	#if show_plot:
	#	plt.show()


def cov_filter(cov_path, taxonomy_string, ANI, TAXONOMY):
	COV = pd.read_csv(cov_path, sep="\t")
	TAXONOMY = ANI[["Genome1", taxonomy_string]]
	COV = pd.merge(TAXONOMY, COV, left_on='Genome1', right_on='#ID').drop_duplicates()
	COV = COV[[taxonomy_string, "Covered_percent"]]
	TOTAL_COUNTS = pd.merge(TOTAL_COUNTS, COV, how="outer", left_on='Species', right_on=taxonomy_string)
	return(TOTAL_COUNTS)

def partitioning(G):
	partition = community.best_partition(G)
	rev_partition = {}
	for key, value in partition.items():
		rev_partition.setdefault(value, set()).add(key)
	return(partition, rev_partition)


# def collapse_centroids(COUNTS, TOTAL_COUNTS, rev_partition, taxonomy_string):
# 	clusters = pd.DataFrame(columns=['Species', 'total_counts', 'Covered_percent'])
# 	for key in rev_partition:
# 		cluster = TOTAL_COUNTS[TOTAL_COUNTS["Species"].isin(rev_partition[key])]
# 		genome_max_counts = cluster.loc[cluster['total_counts'].idxmax()][taxonomy_string]
# 		SUM_COUNTS = COUNTS[COUNTS["Species"].isin(rev_partition[key])].sum(axis=0)
# 		SUM_COUNTS["Species"] = genome_max_counts

# 		COUNTS = COUNTS[~COUNTS["Species"].isin(rev_partition[key])] # ~ not in
# 		COUNTS = COUNTS.append(SUM_COUNTS, ignore_index=True)
# 		#output
# 		return(COUNTS)

# To Do: Change output from Species to Genome.
#Run Script ###########################################################
if __name__ == '__main__':
	main()
