#!/usr/bin/env python

#####################################################################
#
# Welcome the Meiji Snakemake pipeline. This suite of scripts
# is designed around the idea that by aligning short reads to a
# comprehensive set of species-level genomes, we can quickly and
# accurately quantify the taxonomic and functional diversity of the
# microbiome, regardless of environmental origin. Species-level (ANI >
# 95%) genomes and their revised taxonomic classifications were sourced
# from the Genome Taxonomy Database (GTDB) https://gtdb.ecogenomic.org/
# The bbmap suite of tools developed by the Joint Genome Institute
# https://jgi.doe.gov/data-and-tools/bbtools/ is used for ultra-fast
# and accurate read mapping to the GTDB genome database, and for
# coverage calculation.
#
# Mapping to sequenced genomes also provides the advantage of having
# annotated functional data at our fingertips. Here, the functional
# profiling module from the shallow-shotgun metagenome profiling tool
# SHOGUN https://github.com/knights-lab/SHOGUN is used (hence the name
# Meiji). Functional annotations for each genome are pulled from
# ANNOTREE http://annotree.uwaterloo.ca/, which are created using
# UniRef100 clusters  https://www.uniprot.org/, KEGG https://www.kegg.jp/
# KOs and PFAM http://pfam.xfam.org/ annotations.
#
# This work is in its earliest developmental iteration (0.1). Should you
# have any issues or concerns, please contact David Levy-Booth,
# dlevyboo@mail.ubc.ca -- And thank you for choosing Meiji for your
# metagenomic analysis.
#
#
#####################################################################

__authors__ = "David Levy-Booth"
__copyright__ = "Copyright 2020, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"

import os
from pathlib import Path

THREADS = config["threads"]
OUTPUT_DIR = config["OUTPUT_DIR"]
SAMPLE_DIR = config["SAMPLE_DIR"]


# definitions ###########################################################

#Move this to central file parsing module
def read_format(SAMPLE_DIR):
    fq_list = os.listdir(SAMPLE_DIR)
    fq_list.sort()
    fq_list = [f.replace(".fastq", "") for f in fq_list]
    return(fq_list)

# file parsing ##############################################################

input_samples = read_format(SAMPLE_DIR)

# rule all ##################################################################

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        "{OUTPUT_DIR}/reads.qc.lin.fa".format(OUTPUT_DIR = OUTPUT_DIR)
        #expand("{OUTPUT_DIR}/{filtered_samples}", filtered_samples = input_samples, OUTPUT_DIR = OUTPUT_DIR)
        # "{OUTDIR}/merged.sam.covfilt.taxonomy.count.tsv",
        # "{OUTDIR}/clustered/merged.sam.covfilt.taxonomy.count.clustered.tsv",
        # "{OUTDIR}/merged.sam.covfilt.taxonomy.count.species.kegg.txt"

# rules ##################################################################

include: "rules/prepare_reads.smk"
# include: "rules/align_reads.smk"
# include: "rules/process_alignment.smk"
# include: "rules/filter_counts.smk"
# include: "rules/assign_function.smk"