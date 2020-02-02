#!/usr/bin/env python

#####################################################################
#
# align_reads 0.1
# Reads alignment: bbmap
# Database: GTDB 89
#
# Built for use in Meiji Snakemake 0.1
#
######################################################################

__authors__ = "David Levy-Booth"
__copyright__ = "Copyright 2020, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"

THREADS = config["threads"]
OUTPUT_DIR = config["OUTPUT_DIR"]
INDEX_DIR = config["INDEX_DIR"]

# definitions #########################################################

def database_collect(INDEX_DIR):
    index_list = os.listdir(INDEX_DIR+"/ref/index")
    index_list.sort()
    return(index_list)

# parse index #########################################################

index_nodes = database_collect(INDEX_DIR)

# rules ###############################################################

rule bbmap:
    input: "{OUTPUT_DIR}reads/reads.qc.lin.fa".format(OUTPUT_DIR = OUTPUT_DIR)
    output: expand("{OUTPUT_DIR}output/bbmap.{{build}}.sam", OUTPUT_DIR = OUTPUT_DIR)
    params:
        fastareadlen=150,
        minidentity=0.97,
        index=index_nodes,
        index_dir=INDEX_DIR
    threads: THREADS
    run:
        for i in params[2]:
            shell("bbmap.sh path={params.index_dir} in={input} out={output} t={threads} fastareadlen={params.fastareadlen} minidentity={params.minidentity} build={i}")
