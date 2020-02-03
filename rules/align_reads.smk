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

import time

THREADS = config["threads"]
OUTPUT_DIR = config["OUTPUT_DIR"]
INDEX_DIR = config["INDEX_DIR"]

# definitions #########################################################

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# definitions #########################################################
def database_collect(INDEX_DIR):
    index_list = os.listdir(INDEX_DIR+"/ref/index")
    index_list.sort()
    return(index_list)

# parse index #########################################################

index_nodes = database_collect(INDEX_DIR)
index_nodes.sort()
# print(bcolors.OKBLUE + "Aligning against:\n23458 bacterial genomes\n1248 archaeal genomes\nacross " + str(len(index_nodes)) + " index nodes"+ bcolors.ENDC)
# print(bcolors.OKBLUE + "This is going to take a while.\n"+ bcolors.ENDC)
# time.sleep(2)

# rules ###############################################################

rule bbmap:
    input: "{OUTPUT_DIR}reads/reads.qc.lin.fa".format(OUTPUT_DIR = OUTPUT_DIR)
    output: expand("{OUTPUT_DIR}output/bbmap.{build}.sam", OUTPUT_DIR = OUTPUT_DIR, build = index_nodes)
    params:
        index_dir=INDEX_DIR,
        fastareadlen=150,
        minidentity=0.97
    threads: THREADS
    message: bcolors.OKBLUE + "\nRunning Meiji align_reads module 0.1\n" + bcolors.ENDC
    run:
        #To do: Explicitly check if SAM outputs exist to avoid rewritting
        #Consider using protected files for SAM outputs
        for i in index_nodes:
            alignment = "{OUTPUT_DIR}output/bbmap.{i}.sam".format(OUTPUT_DIR = OUTPUT_DIR, i = i)
            shell("bbmap.sh path={params.index_dir} in={input} outm={alignment} t={threads} fastareadlen={params.fastareadlen} minidentity={params.minidentity} build={i}")
