#!/usr/bin/env python

#####################################################################
#
# process alignments 0.1
# Merge SAM files from bbmap
# Calculate SAM coverage with pileup.sh
# Process output:
#    1. Filter alignments (>0.5% genome coverage)
#    2. Assign GTDB taxonomy by genome accession
#    3. Collapse counts to taxonomy by sample table
# Built for use in Meiji Snakemake 0.1
#
######################################################################

__authors__ = "David Levy-Booth"
__copyright__ = "Copyright 2020, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

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


# rules #################################################################

rule merge_sams:
    input: expand("{OUTPUT_DIR}output/bbmap.{build}.sam", OUTPUT_DIR = OUTPUT_DIR, build = index_nodes)
    output: temp("{OUTPUT_DIR}processed/merged.sam".format(OUTPUT_DIR = OUTPUT_DIR))
    message: bcolors.OKBLUE + "\nRunning Meiji process_alignments module 0.1\n" + bcolors.ENDC
    params:
        input_directory="{OUTPUT_DIR}output/".format(OUTPUT_DIR = OUTPUT_DIR)
    script:
        "scripts/sam.merge.py"

rule pileup:
    input:  "{OUTPUT_DIR}processed/merged.sam".format(OUTPUT_DIR = OUTPUT_DIR)
    output: "{OUTPUT_DIR}processed/merged.cov".format(OUTPUT_DIR = OUTPUT_DIR)
    shell:
        """
        pileup.sh in={input} out={output}
        """

rule sam_coverage:
    input:
        "{OUTPUT_DIR}processed/merged.cov".format(OUTPUT_DIR = OUTPUT_DIR),
        "{OUTPUT_DIR}processed/merged.sam".format(OUTPUT_DIR = OUTPUT_DIR)
    output: temp("{OUTPUT_DIR}processed/merged.sam.covfilt".format(OUTPUT_DIR = OUTPUT_DIR))
    params:
        first_cov_cutoff=0.5
    script:
        "scripts/sam.coverage.py"

rule sam_taxonomy:
    input:
        "{OUTPUT_DIR}processed/merged.sam.covfilt".format(OUTPUT_DIR = OUTPUT_DIR),
        "taxonomy/TAXONOMY"
    output: temp("{OUTPUT_DIR}processed/merged.sam.covfilt.taxonomy".format(OUTPUT_DIR = OUTPUT_DIR))
    script:
        "scripts/sam.taxonomy.py"

rule compile_taxonomy:
    input: "{OUTPUT_DIR}processed/merged.sam.covfilt.taxonomy".format(OUTPUT_DIR = OUTPUT_DIR)
    output: "{OUTPUT_DIR}processed/merged.sam.covfilt.taxonomy.count".format(OUTPUT_DIR = OUTPUT_DIR)
    script:
        "scripts/compile.taxonomy.py"

rule cluster_taxonomy:
    input:
        "{OUTPUT_DIR}processed/merged.sam.covfilt.taxonomy.count".format(OUTPUT_DIR = OUTPUT_DIR),
        "{OUTPUT_DIR}processed/merged.cov".format(OUTPUT_DIR = OUTPUT_DIR),
        "taxonomy/ANI"
    output:
        "{OUTPUT_DIR}processed/merged.sam.covfilt.taxonomy.count.clustered".format(OUTPUT_DIR = OUTPUT_DIR)
    params:
        ANI_cutoff=94,
        count_cutoff=100,
        write_cluster_plot=False
    script:
        "scripts/cluster.taxonomy.py"
