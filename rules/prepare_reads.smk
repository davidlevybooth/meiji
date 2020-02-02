#!/usr/bin/env python

#####################################################################
#
# prepare_reads.smk 0.1
# Reads filtering: bbduk
# Read formatting: reshape.sh
# Read linearization: seqtk
#
# Built for use in Meiji Snakemake 0.1
#
######################################################################

__authors__ = "David Levy-Booth"
__copyright__ = "Copyright 2020, David Levy-Booth"
__email__ = "dlevyboo@mail.ubc.ca"
__license__ = "GPL3"

import os
from pathlib import Path
import shutil

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
SAMPLE_DIR = config["SAMPLE_DIR"]
OUTPUT_DIR = config["OUTPUT_DIR"]

print(bcolors.OKBLUE + "\nRunning Meiji prepare_reads module 0.1\n" + bcolors.ENDC)

# file parsing ############################################################

samples = os.listdir(SAMPLE_DIR)
samples = [f.replace(".fastq", "") for f in samples]


# rules ##################################################################

rule bbduk:
    input: "{SAMPLE_DIR}/{{samples}}.fastq".format(SAMPLE_DIR = SAMPLE_DIR)
    output: "{OUTPUT_DIR}/{{samples}}.fasta".format(OUTPUT_DIR = OUTPUT_DIR)
    threads: THREADS
    shell:
        """
        bbduk.sh -Xmx20g t={threads} in={input} out={output} qtrim=r trimq=30 minlen=130 ftl=10 ftr=160
        """
rule format_headers:
    input: expand("{OUTPUT_DIR}/{samples}.fasta", OUTPUT_DIR = OUTPUT_DIR, samples = samples)
    output: expand("{OUTPUT_DIR}/{samples}.head.fasta", OUTPUT_DIR = OUTPUT_DIR, samples = samples)
    run:
        for input_fasta in input:
            print("Refomatted Header: " + input_fasta)
            output_fasta = input_fasta.replace(".fasta", ".head.fasta")
            filename = os.path.splitext(input_fasta)[0]
            filename = os.path.basename(filename)
            rn = 1
            fo = open(output_fasta, "w")
            with open(input_fasta,"r") as fi:
                for ln in fi:
                    if ln.startswith('>'):
                        ln = ">" + filename + "_" + str(rn)+"\n"
                        rn += 1
                    fo.write(ln)
            os.remove(input_fasta)


rule merge_reads:
    input: expand("{OUTPUT_DIR}/{samples}.head.fasta", OUTPUT_DIR = OUTPUT_DIR, samples = samples)
    output: temp("{OUTPUT_DIR}/reads.qc.fa".format(OUTPUT_DIR = OUTPUT_DIR))
    run:
        shell("cat {input} > {output}")
        for input_head_fasta in input:
            os.remove(input_head_fasta)


rule linearize_reads:
    input: "{OUTPUT_DIR}/reads.qc.fa".format(OUTPUT_DIR = OUTPUT_DIR)
    output: "{OUTPUT_DIR}/reads.qc.lin.fa".format(OUTPUT_DIR = OUTPUT_DIR)
    shell: "seqkit seq -w 0 {input} > {output}"
