#!/usr/bin/env python

#####################################################################
#
# prepare_reads.smk 0.1
# Reads filtering: bbduk
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


# definitions ###########################################################

#Move this to central file parsing module
def read_format(SAMPLE_DIR):
    fq_list = os.listdir(SAMPLE_DIR)
    fq_list.sort()
    if ".gz" in fq_list[0]:
        raise Exception(bcolors.FAIL + 'gzipped input files detected.\nProvide directory with gunzipped fastq files\n' + bcolors.ENDC)
    if ".fasta" in fq_list[0]:
        raise Exception(bcolors.FAIL + 'fasta input files detected.\nProvide directory with gunzipped fastq files\n' + bcolors.ENDC)
    if ".fastq" or ".fq" in fq_list[0]:
        print(bcolors.OKGREEN + 'fastq input files detected.\nContinuing with Meiji pipeline\n' + bcolors.ENDC)
    fq_list = [f.replace(".fastq", "") for f in fq_list]
    fq_list = [f.replace(".fq", "") for f in fq_list]
    return(fq_list)

# file parsing ############################################################

samples = read_format(SAMPLE_DIR)

# rules ##################################################################

rule bbduk:
    input: "{SAMPLE_DIR}/{{samples}}.fastq".format(SAMPLE_DIR = SAMPLE_DIR)
    output: "{OUTPUT_DIR}/reads/{{samples}}.fasta".format(OUTPUT_DIR = OUTPUT_DIR)
    threads: THREADS
    message: bcolors.OKBLUE + "\nRunning Meiji prepare_reads module 0.1\n" + bcolors.ENDC
    params:
        trimq=30,
        minlen=130,
        ftr=160,
        ftl=10
    run:
        shell("bbduk.sh -Xmx20g t={threads} in={input} out={output} qtrim=r trimq={params.trimq} minlen={params.minlen} ftl={params.ftl} ftr={params.ftr}")

rule format_headers:
    input: expand("{OUTPUT_DIR}reads/{samples}.fasta", OUTPUT_DIR = OUTPUT_DIR, samples = samples)
    output: expand("{OUTPUT_DIR}reads/{samples}.head.fasta", OUTPUT_DIR = OUTPUT_DIR, samples = samples)
    run:
        try:
            for input_fasta in input:
                print(bcolors.OKGREEN + "Refomatted Header: " + input_fasta + bcolors.ENDC)
                output_fasta = input_fasta.replace(".fasta", ".head.fasta")
                filename = os.path.splitext(input_fasta)[0]
                filename = os.path.basename(filename)
                if "_" in filename:
                    raise Exception(bcolors.FAIL + 'filenames should not contain underscores "_"\n' + bcolors.ENDC)
                rn = 1
                fo = open(output_fasta, "w")
                with open(input_fasta,"r") as fi:
                    for ln in fi:
                        if ln.startswith('>'):
                            ln = ">" + filename + "_" + str(rn)+"\n"
                            rn += 1
                        fo.write(ln)
                os.remove(input_fasta)
        except:
            print(bcolors.FAIL + 'Header formatting failed.' + bcolors.ENDC)


rule merge_reads:
    input: expand("{OUTPUT_DIR}reads/{samples}.head.fasta", OUTPUT_DIR = OUTPUT_DIR, samples = samples)
    output: temp("{OUTPUT_DIR}reads/reads.qc.fa".format(OUTPUT_DIR = OUTPUT_DIR))
    run:
        shell("cat {input} > {output}")
        for input_head_fasta in input:
            os.remove(input_head_fasta)

rule linearize_reads:
    input: "{OUTPUT_DIR}reads/reads.qc.fa".format(OUTPUT_DIR = OUTPUT_DIR)
    output: "{OUTPUT_DIR}reads/reads.qc.lin.fa".format(OUTPUT_DIR = OUTPUT_DIR)
    shell: "seqkit seq -w 0 {input} > {output}"
