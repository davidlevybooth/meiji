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
import os.path
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
READ_LIMIT = config["READ_LIMIT"]

# definitions ###########################################################
#
#Move this to central file parsing module
def read_format(SAMPLE_DIR):
    ext=""
    fq_list = os.listdir(SAMPLE_DIR)
    fq_list.sort()
    test_file = os.path.join(SAMPLE_DIR,fq_list[0])
    extentions = Path(test_file).suffixes
    if ".gz" in extentions:
        if get_compression_type(test_file) == 'gz':
            pass
        else:
            raise Exception(bcolors.FAIL + 'ExceptCompression: .gz files found but compression is invalid\n' + bcolors.ENDC)
    if any(x in [".fasta", ".fa"] for x in extentions):
        raise Exception(bcolors.FAIL + 'ExceptFormat: Fasta input files detected.\nProvide directory with fastq files (can be gzipped)\n' + bcolors.ENDC)
    if any(x in [".fastq", ".fq"] for x in extentions):
        fq_list = os.listdir(SAMPLE_DIR)
        fq_list.sort()
        ext = "".join(extentions)
        fq_list = [f.replace(ext, "") for f in fq_list]
    else:
        raise Exception(bcolors.FAIL + 'ExceptFormat: Unknown input file format.\nProvide directory with fastq files (can be gzipped)\n' + bcolors.ENDC)
    return(ext, fq_list)

def get_compression_type(filename):
	"""
	Attempts to guess the compression (if any) on a file using the first few bytes.
	http://stackoverflow.com/questions/13044562
	"""
	magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
				  'bz2': (b'\x42', b'\x5a', b'\x68'),
				  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
	max_len = max(len(x) for x in magic_dict)

	unknown_file = open(filename, 'rb')
	file_start = unknown_file.read(max_len)
	unknown_file.close()
	compression_type = 'plain'
	for file_type, magic_bytes in magic_dict.items():
		if file_start.startswith(magic_bytes):
			compression_type = file_type
	if compression_type == 'bz2':
		sys.exit('Error: cannot use bzip2 format - use gzip instead')
	if compression_type == 'zip':
		sys.exit('Error: cannot use zip format - use gzip instead')
	return compression_type


# file parsing ############################################################

ext, samples = read_format(SAMPLE_DIR)

# rules ##################################################################

rule bbduk:
    input: expand("{SAMPLE_DIR}{samples}{ext}", SAMPLE_DIR = SAMPLE_DIR, samples = samples, ext=ext)
    output: expand("{OUTPUT_DIR}reads/{samples}.fasta", OUTPUT_DIR = OUTPUT_DIR, samples = samples)
    threads: THREADS
    message: bcolors.OKBLUE + "\nRunning Meiji prepare_reads module 0.1\n" + bcolors.ENDC
    params:
        trimq=30,
        minlen=130,
        ftr=160,
        ftl=10
    run:
        for sample in samples:
            input_fastq="{SAMPLE_DIR}{sample}{ext}".format(SAMPLE_DIR=SAMPLE_DIR, sample = sample, ext=ext)
            output_fasta="{OUTPUT_DIR}reads/{sample}.fasta".format(OUTPUT_DIR=OUTPUT_DIR, sample = sample)
            shell("bbduk.sh -Xmx20g t={threads} in={input_fastq} out={output_fasta} qtrim=r trimq={params.trimq} minlen={params.minlen} ftl={params.ftl} ftr={params.ftr}")
            print("{params.readlimit}")
            if READ_LIMIT:
                temp_fasta=output_fasta.replace(".fasta", "_temp.fasta")
                shell("reformat.sh in={output_fasta} out={temp_fasta} samplereadstarget={READ_LIMIT}")
                if os.path.exists("{temp_fasta}".format(temp_fasta = temp_fasta)):
                    shell("rm {output_fasta}")
                    shell("mv {temp_fasta} {output_fasta}")
                else:
                    raise Exception(bcolors.FAIL + 'Read subsetting failed.\n' + bcolors.ENDC)


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
