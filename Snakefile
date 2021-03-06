#Snakefile for the pipeline to perform alignment ad variant calling on WES or WGS,
#
#
# 27/01/2021
#
# Author: massimiliano [dot] Cocca [at] burlo [dot] trieste [dot] it
#import libraries
import pandas as pd
import pathlib
import io
import os
import re
from snakemake.exceptions import print_exception, WorkflowError
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.31.0")

##### load config file #####
# configfile: "config.yaml"
#read samplesheet
# samples_df = pd.read_table(samples, sep=" ", header=0, dtype='object')
samples_df = pd.read_table(config["samples"], sep=" ", header=0, dtype='object')
#get samples ids from tablesheet
sample_names = list(samples_df.SAMPLE_ID)
# we need file names for R1 and R2
R1 = [ os.path.splitext(os.path.splitext(os.path.basename(fq1))[0])[0] for fq1 in samples_df.fq1]
R2 = [ os.path.splitext(os.path.splitext(os.path.basename(fq2))[0])[0] for fq2 in samples_df.fq2]
# print(sample_names)
# print(R1)
# print(R2)
#get WGS calling intervals from GATK bundle
# samples_df = pd.read_table(config["wgs_gatk"], sep=" ", header=0, dtype='object')
# Define some variables
PROJ= config["proj_name"]
BASE_OUT=config["files_path"]["base_out"]

# call_intervals=expand(config["callable_intervals"]+"/wgs_calling_regions_{chr}.GRCh38.p13.interval_list", chr=config["call_chr"])
call_intervals=expand("wgs_calling_regions_{chr}.GRCh38.p13.interval_list", chr=config["call_chr"])
##### local rules #####
localrules: all

#set minimum snakemake version#
min_version("5.32.0")

# Definition of helper functions
# def get_vcfs():
# vcf_files = list(pd.read_table(config["vcf_files"]).vcf_file)
# vcf_names = dict([(ntpath.basename(vcf_file), vcf_file) for vcf_file in vcf_files])



##### Target rules #####

rule all:
    input:
        # [(BASE_OUT +"/"+ config["fastqc_pre_dir"] + "/{sample_fs1}_fastqc.html").format(sample_fs1=r1_strand) for r1_strand in R1]
        # expand(BASE_OUT + "/"+ config["rules"]["ubam_gen"]["out_dir"]+"/" + "{sample}_unmap.bam", sample=sample_names),
        # expand(BASE_OUT + "/"+ config["rules"]["bwa_mem"]["out_dir"]+"/" + "{sample}_map.bam", sample=sample_names)
        expand(BASE_OUT + "/" + config["rules"]["gatk_hap_caller"]["out_dir"] + "/{sample}/{sample}_{interval_name}_g.vcf.gz",sample=sample_names,interval_name=call_intervals),
        expand(BASE_OUT + "/"+ config["rules"]["stats"]["out_dir"] + "/{sample}_validate.txt", sample=sample_names),
        expand(BASE_OUT + "/"+ config["rules"]["stats"]["out_dir"] + "/{sample}_flagstat.txt", sample=sample_names),
        expand(BASE_OUT +"/" +config["rules"]["stats"]["out_dir"] + "/{sample}_ismetrics.txt", sample=sample_names),
        expand(BASE_OUT +"/" +config["rules"]["stats"]["out_dir"] + "/{sample}_stats.txt", sample=sample_names),
        expand(BASE_OUT +"/" +config["rules"]["stats"]["out_dir"] + "/{sample}_wgsmetrics.txt", sample=sample_names)

        # [(BASE_OUT + "/"+ config["rules"]["ubam_gen"]["out_dir"]+"/" + "{sample}_unmap.bam").format(sample=sample_id) for sample_id in sample_names]
        # BASE_OUT + config["rules"]["bwa_mem"]["out_dir"] + "{sample}_map.bam"

##### Modules #####
include_prefix="rules"
include:
    include_prefix + "/functions.py"
include:
    include_prefix + "/alignment.smk"
include:
    include_prefix + "/bqsr.smk"
include:
    include_prefix + "/stats.smk"
include:
    include_prefix + "/var_call.smk"
# include:
#     include_prefix + "/picard.smk"
# include:
#     include_prefix + "/picard_stats.smk"
# include:
#     include_prefix + "/call_variants.smk"
# include:
#     include_prefix + "/joint_call.smk"
# include:
#     include_prefix + "/qc.smk"
# include:
#     include_prefix + "/vqsr.smk"
# include:
#     include_prefix + "/identity_check.smk"
# include:
#     include_prefix + "/coverage.smk"
# include:
#     include_prefix + "/mtdna.smk"
