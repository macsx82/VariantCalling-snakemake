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
min_version("5.32.0")

##### load config file #####
# configfile: "config.yaml"

# define scatter parameter
SCATTER_COUNT = config['rules']['split_intervals']['scatter_count']
scattergather:
    split=SCATTER_COUNT

#read samplesheet
# samples_df = pd.read_table(samples, sep=" ", header=0, dtype='object')
samples_df = pd.read_table(config["samples"], sep=" ", header=0, dtype='object')
#get samples ids from tablesheet
sample_names = list(samples_df.SAMPLE_ID)
# we need file names for R1 and R2
R1 = [ os.path.splitext(os.path.splitext(os.path.basename(fq1))[0])[0] for fq1 in samples_df.fq1 ]
R2 = [ os.path.splitext(os.path.splitext(os.path.basename(fq2))[0])[0] for fq2 in samples_df.fq2 ]
# print(sample_names)
# print(R1)
# print(R2)
#get WGS calling intervals from GATK bundle
# samples_df = pd.read_table(config["wgs_gatk"], sep=" ", header=0, dtype='object')
# Define some variables
PROJ= config["proj_name"]
BASE_OUT=config["files_path"]["base_out"]

# call_intervals=expand(config["callable_intervals"]+"/wgs_calling_regions_{chr}.GRCh38.p13.interval_list", chr=config["call_chr"])
#this list of intervals is the complete liste we have to use when we work with all samples. If we are performing variant caling, we have to use the smartest get_intervals_by_sex function,
#which allows us to define the interval list based on sample's sex
call_intervals=expand("wgs_calling_regions_{chr}.GRCh38.p13.interval_list", chr=config["call_chr"])

##### 
# since we want to use this pipeline in batch mode but also in joint call mode, we will try to use the value 
# of the pipe_mode declaration in the config file. In BATCH mode, we will execute only some rules
# in JOINT_CALL mode, we will execute only the relevant rules

pipeline_mode = config["pipe_mode"]

##### functions #####
include_prefix="rules"
include:
    include_prefix + "/functions.py"

if pipeline_mode == "BATCH":
    print(pipeline_mode)
    ##### local rules #####
    localrules: all

    ##### Target rules #####
    rule all:
        input:
            call_variants_by_sex(BASE_OUT + "/" + config["rules"]["gatk_hap_caller"]["out_dir"]),
            expand(BASE_OUT + "/"+ config["rules"]["stats"]["out_dir"] + "/{sample}_validate.txt", sample=sample_names),
            expand(BASE_OUT + "/"+ config["rules"]["stats"]["out_dir"] + "/{sample}_flagstat.txt", sample=sample_names),
            expand(BASE_OUT +"/" +config["rules"]["stats"]["out_dir"] + "/{sample}_ismetrics.txt", sample=sample_names),
            expand(BASE_OUT +"/" +config["rules"]["stats"]["out_dir"] + "/{sample}_stats.txt", sample=sample_names),
            expand(BASE_OUT +"/" +config["rules"]["stats"]["out_dir"] + "/{sample}_wgsmetrics.txt", sample=sample_names),
            call_variants_by_sex(config["files_path"]["base_joint_call_path"])

else :
    print(pipeline_mode)
    ##### local rules #####
    localrules: all

    ##### Target rules #####

    rule all:
        input:
            #stuff needed for joint call
            expand(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genotype_gvcfs").get("out_dir"),"all.{interval_name}.vcf.gz"),interval_name=call_intervals),
            expand(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genotype_gvcfs").get("out_dir"),"all.{interval_name}.vcf.gz.tbi"),interval_name=call_intervals),
            expand(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genotype_gvcfs").get("out_dir"),"all.{interval_name}.stats"),interval_name=call_intervals),
            os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("pre_vqsr_rules").get("out_dir"),"ALL.CLEAN.SITES_ONLY.vcf.gz"),
            os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("pre_vqsr_rules").get("out_dir"),"ALL.CLEAN.SITES_ONLY.vcf.gz.tbi"),
            # expand(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("rsid_annotation").get("out_dir"),"/{chr}.PASS_rsID.vcf.gz"),chr=config["call_chr"]),
            # expand(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("rsid_annotation").get("out_dir"),"/{chr}.PASS_rsID.vcf.gz.tbi"),chr=config["call_chr"])
            expand(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("phase").get("out_dir"),"{interval_name}.phased.vcf.gz"),interval_name=call_intervals),
            expand(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("phase").get("out_dir"),"{interval_name}.phased.vcf.gz.tbi"),interval_name=call_intervals),
            expand(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("phase").get("out_dir"),"all.{interval_name}.stats"),interval_name=call_intervals),
            #point to the final merged output file. This final file is the one that has to go through the beautify/last qc pipeline
            os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("AllChrConcat").get("out_dir"),"all."+ PROJ+".vcf.gz"),
            os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("AllChrConcat").get("out_dir"),"all."+ PROJ+".vcf.gz.tbi")


##### Modules #####
if pipeline_mode == "BATCH" :
    include:
        include_prefix + "/alignment.smk"
    include:
        include_prefix + "/bqsr.smk"
    include:
        include_prefix + "/var_call.smk"
    include:
        include_prefix + "/stats.smk"
else :
    include:
        include_prefix + "/joint_call.smk"
    include:
        include_prefix + "/pre_vqsr_clean.smk"
    include:
        include_prefix + "/vqsr.smk"    
    include:
        include_prefix + "/variant_annotation.smk"    
    include:
        include_prefix + "/vcf_stats.smk"

onsuccess:
    print("The workflow finished without errors!")

onerror:
    print("An error occurred in the current workflow execution!!")