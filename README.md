# Variant Calling pipeline

A Snakemake implementation of the variant calling pipeline

This pipeline is designed to perform Variant calling on WGS data, starting from preprocessed fastq files. A pipeline to perform the necessary preprocessing steps is available [here](https://gitlab.burlo.trieste.it/max/SeqPreproc-snakemake/-/tree/main)


---

## Setting things up

In order to run the pipeline, there are some requirements to fullfill and some set up needs to be perfomed.
In this current version, the pipeline is tested and configured to be run on the [ORFEO cluster](https://orfeo-documentation.readthedocs.io/en/latest/) .
It is possible to run the pipeline on the Apollo cluster, but it will require to manually specify the location of all software binaries in the provided config file.

### Required Software

The following software has to be installed system-wide, in a user-defined Conda environment or using the modules architecture (ORFEO cluster).

+ awk
+ sed
+ python3
+ [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
+ [bcftools](http://www.htslib.org/doc/)
+ [samtools](http://www.htslib.org/doc/)
+ [htslib](http://www.htslib.org/doc/)
+ [sambamba](https://lomereiter.github.io/sambamba/index.html)
+ [GATK](https://gatk.broadinstitute.org/hc/en-us)
+ [Picard](https://gatk.broadinstitute.org/hc/en-us)
+ R
+ git
+ snakemake

At the moment, the load directive of each module on the ORFEO cluster is hard-coded in the pipeline code, with the form "module-name/version". Modules and versions used are:

+ gatk/4.1.9.0
+ gatk/4.2.2.0
+ picard/2.24.0
+ bcftools/1.14
+ samtools/1.14
+ sambamba/0.8.0
+ R/4.0.3
+ bwa-mem2/2.1/gnu/9.3.0

**Before switching to a new version of each software/module, a test run should be performed to check that the expected output files are generated and that they are consistent with the previous production version.**

### Required python packages

In order to run the pipeline, the following python packages have to be installed in your conda environment:

+ pandas
+ pathlib
+ io
+ os
+ re
+ snakemake
+ gzip 
+ re
+ sys
+ errno
+ multiprocessing
+ psutil


### ORFEO/general set up
1. Install Snakemake via conda ([link](https://snakemake.readthedocs.io/en/stable/getting\_started/installation.html));
    ```bash
    conda create -c conda-forge -c bioconda -n snakemake snakemake pandas
    ```
2. Activate the environment you created
    ```bash
    conda activate snakemake
    ```

### Apollo set up

1. Add the global snakemake environment to your environment list:
    ```bash
    conda config --append envs_dirs /shared/software/conda/envs
    conda config --prepend envs_dirs ~/.conda/envs
    ```

2. Check that the environment is available to you (you should see an entry "snakemake_g" in the list)
    ```bash
    conda env list
    ```
3. Load the environment
    ```bash
    conda activate snakemake_g
    ```
---

## Resources setup

In order to run the pipeline on a new system, there are some preparatory steps to perform, in order to retrieve or generate the resources needed.

At the moment, all resources needed are already available on the APOLLO cluster, at the following locations:

```
references:
  basepath: "/shared/resources"
  provider: "hgRef"
  release: "GRCh38.p13"

genome_fasta: "GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna"
callable_intervals: "wgs_intervals"

dbsnp_latest: "/shared/resources/dbSNP/human_9606_b154_GRCh38p12/GCF_000001405.38.vcf.gz"
dbsnp: "/shared/resources/gatk4hg38db/dbsnp_151.hg38.vcf.gz"
hapmap: "/shared/resources/gatk4hg38db/hapmap_3.3.hg38.vcf.gz"
g1k: "/shared/resources/gatk4hg38db/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
omni: "/shared/resources/gatk4hg38db/1000G_omni2.5.hg38.vcf.gz"
mills: "/shared/resources/gatk4hg38db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
indels: "/shared/resources/gatk4hg38db/Homo_sapiens_assembly38.known_indels.vcf.gz"
axiom: "/shared/resources/gatk4hg38db/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"

```

Below are some indication on how to retrieve the data from external resources.

### GATK4 data bundle

All the files in the folder */shared/resources/gatk4hg38db/* are available to download from the gatk4 [gcp-public-data--broad-references/hg38/v0](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0)


### Genome reference

We used [GRCh38.p13](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/), as default reference for alignment and variant calling.
Since the aligned data will be saved using the **CRAM** format, it is necessary to have a backup copy of the **exact file** of the reference used.


Since the pipeline uses the same interval list approach used with the standard GATK bundle, we re-generated the callable regions interval list.

We followed this workflow:

1- Use *picard ScatterIntervalsByNs* to extract the non-poly N regions.
2- Convert this interval list to bed with awk.
3- Run *bedtools makewindows* to generate your regions.
	(from https://www.biostars.org/p/473462/#473466 )


The example command below are executed on the ORFEO cluster. Please modify accordingly the job submission and the relevant paths.

```bash
base_out=/storage/burlo/cocca/resources/hgRef/GRCh38.p13/wgs_intervals
ref_seq=/storage/burlo/cocca/resources/hgRef/GRCh38.p13/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna
log_folder=/large/___HOME___/burlo/cocca/analyses/Logs

echo "module load picard/2.24.0;java -jar $picard ScatterIntervalsByNs R=${ref_seq} OT=ACGT O=${base_out}/wgs_calling_regions.GRCh38.p13.interval_list" | qsub -q thin -V -k eod -o ${log_folder}/picard_split_interval.o -e ${log_folder}/picard_split_interval.e -l nodes=1:ppn=1 -l walltime=48:00:00 -l mem=10g

```

The Variant calling with GATK can be performed adding the ploidy information for each sample, to perform the variant calling correctly on the X and Y chromosome.
For this purpose, the PAR and NON PAR regions on chrX have to be treated separately. 
We need to create a bed file with PAR1 and PAR2 regions cohordinates. The file has to be TAB separated:

```
chrX	10001	2781479	PAR1                
chrX	155701383	156030895	PAR2
```

Using the par regions bed file, we can generate intervals for PAR regions on X chromosome:

```bash
(fgrep "@" /storage/burlo/cocca/resources/hgRef/GRCh38.p13/wgs_intervals/wgs_calling_regions_chrX.GRCh38.p13.interval_list;bedtools subtract -a <(fgrep -v "@" /storage/burlo/cocca/resources/hgRef/GRCh38.p13/wgs_intervals/wgs_calling_regions_chrX.GRCh38.p13.interval_list) -b /storage/burlo/cocca/resources/hgRef/GRCh38.p13/wgs_intervals/par_regions.bed )> /storage/burlo/cocca/resources/hgRef/GRCh38.p13/wgs_intervals/wgs_calling_regions_chrX_NONPAR.GRCh38.p13.interval_list
(fgrep "@" /storage/burlo/cocca/resources/hgRef/GRCh38.p13/wgs_intervals/wgs_calling_regions_chrX.GRCh38.p13.interval_list;bedtools intersect -a <(fgrep -v "@" /storage/burlo/cocca/resources/hgRef/GRCh38.p13/wgs_intervals/wgs_calling_regions_chrX.GRCh38.p13.interval_list) -b /storage/burlo/cocca/resources/hgRef/GRCh38.p13/wgs_intervals/par_regions.bed -wb|cut -f -5,9) > /storage/burlo/cocca/resources/hgRef/GRCh38.p13/wgs_intervals/wgs_calling_regions_chrX_PAR.GRCh38.p13.interval_list

```

**We need to manually modify the NON PAR regions start and end, since the bedtools operation leave the data with a 0-based position.**

i.e. the NON PAR file will have the following cohordinates for the first interval:

```
chrX    2781479 37099262        +       ACGTmer
```
Since the START position (2781479) is overlapping with the END position of the PAR1 region (2781479), we need to shift by 1 base, that cohordinate.

Same for the last interval:

```
chrX    144475607       155701383       +       ACGTmer
```

We need to fix the overlap between the END position (155701383) and the START position of the PAR2 region (155701383), shifting by 1 base.


### Other requirements

Absolute paths of the fastq files (saparate files for R1 and R2 strand) to be processed, produced by the [SeqPreprocessing](https://gitlab.burlo.trieste.it/max/SeqPreproc-snakemake/-/tree/main) pipeline.

Manifest file following the template provided in the **resource** folder, space separated and with the following header line:

```
SAMPLE_ID fq1 fq2 sex
```


At this point, we can proceed with the pipeline setup, the same we did for the SeqPreprocessing pipeline.


---

# Usage

It is advised to clone the pipeline to a personal folder from the git repository, in order to be able to correctly execute the workflow.
The clone command to use is:

```bash
	git clone https://gitlab.burlo.trieste.it/max/VariantCalling-snakemake.git
```

This command will create a new folder, in the current location, named **VariantCalling-snakemake**


Before running the pipeline, it is necessary to provide some information using a config file. A config file template, with default values is provided and it is named **config.yaml**

It is advisable not to execute the pipeline in the repository folder, but to create a separate folder and a copy of the config file to fill the relevant parameters.


## Config file setup

A config file template, with default values is provided and it is named **config.yaml** . The final user has to provide some information to the pipeline, in order to retrieve input files and to define output location.

### Execution modes

The pipeline can be executed in different modes:

- **BATCH**
- **JOINT_CALL**
- **JOINT_CALL_UPDATE**


**BATCH mode** is meant to be used to run the pipeline splitting the samples in different batches and will process each sample separately from alignment to HaplotypeCaller rules.
 **On new samples, this is the running mode that has to be used to generate the initial variant call data for each sample.**


**JOINT_CALL mode** represents the second step, where we collect all gVCF and perform :

- joint variant calling
- VQSR filtering
- variant annotation

 
**JOINT_CALL_UPDATE** mode is meant to be used to add new samples to an existing callset, **after performing a BATCH mode run**, using the 'ImportDB' update mode and performing :

- joint variant calling
- VQSR filtering
- variant annotation

To define the execution mode, we can set the parameter **pipe_mode** in the config file:

```bash
##############################
pipe_mode : "BATCH" #change to: "JOINT_CALL" to work in joint call mode, "JOINT_CALL_UPDATE" to add new samples to an existing callset.

```

### Samples manifest and project name

The user has to provide the path to a manifest file, generated as explained in the section **Other requirements** .
It is advisable to provide a meaningful project name, in order to distinguish different pipeline runs. In joint_call mode and joint_call_update mode, the project name will be used to generate the name of the final vcf file.

```bash
########### INPUT DEFINITION ####################
samples: "VARCALL_SAMPLE_FILE_PATH"   # Manifest file with the list of sample ids, relative fastq files to process and sex for each sample
proj_name : "WGS_VarCalling" #name of the project to be used also as suffix for file naming in the final concat rule

```

### Read Group parameters

After defining the manifest file and project name, the parameters **lb**, **pl** and **cl** are set up. Those parameters are used to fill the Read Group (RG) tag in the aligned bam/cram files. In the config file are provided default values that can be used. In case of different platform or in case that the library used is known, the user should fill the values accordingly.

```bash
#########################
#Variables to be used for RG creation and compression level for bam files.
#These variables can be left to this defaults values or set to the correct ones, if known.
#
lb : "L001"      #library name
pl : "Illumina"  #sequencing platform
cl : 5           #bam compression level
```

### References and intervals

This section of the config file contains default values for the execution on the ORFEO cluster. At the moment, the paths provided point to a shared folder under the user "cocca", but they should be changed to a group shared path. The alternative, on ORFEO cluster, is to generate a resource folder for each user.

```bash
################### REFERENCES AND INTERVALS #####################################
#paths and subfolders containing data we need for the pipeline to run
references:
  basepath: "/storage/burlo/cocca/resources"
  provider: "hgRef"
  release: "GRCh38.p13"

genome_fasta: "GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna"
callable_intervals: "wgs_intervals"
#################################################################################
known_variants:
  dbsnp_latest: "/storage/burlo/cocca/resources/dbSNP/human_9606_b154_GRCh38p12/GCF_000001405.38.vcf.gz"
  dbsnp: "/storage/burlo/cocca/resources/gatk4hg38db/dbsnp_151.hg38.vcf.gz"
  hapmap: "/storage/burlo/cocca/resources/gatk4hg38db/hapmap_3.3.hg38.vcf.gz"
  g1k: "/storage/burlo/cocca/resources/gatk4hg38db/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  omni: "/storage/burlo/cocca/resources/gatk4hg38db/1000G_omni2.5.hg38.vcf.gz"
  mills: "/storage/burlo/cocca/resources/gatk4hg38db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  indels: "/storage/burlo/cocca/resources/gatk4hg38db/Homo_sapiens_assembly38.known_indels.vcf.gz"
  axiom: "/storage/burlo/cocca/resources/gatk4hg38db/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"


```

### Chromosome definition

The user can specify the chromosomes to be processed. Regarding chrX , to correctly perform the variant calling, the user has to use the notation **"chrX_PAR"**,**"chrX_NONPAR"**.

```bash
#################################################################################
# Subset of chromosomes to call
call_chr: ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX_PAR","chrX_NONPAR","chrY"]
############################################

```

### Rules parameters

This section of the config file contains parameters that should be modified only if the relevant rules o tools change because of new software version.

### Paths definition

The bottom section of the config file allow the user to specify the output paths and the paths of each tool used.

Regarding output paths, if a relative path is provided, all folders will be generated in a path relative to the current pipeline execution folder.
The **base_joint_call_path** parameter will be used to generate some symlink to allow for the joint variant calling mode to correctly retrieve the input files.

```bash
############# PATHS #############
### - OUTPUT PATHS - ###
files_path:
  base_out: "VARCALL_OUT_FOLDER"
  tmp: "localtemp"
  log_dir: "Log"
  benchmark: "benchmarks"
  base_joint_call_path : "" #this folder should be located in a different path, so that can be accessed by all concurrent batch instances of the pipeline, to create symlink to gvcf files. Than, it could be read in joint call mode and used to perform DbImport and all other steps
```

The PATH TOOL section contains path of the binary used in the pipeline. This path should be changed only if the software used are not available in the user **$PATH** , following a custom installation, for example.

```bash
### - PATH TOOL - ###
ALIGN_TOOL: "bwa-mem2"
PHASE_TOOL: "eagle"
PICARD_TOOL: "picard"
SAMBAMBA_TOOL: "sambamba"
SAMTOOLS: "samtools"
BCFTOOLS: "bcftools"
GATK_TOOL: "gatk"
BEDTOOLS: "bedtools"
```

---
## Details on Execution modes

As explained above, the pipeline can be executed in different modes:

- **BATCH**
- **JOINT_CALL**
- **JOINT_CALL_UPDATE**


**BATCH mode** is always the mode to choose if you are processing samples from scratch, to perform alignment and variant calling using GATK HaplotypeCaller. With this mode the pipeline will work by sample.
If you plan to perform joint variant calling, after this first step, you have to specify a working location for the pipeline, using the parameter "base_joint_call_path" in the config file.
The value of this parameter should be always the same, for each batch mode pipeline run, if you plan to perform joint calling on all your processed samples.
For example, if you have 30 samples splitted in 3 batches of 10 samples each and you process them with 3 different runs of the pipeline for alignment and variant calling, you should specify the same location in the "base_joint_call_path" parameter.
**The parameter "base_joint_call_path" cannot be left with the default setting (empty string), otherwise the pipeline will produce an error.**


**JOINT_CALL mode** is ideally the second step in a population-based callset.
To activate this mode, you should set the "pipe_mode" parameter accordingly. It is advisable to run the pipeline in a different working folder with respect to the BATCH mode runs.
For example:
We run 3 batches of samples in BATCH mode in the following paths:

```bash
${HOME}/WGS_HC/20220216_VARCALL/BATCHES/BATCH000
${HOME}/WGS_HC/20220216_VARCALL/BATCHES/BATCH001
${HOME}/WGS_HC/20220216_VARCALL/BATCHES/BATCH002
```

for each batch we have a config file with all the relevant parameters set and, for all batches, the "base_joint_call_path" parameter points at:

```bash
  base_joint_call_path : "/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/20220216"

```

This way all the processed samples will have a symlink to the joint variant calling folder, so that the pipeline can perform all the subsequent steps.
This folder will be the base folder for the creation of the DBImport database, in the default folder named as "6.DbImport".
The "6.DbImport" folder will be referenced also when using the JOINT_CALL_UPDATE mode.

To start a JOINT_CALL run, we will create a new folder, for example in:

```bash
${HOME}/WGS_HC/20220216_JOINT_CALL
```

And we will need to create a new config file, using as template one of the config files generated for the BATCH mode. This way we will have all the parameters correctly set, even the mandatory "base_joint_call_path".
The user should also provide a new manifest file, containing ALL the processed samples information. This could be done merging together the manifest files of different batches, i.e.:

```bash

(echo "SAMPLE_ID fq1 fq2 sex";cat ${HOME}/WGS_HC/20220216_VARCALL/BATCHES/BATCH000/manifest_file_BATCH000.txt ${HOME}/WGS_HC/20220216_VARCALL/BATCHES/BATCH001/manifest_file_BATCH001.txt ${HOME}/WGS_HC/20220216_VARCALL/BATCHES/BATCH002/manifest_file_BATCH002.txt | fgrep -v "SAMPLE_ID" ) > ${HOME}/WGS_HC/20220216_JOINT_CALL/manifest_file_JOINT.txt
```

That will include all 30 samples already processed that will undergo joint calling procedures and will be the new value of the "samples" parameter in the new config file.


In **JOINT_CALL_UPDATE** mode, the pipeline will exploit a feature introduced by GATK v4.2+ that allow to UPDATE the DBImport object. This way, you will save some time, since you will perform only the GVCFGenotype step, on the whole callset.
This mode should be chosen if you plan to add samples on a population-based callset. As a first step you should perform a standard BATCH mode run, to align and generate the GVCF files for each new sample.
After the batch mode run, similarly to the JOINT_CALL mode, it is advisable to create a new working folder for the pipeline to run, and a new config file, specifying the new pipeline mode.
In addition, you have to add the new samples to the manifest file with all the existing samples in your call set. If we want to add a 4th batch to our callset, for example, we will run:

```bash
mkdir -p ${HOME}/WGS_HC/20220216_JOINT_CALL_UPDATE
cp ${HOME}/WGS_HC/20220216_JOINT_CALL/config.yaml ${HOME}/WGS_HC/20220216_JOINT_CALL_UPDATE/config.yaml

(echo "SAMPLE_ID fq1 fq2 sex";cat ${HOME}/WGS_HC/20220216_JOINT_CALL/manifest_file_JOINT.txt ${HOME}/WGS_HC/20220216_VARCALL/BATCHES/BATCH003/manifest_file_BATCH003.txt | fgrep -v "SAMPLE_ID" ) > ${HOME}/WGS_HC/20220216_JOINT_CALL_UPDATE/manifest_file_JOINT_UPDATED.txt
```

Then we need to point the "samples" parameter to the new manifest file. In this case, there will be no need to modify anything else, since all the parameter needed have been already set by the previous joint calling run.

---
## Output

The pipeline produces different outputs by sample, when in BATCH mode, and by cohort, when in JOINT_CALL mode.

For a BATCH mode run, the output folder structure will look as follow:

```
.
├── 1.ALIGNMENT
├── 2.BQSR
├── 3.VarCall
├── 4.Stats
└── localtemp

```

The folder "3.VarCall" will contain a subfolder for each sample, with the GVCF files generated by GATK HaplotypeCaller. Those files will be referenced to when in JOINT_CALL mode.

When in JOINT_CALL mode, the folder structure will be as follow:

```
.
├── [SAMPLE_ID1]
├── [SAMPLE_ID2]
├── ....
├── [SAMPLE_IDn]
├── 5.SplittedIntervals
├── 6.DbImport
├── 7.GenotypeGVCFs
├── 8.VQSR_input
├── 9.VQSR
├── 10.Apply_VQSR
├── 11.rsID_annotation
├── 12.Concat_vcf

```

There will be a folder for each sample processed in batch mode and included in the call set.
The "6.DbImport" folder will contain the variant DB generated to perform the joint calling.
The "7.GenotypeGVCFs" folder will contain the GVCF actual calls for all samples
The "8.VQSR_input" and "9.VQSR" folders will contain VCF format files used to perform the VQSR filtering by GATK and the recalibration tables and outputs.
The "10.Apply_VQSR" folder will contain the VCF files after the application of the VQSR filter.
The "11.rsID_annotation" folder will contain rsID annotated VCF files (splitted by chromosome)
The "12.Concat_vcf" folder will contain a single file with chromosomes from 1 to X. Chromosome Y is not included since it is called only for males samples. ChrY data can be still found in the previous folder.

---

## Running the pipeline

There are different ways to run the pipeline: **Local mode**, **Cluster mode** or **Single node mode**

### Local mode

In Local mode, the pipeline is execute in an interactive shell session (locally or on a cluster) and all the rules are treated as processes that can be run sequentially or in parallel, depending on the resources provided. One example of a Local execution is:

```bash
conda activate snakemake

base_cwd=/<USER_DEFINED_PATH>/VAR_CALL
log_folder=${base_cwd}/Log
mkdir -p ${log_folder}
cd ${base_cwd}

snakefile=/<USER_DEFINED_PATH>/VariantCalling-snakemake/Snakefile
configfile=/<USER_DEFINED_PATH>/VAR_CALL/variantCalling_pipeline.yaml
cores=24

snakemake -p -r -s ${snakefile} --configfile ${configfile} --use-envmodules --keep-going --cores ${cores}

```

In this example we assumed we had 24 CPU available for our calculation



### Cluster mode

In cluster mode, the pipeline runs on a interactive shell (**screen or tmux**) and each rule is submitted as a job on the cluster.
One example of a Cluster execution, on the **ORFEO cluster**, is:

```bash
module load conda
conda activate snakemake

base_cwd=/<USER_DEFINED_PATH>/VAR_CALL
log_folder=${base_cwd}/Log
mkdir -p ${log_folder}
cd ${base_cwd}

snakefile=/<USER_DEFINED_PATH>/VariantCalling-snakemake/Snakefile
configfile=/<USER_DEFINED_PATH>/VAR_CALL/variantCalling_pipeline.yaml

log_name=variantCalling_pipeline.log
stderr_name=variantCalling_pipeline.err
cores=24
queue=thin

snakemake -p -r -s ${snakefile} --configfile ${configfile} --use-envmodules --keep-going --cluster "qsub -q ${queue} -V -k eod -l select=1:ncpus={threads}:mem={resources.mem_mb}mb -l walltime=96:00:00" -j ${cores} 1> ${log_name} 2> ${stderr_name}
```

In this example we defined also the name for two additional log files, which will help to keep track of the pipeline execution. In this case, the **-j** option will define how many concurrent jobs are submitted on the cluster.



### Single node mode

In Single node mode, the pipeline runs as a job on the cluster and all rules are treated as processes that can be run sequentially or in parallel, depending on the resources provided. Similar to the Local execution mode.
One example of a single node mode execution, on the **ORFEO cluster**, is:

```bash
module load conda
conda activate snakemake

base_cwd=/<USER_DEFINED_PATH>/VAR_CALL
log_folder=${base_cwd}/Log
mkdir -p ${log_folder}
cd ${base_cwd}

snakefile=/<USER_DEFINED_PATH>/VariantCalling-snakemake/Snakefile
configfile=/<USER_DEFINED_PATH>/VAR_CALL/variantCalling_pipeline.yaml

log_name=variantCalling_pipeline.log
stderr_name=variantCalling_pipeline.err
cores=24
queue=thin
mem=750g

echo "cd ${base_cwd};module load conda;conda activate snakemake; snakemake -p -r -s ${snakefile} --configfile ${configfile} --cores ${cores} --use-envmodules --keep-going" | qsub -N snake_preprocessing -q ${queue} -V -k eod -o ${log_folder}/${log_name} -e ${log_folder}/${stderr_name} -l select=1:ncpus=${cores}:mem=${mem} -l walltime=96:00:00
done
```

In this example we selected an entire cluster node on the "thin" queue of the ORFEO cluster, defining the number of CPU (24) and the total amount of RAM required to run the pipeline (750g). We defined also the name for the two additional log files, to keep track of the pipeline execution.

