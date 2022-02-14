# Variant Calling pipeline

A Snakemake implementation of the variant calling pipeline

This pipeline is designed to perform Variant calling on WGS data, starting from preprocessed fastq files. A pipeline to perform the necessary preprocessing steps is available [here](https://gitlab.burlo.trieste.it/max/SeqPreproc-snakemake/-/tree/main)

The pipeline can be executed in different modes:

- BATCH
- JOINT_CALL
- JOINT_CALL_UPDATE


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

### Other requirements

Absolute paths of the fastq files (saparate files for R1 and R2 strand) to be processed, produced by the [SeqPreprocessing](https://gitlab.burlo.trieste.it/max/SeqPreproc-snakemake/-/tree/main) pipeline.

Manifest file following the template provided in the **resource** folder, space separated and with the following header line:

```
SAMPLE_ID fq1 fq2 sex
```


At this point, we can proceed with the system environment setup, the same we did for the SeqPreprocessing pipeline.


---

## Usage

It is advised to clone the pipeline to a personal folder from the git repository, in order to be able to correctly execute the workflow.
The clone command to use is:

```bash
	git clone https://gitlab.burlo.trieste.it/max/VariantCalling-snakemake.git
```

This command will create a new folder, in the current location, named **VariantCalling-snakemake**


Before running the pipeline, it is necessary to provide some information using a config file. A config file template, with default values is provided and it is named **config.yaml**

It is advisable not to execute the pipeline in the repository folder, but to create a separate folder and a copy a copy of the config file to fill the relevant parameters.




### Running the pipeline

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

