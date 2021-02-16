#27/01/2021

##Snakemake implementation of the variant calling pipeline 

We will setup the variant calling pipeline, from fastq alignment to vcf file generation, using snakemake, to make it portable and deployable on different cluster.

---
#Genome reference

We used GRCh38.p13, as default reference for alignment and variant calling.
Since we want to use the same interval list approach used with the standard GATK bundle, we will generate the
callable regions interval list.

We'll the following pipeline:

1- Use picard ScatterIntervalsByNs to extract the non-poly N regions.
2- Convert this interval list to bed with awk.
3- Run bedtools makewindows to generate your regions.
	(from https://www.biostars.org/p/473462/#473466 )


```bash
base_out=/storage/burlo/cocca/resources/hgRef/GRCh38.p13
ref_seq=/storage/burlo/cocca/resources/hgRef/GRCh38.p13/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna
log_folder=/large/___HOME___/burlo/cocca/analyses/Logs


# echo "module load picard/2.24.0;java -jar $picard CreateSequenceDictionary R=${ref_seq} O=/storage/burlo/cocca/resources/hgRef/GRCh38.p13/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.dict" | qsub -q thin -V -k eod -o ${log_folder}/picard_index_ref.o -e ${log_folder}/picard_index_ref.e -l nodes=1:ppn=1 -l walltime=48:00:00 -l mem=10g
echo "module load picard/2.24.0;java -jar $picard ScatterIntervalsByNs R=${ref_seq} OT=ACGT O=${base_out}/wgs_calling_regions.GRCh38.p13.interval_list" | qsub -q thin -V -k eod -o ${log_folder}/picard_split_interval.o -e ${log_folder}/picard_split_interval.e -l nodes=1:ppn=1 -l walltime=48:00:00 -l mem=10g

```

Generate intervals for PAR regions on X chromosome

```bash
(fgrep "@" wgs_calling_regions_chrX.GRCh38.p13.interval_list;bedtools subtract -a <(fgrep -v "@" wgs_calling_regions_chrX.GRCh38.p13.interval_list) -b par_regions.bed )> wgs_calling_regions_chrX_NONPAR.GRCh38.p13.interval_list

(fgrep "@" wgs_calling_regions_chrX.GRCh38.p13.interval_list;bedtools intersect -a <(fgrep -v "@" wgs_calling_regions_chrX.GRCh38.p13.interval_list) -b par_regions.bed -wb|cut -f -5,9) > wgs_calling_regions_chrX_PAR.GRCh38.p13.interval_list

```

We need to manually modify the non par regions start and end, since the bedtools operation leave the data with a 0-based position
