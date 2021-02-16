#!/usr/bin/env bash

#function file for the pipeline, to include in the header, after the CONFIG FILE
if [[ ${cluster_man}=="CINECA" ]]; then
    source ${HOME}/.startup_modules 
    # conda init bash
    # conda activate py37
    source activate py37
fi
#extract stats from the file
function sam_stats(){
    infile=$1

    ${SAMTOOLS} flagstat ${infile}
    echo
    ${SAMTOOLS} view -H ${infile} | grep -e '@RG' -e '@PG'
    echo "- END -"
    # ${SAMTOOLS} stats -r ${GNMhg38} ${infile} > ${infile}.stats
}

#validate the bam file
function sam_validate(){
    in_bam=$1
    java -XX:+UseSerialGC -jar ${PICARD} ValidateSamFile I=${in_bam} MODE=VERBOSE R=${GNMhg38} TMP_DIR=${tmp}/
    #get ready for new Picard syntax: 23/10/2019
    # java -XX:+UseSerialGC -jar ${PICARD} ValidateSamFile -I ${in_bam} -MODE VERBOSE -R ${GNMhg38} -TMP_DIR ${tmp}/
    echo "- END -"
}


#extract stats from the vcf file
function vcf_stats(){
    in_vcf=$1
    out_file=$2
    ${BCFTOOLS} stats -s - -d 0,1000,10 -F ${GNMhg38} ${in_vcf} > ${out_file}

    stats_plots=`dirname ${out_file}`

    /share/apps/bio/bin/plot-vcfstats -s -p ${stats_plots}/stats ${out_file}
    echo "- END -"
}

