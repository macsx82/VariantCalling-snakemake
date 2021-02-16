#Snakemake pipeline to perform variant annotation using VEP
# 
#
# 18/11/2019
#
# Author: massimiliano [dot] Cocca [at] burlo [dot] trieste [dot] it
# configfile: "annotation.yaml"

# vcf_files = list(pd.read_table(config["vcf_files"]).vcf_file)
# vcf_names = [ntpath.basename(vcf_file) for vcf_file in vcf_files]
# Definition of helper functions
# def get_vcfs():



rule var_annotation:
    input:
        lambda wildcards: vcf_names[wildcards.vcf_name]
        # "{vcf_files}"
    output:
        "annotated/{vcf_name}"
    # conda:
    message: "Executing one using {input}"

    shell:
        "bcftools view {input} | "
        "bgzip -c > {output}"
    # infile=${inputfile}
    # outfile=$3
    # chr=$4
    # annot_opt=$5
    
    # cluster=$6
    # #add 3 different modes:
    # # NEW: brand new annotation of the input file
    # # STRIP: stripo old CSQ annotation and reannotate
    # # REGION: select a region for annotation
    # #define vep cache basefolder
    # # maincache=/tmp/VEP_cache
    # # main_sw=z/home/cocca/softwares
    # # ref_fasta=${maincache}/90/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
    # vep_version=$7
    # assembly_build=$8
    # ref_fasta=$9
    
    # if [[ ${cluster} == "CBM" ]]; then
    #     vep_bin=/home/cocca/softwares/VEP/ensembl-vep/vep
    #     maincache=/home/cocca/resources/VEP_cache
    # else
    #     maincache=/shared/resources/VEP_cache
    #     case ${vep_version} in
    #         90)
    #             vep_bin=/shared/softwares/VEP/ensembl-vep/vep
    #             loftee_path=/shared/softwares/VEP/loftee
    #         ;;
    #         94)
    #             vep_bin=/share/apps/bio/miniconda2/bin/vep
    #             loftee_path=/share/apps/bio/miniconda2/share/ensembl-vep-94.5-0
    #         ;;
    #     esac
    # fi



    # echo -e "infile: ${infile}
    # outfile: ${outfile}
    # chr: ${chr}
    # annot_opt: ${annot_opt}
    # cluster: ${cluster}
    # vep_version: ${vep_version}
    # assembly_build: ${assembly_build}
    # ref_fasta: ${ref_fasta}
    # vep binary: ${vep_bin}"

    # case $annot_opt in
    #     NEW)
    #         # ${main_sw}/ensembl-vep/vep --i ${infile} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${maincache}/90/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/90/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,loftee_path:${main_sw}/loftee --symbol --vcf --cache --merged --dir ${maincache}/90 --assembly GRCh37
    #         # ${main_sw}/ensembl-vep/vep --i ${infile} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${maincache}/90/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/90/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,loftee_path:${main_sw}/loftee --plugin CADD,${maincache}/90/cadd13/1000G_phase3.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/90 --assembly GRCh37
    #         ${vep_bin} --i ${infile} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd13/1000G_phase3.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly ${assembly_build}
    #         # ${main_sw}/ensembl-vep/vep --i ${infile} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf.sql,loftee_path:${main_sw}/loftee --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly GRCh37
    #     ;;
    #     STRIP)
    #         # bcftools annotate -x INFO/CSQ ${infile} | ${main_sw}/ensembl-vep/vep --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,loftee_path:${main_sw}/loftee --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly GRCh37
    #         # bcftools annotate -x INFO/CSQ ${infile} | ${main_sw}/ensembl-vep/vep --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,loftee_path:${main_sw}/loftee --plugin CADD,${maincache}/${vep_version}/cadd13/1000G_phase3.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly GRCh37
    #         bcftools annotate -x INFO/CSQ ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd13/1000G_phase3.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly ${assembly_build}
    #         # bcftools annotate -x INFO/CSQ ${infile} | ${main_sw}/ensembl-vep/vep --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${maincache}/${vep_version}/homo_sapiens/${vep_version}_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf.sql,loftee_path:${main_sw}/loftee --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly GRCh37
    #     ;;
    #     REGION)
    #         start=$7
    #         end=$8
    #         # bcftools view -r ${chr}:${start}-${end} ${infile} | ${main_sw}/ensembl-vep/vep --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,loftee_path:${main_sw}/loftee --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly GRCh37
    #         # bcftools view -r ${chr}:${start}-${end} ${infile} | ${main_sw}/ensembl-vep/vep --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,loftee_path:${main_sw}/loftee --plugin CADD,${maincache}/${vep_version}/cadd13/1000G_phase3.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly GRCh37
    #         bcftools view -r ${chr}:${start}-${end} ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd13/1000G_phase3.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly ${assembly_build}
    #         # bcftools view -r ${chr}:${start}-${end} ${infile} | ${main_sw}/ensembl-vep/vep --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${maincache}/90/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/90/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/90/loftee_data/phylocsf.sql,loftee_path:${main_sw}/loftee --symbol --vcf --cache --merged --dir ${maincache}/90 --assembly GRCh37
    #     ;;
    # esac

    # # annotate vcf files with consequences
    # bgzip ${outfile}
    # tabix -p vcf ${outfile}.gz