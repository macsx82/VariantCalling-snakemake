#rule for calling on autosomal chromosomes
rule gatk_hap_caller:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    output:
        BASE_OUT + "/" + config["rules"]["gatk_hap_caller"]["out_dir"] + "/{sample}/{sample}_{interval_name}_g.vcf.gz"
    input:
        cram=rules.apply_bqsr.output[0]
        # crai=rules.apply_bqsr.output[1]
        # interval=config["callable_intervals"]+"/{interval_name}.interval_list"
    params:
        gatk=config['GATK_TOOL'],
        tmp=os.path.join(BASE_OUT,config.get("files_path").get("tmp")),
        java_opt=config['java_opts']['opt2x'],
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        ip1=config['ip1'],
        maa=config['maa'],
        sample_sex=get_sample_sex,
        current_chr=get_chr_from_interval_list,
        intervals_file = get_interval_file,
        chr_mode=get_sex_chr_from_interval
    log:
        config["files_path"]["log_dir"] + "/{sample}-{interval_name}.log",
        config["files_path"]["log_dir"] + "/{sample}-{interval_name}.e"
    threads: 2
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt2x'])
    benchmark:
        config["files_path"]["benchmark"] + "{sample}_{interval_name}.tsv"
    envmodules:
        "gatk/4.1.9.0"
    message: """ HaplotypeCaller """
    shell:
        """
        sample_sex={params.sample_sex}
        current_chr={params.current_chr}
        chr_mode={params.chr_mode}

        if [[ ${{sample_sex}} -eq 2 ]];then
        #we are working with female samples: no need to worry about ploidy
        #but we need to skip chrY calls
            if [[ ${{current_chr}} != "chrY" ]];then
                echo "Job submitted for female sample on non Y chromosome."
                {params.gatk} --java-options {params.java_opt} HaplotypeCaller -R {params.ref_genome} -I {input.cram} -O ${output} -L {params.intervals_file} -ip {params.ip1} --max-alternate-alleles {params.maa} -ERC GVCF --native-pair-hmm-threads 1 --output-mode EMIT_ALL_CONFIDENT_SITES -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation 2> {log[1]} 1> {log[0]}
            fi
        else
            #we are working with a male sample
            if [[ ${{chr_mode}} == "SEXUAL" ]];then
                echo "Job submitted for male sample on sexual chromosome."
                {params.gatk} --java-options {params.java_opt} HaplotypeCaller -R {params.ref_genome} -I {input.cram} -O ${output} -L {params.intervals_file} -ip {params.ip1} -ploidy {params.sample_sex} --max-alternate-alleles {params.maa} -ERC GVCF --native-pair-hmm-threads 1 --output-mode EMIT_ALL_CONFIDENT_SITES -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation 2> {log[1]} 1> {log[0]}
            else 
                echo "Job submitted for male sample on non sexual chromosome."
                {params.gatk} --java-options {params.java_opt} HaplotypeCaller -R {params.ref_genome} -I {input.cram} -O ${output} -L {params.intervals_file} -ip {params.ip1} --max-alternate-alleles {params.maa} -ERC GVCF --native-pair-hmm-threads 1 --output-mode EMIT_ALL_CONFIDENT_SITES -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation 2> {log[1]} 1> {log[0]}
            fi
        fi
        """
# run:
#         if params.sample_sex == "2" :
#             if params.current_chr != "chrY" :
#                 print("Job submitted for female sample on non Y chromosome.")
#                 shell("{params.gatk} --java-options {params.java_opt} HaplotypeCaller -R {params.ref_genome} -I {input.cram} -O ${output} -L {params.intervals_file} -ip {params.ip1} --max-alternate-alleles {params.maa} -ERC GVCF --native-pair-hmm-threads 1 --output-mode EMIT_ALL_CONFIDENT_SITES -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation 2> {log[1]}")
#         else :
#             #we are working with a male sample
#             if params.chr_mode == "SEXUAL":
#                 print("Job submitted for male sample on sexual chromosome.")
#                 shell("{params.gatk} --java-options {params.java_opt} HaplotypeCaller -R {params.ref_genome} -I {input.cram} -O ${output} -L {params.intervals_file} -ip {params.ip1} -ploidy {params.sample_sex} --max-alternate-alleles {params.maa} -ERC GVCF --native-pair-hmm-threads 1 --output-mode EMIT_ALL_CONFIDENT_SITES -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation 2> {log[1]}")
#             elif params.chr_mode == "AUTOSOMAL":
#                 print("Job submitted for male sample on non sexual chromosome.")
#                 shell("{params.gatk} --java-options {params.java_opt} HaplotypeCaller -R {params.ref_genome} -I {input.cram} -O ${output} -L {params.intervals_file} -ip {params.ip1} --max-alternate-alleles {params.maa} -ERC GVCF --native-pair-hmm-threads 1 --output-mode EMIT_ALL_CONFIDENT_SITES -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation 2> {log[1]}")
#rule for calling on sex chromosomes
# rule gatk_hap_caller_sex:
#     output:
#         BASE_OUT + "/" + config["rules"]["bqsr_proc"]["out_dir"] + "/{sample}_recaldata.csv"
#     input:
#         rules.sort_bam.output[0]
#     params:
#         gatk=config['GATK_TOOL'],
#         tmp=config['files_path']['tmp'],
#         java_opt=config['java_opts']['opt2x'],
#         ref_genome=config['genome_fasta'],
#         DBSNP_latest=config['known_variants']['dbsnp_latest'],
#         INDELS=config['known_variants']['indels'],
#         mills=config['known_variants']['mills']
#     log:
#         config["files_path"]["log_dir"] + "/{sample}-bqsr.log",
#         config["files_path"]["log_dir"] + "/{sample}-bqsr.e"
#     threads: 1
#     resources:
#         mem_mb=get_resources_from_jvm(config['java_opts']['opt2x'])
#     benchmark:
#         config["files_path"]["benchmark"] + "{sample}_bqsr.tsv"
#     message: """ BaseRecalibrator """
#     shell:
#         """
#         {params.gatk} --java-options {params.java_opt} BaseRecalibrator -R {params.ref_genome} -I {input} --use-original-qualities --showHidden --tmp-dir {params.tmp}/ -O {output} --known-sites {params.DBSNP_latest} --known-sites {params.INDELS} --known-sites {params.mills}
#         """

# echo "> HaplotypeCaller"
#     # ${GATK4} --java-options ${java_opt2x} ${java_XX1} ${java_XX2} HaplotypeCaller -R ${GNMhg38} -I ${fol4}/${applybqsr} -O ${fol5}/${c_gv} -L "${f1}" -ip ${ip1} --max-alternate-alleles ${maa} -ERC GVCF
#     #We need to manage also the ploidy based on the sex and the interval
#     if [[ ${sex} -eq 2 ]]; then
#         #we are working with female samples: no need to worry about ploidy
#         #but we need to skip chrY calls
#         if [[ -f "${f1}" ]]; then
#             #if the interval is a file, we need to check if its a chrY file 
#             #and skip the submission
#             chrY_content=$(awk '$1~"Y" {print $1}' ${f1} | sort|uniq| wc -l| cut -f 1 -d " ")
#         else
#             #If the interval is submitted in the form of
#             # f1="chrY:2781480-9046914"
#             chrY_content=$(echo ${f1} | awk '$0~"Y"'| wc -l | cut -f 1 -d " ")
#         fi

#         if [[ ${chrY_content} -eq 0 ]]; then
#             #we will submit the job
#             echo "Job submitted for female sample on non Y chromosome."
#             ${GATK4} --java-options "${java_opt2x} ${java_XX1} ${java_XX2}" HaplotypeCaller -R ${GNMhg38} -I ${fol4}/${applybqsr} -O ${fol5}/${c_gv} -L "${f1}" -ip ${ip1} --max-alternate-alleles ${maa} -ERC GVCF --native-pair-hmm-threads 1 --output-mode EMIT_ALL_CONFIDENT_SITES
#         fi
#     else
#         #we are working with male samples: we need to worry about ploidy
#         #to correctly handle non par regions and chrY calls
#         if [[ -f "${f1}" ]]; then
#             #if the interval is a file, we need to check if its a chrY file 
#             #and skip the submission
#             chrY_content=$(awk '$1~"Y" {print $1}' ${f1} | sort|uniq| wc -l| cut -f 1 -d " ")
#             chrX_content=$(awk '$0~"NON_PAR"' ${f1}| sort|uniq| wc -l| cut -f 1 -d " ")
#         else
#             #If the interval is submitted in the form of
#             # f1="chrY:2781480-9046914"
#             chrY_content=$(echo ${f1} | awk '$0~"Y"'| wc -l | cut -f 1 -d " ")
#             #for chrX interval, we need to check if it is in the par or non par regions
#             # f1="chrX:10001-2781479"
#             #X: 10001 - 2781479 (PAR1)
#             #X: 155701383 - 156030895 (PAR2)
#             #we are checking if we are outside the par regions
#             chrX_content=$( echo ${f1} | awk '$0~"X"' | awk '{split($1,b,":");split(b[2],c,"-");if( (c[1] <= 10001 || (c[1] >= 2781479 && c[1] <= 155701383) || c[1] >= 156030895) && (c[2] <= 10001 ||( c[2] >= 2781479 && c[2] <= 155701383)|| c[2] >= 156030895) ) print $0}' | wc -l | cut -f 1 -d " ")
#         fi
#         if [[ ${chrY_content} != "0" ]]; then
#             chrom_status="SEXC"
#         fi
#         if [[ ${chrX_content} != "0" ]]; then
#             chrom_status="SEXC"
#         fi
#         if [[ ${chrX_content} -eq 0 && ${chrY_content} -eq 0 ]]; then
#             chrom_status="AUTOSOMAL"
#         fi
#         # echo ${chrom_status}
#         case ${chrom_status} in
#             AUTOSOMAL )
#                 echo "Job submitted for male sample on autosomal chromosome or chrX PAR regions."
#                 ${GATK4} --java-options "${java_opt2x} ${java_XX1} ${java_XX2}" HaplotypeCaller -R ${GNMhg38} -I ${fol4}/${applybqsr} -O ${fol5}/${c_gv} -L "${f1}" -ip ${ip1} --max-alternate-alleles ${maa} -ERC GVCF --native-pair-hmm-threads 1 --output-mode EMIT_ALL_CONFIDENT_SITES
#                 ;;
#             SEXC )
#                 echo "Job submitted for male sample on Y chromosome or chrX NON_PAR regions."
#                 ${GATK4} --java-options "${java_opt2x} ${java_XX1} ${java_XX2}" HaplotypeCaller -R ${GNMhg38} -I ${fol4}/${applybqsr} -O ${fol5}/${c_gv} -L "${f1}" -ip ${ip1} -ploidy ${sex} --max-alternate-alleles ${maa} -ERC GVCF --native-pair-hmm-threads 1 --output-mode EMIT_ALL_CONFIDENT_SITES
#                 ;;
#         esac
#         # ${GATK4} --java-options "${java_opt2x} ${java_XX1} ${java_XX2}" HaplotypeCaller -R ${GNMhg38} -I ${fol4}/${applybqsr} -O ${fol5}/${c_gv} -L "${f1}" -ip ${ip1} --max-alternate-alleles ${maa} -ERC GVCF --output-mode EMIT_ALL_CONFIDENT_SITES
#     fi
#     