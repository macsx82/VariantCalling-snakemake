#rule for calling on autosomal chromosomes
rule gatk_hap_caller:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    output:
        BASE_OUT + "/" + config["rules"]["gatk_hap_caller"]["out_dir"] + "/{sample}/{sample}_{interval_name}_g.vcf.gz"
    input:
        cram=rules.apply_bqsr.output[0]
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
        config["files_path"]["benchmark"] + "/{sample}_{interval_name}.tsv"
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
                {params.gatk} --java-options "{params.java_opt}" HaplotypeCaller -R {params.ref_genome} -I {input.cram} -O {output} -L {params.intervals_file} -ip {params.ip1} --max-alternate-alleles {params.maa} -ERC GVCF --native-pair-hmm-threads 1 --output-mode EMIT_ALL_CONFIDENT_SITES -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation 2> {log[1]} 1> {log[0]}
            fi
        else
            #we are working with a male sample
            if [[ ${{chr_mode}} == "SEXUAL" ]];then
                echo "Job submitted for male sample on sexual chromosome."
                {params.gatk} --java-options "{params.java_opt}" HaplotypeCaller -R {params.ref_genome} -I {input.cram} -O {output} -L {params.intervals_file} -ip {params.ip1} -ploidy {params.sample_sex} --max-alternate-alleles {params.maa} -ERC GVCF --native-pair-hmm-threads 1 --output-mode EMIT_ALL_CONFIDENT_SITES -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation 2> {log[1]} 1> {log[0]}
            else 
                echo "Job submitted for male sample on non sexual chromosome."
                {params.gatk} --java-options "{params.java_opt}" HaplotypeCaller -R {params.ref_genome} -I {input.cram} -O {output} -L {params.intervals_file} -ip {params.ip1} --max-alternate-alleles {params.maa} -ERC GVCF --native-pair-hmm-threads 1 --output-mode EMIT_ALL_CONFIDENT_SITES -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation 2> {log[1]} 1> {log[0]}
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
