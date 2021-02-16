rule bqsr_proc:
    output:
        BASE_OUT + "/" + config["rules"]["bqsr_proc"]["out_dir"] + "/{sample}_recaldata.csv"
    input:
        rules.sort_bam.output[0]
    params:
        gatk=config['GATK_TOOL'],
        tmp=config['files_path']['tmp'],
        java_opt=config['java_opts']['opt2x'],
        ref_genome=config['genome_fasta'],
        DBSNP_latest=config['known_variants']['dbsnp_latest'],
        INDELS=config['known_variants']['indels'],
        mills=config['known_variants']['mills']
    log:
        config["files_path"]["log_dir"] + "/{sample}-bqsr.log",
        config["files_path"]["log_dir"] + "/{sample}-bqsr.e"
    threads: 1
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt2x'])
    benchmark:
        config["files_path"]["benchmark"] + "{sample}_bqsr.tsv"
    message: """ BaseRecalibrator """
    shell:
        """
        {params.gatk} --java-options {params.java_opt} BaseRecalibrator -R {params.ref_genome} -I {input} --use-original-qualities --showHidden --tmp-dir {params.tmp}/ -O {output} --known-sites {params.DBSNP_latest} --known-sites {params.INDELS} --known-sites {params.mills}
        """

rule apply_bqsr:
    output:
        BASE_OUT +"/" + config["rules"]["apply_bqsr"]["out_dir"] + "/{sample}_bqsr.cram"
    input:
        rules.sort_bam.output[0],rules.bqsr_proc.output
    params:
        gatk=config['GATK_TOOL'],
        java_opt=config['java_opts']['opt2x'],
        tmp=config['files_path']['tmp'],
        ref_genome=config['genome_fasta'],
    log:
        config["files_path"]["log_dir"] + "/{sample}-bqsr.log",
        config["files_path"]["log_dir"] + "/{sample}-bqsr.e"
    benchmark:
        config["files_path"]["benchmark"] + "{sample}_bqsr.tsv"
    threads: 1
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt2x'])
    message: """ ApplyBQSR """
    shell:
        """
        {params.gatk} --java-options {params.java_opt} ApplyBQSR -R {params.ref_genome} -I {input[0]} -O {output} -bqsr {input[1]} --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --add-output-sam-program-record --create-output-bam-md5 --use-original-qualities
        """

