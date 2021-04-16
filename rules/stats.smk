rule sam_validate:
    output:
        BASE_OUT +"/" + config["rules"]["stats"]["out_dir"] + "/{sample}_validate.txt"
    input:
        rules.apply_bqsr.output
    params:
        samtools=config['SAMTOOLS']
    log:
        config["files_path"]["log_dir"] + "/{sample}-validate.log",
        config["files_path"]["log_dir"] + "/{sample}-validate.e"
    threads: 1
    resources:
        mem_mb=config['rules']['stats']['mem']
    benchmark:
        config["files_path"]["benchmark"] + "/{sample}_validate.tsv"
    envmodules:
        "samtools/1.11"
    message: """ Samtools quickcheck validation for bqsr'ed files """
    shell:
        """
        {params.samtools} quickcheck -vvv {input} > {output} 2> {log[1]}
        """

rule sam_stats:
    output:
        BASE_OUT +"/" + config["rules"]["stats"]["out_dir"] + "/{sample}_stats.txt"
    input:
        rules.apply_bqsr.output
    params:
        samtools=config['SAMTOOLS'],
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    log:
        config["files_path"]["log_dir"] + "/{sample}-stats.log",
        config["files_path"]["log_dir"] + "/{sample}-stats.e"
    threads: 1
    resources:
        mem_mb=config['rules']['stats']['mem']
    benchmark:
        config["files_path"]["benchmark"] + "/{sample}_stats.tsv"
    envmodules:
        "samtools/1.11"
    message: """ Samtools stats for bqsr'ed files """
    shell:
        """
        {params.samtools} stats -r {params.ref_genome} {input} > {output} 2> {log[1]}
        """


rule sam_flagstats:
    output:
        BASE_OUT +"/" +config["rules"]["stats"]["out_dir"] + "/{sample}_flagstat.txt"
    input:
        rules.apply_bqsr.output
    params:
        samtools=config['SAMTOOLS']
    log:
        config["files_path"]["log_dir"] + "/{sample}-flagstat.log",
        config["files_path"]["log_dir"] + "/{sample}-flagstat.e"
    threads: 1
    resources:
        mem_mb=config['rules']['stats']['mem']
    benchmark:
        config["files_path"]["benchmark"] + "/{sample}_flagstat.tsv"
    envmodules:
        "samtools/1.11"
    message: """ Samtools flagstat validation for bqsr'ed files """
    shell:
        """
        {params.samtools} flagstat {input} > {output} 2> {log[1]}
        """

rule picard_InsertSizeMetrics:
    output:
        BASE_OUT +"/" +config["rules"]["stats"]["out_dir"] + "/{sample}_ismetrics.txt",
        BASE_OUT +"/" +config["rules"]["stats"]["out_dir"] + "/{sample}_ismetrics.pdf"
    input:
        rules.apply_bqsr.output
    params:
        custom_jvm=config['java_opts']['opt2x'] + " -Djava.io.tmpdir=" + os.path.join(BASE_OUT,config.get("files_path").get("tmp")),
        picard=config.get("PICARD_TOOL"),
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    log:
        config["files_path"]["log_dir"] + "/{sample}-ismetrics.log",
        config["files_path"]["log_dir"] + "/{sample}-ismetrics.e"
    benchmark:
        config["files_path"]["benchmark"] + "/{sample}_ismetrics.tsv"
    envmodules:
        "picard/2.24.0",
        "R/4.0.3"
    threads: 2
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt2x'])
    message: """ extract insert size metrics """
    shell:
        """
        java {params.custom_jvm} -jar ${params.picard} CollectInsertSizeMetrics -R {params.ref_genome} -I {input} -O {output[0]} -H {output[1]} 2> {log[1]} 1> {log[0]}
        """

rule picard_WGSMetrics:
    output:
        BASE_OUT +"/" +config["rules"]["stats"]["out_dir"] + "/{sample}_wgsmetrics.txt"
    input:
        rules.apply_bqsr.output
    params:
        arguments=config.get("rules").get("picard_WGSMetrics").get("arguments"),
        custom_jvm=config['java_opts']['opt2x'] + " -Djava.io.tmpdir=" + os.path.join(BASE_OUT,config.get("files_path").get("tmp")),
        picard=config.get("PICARD_TOOL"),
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    log:
        config["files_path"]["log_dir"] + "/{sample}-wgsmetrics.log",
        config["files_path"]["log_dir"] + "/{sample}-wgsmetrics.e"
    benchmark:
        config["files_path"]["benchmark"] + "/{sample}_wgsmetrics.tsv"
    envmodules:
        "picard/2.24.0",
        "R/4.0.3"
    threads: 2
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt2x'])
    message: """ extract WGSMetrics """
    shell:
        """
        java {params.custom_jvm} -jar ${params.picard} CollectWgsMetrics {params.arguments} I={input} O={output} R={params.ref_genome} 2> {log[1]} 1> {log[0]}
        """
