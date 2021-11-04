#generate stats for vcf files before VQSR
rule vcf_stats:
    output:
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genotype_gvcfs").get("out_dir"),"all.{interval_name}.stats")
    input:
        rules.chrom_intervals_gather.output[0]
    params:
        bcftools=config['BCFTOOLS']
    log:
        config["files_path"]["log_dir"] + "/{interval_name}-vcf_stats.log",
        config["files_path"]["log_dir"] + "/{interval_name}-vcf_stats.e"
    benchmark:
        config["files_path"]["benchmark"] + "/{interval_name}_vcf_stats.tsv"
    resources:
        mem_mb=5000
    envmodules:
        "bcftools/1.11"
    message: """ VCF stats before VQSR """
    shell:
        """
        {params.bcftools} stats -v {input} > {output} 2> {log[1]}
        """