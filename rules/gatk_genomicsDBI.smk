
# Start DBImport: the best would be to generate 5mb chunks on each chromosome, but this could end up in a lot of jobs
rule gatk_genomics_db_import:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    input:
        # gvcfs=expand("variant_calling/{sample.sample}.{{interval}}.g.vcf.gz",
        #              sample=samples.reset_index().itertuples())
        gvcfs=lambda wildcards: expand(config["files_path"]["base_joint_call_path"] + "/{sample}/{sample}_{{interval_name}}_g.vcf.gz", sample=sample_names_from_interval(wildcards.interval_name)),
        import_interval=os.path.join(config.get('files_path').get('base_joint_call_path'),config.get('rules').get('split_intervals').get('out_dir')) + '/{scatteritem}_{interval_name}'
        # gvcfs=expand(config["files_path"]["base_joint_call_path"] + "/{sample}/{sample}_wgs_calling_regions_chr{chr}.GRCh38.p13.interval_list_g.vcf.gz")
    output:
        # directory(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genomics_db_import").get("out_dir"),"{scatteritem}_{interval_name}"))
        protected(directory(config.get("files_path").get("base_joint_call_path")+ "/"+ config.get("rules").get("gatk_genomics_db_import").get("out_dir") + "/{scatteritem}_{interval_name}"))
        # os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genomics_db_import").get("out_dir"),"{scatteritem}_{interval_name}/pippo.txt")
    params:
        gatk=config['GATK_TOOL'],
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        java_opt=config['java_opts']['opt3x'],
        fixed_args=config.get("rules").get("gatk_genomics_db_import").get("arguments"),
        genomics_db_option=config.get("rules").get("gatk_genomics_db_import").get("genomics_db_option"),
        tmp=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("files_path").get("tmp")),
        gvcf_args=lambda wildcards : " -V ".join(expand(config["files_path"]["base_joint_call_path"] + "/{sample}/{sample}_{interval_name}_g.vcf.gz", sample=sample_names_from_interval(wildcards.interval_name), interval_name=wildcards.interval_name))
    log:
        config["files_path"]["log_dir"] + "/{interval_name}-{scatteritem}-genomics_db_import.log",
        config["files_path"]["log_dir"] + "/{interval_name}-{scatteritem}-genomics_db_import.e"
    threads: 4
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt3x'])
    priority: 50
    benchmark:
        config["files_path"]["benchmark"] + "/{interval_name}_{scatteritem}_genomics_db_import.tsv"
    envmodules:
        "gatk/4.2.2.0"
    message: """ GenomicsDBImport """
    shell:
        """
        echo "Let\'s do stuff with {input.import_interval}...."
        mkdir -p {params.tmp}
        # {params.gatk} --java-options "{params.java_opt}" GenomicsDBImport --genomicsdb-workspace-path {output[0]} {params.fixed_args} -L {input.import_interval} -V {params.gvcf_args} --tmp-dir {params.tmp} > {log[0]} 2> {log[1]}
        {params.gatk} --java-options "{params.java_opt}" GenomicsDBImport {params.genomics_db_option} {output[0]} {params.fixed_args} -L {input.import_interval} -V {params.gvcf_args} --tmp-dir {params.tmp} > {log[0]} 2> {log[1]}
        """
