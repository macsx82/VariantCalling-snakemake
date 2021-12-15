# Start DBImport: we need to define the rule in a differnet way, so to avoid the removal of already existing data
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
        # touch(directory(config.get("files_path").get("base_joint_call_path")+ "/"+ config.get("rules").get("gatk_genomics_db_import").get("out_dir") + "/{scatteritem}_{interval_name}"))
        config.get("files_path").get("base_joint_call_path")+ "/"+ config.get("rules").get("gatk_genomics_db_import").get("out_dir") + "/{scatteritem}_{interval_name}_update.txt"
    params:
        gatk=config['GATK_TOOL'],
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        java_opt=config['java_opts']['opt3x'],
        fixed_args=config.get("rules").get("gatk_genomics_db_import").get("arguments"),
        genomics_db_option=config.get("rules").get("gatk_genomics_db_import").get("genomics_db_option"),
        genomics_db_location=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genomics_db_import").get("out_dir"),"{scatteritem}_{interval_name}")
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
        {params.gatk} --java-options "{params.java_opt}" GenomicsDBImport {params.genomics_db_option} {params.genomics_db_location} {params.fixed_args} -L {input.import_interval} -V {params.gvcf_args} --tmp-dir {params.tmp} > {log[0]} 2> {log[1]}
        echo "{params.genomics_db_location}" > {output[0]}
        """

#perform joint genotyping after DB Import UPDATE
rule gatk_genotype_gvcfs:
    input:
        import_db=rules.gatk_genomics_db_import.output,
        import_interval=os.path.join(config.get('files_path').get('base_joint_call_path'),config.get('rules').get('split_intervals').get('out_dir')) + '/{scatteritem}_{interval_name}'
        # "db/imports/{interval}"
    output:
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genotype_gvcfs").get("out_dir"),"{scatteritem}_{interval_name}/all.{scatteritem}_{interval_name}.vcf.gz"),
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genotype_gvcfs").get("out_dir"),"{scatteritem}_{interval_name}/all.{scatteritem}_{interval_name}.vcf.gz.tbi")
        # protected("variant_calling/all.{scatteritem}_{interval_name}.vcf.gz")
    params:
        gatk=config['GATK_TOOL'],
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        java_opt=config['java_opts']['opt3x'],
        fixed_args=config.get("rules").get("gatk_genotype_gvcfs").get("arguments"),
        updated_import_db=get_db_path(input.import_db),
        # tmp=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("files_path").get("tmp"))
        tmp=config.get("rules").get("gatk_genotype_gvcfs").get("temp_folder")
    log:
        config["files_path"]["log_dir"] + "/{interval_name}-{scatteritem}-genotype_gvcfs.log",
        config["files_path"]["log_dir"] + "/{interval_name}-{scatteritem}-genotype_gvcfs.e"
    threads: 4
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt3x'])
    priority: 49
    benchmark:
        config["files_path"]["benchmark"] + "/{interval_name}-{scatteritem}_genotype_gvcfs.tsv"
    envmodules:
        "gatk/4.1.9.0"
    message: """ GenotypeGVCFs """
    shell:
        """
        mkdir -p {params.tmp}

        {params.gatk} --java-options "{params.java_opt}" GenotypeGVCFs -R {params.ref_genome} -L {input.import_interval} -V gendb://{params.updated_import_db} -O {output[0]} {params.fixed_args} --tmp-dir {params.tmp} > {log[0]} 2> {log[1]}
        """
