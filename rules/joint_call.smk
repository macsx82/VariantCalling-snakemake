#rule for joint call
# Start DBImport: the best would be to generate 5mb chunks on each chromosome, but this could end up in a lot of jobs
rule gatk_genomics_db_import:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    input:
        # gvcfs=expand("variant_calling/{sample.sample}.{{interval}}.g.vcf.gz",
        #              sample=samples.reset_index().itertuples())
        gvcfs=expand(config["files_path"]["base_joint_call_path"] + "/{sample}/{sample}_{{interval_name}}_g.vcf.gz", sample=sample_names )
        # gvcfs=expand(config["files_path"]["base_joint_call_path"] + "/{sample}/{sample}_wgs_calling_regions_chr{chr}.GRCh38.p13.interval_list_g.vcf.gz")
    output:
        touch(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genomics_db_import").get("out_dir"),"{interval_name}"))
    params:
        # custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=2),
        # genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        # gvcfs=_multi_flag_dbi("-V", expand("variant_calling/{sample.sample}.{{interval}}.g.vcf.gz",sample=samples.reset_index().itertuples()))
    log:
        config["files_path"]["log_dir"] + "/{interval}-genomics_db_import.log",
        config["files_path"]["log_dir"] + "/{interval}-genomics_db_import.e"
    threads: 2
    # resources:
    #     mem_mb=get_resources_from_jvm(config['java_opts']['opt2x'])
    # benchmark:
    #     config["files_path"]["benchmark"] + "/{interval}_genomics_db_import.tsv"
    envmodules:
        "gatk/4.2.2.0"
        # "gatk/4.1.9.0"
    message: """ GenomicsDBImport """
    shell:
        """
        echo {input.gvcfs}
        """



# rule gatk_genotype_gvcfs:
#     input:
#         "db/imports/{interval}"
#     output:
#         protected("variant_calling/all.{interval}.vcf.gz")
#     envmodules:
#         "gatk/4.2.2.0"
#         # "gatk/4.1.9.0"
#     params:
#         custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=2),
#         genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
#     log:
#         config["files_path"]["log_dir"] + "/{interval}-genotype_gvcfs.log",
#         config["files_path"]["log_dir"] + "/{interval}-genotype_gvcfs.e"
#     threads: 2
#     resources:
#         mem_mb=get_resources_from_jvm(config['java_opts']['opt2x'])
#     benchmark:
#         config["files_path"]["benchmark"] + "/{interval}_genotype_gvcfs.tsv"
#     shell:
#         "gatk GenotypeGVCFs --java-options {params.custom} "
#         "-R {params.genome} "
#         "-V gendb://db/{wildcards.interval} "
#         "-G StandardAnnotation "
#         "-O {output} "
#         ">& {log} "
