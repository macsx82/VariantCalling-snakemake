#rules for joint call

#split the callable intervals in different chunks, so we can better parallelize
rule split_intervals:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    output:
        scatter.split(os.path.join(config.get('files_path').get('base_joint_call_path'),config.get('rules').get('split_intervals').get('out_dir')) + '/{scatteritem}_{{interval_name}}')
        # scatter.split(config.get("files_path").get("base_joint_call_path") + "/" + config.get("rules").get("split_intervals").get("out_dir") + "/{interval_name}_{scatteritem}")
        # scatter.split(config["files_path"]["base_joint_call_path"] + "/" + config["rules"]["split_intervals"]["out_dir"] + "/{interval_name}_{scatteritem}")
        # scatter.split("/large/___SCRATCH___/burlo/cocca/JOINT_CALL_TEST_20210930/JOINT_CALL/6.SplittedIntervals/{interval_name}")
        # scatter.split(config["files_path"]["base_joint_call_path"])
    input:
        get_interval_file
    params:
        gatk=config['GATK_TOOL'],
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        java_opt=config['java_opts']['opt1x'],
        scatter_count=config.get("rules").get("split_intervals").get("scatter_count"),
        output_folder=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("split_intervals").get("out_dir"))
        # intervals_file = get_interval_file,
        # custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=2),
        # genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        # gvcfs=_multi_flag_dbi("-V", expand("variant_calling/{sample.sample}.{{interval}}.g.vcf.gz",sample=samples.reset_index().itertuples()))
    log:
        config["files_path"]["log_dir"] + "/{interval_name}-split_intervals.log",
        config["files_path"]["log_dir"] + "/{interval_name}-split_intervals.e"
    threads: 2
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt2x'])
    benchmark:
        config["files_path"]["benchmark"] + "/{interval_name}_split_intervals.tsv"
    # envmodules:
    #     "gatk/4.2.2.0"
    message: """ Split intervals """
    run:
        gatk_split_intervals_cmd="module load conda;module load gatk/4.2.2.0;%s SplitIntervals --java-options \"%s\" -R %s -L %s --scatter-count %s -O %s --extension \"_%s\" > %s 2> %s" %(params.gatk,params.java_opt,params.ref_genome[0],input,params.scatter_count,params.output_folder,wildcards.interval_name,log[0], log[1])
        shell(gatk_split_intervals_cmd)
        for scattered_item in output:
            print(scattered_item)
            #need to get the output of the gatk command using the numbering of the scattered items!
            scattered_index=os.path.basename(scattered_item).split("_")[0].split("-")[0]
            # now rename the files generated by GATK
            gatk_out="%s/%s_%s" %(params.output_folder,'{num:04d}'.format(num=int(scattered_index)-1),wildcards.interval_name)
            print(gatk_out)
            os.rename(gatk_out,scattered_item)

    # shell:
    #     """
    #     {params.gatk} SplitIntervals --java-options "{params.java_opt}" -R {params.ref_genome} -L {input} --scatter-count {params.scatter_count} -O {params.output_folder} --extension "_{wildcards.interval_name}" > {log[0]} 2> {log[1]}
    #     """
    

# Start DBImport: the best would be to generate 5mb chunks on each chromosome, but this could end up in a lot of jobs
rule gatk_genomics_db_import:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    input:
        # gvcfs=expand("variant_calling/{sample.sample}.{{interval}}.g.vcf.gz",
        #              sample=samples.reset_index().itertuples())
        gvcfs=expand(config["files_path"]["base_joint_call_path"] + "/{sample}/{sample}_{{interval_name}}_g.vcf.gz", sample=sample_names_from_interval),
        import_interval=os.path.join(config.get('files_path').get('base_joint_call_path'),config.get('rules').get('split_intervals').get('out_dir')) + '/{scatteritem}_{interval_name}'
        # gvcfs=expand(config["files_path"]["base_joint_call_path"] + "/{sample}/{sample}_wgs_calling_regions_chr{chr}.GRCh38.p13.interval_list_g.vcf.gz")
    output:
        directory(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genomics_db_import").get("out_dir"),"{scatteritem}_{interval_name}"))
        # os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genomics_db_import").get("out_dir"),"{scatteritem}_{interval_name}/pippo.txt")
    params:
        gatk=config['GATK_TOOL'],
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        java_opt=config['java_opts']['opt3x'],
        fixed_args=config.get("rules").get("gatk_genomics_db_import").get("arguments"),
        tmp=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("files_path").get("tmp")),
        gvcf_args=" -V ".join(expand(config["files_path"]["base_joint_call_path"] + "/{sample}/{sample}_{{interval_name}}_g.vcf.gz", sample=sample_names_from_interval)
    log:
        config["files_path"]["log_dir"] + "/{interval_name}-{scatteritem}-genomics_db_import.log",
        config["files_path"]["log_dir"] + "/{interval_name}-{scatteritem}-genomics_db_import.e"
    threads: 4
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt3x'])
    benchmark:
        config["files_path"]["benchmark"] + "/{interval_name}_{scatteritem}_genomics_db_import.tsv"
    envmodules:
        "gatk/4.2.2.0"
    message: """ GenomicsDBImport """
    shell:
        """
        echo "Let\'s do stuff with {input.import_interval}...."
        {params.gatk} --java-options "{params.java_opt}" GenomicsDBImport --genomicsdb-workspace-path {output[0]} {params.fixed_args} -L {input.import_interval} -V {params.gvcf_args} --tmp-dir {params.tmp} > {log[0]} 2> {log[1]}
        """

#perform joint genotyping
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
        tmp=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("files_path").get("tmp"))
    log:
        config["files_path"]["log_dir"] + "/{interval_name}-{scatteritem}-genotype_gvcfs.log",
        config["files_path"]["log_dir"] + "/{interval_name}-{scatteritem}-genotype_gvcfs.e"
    threads: 4
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt3x'])
    benchmark:
        config["files_path"]["benchmark"] + "/{interval_name}-{scatteritem}_genotype_gvcfs.tsv"
    envmodules:
        "gatk/4.1.9.0"
    message: """ GenotypeGVCFs """
    shell:
        """
        {params.gatk} --java-options "{params.java_opt}" GenotypeGVCFs -R {params.ref_genome} -L {input.import_interval} -V gendb://{input.import_db} -O {output[0]} {params.fixed_args} > {log[0]} 2> {log[1]}
        """
        # {params.gatk} --java-options "{params.java_opt}" GenotypeGVCFs -R {params.ref_genome} -V gendb://{input.import_db} -O {output[0]} {params.fixed_args} > {log[0]} 2> {log[1]}


# this rule should be then expanded to be the rule that concat the vcf back together, after genotyping
rule chrom_intervals_gather:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    output:
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genotype_gvcfs").get("out_dir"),"all.{interval_name}.vcf.gz"),
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genotype_gvcfs").get("out_dir"),"all.{interval_name}.vcf.gz.tbi")
    input:
        gather.split(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_genotype_gvcfs").get("out_dir"),"{scatteritem}_{{interval_name}}/all.{scatteritem}_{{interval_name}}.vcf.gz"))
    params:
        bcftools=config['BCFTOOLS'],
        tmp=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("files_path").get("tmp")),
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    log:
        config["files_path"]["log_dir"] + "/{interval_name}-chrom_intervals_gather.log",
        config["files_path"]["log_dir"] + "/{interval_name}-chrom_intervals_gather.e"
    threads: 3
    benchmark:
        config["files_path"]["benchmark"] + "/{interval_name}_chrom_intervals_gather.tsv"
    envmodules:
        "bcftools/1.11"
    message: """Let\'s gather things together, by chromosome, basically!"""
    shell:
        """
        temp=$(mktemp -u -d -p {params.tmp})

        {params.bcftools} concat {input} | {params.bcftools} sort -T ${{temp}} | {params.bcftools} norm -f {params.ref_genome} -O z -o {output[0]} > {log[0]} 2> {log[1]}
        {params.bcftools} index -t {output[0]}
        """
