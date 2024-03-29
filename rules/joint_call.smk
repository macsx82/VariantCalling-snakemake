#rules for joint call

#split the callable intervals in different chunks, so we can better parallelize
rule split_intervals:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    output:
        scatter.split(os.path.join(config.get('files_path').get('base_joint_call_path'),config.get('rules').get('split_intervals').get('out_dir')) + '/{scatteritem}_{{interval_name}}')
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
    
if pipeline_mode == "JOINT_CALL" :
    include:
        # include_prefix + "/gatk_genomicsDBI.smk"
        "gatk_genomicsDBI.smk"
else :
    include:
        "gatk_genomicsDBI_update.smk"


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
    resources:
        mem_mb=5000
    benchmark:
        config["files_path"]["benchmark"] + "/{interval_name}_chrom_intervals_gather.tsv"
    envmodules:
        "bcftools/1.14"
    message: """Let\'s gather things together, by chromosome, basically!"""
    shell:
        """
        temp=$(mktemp -u -d -p {params.tmp})

        {params.bcftools} concat {input} | {params.bcftools} sort -T ${{temp}} | {params.bcftools} norm -f {params.ref_genome} -O z -o {output[0]} > {log[0]} 2> {log[1]}
        {params.bcftools} index -t {output[0]}
        """
