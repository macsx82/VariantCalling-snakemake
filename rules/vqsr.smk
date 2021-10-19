
rule gatk_variant_recalibrator:
    output:
        recal=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_variant_recalibrator").get("out_dir"),"/ALL.{prefix}.{type,(snp|indel)}.recal"),
        tranches=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_variant_recalibrator").get("out_dir"),"/ALL.{prefix}.{type,(snp|indel)}.tranches"),
        plotting=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_variant_recalibrator").get("out_dir"),"/ALL.{prefix}.{type,(snp|indel)}.plotting.R")
    input:
        resolve_multi_filepath(*references_abs_path(), config["known_variants"]).values(),
        # vcf=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("pre_vqsr_rules").get("out_dir"),"ALL.{prefix}.vcf.gz")
        vcf=rules.concat_vcfs.output[0]
    params:
        gatk=config['GATK_TOOL'],
        java_opt=config['java_opts']['opt3x'],
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        recal=_get_recal_params
    log:
        config["files_path"]["log_dir"] + "/ALL.{prefix}.{type}-gatk_variant_recalibrator.log",
        config["files_path"]["log_dir"] + "/ALL.{prefix}.{type}-gatk_variant_recalibrator.e"
    threads: config.get("rules").get("gatk_variant_recalibrator").get("threads")
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt3x'])
    benchmark:
        config["files_path"]["benchmark"] + "/ALL.{prefix}.{type}_gatk_variant_recalibrator.tsv"
    envmodules:
        "gatk/4.2.2.0"
    message: """ VariantRecalibrator """
    shell:
        """
        {params.gatk} --java-options "{params.java_opt}" VariantRecalibrator -R {params.genome} -V {input.vcf} {params.recal} -O {output.recal} --tranches-file {output.tranches} --rscript-file {output.plotting} > {log[0]} 2> {log[1]}
        """

#this rule will work again by chromosome and variant type. At the moment this is written in a dumb way, to be optimized later.
#we are not generating separate files for indel and snp, but applying the filter to the same original vcf.
rule gatk_apply_VQSR:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    output:
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_apply_VQSR").get("out_dir"),"/all.{interval_name}.indel_recalibrated.vcf.gz"),
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_apply_VQSR").get("out_dir"),"/all.{interval_name}.indel_recalibrated.snp_recalibrated.vcf.gz")
    input:
        vcf=rules.clean_and_excess_het_filter.output,
        recal_s=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_variant_recalibrator").get("out_dir"),"variant_calling/{prefix}.snp.recal"),
        recal_i=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_variant_recalibrator").get("out_dir"),"variant_calling/{prefix}.indel.recal"),
        tranches_s=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_variant_recalibrator").get("out_dir"),"variant_calling/{prefix}.snp.tranches"),
        tranches_i=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("gatk_variant_recalibrator").get("out_dir"),"variant_calling/{prefix}.indel.tranches")
    params:
        gatk=config['GATK_TOOL'],
        java_opt=config['java_opts']['opt3x'],
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        args_snp= config.get("rules").get("gatk_apply_VQSR").get("SNP").get("arguments"),
        args_indel= config.get("rules").get("gatk_apply_VQSR").get("INDEL").get("arguments")
        # mode=lambda wildcards: wildcards.type.upper()
    log:
        config["files_path"]["log_dir"] + "/all.{interval_name}-gatk_variant_recalibrator.log",
        config["files_path"]["log_dir"] + "/all.{interval_name}-gatk_variant_recalibrator.e"
    threads: config.get("rules").get("gatk_apply_VQSR").get("threads")
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt3x'])
    benchmark:
        config["files_path"]["benchmark"] + "/all.{interval_name}_gatk_variant_recalibrator.tsv"
    envmodules:
        "gatk/4.2.2.0"
    message: """ ApplyVQSR """
    shell:
        """
        echo "INDEL first..."
        {params.gatk} --java-options "{params.java_opt}" ApplyVQSR -R {params.genome} -V {input.vcf} -mode INDEL --recal-file {input.recal_i} --tranches-file {input.tranches_i} {params.args_indel} -O {output[0]} 1> {log[0]} 2> {log[1]}
        
        echo "SNP now..."
        {params.gatk} --java-options "{params.java_opt}" ApplyVQSR -R {params.genome} -V {output[0]} -mode SNP --recal-file {input.recal_s} --tranches-file {input.tranches_s} {params.args_snp} -O {output[1]} 1>> {log[0]} 2>> {log[1]}
        """

