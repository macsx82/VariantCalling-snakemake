rule recal_pass_filter:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    output:
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("recal_pass_filter").get("out_dir"),"/all.{interval_name}.PASS.vcf.gz"),
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("recal_pass_filter").get("out_dir"),"/all.{interval_name}.PASS.vcf.gz.tbi")
    input:
        rules.gatk_apply_VQSR.output[1]
    params:
        bcftools=config["BCFTOOLS"],
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    log:
        config["files_path"]["log_dir"] + "/all.{interval_name}-recal_pass_filter.log",
        config["files_path"]["log_dir"] + "/all.{interval_name}-recal_pass_filter.e"
    benchmark:
        config["files_path"]["benchmark"] + "/all.{interval_name}_recal_pass_filter.tsv"
    envmodules:
        "bcftools/1.11"
    threads: 2
    message: """ Filter only PASS variants after VQSR"""
    shell:
        """
        {params.bcftools} view -i "FILTER=='PASS'" {input} | {params.bcftools} norm -f {params.ref_genome} -O z -o {output[0]}
        {params.bcftools} intex -t {output[0]}
        """


rule rsid_annotation:
    output:
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("rsid_annotation").get("out_dir"),"/{current_chr}.PASS_rsID.vcf.gz"),
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("rsid_annotation").get("out_dir"),"/{current_chr}.PASS_rsID.vcf.gz.tbi")
    input:
        rules.recal_pass_filter.output[0]
        # os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("recal_pass_filter").get("out_dir"),"/all.{interval_name}.PASS.vcf.gz"),
    params:
        bcftools=config["BCFTOOLS"],
        dbsnp_latest=config.get("known_variants").get("dbsnp_latest"),
        current_chr=lambda wildcards, input : get_chr_from_vcf(input),
        out_folder=os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("rsid_annotation").get("out_dir"))
    log:
        config["files_path"]["log_dir"] + "/all.{current_chr}-rsid_annotation.log",
        config["files_path"]["log_dir"] + "/all.{current_chr}-rsid_annotation.e"
    benchmark:
        config["files_path"]["benchmark"] + "/all.{current_chr}_rsid_annotation.tsv"
    envmodules:
        "bcftools/1.11"
    threads: 2
    message: """ Add rsID using latest dbSNP information """
    shell:
        """
        echo "working on chromosome {params.current_chr}"
        
        {params.bcftools} annotate -a {params.dbsp_latest} -c ID {input} | {params.bcftools} +fill-tags -O z -o {params.outfolder}/{params.current_chr}.PASS_rsID.vcf.gz
        {params.bcftools} index -t {params.outfolder}/{params.current_chr}.PASS_rsID.vcf.gz
        """