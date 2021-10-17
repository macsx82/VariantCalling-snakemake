#set of rules to perform cleaning and to get ready for vqsr.We are following the best practices detailed here:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#1

#get only sites with AN>0, in order to remove useless stuff
#and apply the ExcessHet filter to our data
rule clean_and_excess_het_filter:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    output:
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("pre_vqsr_rules").get("out_dir"),"all.{interval_name}.CLEAN.vcf.gz")
    input:
        rules.chrom_intervals_gather.output[0]
    params:
        bcftools=config['BCFTOOLS']
    log:
        config["files_path"]["log_dir"] + "/{interval_name}-clean_vcf.log",
        config["files_path"]["log_dir"] + "/{interval_name}-clean_vcf.e"
    threads: 3
    benchmark:
        config["files_path"]["benchmark"] + "/{interval_name}_clean_vcf.tsv"
    envmodules:
        "bcftools/1.11"
    message: """Let\'s clean a little bit, by chromosome!"""
    shell:
        """
        {params.bcftools} view -i "ALT!='.'" {input} | {params.bcftools} filter -s ExcessHet -e "ExcessHet > 54.69"  -O z -o {output} > {log[0]} 2> {log[1]}
        """

#generate site only vcf files to pass to vqsr
rule sites_only_vcf:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    output:
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("pre_vqsr_rules").get("out_dir"),"all.{interval_name}.CLEAN.SITES_ONLY.vcf.gz"),
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("pre_vqsr_rules").get("out_dir"),"all.{interval_name}.CLEAN.SITES_ONLY.vcf.gz.tbi")
    input:
        rules.clean_and_excess_het_filter.output
    params:
        bcftools=config['BCFTOOLS']
    log:
        config["files_path"]["log_dir"] + "/{interval_name}-sites_only_vcf.log",
        config["files_path"]["log_dir"] + "/{interval_name}-sites_only_vcf.e"
    benchmark:
        config["files_path"]["benchmark"] + "/{interval_name}_sites_only_vcf.tsv"
    envmodules:
        "bcftools/1.11"
    message: """Create sites only files, by chromosome!"""
    shell:
        """
        {params.bcftools} view -G {input} -O z -o {output[0]} > {log[0]} 2> {log[1]}
        {params.bcftools} index -t {output[0]}

        """

#concat the site only vcfs to be passed to vqsr steps
rule concat_vcfs:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    input:
        expand(os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("pre_vqsr_rules").get("out_dir"),"all.{interval_name}.CLEAN.SITES_ONLY.vcf.gz"), interval_name=call_intervals)
    output:
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("pre_vqsr_rules").get("out_dir"),"ALL.CLEAN.SITES_ONLY.vcf.gz"),
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("pre_vqsr_rules").get("out_dir"),"ALL.CLEAN.SITES_ONLY.vcf.gz.tbi")
    params:
        bcftools=config['BCFTOOLS'],
        tmp=os.path.join(BASE_OUT,config.get("files_path").get("tmp"))
    log:
        config["files_path"]["log_dir"] + "/ALL-concat_vcfs.log",
        config["files_path"]["log_dir"] + "/ALL-concat_vcfs.e"
    threads: 2
    benchmark:
        config["files_path"]["benchmark"] + "/ALL_concat_vcfs.tsv"
    envmodules:
        "bcftools/1.11"
    message: """Concat sites only files!"""
    shell:
        """
        {params.bcftools} concat {input} | {params.bcftools} sort -T {params.tmp} -O z -o {output[0]} > {log[0]} 2> {log[1]}
        {params.bcftools} index -t {output[0]}
        """