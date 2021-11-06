#phasing rule using Eagle (so we can phase also non par chrX!!)
    rule phase:
    wildcard_constraints:
        interval_name='wgs_calling_regions_.+.interval_list'
    output:
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("phase").get("out_dir"),"{interval_name}.phased.vcf.gz"),
        os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("phase").get("out_dir"),"{interval_name}.phased.vcf.gz.tbi")
    input:
        rules.rsid_annotation.output[0],
        rules.rsid_annotation.output[1]
    params:
        phase_bin=config["PHASE_TOOL"],
        bcftools=config["BCFTOOLS"],
        genetic_map=config.get("rules").get("phase").get("genetic_map")
    log:
        config["files_path"]["log_dir"] + "/all.{interval_name}-phase.log",
        config["files_path"]["log_dir"] + "/all.{interval_name}-phase.e"
    benchmark:
        config["files_path"]["benchmark"] + "/all.{interval_name}_phase.tsv"
    envmodules:
        "bcftools/1.11"
    threads: 8
    resources:
        mem_mb=50000
    message: """ Phase vcf files """
    shell:
        """
        {params.phase_bin} --geneticMapFile {params.genetic_map} --vcf {input[0]} --outPrefix {output[0]} --vcfOutFormat z --numThreads {threads} --noImpMissing 1> {log[0]} 2> {log[1]} 
        {params.bcftools} index -t {output[0]}
        """