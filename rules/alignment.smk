#alignment rule
rule bwa_mem:
    output:
        temp(BASE_OUT +"/" + config["rules"]["bwa_mem"]["out_dir"] + "/{sample}/{sample}_map.bam"),
        temp(BASE_OUT +"/" + config["rules"]["bwa_mem"]["out_dir"] + "/{sample}/{sample}_map_sorted.bam")
    input:
        r1 = lambda wc: samples_df[samples_df.SAMPLE_ID == (wc.sample).split(sep="_")[0]].fq1,
        r2 = lambda wc: samples_df[samples_df.SAMPLE_ID == (wc.sample).split(sep="_")[0]].fq2
    params:
        bwa=config['ALIGN_TOOL'],
        sambamba=config['SAMBAMBA_TOOL'],
        samtools=config['SAMTOOLS'],
        tmp=os.path.join(BASE_OUT,config.get("files_path").get("tmp")),
        LB=config['lb'],
        PL=config['pl'],
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        rg_tag=get_rg_id_data
    log:
        config["files_path"]["log_dir"] + "/{sample}-mapping.log",
        config["files_path"]["log_dir"] + "/{sample}-mapping.e"
    benchmark:
        config["files_path"]["benchmark"] + "/{sample}_mapping.tsv"
    envmodules:
        "bwa-mem2/2.1/gnu/9.3.0",
        "sambamba/0.8.0",
        "samtools/1.14"
    threads:
        config["rules"]["bwa_mem"]["threads"]
    resources:
        mem_mb=config["rules"]["bwa_mem"]["mem"]
    priority: 50
    message: """ Generating a mapped bam file to be merged with the unmapped one, to not lose information."""
    shell:
        """
        ({params.bwa} mem -R \"@RG\tID:{params.rg_tag[PU1]}.{params.rg_tag[PU2]}\tSM:{wildcards.sample}\tLB:{params.LB}\tPL:{params.PL}\" -K 10000000 -v 3 -t {threads} -Y {params.ref_genome} {input.r1} {input.r2} | {params.samtools} view -1 -o {output[0]}) 2> {log[1]} 1> {log[0]}
        {params.sambamba} sort --tmpdir={params.tmp} -t {threads} --sort-picard --out={output[1]} {output[0]} 2>> {log[1]} 1>> {log[0]}
        """

#mark duplicate rule
rule mark_dup:
    output:
        BASE_OUT +"/" + config["rules"]["mark_dup"]["out_dir"] + "/{sample}/{sample}_mapped_md.bam"
    input:
        rules.bwa_mem.output[1]
        # rules.map_unmap_merge.output
    params:
        sambamba=config['SAMBAMBA_TOOL'],
        tmp=os.path.join(BASE_OUT,config.get("files_path").get("tmp")),
        cl=config['cl']
    log:
        config["files_path"]["log_dir"] + "/{sample}-mark_dup.log",
        config["files_path"]["log_dir"] + "/{sample}-mark_dup.e"
    benchmark:
        config["files_path"]["benchmark"] + "/{sample}_mark_dup.tsv"
    envmodules:
        "sambamba/0.8.0"
    threads:
        config["rules"]["mark_dup"]["threads"]
    resources:
        mem_mb=config["rules"]["mark_dup"]["mem"]
    priority: 49
    message: """ Mark duplicate reads to avoid counting non-independent observations"""
    shell:
        """
        {params.sambamba} markdup -t {threads} --compression-level={params.cl} --tmpdir={params.tmp}/ {input} {output} 2> {log[1]} 1>{log[0]}
        """

rule sort_bam:
    output:
        BASE_OUT +"/" + config["rules"]["sort_bam"]["out_dir"] + "/{sample}/{sample}_mapped_md_sorted.cram",
        BASE_OUT +"/" + config["rules"]["sort_bam"]["out_dir"] + "/{sample}/{sample}_mapped_md_sorted.crai"
    input:
        rules.mark_dup.output
    params:
        sambamba=config['SAMBAMBA_TOOL'],
        samtools=config['SAMTOOLS'],
        tmp=os.path.join(BASE_OUT,config.get("files_path").get("tmp")),
        cl=config['cl'],
        java_opt=config['java_opts']['opt2x'],
        ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    log:
        config["files_path"]["log_dir"] + "/{sample}-mark_dup.log",
        config["files_path"]["log_dir"] + "/{sample}-mark_dup.e"
    benchmark:
        config["files_path"]["benchmark"] + "/{sample}_mark_dup.tsv"
    envmodules:
        "sambamba/0.8.0",
        "samtools/1.14"
    threads:
        config["rules"]["sort_bam"]["threads"]
    resources:
        mem_mb=get_resources_from_jvm(config['java_opts']['opt2x'])
    priority: 48
    message: """ Sort BAM file by coordinate order and generate the CRAM file: CRAM recalculates MD and NM tags on the fly, so we don't need to calculate them. """
    shell:
        """
        ({params.sambamba} sort -t {threads} -m 1G --tmpdir={params.tmp}/ --out=/dev/stdout {input} | {params.samtools} view -@ {threads} -h -T {params.ref_genome} -C -o {output[0]}) 2> {log[1]} 1> {log[0]}
        {params.samtools} index {output[0]} {output[1]} 2>> {log[1]} 1>> {log[0]}
        """

