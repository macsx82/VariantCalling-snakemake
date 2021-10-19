import os.path
import gzip 
import re
import sys
import errno
import multiprocessing
import psutil

def get_rg_id_data(wildcards):
    fastq= list(samples_df[samples_df.SAMPLE_ID == (wildcards.sample).split(sep="_")[0]].fq1)[0]
    fastq_file= gzip.open('%s' %(fastq) ,'r')
    header = fastq_file.readline().strip()
    fastq_file.close()
    #we need to split by the space: data should be in Casava 1.8 or above format
    header_sections=[rg_id for section in header.decode('utf-8').split(' ') for rg_id in section.split(":") ]
    #define the RG tag dictionary to be used 
    RG_TAG={}
    RG_TAG["PU1"]=header_sections[0]
    RG_TAG["PU2"]=header_sections[1]
    RG_TAG["ID1"]=header_sections[2]
    RG_TAG["FL"]=header_sections[3]
    RG_TAG["TNFL"]=header_sections[4]
    RG_TAG["XX"]=header_sections[5]
    RG_TAG["YY"]=header_sections[6]
    RG_TAG["PAIR"]=header_sections[7]
    RG_TAG["FIL"]=header_sections[8]
    RG_TAG["BITS"]=header_sections[9]
    RG_TAG["ID2"]=header_sections[10]
    return RG_TAG

def get_input_fasta(wildcards):
    r1 = samples_df[samples_df.SAMPLE_ID == (wildcards.sample).split(sep="_")[0]].fq1
    r2 = samples_df[samples_df.SAMPLE_ID == (wildcards.sample).split(sep="_")[0]].fq2
    return r1, r2
    # @ST-E00233 	the unique instrument name
    # 168 	the run id
    # HMKNCCCXX 	the flowcell id
    # 7 	flowcell lane
    # 1101 	tile number within the flowcell lane
    # 6593 	'x'-coordinate of the cluster within the tile
    # 1344 	'y'-coordinate of the cluster within the tile
    # 1 	the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
    # N 	Y if the read is filtered, N otherwise
    # 0 	0 when none of the control bits are on, otherwise it is an even number
    # NGTACG 	index sequence 
def cpu_count():
    return multiprocessing.cpu_count()

def conservative_cpu_count(reserve_cores=1, max_cores=5):
    cores = max_cores if cpu_count() > max_cores else cpu_count()
    return max(cores - reserve_cores, 1)    

#get resources in mb from java options
def get_resources_from_jvm(jvm_option):
    jvm_Xms_mb=int(re.search(r"\d+",jvm_option.split(" ")[0]).group(0))*1024*1.5
    return int(jvm_Xms_mb)

#get sex from sample file
def get_sample_sex(wildcards):
    sex=list(samples_df[samples_df.SAMPLE_ID == (wildcards.sample).split(sep="_")[0]].sex)[0]
    return sex

def get_chr_from_interval_list(wildcards):
    interval_filename=wildcards.interval_name
    interval_file=os.path.join(*references_abs_path(), config.get("callable_intervals"),interval_filename)
    for line in open(interval_file, 'r'):
        if not(re.match('@', line.strip())):
            interval_chr=line.strip().split("\t")[0]
            break

    return interval_chr

#get the chromosome name from the input vcf file
def get_chr_from_vcf(vcf):
    vcf_filename=vcf
    # for line in open(vcf_filename, 'r'):
    for line in gzip.open(vcf_filename, 'r') if vcf_filename.endswith('.gz') else open(vcf_filename, 'r'):
        if not(re.match('#', line.strip())):
            vcf_chr=line.strip().split("\t")[0]
            break

    return vcf_chr


def references_abs_path(ref='references'):
    references = config.get(ref)
    basepath = references['basepath']
    provider = references['provider']
    release = references['release']

    return [os.path.join(basepath, provider, release)]

def resolve_single_filepath(basepath, filename):
    return [os.path.join(basepath, filename)]

def resolve_multi_filepath(basepath, dictionary):
    for k, v in dictionary.items():
        dictionary[k] = os.path.join(basepath, v)
    return dictionary

def get_interval_file(wildcards):
    interval_file=os.path.join(*references_abs_path(), config.get("callable_intervals"),wildcards.interval_name)
    return interval_file

def get_sex_chr_from_interval(wildcards):
    interval_filename=wildcards.interval_name
    interval_file=os.path.join(*references_abs_path(), config.get("callable_intervals"),interval_filename)
    for line in open(interval_file, 'r'):
        if not(re.match('@', line.strip())):
            interval_chr=line.strip().split("\t")[0]
            if not(re.match('chrY',interval_chr)):
                if not(re.match('chrX',interval_chr)):
                    chr_mode="AUTOSOMAL"
                    break
                else:
                    #need to check if we are working with PAR regions
                    # chrX    10001   2781479 PAR1
                    # chrX    155701383       156030895       PAR2
                    start=int(line.strip().split("\t")[1])
                    end=int(line.strip().split("\t")[2])
                    #we are checking if we are outside the par regions
#                   (c[1] <= 10001 || (c[1] >= 2781479 && c[1] <= 155701383) || c[1] >= 156030895) && (c[2] <= 10001 ||( c[2] >= 2781479 && c[2] <= 155701383)|| c[2] >= 156030895) ) print $0}' | wc -l | cut -f 1 -d " ")
                    if (start < 10001 or start > 2781479) and (start < 155701383 or start > 156030895) and (end < 10001 or end > 2781479) and (end < 155701383 or end > 156030895) :
                        chr_mode="SEXUAL"
                        break
                    else:
                        chr_mode="AUTOSOMAL"
                        break
            else:
                chr_mode="SEXUAL"
                break
    return chr_mode


#function to get parameters needed to perform VQSR step with GATK
def _get_recal_params(wildcards):
    known_variants = resolve_multi_filepath("", config["known_variants"])
    if wildcards.type == "snp":
        base_params=("--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} --resource:omni,known=false,training=true,truth=true,prior=12.0 {omni} --resource:1000G,known=false,training=true,truth=false,prior=10.0 {g1k} --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 {dbsnp}").format(**known_variants)
        tranches=config["rules"]["gatk_variant_recalibrator"]["SNP"]["tranches"]
        annotations=config["rules"]["gatk_variant_recalibrator"]["SNP"]["annotations"]
        arguments=config["rules"]["gatk_variant_recalibrator"]["SNP"]["arguments"]
    else:
        base_params=("--resource:mills,known=false,training=true,truth=true,prior=12.0 {mills} --resource:axiomPoly,known=false,training=true,truth=false,prior=10 {axiom} --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp}").format(**known_variants)
        tranches=config["rules"]["gatk_variant_recalibrator"]["INDEL"]["tranches"]
        annotations=config["rules"]["gatk_variant_recalibrator"]["INDEL"]["annotations"]
        arguments=config["rules"]["gatk_variant_recalibrator"]["INDEL"]["arguments"]
    return "%s %s %s %s" %(base_params, tranches,annotations,arguments)




