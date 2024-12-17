import os
import pandas as pd
from collections import Counter


import yaml
with open('config.yaml') as s:
    d=yaml.safe_load(s)


sample_sheet=pd.read_csv(d['sample_sheet'], sep="\t",dtype=object)
genome = d['genome']
bed= d['bed']


DIRECTION=["1","2"]
GROUPS=set()


def check_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except FileExistsError:
        pass
        #print("Link for " + file1 + " is already present in 01_raw")

SAMPLES=list(sample_sheet["ID"])


for b in set(SAMPLES):
    os.makedirs("../01_raw/" +b+ "/fastqc", exist_ok=True)
    check_symlink(sample_sheet[sample_sheet["ID"]==b]["forward reads"].values[0], "../01_raw/"+b +"/"+ b +"_1P.fastq.gz")
    check_symlink(sample_sheet[sample_sheet["ID"]==b]["reverse reads"].values[0], "../01_raw/"+b +"/"+ b +"_2P.fastq.gz")
    os.makedirs("../02_trimmed/" + b +"/fastqc", exist_ok=True)
    os.makedirs("../03_mapped/"+b, exist_ok=True)


rule all:
    input:
        expand("../04_MethylDackel/{sample}/{sample}_CpG.bedGraph",sample=SAMPLES),
        expand("../03_mapped/{sample}/{sample}.mosdepth.summary.txt",sample=SAMPLES)



rule trimm_galore:
    input:
        r1= '../01_raw/{sample}/{sample}_1P.fastq.gz',
        r2= '../01_raw/{sample}/{sample}_2P.fastq.gz',
    output:
        p1="../02_trimmed/{sample}/{sample}_1P_val_1.fq.gz",
        p2="../02_trimmed/{sample}/{sample}_2P_val_2.fq.gz",
    conda:
        "envs/nf.yaml"
    threads: 4
    priority: 50
    shell:
        'trim_galore --fastqc --gzip --paired {input.r1} {input.r2} -o ../02_trimmed/{wildcards.sample} --clip_r1 8 --clip_r2 8 --three_prime_clip_r1 8 --three_prime_clip_r2 8 --cores {threads}'



rule bwameth:
    input:
        p1="../02_trimmed/{sample}/{sample}_1P_val_1.fq.gz",
        p2="../02_trimmed/{sample}/{sample}_2P_val_2.fq.gz" 
    output:
        "../03_mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam"
    conda:
        "envs/nf.yaml"
    threads: 8
    priority: 50
    shell:
        'bwameth.py --threads {threads} --reference {genome} {input.p1} {input.p2} |samtools view -@ 6 -m 4G -b | samtools sort -@ 6 -m 4G -o ../03_mapped/{wildcards.sample}/{wildcards.sample}Aligned.sortedByCoord.out.bam'



rule dedup_picard:
    input:
        bam = "../03_mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam"

    output:
        mkdupbam = "../03_mapped/{sample}/{sample}Aligned.sortedByCoord.markdup.out.bam",
        mkdupbai = "../03_mapped/{sample}/{sample}Aligned.sortedByCoord.markdup.out.bam.bai",
        picard_metrics = "../03_mapped/{sample}/{sample}_picard.txt"

    conda:
        "envs/nf.yaml"
    threads: 4
    priority: 50
    shell:
        "picard MarkDuplicates INPUT={input.bam} OUTPUT={output.mkdupbam} METRICS_FILE={output.picard_metrics} REMOVE_DUPLICATES=false ASSUME_SORTED=true PROGRAM_RECORD_ID='null'  VALIDATION_STRINGENCY=LENIENT ;\
        samtools index {output.mkdupbam}"

rule methyldackel:
    input:
        mkdupbam = "../03_mapped/{sample}/{sample}Aligned.sortedByCoord.markdup.out.bam",

    conda: 'envs/nf.yaml'
    output:
        cov1="../04_MethylDackel/{sample}/{sample}_CpG.bedGraph",
    threads:4
    shell: 'MethylDackel extract   -o ../04_MethylDackel/{wildcards.sample}/{wildcards.sample} {genome} {input.mkdupbam};\
             MethylDackel mbias {genome} {input.mkdupbam} ../04_MethylDackel/{wildcards.sample}/{wildcards.sample}  --txt > ../04_MethylDackel/{wildcards.sample}/{wildcards.sample}_methyldackel.txt'


rule TargetCov:
    input:
        mkdupbam = "../03_mapped/{sample}/{sample}Aligned.sortedByCoord.markdup.out.bam",

    conda: 'envs/mosdepth.yaml'
    output:
        cov1="../03_mapped/{sample}/{sample}.mosdepth.summary.txt",
    threads:4
    shell: 'mosdepth -t {threads} -b {bed} ../03_mapped/{wildcards.sample}/{wildcards.sample} {input.mkdupbam}'
