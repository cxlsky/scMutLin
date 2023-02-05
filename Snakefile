"""
  scGLAS: single-cell Genetic Lineage Analysis Suite
  author: Xiaolong Cheng from Dr. Wei Li lab
  email: xcheng@childrensnational.org
"""
import glob
import os

configfile: "config.yaml"

def get_filelist(fpath):
    f = []
    for (dirpath, dirnames, filenames) in os.walk(fpath):
        f.extend(filenames)
        break
    return [i[0:-4] for i in f if i.endswith('.bam')]

SAMPLES = get_filelist(config["input"]["path"])


rule all:
    input: 
        expand("results/mutation/{sample}_raw.vcf.gz", sample=SAMPLES)

filepath = config["input"]["path"]

### preprocessing
if config["preprocessing"]["gatk"]:
    rule MarkingDuplicates:
        input:
            filepath+"/{sample}.bam"
        output:
            bam = temp("results/preprocess/{sample}_sorted_dupMarked.bam"),
            mtrx= temp("results/preprocess/{sample}_dupMetrics")
        shell:
            "gatk MarkDuplicates --INPUT {input}  --OUTPUT {output.bam}"
            " --METRICS_FILE {output.mtrx} > /dev/null 2>&1"
    rule RepairingReadgroups:
        input:
            "results/preprocess/{sample}_sorted_dupMarked.bam"
        output:
            temp("results/preprocess/{sample}_sorted_DM_RG.bam")
        shell:
            "gatk AddOrReplaceReadGroups --INPUT {input}"
            " --OUTPUT {output} --RGLB illumina_{wildcards.sample}"
            " --RGPL illumina --RGPU JR001 --RGSM {wildcards.sample} > /dev/null 2>&1"
    rule BaseRecalibrator:
        input:
            "results/preprocess/{sample}_sorted_DM_RG.bam"
        output:
            temp("results/preprocess/{sample}_recal.table")
        params:
            ref = config["input"]["ref"],
            dbsnp = config["preprocessing"]["dbsnp"]
        shell:
            "gatk BaseRecalibrator --input {input} --output {output}"
            " --reference {params.ref} --known-sites {params.dbsnp} > /dev/null 2>&1"
    rule ApplyBQSR:
        input:
            bam = "results/preprocess/{sample}_sorted_DM_RG.bam",
            tab = "results/preprocess/{sample}_recal.table"
        output:
            "results/preprocess/{sample}_sorted.bam"
        params:
            ref = config["input"]["ref"]
        shell:
            "gatk ApplyBQSR -R {params.ref} -I {input.bam}"
            " --bqsr-recal-file {input.tab} -O {output} > /dev/null 2>&1"
else:
    rule sortbam:
        input:
            filepath+"/{sample}.bam"
        output:
            bam = "results/preprocess/{sample}_sorted.bam",
            idx = "results/preprocess/{sample}_sorted.bam.bai"
        threads:
            config["mutation"]["threads"]
        run:
            shell("samtools sort -@ {threads} -o {output.bam} {input} > /dev/null 2>&1")
            shell("samtools index -@ {threads} {output.bam} {output.idx} > /dev/null 2>&1")

### Call Mutation
lst_call = list(config["mutation"]["caller"].split(','))
for caller in lst_call:
    if caller == 'gatk':
        rule gatkMutect2:
            input:
                "results/preprocess/{sample}_sorted.bam"
            output:
                "results/mutation/{sample}_raw.vcf.gz"
            params:
                ref = config["input"]["ref"]
            shell:
                "gatk Mutect2 -R {params.ref} -I {input}"
                " -O {output} > /dev/null 2>&1"
    elif caller == 'bcftools':
        rule bcftoolsSNP:
            input:
                "results/preprocess/{sample}_sorted.bam"
            output:
                "results/mutation/{sample}_raw.vcf.gz"
            params:
                ref = config["input"]["ref"]
            threads:
                config["mutation"]["threads"]
            shell:
                "bcftools mpileup --threads {threads} -Ou -f {params.ref} {input} |"
                "  bcftools call --threads {threads} -vm -Oz  > {output}"
    elif caller == 'strelka2':
        rule strelkaConfig:
            input:
                "results/preprocess/{sample}_sorted.bam"
            output:
                "results/strelka2/{sample}/runWorkflow.py"
            conda:
                "src/strelka.yaml"
            params:
                ref = config["input"]["ref"],
                path = "results/strelka2/{sample}"
            shell:
                "configureStrelkaGermlineWorkflow.py --bam {input}"
                " --referenceFasta {params.ref} --runDir {params.path}"
        rule runStrelka2:
            input:
                "results/strelka2/{sample}/runWorkflow.py"
            output:
                "results/strelka2/{sample}/results/variants/variants.vcf.gz" 
            threads:
                config["mutation"]["threads"]
            conda:
                "src/strelka.yaml"
            shell:
                "python {input} -m local -j {threads} > /dev/null 2>&1"
        rule moveVCF:
            input:
                "results/strelka2/{sample}/results/variants/variants.vcf.gz"
            output:
                "results/mutation/{sample}_raw.vcf.gz"
            shell:
                "cp -rf {input} {output}"




