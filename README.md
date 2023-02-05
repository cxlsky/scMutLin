# scMutLin: single-cell Mutation calling and genetic Lineage analysis #

scMutLin is designed to provide a stable and scalable tool/pipeline for calling somatic mutations and analyzing genetic lineage relationships between cells and evaluate the performance of somatic-mutation-based lineage analysis. We will keep to update new functions and explore the application of genetic lineage tracing.

### Installation ###

clone the repository to your work folder and run this command:

```
    conda env create --file environment.yml
```
	
Activate the environment:

```
    conda activate scMutLin
```

### Usage ###

* scMutLin starts from aligned bam files. If you only have one merged bam file from 10x platform, we provide a python script to split it based on the *CB* tag. See the ***Demo***.

* All the configurations can be set in the config.yaml file

* Once the parameters are configured, run the pipeline by command:

```
    snakemake --cores 8 --use-conda
```

### How to set config.yaml file ###

##### Configurartion of INPUT ####

* We should specify the path for .bam files

* We should specify the path for reference genome fasta file

* Other parameters are used for testing currently

##### Configurartion of PREPROCESSING  ####

* If you'd like to preprocess the bam files by [GATK Best Practices Workflows](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912), please set *gatk* to **true**. Our pipeline will Mark Duplicates and Recalibrate Base Quality Score for each bam file.

* You can also skip the preprocessing step by setting *gatk* to **false**. Then our pipeline will just sort and index the bam files.

* Note: 
    1. A dbSNP file is required for the PREPROCESSING step (you can download the file in the **Resource** section)
    2. The PREPROCESSING step will cost huge memory (up to 35GB)


##### Configurartion of MUTATION CALLING ####

* Please specify the mutation calling method, you can choose gatk, bcftools or strelka2. 

* Set threads for mutation caller (not suitable for GATK)


### Demo:  Split 10X bam file and call mutation by Strelka2###

We provide a [demo](https://bitbucket.org/weililab/scmutlin/downloads/demo.zip) to show how to use this pipeline to split bam file and call mutation by Strelka2. You can run the demo by following steps:

* activate the env

```
    conda activate scMutLin
```

* download the demo file

```
    wget https://bitbucket.org/weililab/scmutlin/downloads/demo.zip
    
    unzip demo.zip
    
    cd demo
```

* split the merged bam file

```
    mkdir input/bam
    
    python src/split10xbam.py -i input/scRNA_tiny.bam -b input/barcodes.tsv -o input/bam
```

* edit the config.yaml file and configure the parameters (in this demo, you just need to set path for the mm10 fasta file)

* run the pipeline by command:

```
    snakemake --cores 8 --use-conda
```

* **Note**: the parameter *--use-conda* is just for Strelka2 as is is designed in Python 2. For bcftools or GATK, you can run this command:

```
    snakemake --cores 8
```




### Resource ###

* The preprocessing step is based on [GATK tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery) and a dbSNP file is required. You can download the dbSNP file from:
 
  1. [dnSNP for hg38](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/)
  
  2. [dbSNP for GRCm38](ftp://ftp-mouse.sanger.ac.uk/REL-1303-SNPs_Indels-GRCm38/)
