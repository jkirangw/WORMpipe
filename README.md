# WORMpipe

## This is a worm genome analysis pipeline for:

* Quality control of hifi long reads
* Estimation of major genome features such as size, repeatitiveness, heterozygosity 
* De novo genome assembly using either flye or Hifiasm
* Assesment of assembly metrics such as N50 contig
* Interrogation of assembly contamination and decontamination 
* BUSCO assessment of draft genome assemblies

## Run the pipeline with default hifiasm mode

```nextflow run WORMpipe.nf -c WORMpipe.config --reads 312-11_deduplicated.ccs.fastq --assembler hifiasm```

## For help 

```nextflow run WORMpipe.nf -help```

## To restart the run from last succesful steps with flye assembly mode

```nextflow run WORMpipe.nf -c WORMpipe.config --reads 312-11_deduplicated.ccs.fastq --assembler --flye -resume```
