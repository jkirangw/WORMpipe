# WORMpipe

## This is a worm genome analysis pipeline for:

* Quality control of hifi long reads
* Estimation of major genome features such as size, repeatitiveness, heterozygosity 
* De novo genome assembly
* Assesment of assembly metrics
* Assembly contanination interrogation and decontamination
* Assembly gene density assesment 

## Run the pipeline

```nextflow run WORMpipe.nf -c WORMpipe.config --reads 312-11_deduplicated.ccs.fastq --assembler hifiasm```

## For help options

```nextflow run WORMpipe.nf -help```

## For restart the run from last succesful steps

```nextflow run WORMpipe.nf -c WORMpipe.config --reads 312-11_deduplicated.ccs.fastq --assembler --flye -resume```
