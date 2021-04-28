# WORMpipe

## This is a worm genome analysis pipeline for:

<<<<<<< HEAD
* Quality control of hifi long reads
* Estimation of major genome features such as size, repeatitiveness, heterozygosity 
* De novo genome assembly
* Assesment of assembly metrics such as N50 contig
* Assembly contamination interrogation and decontamination
=======
* Quality control of hifi long reads 
* Estimation of major genome features such as size, repeatitiveness, heterozygosity
* De novo genome assembly
* Assesment of assembly metrics such as N50 contig
* Assembly contamination interrogation and decontamination 
>>>>>>> 39c11fb08fcad42c8e59a6a6464cc281e5cae536
* Assembly gene density assessment using BUSCOs

## Run the pipeline

```nextflow run WORMpipe.nf -c WORMpipe.config --reads 312-11_deduplicated.ccs.fastq --assembler hifiasm```
<<<<<<< HEAD

## For help 
=======
## For help
>>>>>>> 39c11fb08fcad42c8e59a6a6464cc281e5cae536

```nextflow run WORMpipe.nf -help```

## To restart the run from last succesful steps

```nextflow run WORMpipe.nf -c WORMpipe.config --reads 312-11_deduplicated.ccs.fastq --assembler --flye -resume```
