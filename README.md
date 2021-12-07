# WORMpipe

## This is a worm genome analysis pipeline for:

* Quality control of hifi long reads using fastp, kmer counting using either jellyfish, kmc, kat
* Estimation of major genome features such as size, repeatitiveness, heterozygosity using genomescope or Bbtools
* De novo genome assembly using either flye,Hifiasm,Hicanu,Nextdenovo
* Assesment of assembly metrics such as N50 contig using Quast
* Interrogation of assembly contamination and decontamination using blobtools2
* BUSCO assessment of draft genome assemblies
* Assembly polishing with Hapo-G,NextPolish, Racon or Pilon
* Purge false duplicates with purge_dups or purge_haplotigs
* Mapping RNA-SEQ reads to draft assembly using GMAP-GSNAP
* Evaluation of mapped reads with qualimap
* Genome annotation with BRAKER2 pipeline


## Pipeline flow chart

![alt text](https://github.com/jkirangw/WORMpipe/blob/main/WORMpipe_flowchart.png)

## Run the pipeline with default hifiasm mode

```nextflow run WORMpipe.nf -c WORMpipe.config --reads 312-11_deduplicated.ccs.fastq --assembler hifiasm```

## For help 

```nextflow run WORMpipe.nf -help```

## To restart the run from last succesful steps with flye assembly mode

```nextflow run WORMpipe.nf -c WORMpipe.config --reads 312-11_deduplicated.ccs.fastq --assembler --flye -resume```
