#!/usr/bin/env nextflow

//This is WORMpipe used for long read genome preprocessing and genome assembly and conatamination assessement

def help_message() {
   log.info """
      This is WORMpipe used for HiFi long read preprocessing, de novo genome assembly and conatamination assessement

      Usage:

      Typical command for running WORMpipe pipeline is as follows:
        nextflow run WORMpipe.nf -c WORMpipe.config --read readsample.fastq

      Mandatory arguments:
        --read		fastq read file that you want to analyse

      Optional arguments:
	--help		This is usage statement
	--outdir	Output directory to place the final results [results]
	--fastp_out	Prefix the name of fastp output file [fast.out]
	--outFileName	Prefix the name of output file [out.file]
	--flye		--flye flag on command line for flye assembly [hifiasm] 
      """
}

// Show help message

if(params.help){

	help_message()
	exit 0
}

// create a queue channel
Channel
   .fromPath(params.read)
   .ifEmpty { error "Cannot find any reads matching: ${params.read}" }
   .into { fastp_ch; jellyfish_ch; assembly_ch; flye_ch; hifiasm_ch; minimap2_ch }
	
// read preprocessing
process run_fastp {
	tag "$read_file"
	publishDir "${params.outdir}/fastp_result", mode: 'copy'

	input:
	path(read_file) from fastp_ch

	output:
	file "fastp*"	

	script:
	fastp_out = read_file.baseName
	"""
	fastp -i $read_file -o ${fastp_out} 

	"""

}

// Estimate genome size, heterozygosity and repeatitiveness

process run_jellyfish {
	tag "$read_file"
	publishDir "${params.outdir}/jellyfish_result", mode: 'copy'

	input:
	path read_file from jellyfish_ch


	output:
	file "*"

	
	script:
	jelly_out = read_file.baseName
	"""
	jellyfish count -C -m 21 -s 100M -t $params.threads $read_file -o ${jelly_out}.jf
	jellyfish histo -t 10 ${jelly_out}.jf > ${jelly_out}.histo
	"""
}


// Genome assembly


process run_genome_assembly {
	tag "$read_file"
	publishDir "${params.outdir}/assembly_result", mode: 'copy'
	
	input:
	path read_file from assembly_ch
	output:
	file "*"
	file('*hifi_out.p_ctg.gfa') into fasta_ch
 
	script:
	file_out = read_file.baseName
	if( params.assembler == 'flye' )
	   """
	    flye --pacbio-hifi $read_file --out-dir ${file_out}.flye_out --threads $params.threads
	   """
	else if( params.assembler == 'hifiasm' )
	   """
	    hifiasm -t $params.threads $read_file -o ${file_out}.hifi_out
	   """
	else
	   throw new IllegalArgumentException("Unknown genome_assembler $params.assembler") 

}

// Get primary contigs in fasta

process run_convertFasta {

	tag "$fa"
	publishDir "${params.outdir}/assembly_result/assembly_fasta", mode: 'copy'

	input:
	path fa from fasta_ch

	output:
	file "${fa.baseName}.fa" into blast, minimap2, blobtools2 

	script:
	"""
	awk '/^S/{print ">"\$2;print \$3}' ${fa} > "${fa.baseName}".fa
	"""
}
workflow.onComplete { 
	println ( workflow.success ? "\nDone! -> Thank you for using WORMpipe" : "error!" )
}

