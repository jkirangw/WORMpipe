# assemble long reads
wtdbg2 -x rs -g 340m edge-min 2 --rescue-low-cov-edges --tidy-reads 5000 -i <(zcat 312-4_deduplicated.ccs.fastq.gz 312-4_2.fastq.gz 312-4_cell3.hifi_reads.fastq.gz 312-4_cell4.hifi_reads.fastq.gz) -t 32 -fo 314-2-Rc 2>wtdbg2.log

# derive consensus
wtpoa-cns -t 32 -i 314-2-Rc.ctg.lay.gz -fo 314-2-Rc.raw.fa

#polish consensus, not necessary if you want to polish the assemblies using other tools
minimap2 -t 32 -ax map-pb -r2k 314-2-Rc.raw.fa 312-4_deduplicated.ccs.fastq.gz 312-4_2.fastq.gz 312-4_cell3.hifi_reads.fastq.gz 312-4_cell4.hifi_reads.fastq.gz | samtools sort -@4 >314-2-Rc.bam
samtools view -F0x900 314-2-Rc.bam | wtpoa-cns -t 32 -d 314-2-Rc.raw.fa -i - -fo 314-2-Rc.cns.fa

# Addtional polishment using short reads
bwa index 314-2-Rc.cns.fa
bwa mem -t 32 314-2-Rc.cns.fa ERR218336_1.fastq.gz ERR218336_2.fastq.gz | samtools sort -O SAM | wtpoa-cns -t 32 -x sam-sr -d 314-2-Rc.cns.fa -i - -fo 314-2-Rc.srp.polished.fa
