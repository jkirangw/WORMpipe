
#reads in FASTA/FASTQ format and overlaps/alignments between the reads and the contigs in MHAP/PAF/SAM format.
pri_asm=314-2-Rc.raw.fa 
#minimap2 -xmap-pb $pri_asm <(zcat 312-4_deduplicated.ccs.fastq.gz 312-4_2.fastq.gz 312-4_cell3.hifi_reads.fastq.gz 312-4_cell4.hifi_reads.fastq.gz) | gzip -c - > $i.al.gz

minimap2 -ax map-pb 314-2-Rc.raw.fa  <(zcat 312-4_deduplicated.ccs.fastq.gz 312-4_2.fastq.gz 312-4_cell3.hifi_reads.fastq.gz 312-4_cell4.hifi_reads.fastq.gz) | gzip -c - > 314-2-Rc.raw.fa.aln.sam.gz


 hapog \
  --genome 312-4-Rc.raw.fa \   	          # Fasta file of the genome to polish
  -b  312-4-Rc.raw.fa.bam.sorted.bam \    # Sorted BAM file
  --output 312-4_hypoG-polishing \        # Output directory
  -t 32 \                                 # Number of threads to use
  -u                                      # Include unpolished sequences in the output
