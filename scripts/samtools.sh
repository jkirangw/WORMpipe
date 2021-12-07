#convert sam to bam
samtools view -S -b 314-2-Rc.raw.fa.aln.sam > 314-2-Rc.raw.fa.bam

#sort
samtools sort 314-2-Rc.raw.fa.bam -o 314-2-Rc.raw.fa.bam.sorted.bam

#index
samtools index 314-2-Rc.raw.fa.bam.sorted.bam

