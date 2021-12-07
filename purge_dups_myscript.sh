#Add  purge_dups scripts in your PATH environment
export PATH=$PATH:/home/jkirangw/Scratch/jkirangw/wtgb2-RC/wtdg2/purge_dups/scripts
export PATH=$PATH:/home/jkirangw/Scratch/jkirangw/wtgb2-RC/wtdg2/purge_dups/bin
export PATH=$PATH:/home/jkirangw/Scratch/jkirangw/wtgb2-RC/wtdg2/purge_dups/src
export PATH=$PATH:/home/jkirangw/Scratch/jkirangw/wtgb2-RC/wtdg2/purge_dups

#Run minimap to align  pacbio  data and generate paf files and calculate read depth

pri_asm=312-4-hapog.fasta
for i in $(cat pbfofn)
do
	echo $pri_asm
	minimap2 -xmap-pb $pri_asm $i | gzip -c - > $i.paf.gz
done
#pbcstat *.paf.gz #(produces PB.base.cov and PB.stat files)
calcuts PB.stat > cutoffs 2>calcults.lo

#split an assembly  and do self-self alignment.

split_fa $pri_asm > $pri_asm.split

minimap2 -xasm5 -DP $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz 

#Purge haplotigs and overlaps

purge_dups -2 -T cutoffs -c PB.base.cov 312-4-hapog.fasta.split.self.paf.gz  > dups.bed 2> purge_dups.log

#Get purged primary and haplotig sequences from draft assembly

get_seqs -e dups.bed $pri_asm 
