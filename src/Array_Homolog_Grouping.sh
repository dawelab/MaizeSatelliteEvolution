#!/bin/bash

module load seqtk/1.3-GCC-11.2.0
module load BEDTools/2.29.2-GCC-8.3.0

#grouping all the genomes into one
while read lin
do
cat /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/"$lin"/ragtag.scaffold.fasta | seqtk subseq - chr.list |  sed "s/RagTag/$lin/"  > "$lin"_chr_ragtag.scaffold_Mo17.fasta
done<Mo17_dirs

cat *_chr_ragtag.scaffold_Mo17.fasta > ALL_chr_ragtag.scaffold_Mo17.fasta

mkdir conserved_positions
sed -i 1d grouped_arrays_pos_member
cat grouped_arrays_pos_member |awk '{print $12}' | sort | uniq > labels

#get all the monomers in conserved positions into single fastas
while read g
do
cat grouped_arrays_pos_member | awk -v gen="$g" '{if($12 == gen) print $4"_"$3"\t"$5"\t"$6}' | bedtools sort -i - | uniq > conserved_positions/"$g"_arrays.bed
bedtools getfasta -fi ALL_chr_ragtag.scaffold_Mo17.fasta -bed conserved_positions/"$g"_arrays.bed > conserved_positions/"$g"_arrays.fasta 
done<labels

#AB10 too!
cd /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/
lin="AB10_HALF_AB10"
#sed -i 1d Ab10_arrays_tail_grouped
#cat Ab10_arrays_tail_grouped | awk '{print $12}' | sort | uniq > labels_AB10

while read g
do
cat Ab10_arrays_tail_grouped | awk -v gen="$g" '{if($12 == gen) print $4"\t"$5"\t"$6}' | bedtools sort -i - | uniq > conserved_positions/"$g"_arrays.bed
bedtools getfasta -fi /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/"$lin"/ragtag.scaffold.fasta -bed conserved_positions/"$g"_arrays.bed > conserved_positions/"$g"_arrays.fasta 
done<labels_AB10

####
#dot plot-- each monomer similarity to consensus for rough pattern plotting
makehmmerdb "$g"_arrays.fasta "$g"_arrays.fasta.db 

for element in CentC Cent4 knob180 TR1
do
##finding monomers, comparing to reference mono
nhmmer -o "$g"_"$element".out  /scratch/rdp22327/Dawe/consensus/"$element"_consensus_fix_RC.fasta.afa.sto.hmm  "$g"_arrays.fasta.db   
cat  "$g"_"$element".out  | awk '{print $4"\t"$5"\t"$6}'| grep "chr" |  awk '{if($3>$2) {print $1"\t"$2"\t"$3} else {print $1"\t"$3"\t"$2}}'| bedtools sort -i - > "$g"_"$element".out.bed
bedtools getfasta -fi "$g"_arrays.fasta -bed "$g"_"$element".out.bed > mono_"$g"_"$element".fasta
blat  mono_"$g"_"$element".fasta /scratch/rdp22327/Dawe/consensus/"$element"_main_consensus.fasta -t=dna -q=dna -maxGap=10 -minScore=0 -repMatch=2147483647  consensus_mono_"$g"_"$element".fasta.fa.blat
awk '$1 ~ /^[0-9]*$/' consensus_mono_"$g"_"$element".fasta.fa.blat | awk 'NR>1' | awk '{print $14"\t"($1/(($11/2)+$15-$1))}' > consensus_mono_"$g"_"$element".fasta.fa.blat.sub
awk -v var="$element" '{print $0"\t"var}' consensus_mono_"$g"_"$element".fasta.fa.blat.sub >nam_consensus_mono_"$g"_"$element".fasta.fa.blat.sub
done
cat nam_consensus_mono_"$g"_*.fasta.fa.blat.sub > all_nam_consensus_mono_"$g".fasta.fa.blat.sub

