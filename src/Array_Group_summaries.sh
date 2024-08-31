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

#####

while read i
do
 for j in "$i"/*/*_fin_bins_combined.txt
 do
  grep "chr" $j  > "$j"_sub
  cat "$j"_sub  | awk '{print $1}' | sed 's/.*chr/chr/' | sed 's/_.*//' > chr
  cat "$j"_sub  | awk '{print $1}' | sed 's/.*://' |  sed 's/-.*//'> array_start
  cat "$j"_sub  | awk '{print $1}' |  sed 's/.*-//' | sed 's/_comp.*//' | sed 's/.*_//'> rep
  cat "$j"_sub   | paste chr array_start rep - > all
  cat all | awk -v l="$i" '{print $1"_"l"\t"$2+$17-10000"\t"$2+$18-10000"\t"$3"\t"$5"\t"$6"\t"$20"\t"$1"\t"l}'>> all_Mo17_scaff_structure_bins_2
 done
done<Mo17_dirs_plus
cat  Ab10_arrays_tail_grouped grouped_arrays_pos_member | awk '{print $4"_"$3"\t"$5"\t"$6"\t"$7"\t"$12}' > grouped_arrays_pos_member_ab10_sub


module load BEDTools/2.29.2-GCC-8.3.0
bedtools sort -i all_Mo17_scaff_structure_bins_2 > all_Mo17_scaff_structure_bins.bed
bedtools sort -i grouped_arrays_pos_member_ab10_sub > grouped_arrays_pos_member_ab10_sub.bed
 
bedtools intersect  -nonamecheck  -a grouped_arrays_pos_member_ab10_sub.bed -b all_Mo17_scaff_structure_bins.bed -wa -wb > all_structure_conserved_wMo17_v2
cat all_structure_conserved_wMo17_v2 | awk '{if($4 == $9)print $0}'| awk '{print $14"\t"$4"\t"$13"\t"$2"\t"$7"\t"$8"\t"$5"\t"$10"\t"$11"\t"$12}' > all_structure_conserved_wMo17_v2_mod

