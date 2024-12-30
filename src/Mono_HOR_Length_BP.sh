#!/bin/bash

cd /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/CG108_half_Mo17

module load bioawk/1.0-GCC-11.2.0
module load BEDTools/2.29.2-GCC-8.3.0

for i in mono_*_ragtag.scaffold.fasta_arrays.fasta 
do
cat $i | bioawk -c fastx '{ print $name, length($seq) }' | awk '{print $2}' > "$i"_lens
done


for i in AB10_HALF_AB10 CG119_half_Mo17 AB10_half_Mo17 B73_half_Mo17 K64_half_Mo17 Tx777_half_Mo17 Tx779_half_Mo17 CG44_half_Mo17 CML442_half_Mo17 TIL11.2cell.HiFi_half_Mo17 Zea-mays-ssp-mexicana-TIL25_4cell-hifi_half_Mo17 TIL01.3cell.HiFi_half_Mo17 Mo17_half_Mo17 CG108_half_Mo17
do

cat "$i"/mono_ALL_positions.bed | awk '{print $1}' | sed "s/_.*//" > chr_nam
paste chr_nam "$i"/mono_ALL_positions.bed  | awk '{print $1"\t"$3"\t"$4"\t"$5}' | bedtools sort -i - > sort_mono_ALL_positions.bed

cat "$i"/HOR_re-eval_clusters.csv_intersect.bed | awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2}'  > coord_HOR_re-eval_clusters.csv_intersect.bed

bedtools intersect -a coord_HOR_re-eval_clusters.csv_intersect.bed -b  sort_mono_ALL_positions.bed -wa -wb | awk '{print $5"\t"$6"\t"$7"\t"$4"\t"$8"\t"$7-$6}' | sort -k4 > mono_info_lens.bed

cat  "$i"/string_out.csv | awk -F"," '{print $3}' | sort > binns
 
join -1 4 -2 1 mono_info_lens.bed binns > "$i"_sub_mono_info_lens.bed

done


cat *_sub_mono_info_lens.bed > mono_HOR_coords

