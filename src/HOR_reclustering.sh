#!/bin/bash

module load BEDTools/2.29.2-GCC-8.3.0
module load seqtk/1.3-GCC-11.2.0
module load Clustal-Omega/1.2.4-GCC-11.2.0
module load EMBOSS/6.6.0-foss-2021b

#all the array positions
cat grouped_arrays_pos_member Ab10_arrays_tail_grouped > all_groups
cat  all_groups | awk '{print $12}' | sort | uniq > distinct_groups

#loop to pull all the HOR bins with a conserved array positions, extract all those monos
while read dir #group name
do
cat all_groups | grep $dir | awk '{print $4"\t"$5"\t"$6"\t"$12}' | bedtools sort -i - | sed 's/_RagTag//' > "$dir"/grouped_bins.bed
cat "$dir"/HOR_re-eval_clusters.csv | awk -F"," '{if($6 != "NA") print $3"\t"$4"\t"$5"\t"$6"\t"$2}'  > "$dir"/HOR_re-eval_clusters.csv.bed
bedtools intersect -a "$dir"/grouped_bins.bed -b "$dir"/HOR_re-eval_clusters.csv.bed -wa -wb | \
awk '{print $5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$4}' > "$dir"/HOR_re-eval_clusters.csv_intersect.bed #HOR bins with labels and group names

if [ $dir = "AB10_HALF_AB10" ]; then
bedtools getfasta -fi /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/"$dir"/rename_ragtag.scaffold.fasta -bed  "$dir"/mono_ALL_positions.bed > "$dir"/mono_ALL_positions.bed.fasta
else
bedtools getfasta -fi /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/ALL_chr_ragtag.scaffold_Mo17.fasta -bed  "$dir"/mono_ALL_positions.bed > "$dir"/mono_ALL_positions.bed.fasta
fi

#then, loop through the optimal clustering thresholds
while read lab #group info
do
for clust in .9 .91 .92 .93 .94 .95 .96 .97 .98 .99
do
cat "$dir"/HOR_re-eval_clusters.csv_intersect.bed | awk -v l="$lab"  '{if( $6 == l ) print $0}' | awk -v c="$clust" '{if( $4 == c ) print $5}'| awk -v d="$dir" '{print $1}' > "$dir"/bins_interest

#for each distinct monomer subtype (i.e. all the A's from bin 1 ) -- group monomers of the same subtype and get a consensus sequence
if [[ $(wc -l < "$dir"/bins_interest) -ge 0 ]]; then
while read bin
do
echo $bin> "$dir"/bin_info
grep -Fw -f  "$dir"/bin_info "$dir"/HOR_bed.csv | awk -F"," '{print $3}'| tr --delete '\n'|  sed 's/./& /g' | tr ' ' '\n' | sort | uniq >  "$dir"/let
cat "$dir"/"$bin"*out >  "$dir"/bin_group

while read letter_mon
do
echo $letter_mon >  "$dir"/letter_mon_value
grep -Fw -f  "$dir"/letter_mon_value  "$dir"/bin_group | awk '{print $1}' >  "$dir"/nams_bin_group
seqtk subseq "$dir"/mono_ALL_positions.bed.fasta  "$dir"/nams_bin_group >  "$dir"/letter_mon.fasta
clustalo -i  "$dir"/letter_mon.fasta -o  "$dir"/letter_mon_clust.fasta  --force
cons -sequence  "$dir"/letter_mon_clust.fasta -name "$dir":"$bin":"$letter_mon" -outseq  "$dir"/letter_mon_clust_CONS.fasta
cat   "$dir"/letter_mon_clust_CONS.fasta >> HOR_clust/"$dir"_"$lab"_"$clust".fasta
done< "$dir"/let #letter

done< "$dir"/bins_interest #bin
fi

done #clust

done<distinct_groups


