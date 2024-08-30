#!/bin/bash

#dir is a list of all the directories -- one genetic line per directory
##collecting all the smoothing classification bins with prepped columns
while read i
do
 for j in "$i"/*/*_fin_bins_combined.txt
 do
  grep "chr" $j > "$j"_sub
  cat "$j"_sub  | awk '{print $1}' | sed 's/.*chr/chr/' | sed 's/_.*//' > chr
  cat "$j"_sub  | awk '{print $1}' | sed 's/.*://' |  sed 's/-.*//'> array_start
  cat "$j"_sub  | awk '{print $1}' |  sed 's/.*-//' | sed 's/_comp.*//' | sed 's/.*_//'> rep
  cat "$j"_sub   | paste chr array_start rep - > all
  cat all | awk -v l="$i" '{print $1"_"l"\t"$2+$17"\t"$2+$18"\t"$3"\t"$5"\t"$6"\t"$20"\t"$1"\t"l}'>> all_Mo17_scaff_structure_bins
 done
done< $dir

#######
module load BEDTools/2.29.2-GCC-8.3.0
module load BLAT/3.5-GCC-11.2.0

lin=$(echo $dir | sed 's/_.*//' )
cd $dir

cat all_Mo17_scaff_structure_bins | awk -v l="$dir" '{if($9 == l) print $0}' | awk '{if($7 == "HOR") print $0}' > "$dir"/"$dir"_HOR_bins.bed

#again, extracting monomers
for i in knob180 Cent4 CentC TR1
do
cat filt_"$i"_hmmoutF.out.bed  | awk '{print $1}' | sed 's/.*chr/chr/' | sed 's/_.*//' > chr
cat filt_"$i"_hmmoutF.out.bed | awk '{print $1}' | sed 's/.*://' |  sed 's/-.*//'> array_start
paste chr array_start  filt_"$i"_hmmoutF.out.bed > info_filt_"$i"_hmmoutF.out.bed
cat info_filt_"$i"_hmmoutF.out.bed | awk -v rep="$i"  '{print $1"\t"$2+$4"\t"$2+$5"\t"rep}'  | awk -v d="$dir" '{if($3-$2>= 120) print $1"_"d"\t"$2"\t"$3"\t"$4}' > all_monos_"$i".bed
done

cat all_monos_*.bed | bedtools sort -i - | uniq | grep "chr"> mono_ALL_positions.bed
#for AB10_HALF_AB10, which was separated
#cat  /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/AB10_HALF_AB10/ragtag.scaffold.fasta  | sed 's/_RagTag/_AB10_HALF_AB10/'| sed 's/_corrected//'> /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/AB10_HALF_AB10/rename_ragtag.scaffold.fasta

#for each bin with an HOR classification, re-compare them. Previous data was not saved for space
while read bucket
do
   echo $bucket | awk '{print $1"\t"$2"\t"$3"\t"$4 }'> bucket.bed
   rep_interest=$(cat bucket.bed| awk '{print $4}')
   chr=$(cat bucket.bed| awk '{print $1}')
   ar_start=$(cat bucket.bed| awk '{print $2}')
   descript=$(echo $chr"_""$ar_start""_""$rep_interest")
   
   ##grab the monomers that are within region of interest, subset to the ones that are at least 120 bp (roughly 80% of the consensus length)#
       bedtools intersect -a  bucket.bed  -b mono_ALL_positions.bed | sort | uniq > bucket_mono.bed
       awk '{print $0"\t"$3-$2}' bucket_mono.bed| awk  -v rep="$rep_interest"   '{if($4 == rep) print $0}' | awk '{if($5 >= 120) print $0}' > sub_bucket_mono.bed

        #turn them into fasta
        bedtools getfasta -fi /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/ALL_chr_ragtag.scaffold_Mo17.fasta -bed sub_bucket_mono.bed > sub_bucket_mono.fa
#for AB10 on AB10
#       bedtools getfasta -fi /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/"$dir"/rename_ragtag.scaffold.fasta -bed sub_bucket_mono.bed > sub_bucket_mono.fa 
       blat sub_bucket_mono.fa sub_bucket_mono.fa -t=dna -q=dna -maxGap=10 -minScore=0 -repMatch=2147483647 HOR_"$descript".blat
       grep '^[0-9]' HOR_"$descript".blat > sub_HOR_"$descript".blat
done< "$dir"_HOR_bins.bed

module load R/3.6.2-foss-2019b
Rscript --vanilla /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/network_summary_HOR.R  "$dir"_HOR_bins.bed /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/"$dir"
Rscript --vanilla /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/HOR_pattern_3.R  "$dir"_HOR_bins.bed /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/"$dir"

