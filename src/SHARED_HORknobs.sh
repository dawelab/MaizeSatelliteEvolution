#!/bin/bash
module load BEDTools/2.29.2-GCC-8.3.0
module load seqtk/1.3-GCC-11.2.0
module load Clustal-Omega/1.2.4-GCC-11.2.0
module load EMBOSS/6.6.0-foss-2021b


dir="CG108_half_Mo17"

cat CG108_string | awk '{if($7== "knob180") print $0}' | awk '{if($12>0) print $4"\t"$5}' > HOR_knob180

for clust in .9 .91 .92 .93 .94 .95 .96 .97 .98 .99
do
cat HOR_knob180 | awk -v c="$clust" '{if( $2 == c ) print $1}' > "$dir"/bins_interest

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
cat   "$dir"/letter_mon_clust_CONS.fasta >> HOR_clust/"$dir"_HORknobs_"$clust".fasta
done< "$dir"/let #letter loop

done< "$dir"/bins_interest #bin loop


module load R/4.1.2-foss-2021b
Rscript --vanilla /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/HOR_newpatt_HORknobs.R /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/HOR_newpatt
Rscript --vanilla /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/HOR_newpatt_repatt_HORknobs.R /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/HOR_newpatt


for i in HORknobs_*SHARED_ALL_HOR.bed
do
group=$(echo $i| awk -F"_" '{print $1}')
level=$(echo $i| awk -F"_" '{print $2}')
cat $i| awk -v var="$level" '{print $0"\t"var}'| awk -v var="$group" '{print $0"\t"var}' >> HORknobs_ALL_SHARED_ALL_HOR.bed
done

Rscript --vanilla /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/Shared_purity_HORknobs.R /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/HOR_newpatt

