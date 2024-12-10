#!/bin/bash

module load BEDTools/2.29.2-GCC-8.3.0
module load  BLAT/3.5-GCC-11.2.0

#bin bed files
#==> bins_chr7_1.bed <==
#chr7_CG108_half_Mo17	158500000	159500000

#==> bins_chr7_2.bed <==
#chr7_CG108_half_Mo17	162000000	163000000

#==> bins_chr7_3.bed <==
#chr7_CG108_half_Mo17	173200000	174200000

#==> bins_chr8_1.bed <==
#chr8_CG108_half_Mo17	168000000	169000000

#==> bins_chr8_2.bed <==
#chr8_CG108_half_Mo17	180750000	181750000

#==> bins_chr8_3.bed <==
#chr8_CG108_half_Mo17	187000000	188000000



cd /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/CG108_half_Mo17

#knob monomers 
#all_monos_knob180.bed

bedtools sort -i  all_monos_knob180.bed > s_all_monos_knob180.bed

for c in chr7 chr8
do
for i in bins_"$c"*bed
do
bedtools intersect -a $i -b s_all_monos_knob180.bed > monos_"$i"
sed 's/_CG108_half_Mo17/_RagTag/' monos_"$i" > rename_monos_"$i"
bedtools getfasta -fi ragtag.scaffold.fasta -bed rename_monos_"$i" > monos_"$i".fasta
done
done

blat monos_bins_chr7_1.bed.fasta monos_bins_chr7_2.bed.fasta -t=dna -q=dna -minScore=80 -repMatch=2147483647 chr7_bin_12.blat
blat monos_bins_chr7_1.bed.fasta monos_bins_chr7_3.bed.fasta -t=dna -q=dna -minScore=80 -repMatch=2147483647 chr7_bin_13.blat
blat monos_bins_chr8_1.bed.fasta monos_bins_chr8_2.bed.fasta -t=dna -q=dna -minScore=80 -repMatch=2147483647 chr8_bin_12.blat
blat monos_bins_chr8_1.bed.fasta monos_bins_chr8_3.bed.fasta -t=dna -q=dna -minScore=80 -repMatch=2147483647 chr8_bin_13.blat
blat monos_bins_chr7_1.bed.fasta monos_bins_chr8_1.bed.fasta -t=dna -q=dna -minScore=80 -repMatch=2147483647 chr78_bin_1.blat

for i in chr*blat
do
awk '$1 ~ /^[0-9]*$/' $i| awk 'NR>1' | awk '{print $10"\t"$14"\t"($1/($11+$15-$1))}' | awk '{if($3>=.98) print $0}'> "$i".sub
done


