#!/bin/bash

dir=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Mo17_dir) #directories with Mo17 scaffolding

lin=$(echo $dir | sed 's/_.*//' )
cat "$dir"/ragtag.scaffold.agp | awk '{if($5 == "U") print $1"\t"$2"\t"$3}' > "$dir"/100n.bed
cat "$dir"/ragtag.scaffold.agp | awk -v var="$dir"  '{if($5 == "U") print var"\t"$1"\t"$2"\t"$3}' > "$dir"/labeled_100n.bed

#defining arrays
for rep in CentC Cent4 knob180 TR1
do
grep "$rep" "$dir"/ragtag.scaffold.fasta_repeats.blast | awk '{if($4>=30){print$0}}'|cut -f2,9,10| awk '{if($3<$2){print $1"\t"$3"\t"$2}else{print $1"\t"$2"\t"$3}}' |bedtools sort -i -| bedtools merge -i - -c 1 -o count -d 10000 > "$dir"/ragtag.scaffold.fasta_repeats.blast_"$rep".bed 
awk '{if($4>=10) print $0"\t"$3-$2}' "$dir"/ragtag.scaffold.fasta_repeats.blast_"$rep".bed   | awk  -v rep="$rep"  '{if($5/($3-$2)>.1) print $1"\t"$2"\t"$3"\t"rep"\t"$3-$2"\t"$4}' > "$dir"/ragtag.scaffold.fasta_repeats.blast_"$rep"_arrays.bed 
##monomers merged within 10kb, array must have at least 10% rep content and 10 monomers
cat "$dir"/ragtag.scaffold.fasta_repeats.blast_"$rep".bed  | grep "chr" | awk '{print $0"\t"$3-$2}' | awk -v rep="$rep" '{Total=Total+$5} END{print rep"\t"Total}' >"$dir"/ragtag.scaffold.fasta_repeats.blast_sum_"$rep"_Chr_Assem
bedtools intersect -c -a "$dir"/ragtag.scaffold.fasta_repeats.blast_"$rep"_arrays.bed -b "$dir"/100n.bed > "$dir"/ragtag.scaffold.fasta_repeats.blast_"$rep"_arrays_gaps.bed 
grep "$rep" "$dir"/ragtag.scaffold.fasta_repeats.blast | awk '{if($4>=30){print$0}}'|cut -f2,9,10| awk '{if($3<$2){print $1"\t"$3"\t"$2}else{print $1"\t"$2"\t"$3}}'  |awk '{print $0"\t"$3-$2}' | awk -v rep="$rep" '{Total=Total+$4} END{print rep"\t"Total}' >"$dir"/ragtag.scaffold.fasta_repeats.blast_sum_"$rep"_Assem
grep "$rep" "$dir"/ragtag.scaffold.fasta_repeats.blast | awk '{if($4>=30){print$0}}'|cut -f2,9,10| awk '{if($3<$2){print $1"\t"$3"\t"$2}else{print $1"\t"$2"\t"$3}}'  | grep "chr" | awk '{print $0"\t"$3-$2}' | awk -v rep="$rep" '{Total=Total+$4} END{print rep"\t"Total}' >"$dir"/ragtag.scaffold.fasta_repeats.blast_sum_"$rep"_Chr_Assem
done
cat "$dir"/ragtag.scaffold.fasta_repeats.blast_*_arrays_gaps.bed | awk '{print $1"\t"$2"\t"$3"\t"$3-$2"\t"$5"\t"$7"\t"$4}' |grep "chr" > "$dir"/"$dir"_all_arrays.bed
cat  "$dir"/"$dir"_all_arrays.bed |bedtools sort -i - | uniq >  "$dir"/sort_"$dir"_all_arrays.bed

#conserved gene map
join -1 4 -2 1 "$dir"/"$dir".gff.bed Mo17_core_nams | awk '{print $2"\t"$3"\t"$4"\t"$1}'> "$dir"/sub_"$dir".gff.bed
#Mo17_core_nams from https://github.com/dawelab/ETR_in_HiFi/blob/0590c0c53365594aa8d1373276802e8e71c8d9ee/analysis/core_genes.sh#L45

cat  *half_Mo17/sub_*.gff.bed | awk '{print $4}' | sort | uniq -c | sort -k1 | awk '{if($1==12) print $2}' > gene_counts # genes that mapped well
join -1 4 -2 1 "$dir"/sub_"$dir".gff.bed gene_counts  | awk '{print $2"\t"$3"\t"$4"\t"$1}' >"$dir"/CORE_sub_"$dir".gff.bed  
bedtools sort -i "$dir"/CORE_sub_"$dir".gff.bed |bedtools merge -i - -c 4 -o distinct > "$dir"/sort_CORE_sub_"$dir".gff.bed


#combining gene map with arrays
bedtools closest -id -D a -a "$dir"/sort_"$dir"_all_arrays.bed -b "$dir"/sort_CORE_sub_"$dir".gff.bed > "$dir"/"$dir"_all_arrays.bed_up
bedtools closest -iu -D a -a "$dir"/sort_"$dir"_all_arrays.bed -b "$dir"/sort_CORE_sub_"$dir".gff.bed > "$dir"/"$dir"_all_arrays.bed_down
cat "$dir"/"$dir"_all_arrays.bed_up | sort | uniq > "$dir"/"$dir"_all_arrays.bed_up_uni
cat "$dir"/"$dir"_all_arrays.bed_down | sort | uniq > "$dir"/"$dir"_all_arrays.bed_down_uni

paste "$dir"/"$dir"_all_arrays.bed_up_uni "$dir"/"$dir"_all_arrays.bed_down_uni| awk  -v dir="$dir"  '{print dir"\t"$1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$10"\t"$11"\t"$21"\t"$23}' > "$dir"/array_gen_cords

#cat *_HALF_AB10/array_gen_cords *_out_HALF/array_gen_cords >  B73_scaff_array_gen_cords
#cat *Mo17/array_gen_cords CG119_fix/array_gen_cords > Mo17_scaff_array_gen_cords
# cat */labeled_100n.bed | grep "chr" > all_100_gaps

