#!/bin/bash

module load BEDTools/2.29.2-GCC-8.3.0
module load BLAST+/2.2.31
module load SAMtools/0.1.20-GCC-11.2.0

g=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/labels_AB10) #positions_wMo17) 
##identifying non-repeat space and masking it
cat "$g"_*.out.bed | bedtools sort -i - > M_"$g"_all.out.bed
samtools faidx "$g"_arrays.fasta 
cat "$g"_arrays.fasta.fai | awk '{print $1"\t""1""\t"$2}' > "$g"_arrays.fasta.bed
bedtools subtract -a "$g"_arrays.fasta.bed  -b M_"$g"_all.out.bed > non_mono_"$g"_arrays.fasta.bed  
bedtools maskfasta -fi "$g"_arrays.fasta -bed non_mono_"$g"_arrays.fasta.bed  -fo "$g"_arrays_otherMASK.fasta
makeblastdb -in "$g"_arrays_otherMASK.fasta -out "$g"_arrays_otherMASK.db -dbtype 'nucl' -hash_index
blastn  -db "$g"_arrays_otherMASK.db -query  "$g"_arrays_otherMASK.fasta -outfmt "6 qseqid sseqid qlen slen length nident pident qstart qend sstart send" -num_threads 15 >  "$g"_arrays_otherMASK.blast

cat "$g"_arrays_otherMASK.blast | awk '{if($9>$8){print $1"_"$2"\t"$8"\t"$9"\t"$6"\t"$3"\t"$4}else {print $1"_"$2"\t"$9"\t"$8"\t"$6"\t"$3"\t"$4}}' | \
awk '{if($5>$6){print $1"\t"$2"\t"$3"\t"$4"\t"$6}else {print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' |  bedtools sort -i - | \
bedtools merge -i - -c 4,5 -o max,mean > "$g"_arrays_otherMASK.blast.sub_test2



