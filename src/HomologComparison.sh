#!/bin/bash

module load BEDTools/2.29.2-GCC-8.3.0
module load HMMER/3.3.2-gompi-2022a
module load BLAT/3.5-GCC-11.2.0
module load BLAST+/2.2.31
module load SAMtools/0.1.20-GCC-11.2.0

g=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/labels_AB10) #positions_wMo17)  #positions

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

##identifying non-repeat space and masking it
cat "$g"_*.out.bed | bedtools sort -i - > M_"$g"_all.out.bed
samtools faidx "$g"_arrays.fasta 
cat "$g"_arrays.fasta.fai | awk '{print $1"\t""1""\t"$2}' > "$g"_arrays.fasta.bed
bedtools subtract -a "$g"_arrays.fasta.bed  -b M_"$g"_all.out.bed > non_mono_"$g"_arrays.fasta.bed  
bedtools maskfasta -fi "$g"_arrays.fasta -bed non_mono_"$g"_arrays.fasta.bed  -fo "$g"_arrays_otherMASK.fasta
makeblastdb -in "$g"_arrays_otherMASK.fasta -out "$g"_arrays_otherMASK.db -dbtype 'nucl' -hash_index
blastn  -db "$g"_arrays_otherMASK.db -query  "$g"_arrays_otherMASK.fasta -outfmt "6 qseqid sseqid qlen slen length nident pident qstart qend sstart send" -num_threads 15 >  "$g"_arrays_otherMASK.blast
cat "$g"_arrays_otherMASK.blast | awk '{if($9>$8){print $1"_"$2"\t"$8"\t"$9"\t"$6"\t"$3"\t"$4}else {print $1"_"$2"\t"$9"\t"$8"\t"$6"\t"$3"\t"$4}}' | awk '{if($5>$6){print $1"\t"$2"\t"$3"\t"$4"\t"$6}else {print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' |  bedtools sort -i - | bedtools merge -i - -c 4,5 -o max,mean > "$g"_arrays_otherMASK.blast.sub_test2

##then just comparing the non-repeat space 
bedtools maskfasta -fi "$g"_arrays.fasta -bed M_"$g"_all.out.bed  -fo "$g"_arrays_monoMASK.fasta
makeblastdb -in "$g"_arrays_monoMASK.fasta -out "$g"_arrays_monoMASK.db -dbtype 'nucl' -hash_index
blastn  -db "$g"_arrays_monoMASK.db -query  "$g"_arrays_monoMASK.fasta -outfmt "6 qseqid sseqid qlen slen length nident pident qstart qend sstart send" -num_threads 15 >  "$g"_arrays_monoMASK.blast
cat "$g"_arrays_monoMASK.blast | awk '{if($9>$8){print $1"_"$2"\t"$8"\t"$9"\t"$6"\t"$3"\t"$4}else {print $1"_"$2"\t"$9"\t"$8"\t"$6"\t"$3"\t"$4}}' | awk '{if($5>$6){print $1"\t"$2"\t"$3"\t"$4"\t"$6}else {print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' |  bedtools sort -i - |bedtools merge -i - -c 4,5 -o max,mean > "$g"_arrays_monoMASK.blast.sub_test2



