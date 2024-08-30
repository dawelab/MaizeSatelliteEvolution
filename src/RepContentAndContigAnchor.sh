#!/bin/bash

module load bioawk/1.0-GCC-11.2.0
module load EMBOSS/6.6.0-foss-2021b
module load BLAST+/2.10.1-gompi-2022a
module load BEDTools/2.29.2-GCC-8.3.0
module load hifiasm/0.19.6-GCCcore-11.3.0

#rep_content function in src file
#rep_contigs function in src file

##HiFi raw read rep content
lin=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Data_Lines) #list with additional data line names
rep_content all_"$lin".fastq repeats.fasta HiFi

##HiFi assembly and rep content
hifiasm -u0 -o "$lin"_assem -l0 -t 20 all_"$lin".fastq
awk '/^S/{print ">"$2;print $3}' "$lin"_assem.bp.p_ctg.gfa > "$lin"_assem.bp.p_ctg.gfa.fasta

rep_content "$lin"_assem.bp.p_ctg.gfa.fasta repeats.fasta HiFi_Assem
rep_contigs "$lin"_assem.bp.p_ctg.gfa.fasta_repeats.blast "$lin"_assem.bp.p_ctg.noseq.gfa

cat "$lin"_assem.bp.p_ctg.noseq.gfa_Filt_repvals | sed "s/LN:i://" | sed "s/rd:i://">sub_"$lin"_assem.bp.p_ctg.noseq.gfa_Filt_repvals


##unanchored contigs
anchored_contigs "$lin"_assem.bp.p_ctg.gfa.fasta_repeats.blast "$lin"_assem.bp.p_ctg.noseq.gfa_Filt_repvals 
join -a 1 -1 1 -2 2 -e 0 -o 1.1,1.2,1.3,1.4,1.5,2.1  "$lin"_assem.bp.p_ctg.noseq.gfa_Filt_repvals  anchor_counts_"$lin"_assem.bp.p_ctg.gfa.fasta_repeats.blast> "$lin"_assem.bp.p_ctg.noseq.gfa_Filt_repvals_unanchoredEnds




##repeat content sum
for reads in HiFi HiFi_Assem 
do
if [[ $reads == HiFi ]]; then
  tot_len=$(cat all*"$lin"*_len_"$reads" | awk '{Total=Total+$2} END{print Total}')
    else
  tot_len=$(cat "$lin"*_len_"$reads" | awk '{Total=Total+$2} END{print Total}')
  fi
for rep in CentC Cent4 knob180 TR1
  do
  rep_len=$(cat *"$lin"*fasta_repeats.blast_sum_"$rep"_"$reads" |  awk  '{Total=Total+$2} END{print rep"\t"Total}')
  echo -e $reads"\t"$rep"\t"$((rep_len))"\t"$((tot_len))>> output.vals_"$lin" #vals in Mb

done
done
done<Data_Lines
