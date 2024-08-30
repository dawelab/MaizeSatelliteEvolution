#!/bin/bash

rep_content(){ #repeat_content reads.fasta library.fasta label
  if [[ $1 == *fastq ]]; then
    seqret -sequence  $1 -outseq $1.fasta
    bioawk -c fastx '{ print $name, length($seq) }' $1.fasta> $1.fasta_len_$3
    blastn -subject $1.fasta -query $2 -outfmt 6 -num_threads 10 -max_target_seqs 5000000 > $1_repeats.blast
    else
    blastn -subject $1 -query $2 -outfmt 6 -num_threads 10 -max_target_seqs 5000000 > $1_repeats.blast
    bioawk -c fastx '{ print $name, length($seq) }' $1 > $1_len_$3
  fi

  for rep in CentC Cent4 knob180 TR1
    do
    grep "$rep" $1_repeats.blast | awk '{if($4>=30){print$0}}'| 
    cut -f2,9,10| awk '{if($3<$2){print $1"\t"$3"\t"$2}else{print $1"\t"$2"\t"$3}}' |
    bedtools sort -i -|bedtools merge -i - |awk '{print$0"\t"$3-$2}' |
    awk -v rep="$rep" '{Total=Total+$4} END{print rep"\t"Total}' > $1_repeats.blast_sum_"$rep"_$3
  done
}


##to run:
#module load bioawk/1.0-GCC-11.2.0
#module load EMBOSS/6.6.0-foss-2021b
#module load BLAST+/2.10.1-gompi-2022a
#module load BEDTools/2.29.2-GCC-8.3.0

#source rep_content.sh
#rep_content $file_of_interest.fasta $reference_repeats $label
