#!/bin/bash

module load bioawk/1.0-foss-2019b 
module load  BEDTools/2.30.0-GCC-10.2.0
module load BLAST+/2.2.31


#for K64, individual flow cells: m64041_230111_122536.hifi_reads.fastq, m64041_230112_215536.hifi_reads.fastq

#for CG108, individual flow cells: m64310e_220116_060322.hifi_reads.bam.fastq, m64310e_220201_005109.hifi_reads.bam.fastq

#for Ab10,  individual flow cells: 2_B01.fasta,  3_C01.fasta

read_len_density(){
bioawk -c fastx '{ print $name, length($seq) }' $1 > $1_lens
tot_len=$(cat $1_lens | awk '{Total=Total+$2} END{print Total}')
blastn -subject $1 -query $2 -outfmt 6 -num_threads 10 -max_target_seqs 5000000 > $1_repeats.blast
    
for rep in CentC Cent4 knob180 TR1
    do
    grep "$rep" $1_repeats.blast | awk '{if($4>=30){print$0}}'| 
    cut -f2,9,10| awk '{if($3<$2){print $1"\t"$3"\t"$2}else{print $1"\t"$2"\t"$3}}' |
    bedtools sort -i -|bedtools merge -i - |awk '{print$0"\t"$3-$2}' |
    rep_len=$(awk -v rep="$rep" '{Total=Total+$4} END{print rep"\t"Total}')
     echo -e $4"\t"$rep"\t"$((rep_len))"\t"$((tot_len))>> output.vals_$3_sub
  done
}
}

read_len_density m64041_230111_122536.hifi_reads.fastq repeats.fasta K64 K64_1
read_len_density m64041_230112_215536.hifi_reads.fastq repeats.fasta K64 K64_2

read_len_density m64310e_220116_060322.hifi_reads.bam.fastq repeats.fasta CG108 CG108_1
read_len_density m64310e_220201_005109.hifi_reads.bam.fastq repeats.fasta CG108 CG108_2

read_len_density 2_B01.fasta repeats.fasta Ab10 Ab10_1
read_len_density 3_C01.fasta repeats.fasta Ab10 Ab10_2


