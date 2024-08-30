#!/bin/bash

module load BEDTools/2.29.2-GCC-8.3.0

anchored_contigs(){
  cat $1 | awk '{print $2"\t"$9-100"\t"$9}' | awk '{if ($2<1) {print $1"\t"1"\t"$3} else {print $0}}' | bedtools sort -i - > upstream.bed
  cat $1 | awk '{print $2"\t"$10"\t"$10+100}' | bedtools sort -i - > downstream.bed

  awk '{print $1"\t"1"\t"2}' $2  > $2_1.bed
  cat $2 | sed 's/LN:i://' | awk '{print $1"\t"$2-1"\t"$2}' > $2_2.bed
  cat $2_1.bed $2_2.bed |  sort | uniq  > $2.bed
  bedtools intersect -a $2.bed -b upstream.bed | awk '{print $1}' | sort | uniq > up_$2.nams
  bedtools intersect -a $2.bed -b downstream.bed | awk '{print $1}' | sort | uniq > down_$2.nams

  cat  up_$2.nams  down_$2.nams | sort | uniq -c > anchor_counts_$1 

}

anchored_contigs $blast $repvals
