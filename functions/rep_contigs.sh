#!/bin/bash

module load BEDTools/2.29.2-GCC-8.3.0

#total reps by contig/chrom
rep_contigs(){ #blast and noseq contig file
  for rep in CentC Cent4 knob180 TR1
    do
    grep "$rep" $1 | awk '{if($4>=30){print$0}}'| 
    cut -f2,9,10| awk '{if($3<$2){print $1"\t"$3"\t"$2}else{print $1"\t"$2"\t"$3}}' |
    bedtools sort -i -|bedtools merge -i - |awk '{print$0"\t"$3-$2}' | 
    awk '{a[$1]+=$4;}END{for (i in a)print i"\t"a[i];}' |
    sort | awk -v r="$rep" '{print $0"\t"r}'> conttots_$1_"$rep".bed 
  done
  cat conttots_$1_*.bed   | sort > all_conttots_$1.bed
  awk '{if($1=="S") print $2"\t"$4"\t"$5}' $2  | sort > $2_Filt
  join  -e '0' -j 1 -a 1  -o 1.1,1.2,1.3,2.2,2.3 $2_Filt  all_conttots_$1.bed > $2_Filt_repvals 
}



