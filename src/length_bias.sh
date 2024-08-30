#!/bin/bash
module load BEDTools/2.26.0-GCCcore-8.3.0

read_length_bias(){ #read_lengths, blast output

cat $1 | sort -n -k2 > sorted_$1
rows=$(wc -l sorted_$1| awk '{print $1}')
low_10=$((rows/10))

head sorted_$1 -n $low_10 | sort > shortest_10_$1
tail sorted_$1 -n $low_10 | sort > longest_10_$1

sort $2 -k 2 > sort_$2

join shortest_10_$1 sort_$2 -1 1 -2 2 > shortest_10_$1.blast
join longest_10_$1 sort_$2 -1 1 -2 2 > longest_10_$1.blast

for i in shortest_10_$1.blast longest_10_$1.blast
do
for rep in CentC Cent4 knob180 TR1
    do
    grep "$rep" $i | awk '{if($5>=30){print $1"\t"$10"\t"$11}}'| 
    awk '{if($3<$2){print $1"\t"$3"\t"$2}else{print $1"\t"$2"\t"$3}}' |
    bedtools sort -i -|bedtools merge -i - |awk '{print$0"\t"$3-$2}' |
    awk -v rep="$rep" '{Total=Total+$4} END{print rep"\t"Total}' > "$i"_sum_"$rep"
  done
done
}


