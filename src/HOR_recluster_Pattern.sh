#!/bin/bash

module load BLAT/3.5-GCC-11.2.0
module load SeqKit/0.16.1


mkdir HOR_newpatt

while read lab #group info
do
for clust in .9 .91 .92 .93 .94 .95 .96 .97 .98 .99 #threshold
do

#group all those consensus monomers and compare with blat
cat HOR_clust/*_"$lab"_"$clust".fasta > HOR_newpatt/all_"$lab"_"$clust".fasta
awk '/^>/{f=!d[$1];d[$1]=1}f' HOR_newpatt/all_"$lab"_"$clust".fasta > HOR_newpatt/all_fix_"$lab"_"$clust".fasta
blat HOR_newpatt/all_fix_"$lab"_"$clust".fasta HOR_newpatt/all_fix_"$lab"_"$clust".fasta -t=dna -q=dna -maxGap=10 -minScore=0 -repMatch=2147483647 HOR_newpatt/HOR_all_"$lab"_"$clust".fasta.blat
grep '^[0-9]' HOR_newpatt/HOR_all_"$lab"_"$clust".fasta.blat > HOR_newpatt/sub_HOR_all_"$lab"_"$clust".fasta.blat

done #clust
done<labels_AB10 labels

#send those similarity scores to R scripts that will 
#1. convert to matrix, then network, then label monomers by cluster group (number here)
# Then 2. turn those numbers into letters. Make a character string descirbing each bins. And finally, identify to kmers of interest
module load R/4.1.2-foss-2021b
Rscript --vanilla /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/HOR_newpatt.R /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/HOR_newpatt
Rscript --vanilla /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/HOR_newpatt_repatt.R /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/HOR_newpatt



###################
##grouping shared patterns in/with Mo17 
while read lab #group info
do
for v in 9 91 92 93 94 95 96 97 98 99
do
cat HOR_newpatt/"$lab"_"$v"_SHARED_FILT_HOR.bed | grep "Mo17_half_Mo17" | awk -v a="$v" '{print $0"\t"a}' >> HOR_newpatt/lab2_Mo17_half_Mo17_"$lab"_shared_filt.bed
done
echo $lab
cat HOR_newpatt/lab2_Mo17_half_Mo17_"$lab"_shared_filt.bed | awk '{print $3"\t"$2"_"$6}' | sort | uniq | awk '{print $1}' | uniq -c| awk '{if($1>1) print $0}'

done<labels
