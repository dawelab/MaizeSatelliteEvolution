#!/bin/bash
lin=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Data_Lines)  #all 13 lines

module load BEDTools/2.29.2-GCC-8.3.0
module load HMMER/3.1b2-foss-2022a
module load BLAT/3.7-GCC-11.3.0
module load EMBOSS/6.6.0-foss-2021b
module load seqtk/1.3-GCC-11.2.0
module load SAMtools/0.1.20-GCC-11.2.0
module load R/4.1.2-foss-2021b

#finding all the reads that contain repeat content-- save some resources. blast file is from issues*assmeblies scripts
cat all_"$lin".fastq.fasta_repeats.blast| awk '{if($10>$9) {print $2"\t"$9"\t"$10"\t"$10-$9"\t"$1} else{print $2"\t"$10"\t"$9"\t"$9-$10"\t"$1}}' | bedtools sort -i - | bedtools merge -i - -d 100000 -c 4,4,5 -o sum,count,distinct | sort  | awk '{print $1}' | uniq >m_count_all_"$lin".fastq.fasta_repeats.blast
seqtk subseq all_"$lin".fastq.fasta m_count_all_"$lin".fastq.fasta_repeats.blast > rep_all_"$lin".fastq.fasta
reads=rep_all_"$lin".fastq.fasta


#find monomers and compare all the consensus
makehmmerdb "$reads" "$reads".db
for i in CentC TR1 knob180 Cent4
do
f_hmm_consensus=/scratch/rdp22327/Dawe/consensus/"$i"_consensus_fix.fasta.afa.sto.hmm
nhmmer --tblout hmmoutF.tbl -o "$i"_"$reads"_hmmoutF.out  $f_hmm_consensus  "$reads"  
cat  "$i"_"$reads"_hmmoutF.out  | awk '{print $4"\t"$5"\t"$6}'| grep "ccs" |  awk '{if($3>$2) {print $1"\t"$2"\t"$3} else {print $1"\t"$3"\t"$2}}'| bedtools sort -i - > filt_"$i"_"$reads"_sub_hmmoutF.out.bed        
cat filt_"$i"_"$reads"_sub_hmmoutF.out.bed | bedtools sort -i - | uniq > all_mon_"$i"_"$reads".bed
samtools faidx "$reads" 
bedtools getfasta -fi "$reads" -bed all_mon_"$i"_"$reads".bed  > all_"$i"_"$reads"_monomers.fa
blat all_"$i"_"$reads"_monomers.fa /scratch/rdp22327/Dawe/consensus/"$i"_main_consensus.fasta -t=dna -q=dna -maxGap=10 -minScore=0 -repMatch=2147483647 all_"$i"_"$reads"_monomers.fa.blat
awk '$1 ~ /^[0-9]*$/' all_"$i"_"$reads"_monomers.fa.blat | awk 'NR>1' | awk '{print $14"\t"($1/($11+$15-$1))}' > all_"$i"_"$reads"_monomers.fa.blat.sub
awk '{print $1}' all_"$i"_"$reads"_monomers.fa.blat.sub | sed 's/:.*//' | sort | uniq -c | awk '{if($1>=10) print $2}' > filt_all_"$i"_"$reads"_monomers.fa.blat.sub 
done 

#now compare all to all and predict structure with R script
for rep in CentC knob180 TR1 Cent4
do
reads=rep_all_"$lin".fastq.fasta
mkdir read_out_"$rep"_"$reads"_norep
while read seq
do
seq2=$(echo $seq| sed "s:/:_:g")
grep -wEA1 --no-group-separator "$seq" all_"$rep"_"$reads"_monomers.fa | sed "s:/:_:g"> sub_all_"$rep"_"$reads"_monomers2.fa
blat  sub_all_"$rep"_"$reads"_monomers2.fa  sub_all_"$rep"_"$reads"_monomers2.fa -t=dna -q=dna -maxGap=10 -minScore=0 -repMatch=2147483647  sub_all_"$rep"_"$reads"_monomers2.fa.blat
cat sub_all_"$rep"_"$reads"_monomers2.fa.blat | awk '$1 ~ /^[0-9]*$/'| awk 'NR>1'| awk '{print $1"\t"$10"\t"$11"\t"$14"\t"$15}' > ./read_out_"$rep"_"$reads"_norep/"$seq2".sub
done<filt_all_"$rep"_"$reads"_monomers.fa.blat.sub 
if [ ! -f read_out_"$rep"_"$reads"_norep/*_reads_class.csv ]
then
Rscript --vanilla ../mono/class.R $rep ./read_out_"$rep"_"$reads"_norep
fi
done


##summing up the monomers by structure
for rep in CentC Cent4 TR1 knob180
do
cat ./read_out_"$rep"_"$reads"_norep/"$rep"_reads_class.csv | awk -F"," '{print $12"\t"$3}'|
awk '{
    arr[$1]+=$2
   }
   END {
     for (key in arr) printf("%s\t%s\n", key, arr[key])
   }' \
   | sort -k1,1 | awk -v var1="$dir" -v var2="$rep" '{if($2>0) print var1"\t"var2"\t"$1"\t"$2}' >> All_Raw_"$lin"_sum
done

#######Combining at the end
#cat All_Raw_*_sum > All_Raw
