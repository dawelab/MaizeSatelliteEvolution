#!/bin/bash

reads_dir=#reads dir
lines=#line of interest
assem=#assembly of interest


module load bioawk/1.0-GCC-11.2.0
module load EMBOSS/6.6.0-foss-2021b
module load BLAST+/2.10.1-gompi-2022a
module load BEDTools/2.29.2-GCC-8.3.0

#rep_content function in src file

#raw read repeat content
##ONT
i=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $reads_dir/ont_raw_list ) #ONT raw, split into 15 sub files for analysis
rep_content $i repeats.fasta ONT

#PB
i=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $reads_dir/PB_accessions )
rep_content $i repeats.fasta PB

#illum
i=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $reads_dir/illum_dat )
rep_content $i repeats.fasta IL

##HiFi
i="all_Ab10_HiFiDat.fasta"
rep_content $i repeats.fasta HiFi

#Assem
rep_content $assem_dir/$assem repeats.fasta Assem


##combining info
for i in B73_Ab10_Gapless_fixed.fa Zm-B73-REFERENCE-NAM-5.0.fa Zm-NC350-REFERENCE-NAM-1.0.fa
do
for reads in ONT PB IL HiFi Assem Assem_Chr
do
  tot_len=$(cat "$i"*len_"$reads" | awk '{Total=Total+$2} END{print Total}')
for rep in CentC Cent4 knob180 TR1
  do
  rep_len=$(cat "$i"*_repeats.blast_sum_"$rep"_"$reads" |  awk  '{Total=Total+$2} END{print rep"\t"Total}')
  echo -e $reads"\t"$rep"\t"$((rep_len))"\t"$((tot_len)) >> output.vals_"$i"
done
done
done
#################################

#raw read alignment to assemblies
reads_dir=#reads dir
assem_dir=#assem dir

##tools
module load minimap2/2.22-GCCcore-11.2.0
module load SAMtools/0.1.20-GCC-11.2.0
module load BEDTools/2.29.2-GCC-8.3.0
module load seqtk/1.3-GCC-11.2.0

window_maker_bed(){ # window_maker_bed assem.fa window_size 
  while read lin
  do
  chr=$(echo $lin | awk '{print $1}' )
  max_len=$(echo $lin | awk '{print $3}' )
  start=0
  end=$2
  while [ $end -le $max_len ]
  do
  echo $chr $start $end >> "$1".fai.bed_$2.bed
  ((start+=$2))
  ((end+=$2))
  done
  echo $lin
  done<"$1".fai.bed

  awk '{print $1"\t"$2"\t"$3}' "$1".fai.bed_$2.bed > "$1".fai.bed_$2_fix.bed
}

##filtering the mapped reads
filter_map() { # filter_map out.sam map.bed
  samtools view --threads 10 -bF 2308 -o filt_$1.bam $1
  bedtools bamtobed -i filt_$1.bam | bedtools sort >  filt_$1.bam.bed
  bedtools coverage -a $2  -b filt_$1.bam.bed -hist  | awk '{print $1"\t"$3-2"\t"$3"\t"$4*$7}' | bedtools merge -c 4 -o sum -i - > $2_depth_filt_$1.bam.bed
}

#subset pseudo molecules
seqtk subseq $assem chr_nams > chr_"$assem"

##PB data, split into 22 files##
i=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $reads_dir/PB_accessions) #submitted as an array job
minimap2 -t 20 -ax  map-pb $assem_dir/$assem  $i > "$i".PB_$assem.sam
minimap2 -t 20 -ax  map-pb $assem_dir/chr_"$assem"  $i > "$i".PB_chr_$assem.sam

##hifi data for Ab10##
i="all_Ab10_HiFiDat.fasta"
minimap2 -t 20 -cx map-hifi $assem_dir/$assem  all_Ab10_HiFiDat.fasta > "$i".HiFi_$assem.sam
minimap2 -t 20 -cx map-hifi $assem_dir/chr_"$assem"  all_Ab10_HiFiDat.fasta > "$i".HiFi_chr_$assem.sam
#ONT data, split into 15 files, for Ab10
i=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $reads_dir/ont_raw_list) #submitted as an array job
minimap2 -t 20 -ax map-ont $assem_dir/$assem  $i > "$i".ONT_$assem.sam
minimap2 -t 20 -ax map-ont $assem_dir/chr_"$assem"  $i > "$i".ONT_chr_$assem.sam

for window in 100 10000
do
samtools faidx $assem
grep "chr" "$assem".fai | awk '{print $1"\t"1"\t"$2}' | bedtools sort -i - > "$assem".fai.bed
window_maker_bed $assem $window
for i in *sam
do
filter_map $i ../"$assem".fai.bed_"$window"_fix.bed
done
done
