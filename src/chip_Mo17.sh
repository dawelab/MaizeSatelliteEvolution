#!/bin/bash
ref=" Mo17_half_Mo17_assem.fasta"

#find unique k-mers
module load Meryl/1.4.1
module load BWA/0.7.17-GCCcore-11.2.0

bwa index Mo17_half_Mo17_assem.fasta

meryl count k=51 Mo17_half_Mo17_assem.fasta output Mo17_half_Mo17_assem.fasta.meryl #kmers of 51
meryl equal-to 1 Mo17_half_Mo17_assem.fasta.meryl output unique.meryl #only occur once
meryl-lookup -bed-runs -sequence Mo17_half_Mo17_assem.fasta -mers unique.meryl -output unique.bed  #grab all the single occurrence k-mers, put in bed file


#get raw reads and prep
module load  SRA-Toolkit/3.0.1-centos_linux64
module load fastp/0.23.2-GCC-11.2.0
module load BWA/0.7.17-GCCcore-11.2.0

for i in SRR21509776 SRR21509777 SRR21509778 SRR21509779
do
fasterq-dump $i
fastp -i "$i"_1.fastq -I "$i"_2.fastq -o "$i".R1_Chip.fq.gz -O "$i".R2_Chip.fq.gz --dedup --detect_adapter_for_pe --cut_front --cut_tail
done

#now align reads to reference assembly
module load SAMtools/0.1.19-foss-2019b
module load BEDTools/2.29.2-GCC-8.3.0

for i in SRR21509776 SRR21509777 SRR21509778 SRR21509779
do
bwa mem -k 50 -c 1000000 -t 20 $ref "$i".R1_*.fq.gz "$i".R2_*.fq.gz  > "$i"_Mo17new.sam #align
samtools view -b -S -F 2308 "$i"_Mo17new.sam >  "$i"_Mo17new.bam #filter to unique mapping
bedtools intersect -b unique.bed -abam "$i"_Mo17new.bam -ubam -wa -F 1 > unique_"$i"_Mo17new.bam #only hits that completely overlap a unique k-mer in the assembly
samtools sort "$i"_Mo17new.bam -o sort_"$i"_Mo17new > sort_"$i"_Mo17new.bam
done

#normalize input and chip
module load deepTools/3.5.2-foss-2022a
#rep1
chip="SRR21509777"
in="SRR21509779"
bamCompare -b1 sort_unique_"$chip"_Mo17new.bam -b2 sort_unique_"$in"_Mo17new.bam --operation ratio --binSize 1000 --scaleFactorsMethod None  --normalizeUsing RPKM -o "$chip"_"$in"_1000_rpkmNORM.bw
bamCompare -b1 sort_unique_"$chip"_Mo17new.bam -b2 sort_unique_"$in"_Mo17new.bam --operation ratio --binSize 1000 --scaleFactorsMethod None  --normalizeUsing ratio -o "$chip"_"$in"_1000_ratioNORM.bw

#rep2
chip="SRR21509776"
in="SRR21509778"
bamCompare -b1 sort_unique_"$chip"_Mo17new.bam -b2 sort_unique_"$in"_Mo17new.bam --operation ratio --binSize 1000 --scaleFactorsMethod None  --normalizeUsing RPKM -o "$chip"_"$in"_1000_rpkmNORM.bw
bamCompare -b1 sort_unique_"$chip"_Mo17new.bam -b2 sort_unique_"$in"_Mo17new.bam --operation ratio --binSize 1000 --scaleFactorsMethod None  --normalizeUsing ratio -o "$chip"_"$in"_1000_ratioNORM.bw





