module load BEDTools/2.29.2-GCC-8.3.0
module load HMMER/3.1b2-foss-2019b
module load BLAT/3.5-GCC-8.3.0
module load EMBOSS/6.6.0-GCC-8.3.0-Java-11
module load seqtk/1.3-GCC-10.2.0


# *_repeats.blast is na output from raw_repeat_content.sh

reads=${fastq file of interest}
consensus_dir=${consensus}
cat "$reads"_repeats.blast| awk '{if($10>$9) {print $2"\t"$9"\t"$10"\t"$10-$9"\t"$1} else{print $2"\t"$10"\t"$9"\t"$9-$10"\t"$1}}' | bedtools sort -i - | bedtools merge -i - -d 100000 -c 4,4,5 -o sum,count,distinct | sort  >m_count_"$reads"_repeats.blast

#getting reads that have some ETR
seqtk subseq  "$reads" m_count_"$reads"_repeats.blast > rep_"$reads".fasta

for i in CentC TR1 knob180 Cent4
do
f_hmm_consensus="$consensus_dir"/"$i"_consensus_fix.fasta.afa.sto.hmm #forward consensus

hmmsearch --tblout hmmoutF.tbl -o "$i"_"$reads"_hmmoutF.out $f_hmm_consensus rep_$reads.fasta
awk '($1 ~ /^>>/) || ($2=="!") ' "$i"_"$reads"_hmmoutF.out > filt_"$i"_"$reads"_hmmoutF.out

while read lin
do
if [[ ${lin:0:1} == ">" ]] 
then
read=${lin#">> "}
else
echo $lin $read >> filt_"$i"_"$reads"_hmmoutF.out_mod
fi
done<filt_"$i"_"$reads"_hmmoutF.out

cat filt_"$i"_"$reads"_hmmoutF.out_mod | awk '{print $0"\t"$8-$7}' | awk '{if($18 >= 120 ) print $17"\t"$13"\t"$14}'| bedtools sort -i -| uniq  > all_mon_"$i"_"$reads".bed
bedtools getfasta -fi rep_"$reads".fasta -bed all_mon_"$i"_"$reads".bed  > all_"$i"_"$reads"_monomers.fa
        
blat all_"$i"_"$reads"_monomers.fa /scratch/rdp22327/Dawe/consensus/"$i"_main_consensus.fasta -t=dna -q=dna -maxGap=10 -minScore=0 -repMatch=2147483647 all_"$i"_"$reads"_monomers.fa.blat
awk '$1 ~ /^[0-9]*$/' all_"$i"_"$reads"_monomers.fa.blat | awk 'NR>1' | awk '{print $14"\t"($1/($11+$15-$1))}' > all_"$i"_"$reads"_monomers.fa.blat.sub
awk '{print $1}' all_"$i"_"$reads"_monomers.fa.blat.sub | sed 's/:.*//' | sort | uniq -c | awk '{if($1>=10) print $2}' > filt_all_"$i"_"$reads"_monomers.fa.blat.sub 
done 

####
module load R/4.3.0-foss-2020b

for rep in TR1 knob180 Cent4 CentC
do
mkdir read_out_"$rep"_"$reads"
while read seq #for each read, subset all the monomers and compare
do
seq2=$(echo $seq| sed "s:/:_:g")
grep -wEA1 --no-group-separator "$seq" all_"$rep"_"$reads"_monomers.fa | sed "s:/:_:g"> sub_all_"$rep"_"$reads"_monomers.fa
blat  sub_all_"$rep"_"$reads"_monomers.fa  sub_all_"$rep"_"$reads"_monomers.fa -t=dna -q=dna -maxGap=10 -minScore=0 -repMatch=2147483647  sub_all_"$rep"_"$reads"_monomers.fa.blat
cat sub_all_"$rep"_"$reads"_monomers.fa.blat | awk '$1 ~ /^[0-9]*$/'| awk 'NR>1' > ./read_out_"$rep"_"$reads"_fixed/"$seq2".sub
done<filt_all_"$rep"_"$reads"_monomers.fa.blat.sub 
Rscript --vanilla ../mono/network_sliding_smooth_reads.R $rep ./read_out_"$rep"_"$reads"_fixed #classification for reads with >=5 reads
done
