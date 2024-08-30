#!/bin/bash
module load EMBOSS/6.6.0-foss-2021b
module load BLAT/3.5-GCC-11.2.0
module load seqtk/1.3-GCC-11.2.0
module load BEDTools/2.30.0-GCC-11.2.0
module load HMMER/3.3.2-gompi-2022a
module load BLAST+/2.2.31
module load R/4.1.2-foss-2021b

cd /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff

dir=$(awk "NR==${SLURM_ARRAY_TASK_ID}"  Mo17_dirs) 
#dir="CG119_fix"
lin=$(echo $dir | sed 's/_.*//' )

echo $dir
echo $lin

cat "$dir"/ragtag.scaffold.fasta_repeats.blast_*_arrays.bed > "$dir"/ragtag.scaffold.fasta_repeats.blast_ALL_arrays.bed
bedtools getfasta -fi "$dir"/ragtag.scaffold.fasta -bed "$dir"/ragtag.scaffold.fasta_repeats.blast_ALL_arrays.bed > "$dir"/ragtag.scaffold.fasta_arrays.fasta

makehmmerdb  "$dir"/ragtag.scaffold.fasta_arrays.fasta  "$dir"/ragtag.scaffold.fasta_arrays.fasta.db 

for i in knob180 Cent4 CentC TR1
do
f_hmm_consensus=/scratch/rdp22327/Dawe/consensus/"$i"_consensus_fix.fasta.afa.sto.hmm
nhmmer --tblout "$dir"/"$i"_hmmoutF.tbl -o "$dir"/"$i"_hmmoutF.out $f_hmm_consensus "$dir"/ragtag.scaffold.fasta_arrays.fasta.db
cat  "$dir"/"$i"_hmmoutF.out  | awk '{print $4"\t"$5"\t"$6}' | grep ":"|  awk '{if($3>$2) {print $1"\t"$2"\t"$3} else {print $1"\t"$3"\t"$2}}'| grep 'chr\|ptg' | bedtools sort -i - > "$dir"/filt_"$i"_hmmoutF.out.bed
bedtools getfasta -fi "$dir"/ragtag.scaffold.fasta_arrays.fasta -bed "$dir"/filt_"$i"_hmmoutF.out.bed > "$dir"/mono_"$i"_ragtag.scaffold.fasta_arrays.fasta
blat  "$dir"/mono_"$i"_ragtag.scaffold.fasta_arrays.fasta /scratch/rdp22327/Dawe/consensus/"$i"_main_consensus.fasta -t=dna -q=dna -maxGap=10 -minScore=0 -repMatch=2147483647  "$dir"/consensus_mono_"$i"_ragtag.scaffold.fasta_arrays.fasta.blat 
awk '$1 ~ /^[0-9]*$/' "$dir"/consensus_mono_"$i"_ragtag.scaffold.fasta_arrays.fasta.blat | awk 'NR>1' | awk '{print $14"\t"($1/(($11/2)+$15-$1))}' > "$dir"/consensus_mono_"$i"_ragtag.scaffold.fasta_arrays.fasta.blat.sub
awk -v var="$i" '{print $0"\t"var}' "$dir"/consensus_mono_"$i"_ragtag.scaffold.fasta_arrays.fasta.blat.sub > "$dir"/nam_consensus_mono_"$i".fasta_arrays.fasta.fa.blat.sub
done

cat nam_consensus_mono_*.fasta.fa.blat.sub > all_nam_consensus_mono.fasta.fa.blat.sub

val1_orig=0
val_orig=10000
up=10000
suff=$(echo sliding_"$val_orig")

for i in knob180 Cent4 CentC TR1
do
mkdir "$dir"/structure_"$i"_"$lin"_arrays2
length_mono=$(cat /scratch/rdp22327/Dawe/consensus/"$i"_consensus_fix_RC.fasta.afa.sto.hmm | grep "LENG" | awk '{print $2*.8}')
cat "$dir"/filt_"$i"_hmmoutF.out.bed | awk '{print $1}' | sort | uniq > "$dir"/arrays_"$lin" 

while read ar
do
grep $ar "$dir"/filt_"$i"_hmmoutF.out.bed | bedtools sort -i - | awk '{if($3-$2>= 120) print $0}' > "$dir"/monos_"$lin".bedrows_mono=$(cat "$dir"/monos_"$lin".bed | wc -l)
if [ $rows_mono -gt 5 ]
then
high=$(cat "$dir"/monos_"$lin".bed | tail -n1 | awk '{print $3}' | sort -n  | tail -n1)
low=$(cat "$dir"/monos_"$lin".bed | head  -n1 | awk '{print $2}' | sort -n  | head -n1)
val1=$low
val=$((low+10000))
len=$((high-low))
lenmax=$((high+9999))
 while [ $val -le $lenmax ]
  do
        awk -v v1="$val1" -v v="$val"  '{if( $2 >= v1 && $3 <= v) print $1":"$2"-"$3}'  "$dir"/monos_"$lin".bed > "$dir"/bin_"$lin".bed
        rows=$(cat "$dir"/bin_"$lin".bed  | wc -l)
        if [ $rows -gt 5 ]
        then
        seqtk subseq "$dir"/mono_"$i"_ragtag.scaffold.fasta_arrays.fasta "$dir"/bin_"$lin".bed > "$dir"/all_bin_"$lin".fa
        mon=$(grep ">" "$dir"/all_bin_"$lin".fa| wc -l)
        blat "$dir"/all_bin_"$lin".fa "$dir"/all_bin_"$lin".fa -t=dna -q=dna -maxGap=10 -minScore=0 -repMatch=2147483647 "$dir"/all_bin_"$lin".fa_self.blat
        cat "$dir"/all_bin_"$lin".fa_self.blat | awk '$1 ~ /^[0-9]*$/' | awk 'NR>1'  | awk '{print $1"\t"$10"\t"$11"\t"$14"\t"$15}' > ./"$dir"/structure_"$i"_"$lin"_arrays2/"$val"_"$ar"_"$i"_comp.blat_"$suff".sub
        fi
        ((val1+=up))
        ((val+=up))
  done 
echo $ar
fi
done<"$dir"/arrays_"$lin"
#predict structure
Rscript --vanilla /scratch/rdp22327/Dawe/mono/class.R $i /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/"$dir"/structure_"$i"_"$lin"_arrays2
#smooth bins
Rscript --vanilla /scratch/rdp22327/Dawe/mono/smoothing.R $i /scratch/rdp22327/Dawe/scaffolding/Mo17_scaff/"$dir"/structure_"$i"_"$lin"_arrays2
done

###### Totals by line
for dir in TIL01 TIL11 TIL25 K64 CG119 CML442 CG108 Tx777 Tx779 CG44 AB10 B73 Mo17
do
for rep in CentC Cent4 TR1 knob180
do
cat *"$dir"*_half_Mo17/structure_"$rep"_*"$dir"*arrays2/"$rep"_reads_class.csv | awk -F"," '{print $12"\t"$3}'|
awk '{
    arr[$1]+=$2
   }
   END {
     for (key in arr) printf("%s\t%s\n", key, arr[key])
   }' \
   | sort -k1,1 | awk -v var1="$dir" -v var2="$rep" '{if($2>0) print var1"\t"var2"\t"$1"\t"$2}' >> All_Arrays_Totals
done
done

