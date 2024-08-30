#!/bin/bash

module load RagTag/2.0.1-foss-2022a
module load seqtk/1.3-GCC-11.2.0
module load Liftoff/1.6.3

lin=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Data_Lines)
#subsetting the contigs to contigs that are >1/2 of the expected reads depth -- high confidence contigs.
grep $lin HiFi_ReadCoverage.csv  | awk -F"," '{print $2}' | sed 's/x//' | awk '{print $1/2}' > cov_"$lin"
cov=$( grep $lin HiFi_ReadCoverage.csv  | awk -F"," '{print $2}' | sed 's/x//' | awk '{print $1/2}')
cov_half=$((cov_val/2))
cat sub_"$lin"_assem.bp.p_ctg.noseq.gfa_Filt_repvals  | awk  -v var="$cov" '{if($2>var) print $1}' | sort | uniq >half_cov_contigs_"$lin"
seqtk subseq  "$lin"_HM_assem.bp.p_ctg.gfa.fasta half_cov_contigs_"$lin" > half_"$lin"_assem.bp.p_ctg.gfa.fasta 

#############################################################################################################################
##with Mo17 assembly as reference
ref="/scratch/rdp22327/Dawe/Mo17/Mo17_t2t.fna" #T2T Mo17

#scaffolding all the contigs on T2T Mo17
mkdir "$lin"_half_Mo17
ragtag.py scaffold -o "$lin"_half_Mo17  $ref half_"$lin"_assem.bp.p_ctg.gfa.fasta 

dir="$lin"_half_Mo17
db="/scratch/rdp22327/Dawe/Mo17/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.gff3_db"
ref_gff="/scratch/rdp22327/Dawe/Mo17/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.gff3"

liftoff -o ./"$dir"/"$dir".gff -dir ./"$dir" -db $db ./"$dir"/ragtag.scaffold.fasta $ref

#chr.list is a file with chr1 chr2 etc.
cat /scratch/rdp22327/Dawe/scaffolding/"$dir"/ragtag.scaffold.fasta | seqtk subseq - chr.list |  sed "s/RagTag/$lin/"  > "$lin"_chr_ragtag.scaffold.fasta
rep_content  "$lin"_chr_ragtag.scaffold.fasta repeats.fasta Assem_Chr





##OLDER VERISON WITH ALTERNATE REFERENCE
#############################################################################################################################
##with B73 assembly as reference
ref="Zm-B73-REFERENCE-NAM-5.0.fa" #B73v5

#scaffolding high confidence contigs on B73 v5
mkdir "$lin"_out_HALF
ragtag.py scaffold -w -o "$lin"_out_HALF  $ref half_"$lin"_assem.bp.p_ctg.gfa.fasta

#annotate
dir="$lin"_out_HALF
db="/scratch/rdp22327/Dawe/ref_assemblies/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3_db"
ref_gff="/scratch/rdp22327/Dawe/ref_assemblies/chr_Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
liftoff -o "$dir"/"$dir".gff -dir "$dir" -g "$ref_gff" "$dir"/ragtag.scaffold.fasta "$ref"

#chr.list is a file with chr1 chr2 etc.
cat /scratch/rdp22327/Dawe/scaffolding/"$dir"/ragtag.scaffold.fasta | seqtk subseq - chr.list |  sed "s/RagTag/$lin/"  > "$lin"_chr_ragtag.scaffold.fasta
rep_content  "$lin"_chr_ragtag.scaffold.fasta repeats.fasta Assem_Chr



######### combine
#core genes = B73.all.csv, from Yibing'a paper
# *Sb_subgenomes_complete_longest.bed from NAM paper
awk -F"," '{if($12=="Core Gene") print "chr"$2"\t"$3"\t"$4"\t"$1}' B73.all.csv > B73.all.csv.bed
awk '{print $4}' B73.all.csv.bed | sort > nams_core
i="B73" #ref
cat Zm-"$i"-REFERENCE-NAM-*_Sb_subgenomes_complete_longest.bed | awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -k4 > sub_"$i".bed
join -1 4 -2 1 sub_"$i".bed nams_core | sort -k5 | join -1 5 -2 1 - /scratch/rdp22327/Dawe/all_gff_count.bed_single2 |awk '{print $3"\t"$4"\t"$5"\t"$1}' | bedtools sort -i - |bedtools merge -i - -c 4 -o distinct > nam_NEW_sub_"$i".bed


##NAM lines
while read i
do
for rep in CentC Cent4 knob180 TR-1
do
cat TRANSPOSABLE_ELEMENTS/Zm-"$i"-REFERENCE-NAM-*TE.gff3  | grep -a "$rep" | grep -a  "chr" | awk  -v rep="$rep" '{ print $1"\t"$4"\t"$5"\t"rep}' |bedtools sort -i - | bedtools merge -i - -d 10000 -c 4 -o distinct >"$i"_"$rep"_arrays_NOFILT.bed
done
cat "$i"_*_arrays_NOFILT.bed | bedtools sort -i - > "$i"_all_rep_gaps_NOFILT.bed

bedtools closest -id -D a -a "$i"_all_rep_gaps_NOFILT.bed -b nam_NEW_sub_"$i".bed > "$i"_arrays.bed_new_up
bedtools closest -iu -D a -a "$i"_all_rep_gaps_NOFILT.bed -b nam_NEW_sub_"$i".bed> "$i"_arrays.bed_new_down
paste "$i"_arrays.bed_new_up "$i"_arrays.bed_new_down | awk  -v lin="$i"  '{print lin"\t"$1"\t"$2"\t"$3"\t"$4"\t"$3-$2"\t"$7"\t"$8"\t"$15"\t"$17}' >> array_gen_coords
echo $i
done<Lines_of_Interest
