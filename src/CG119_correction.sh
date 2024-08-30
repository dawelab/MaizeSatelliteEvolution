cat CG119_HM_assem.bp.p_ctg.noseq.gfa_Filt_repvals_unanchoredEnds | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6}' | sed 's/rd:i://'> contig_info

#cat CG119_half_Mo17_ragtag.scaffold.agp | grep "chr7_RagTag"
#cat CG119_out_HALF_ragtag.scaffold.agp | grep "chr7_RagTag"
##ptg000006l is the misplaced contig

#create bedfile to fix chr7 CG119_half_Mo17_ragtag.scaffold.fasta
# CG119_half_Mo17_fix_chr7.bed
#chr7_RagTag 1 22179126
#chr3_RagTag    119322  85947540 #where ptg000006l contig incorrectly is placed
#chr7_RagTag    22179027        22179126 #repeat 100Ngap
#chr7_RagTag    22179127        111878821 #rest of the chromosome

#now bedfile to fic chr3 CG119_half_Mo17_ragtag.scaffold.fasta
# CG119_half_Mo17_fix_chr3.bed
#chr3_RagTag    1       119321
#chr3_RagTag    85947641        307486706

#cat CG119_half_Mo17_fix_chr7.bed | awk '{print $1"\t"$2"\t"$3}' > fix_CG119_half_Mo17_fix_chr7.bed
#cat CG119_half_Mo17_fix_chr3.bed | awk '{print $1"\t"$2"\t"$3}' > fix_CG119_half_Mo17_fix_chr3.bed


#module load BEDTools/2.29.2-GCC-8.3.0
#bedtools getfasta -fi CG119_half_Mo17_ragtag.scaffold.fasta -bed fix_CG119_half_Mo17_fix_chr7.bed > CG119_half_Mo17_fix_chr7.fasta
#bedtools getfasta -fi CG119_half_Mo17_ragtag.scaffold.fasta -bed fix_CG119_half_Mo17_fix_chr3.bed > CG119_half_Mo17_fix_chr3.fasta

#module load EMBOSS/6.6.0-foss-2021b
#union -filter  CG119_half_Mo17_fix_chr7.fasta > merge_CG119_half_Mo17_fix_chr7.fasta
#union -filter  CG119_half_Mo17_fix_chr3.fasta > merge_CG119_half_Mo17_fix_chr3.fasta

#cat CG119_half_Mo17_ragtag.scaffold.fasta merge_CG119_half_Mo17_fix_chr7.fasta merge_CG119_half_Mo17_fix_chr3.fasta > cat_CG119_half_Mo17_ragtag.scaffold.fasta

#grep ">" cat_CG119_half_Mo17_ragtag.scaffold.fasta |  sed 's/>//' > names_cat_CG119_half_Mo17_ragtag.scaffold.fasta
#array=$(cat names_cat_CG119_half_Mo17_ragtag.scaffold.fasta)
#del="chr7_RagTag"
#array=("${array[@]/$del}")
#del="chr3_RagTag"
#array=("${array[@]/$del}")

#echo $array |  tr ' ' '\n' > sub_names_cat_CG119_half_Mo17_ragtag.scaffold.fasta

#module load seqtk/1.3-GCC-11.2.0
#seqtk subseq cat_CG119_half_Mo17_ragtag.scaffold.fasta sub_names_cat_CG119_half_Mo17_ragtag.scaffold.fasta > sub2_cat_CG119_half_Mo17_ragtag.scaffold.fasta
#cat sub2_cat_CG119_half_Mo17_ragtag.scaffold.fasta| sed 's/1-119321/chr3_RagTag/'|sed 's/1-22179126/chr7_RagTag/'  > CG119_half_Mo17_ragtag.scaffold_FIXED.fasta
#cp CG119_half_Mo17_ragtag.scaffold_FIXED.fasta ragtag.scaffold.fasta

module load  SeqKit/0.16.1
#seqkit locate --ignore-case --only-positive-strand --pattern "N" ragtag.scaffold.fasta | grep "chr" | awk '{print $1"\t"$5"\t"$6}' | bedtools merge -i - > 100n.bed


dir="$lin"_half_Mo17
db="/scratch/rdp22327/Dawe/Mo17/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.gff3_db"
ref_gff="/scratch/rdp22327/Dawe/Mo17/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.gff3"

liftoff -o ./"$dir"/"$dir".gff -dir ./"$dir" -db $db ./"$dir"/ragtag.scaffold.fasta $ref

#chr.list is a file with chr1 chr2 etc.
cat /scratch/rdp22327/Dawe/scaffolding/"$dir"/ragtag.scaffold.fasta | seqtk subseq - chr.list |  sed "s/RagTag/$lin/"  > "$lin"_chr_ragtag.scaffold.fasta
rep_content  "$lin"_chr_ragtag.scaffold.fasta repeats.fasta Assem_Chr



##rerun array_positions.sh
