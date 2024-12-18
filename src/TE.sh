#!/bin/bash

dir=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Mo17_dir)
cd $dir
module load BLAST+/2.2.31


blastn -subject /scratch/rdp22327/Dawe/consensus/maizeTE02052020.fasta -query ragtag.scaffold.fasta -outfmt 6 -num_threads 10 > ragtag.scaffold.fasta_arrays.fasta_TE.blast
cat ragtag.scaffold.fasta_arrays.fasta_TE.blast | awk '{print $1}' | awk -F":" '{print $2}' | awk -F"-" '{print $1}'> coord1
cat ragtag.scaffold.fasta_arrays.fasta_TE.blast | awk '{print $1}' | awk -F":" '{print $1}' > chr
cat ragtag.scaffold.fasta_arrays.fasta_TE.blast | awk '{print $2}' | awk -F"#" '{print $2}' > nam
paste nam ragtag.scaffold.fasta_arrays.fasta_TE.blast  | awk '{print $2"\t"$8"\t"$9"\t"$1}' | grep "chr" > "$dir"_TE.bed
