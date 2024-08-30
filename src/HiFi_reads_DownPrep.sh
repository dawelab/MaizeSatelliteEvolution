#!/bin/bash

#HiFi_reads.sh <input_bam> <ccs_nam> #this script uses 
#HiFi_reads.sh /scratch/rdp22327/Dawe/hifi_assem/r64060_20210619_183804/2_B01/m64060_210621_022359.subreads.bam 2_B01

module load CCS/6.3.0_conda
#path to one of my flow cells:/scratch/rdp22327/Dawe/hifi_assem/r64060_20210619_183804/2_B01/m64060_210621_022359.subreads.bam 
#ccs converts subreads to ccs reads. below is how I chunked it to run faster with 20 threads

ccs $1 "$2".ccs.1.bam --chunk 1/10 -j 20
ccs $1 "$2".ccs.2.bam --chunk 2/10 -j 20
ccs $1 "$2".ccs.3.bam --chunk 3/10 -j 20
ccs $1 "$2".ccs.4.bam --chunk 4/10 -j 20
ccs $1 "$2".ccs.5.bam --chunk 5/10 -j 20
ccs $1 "$2".ccs.6.bam --chunk 6/10 -j 20
ccs $1 "$2".ccs.7.bam --chunk 7/10 -j 20
ccs $1 "$2".ccs.8.bam --chunk 8/10 -j 20
ccs $1 "$2".ccs.9.bam --chunk 9/10 -j 20
ccs $1 "$2".ccs.10.bam --chunk 10/10 -j 20

module load SAMtools/0.1.19-foss-2019b
#merge those files together for your final bam
samtools merge -@20  "$2".ccs.bam  "$2".ccs.*.bam 

module load BamTools/2.5.2-GCC-11.2.0
bamtools filter -in "$2".ccs.bam -out "$2"_99.ccs.bam -tag "rq":">=0.99" 

#I also converted to fastq for other checks, and continued to use fastq through my process
module load  BEDTools/2.30.0-GCC-10.2.0
bamToFastq -i "$2".ccs.bam -fq "$2".fastq


