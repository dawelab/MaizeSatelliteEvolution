# MaizeSatelliteEvolution


Repeat Consensus Sequences : https://github.com/dawelab/ETR_in_HiFi/tree/main/Consensus_Repeats

 outputs: https://github.com/dawelab/MaizeSatelliteEvolution/tree/main/out
 
 figures: https://github.com/dawelab/MaizeSatelliteEvolution/tree/main/figures 

 
## Project outline:
1. Issues in OLD Assemblies
**_Repeat Content Estimation._**
	- Repeats in Raw Data, mapping : Long reads from Gapless Ab10 assemly are aligned to the gapless Ab10 assembly using minimap2. Read hits are filtered to a single positions (2308) and averaged in 10kb windows.
https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/OldAssem_ReadsAlign_RepContent.sh

2. Issues in NEW Assemblies
   
   **_Repeat Content Estimation._**
    - HiFi Reads Prep : Raw HiFi are converted to CCS in 10 chunks, then bam files are merged together, filtered to a minimum quality of 99%, and converted to a fastq.
      https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/HiFi_reads_DownPrep.sh

   **_HiFi Assembly and Annotation._**
    - HiFi Contig Generation, summary, scaffolding, annotation : Contigs with > 1/2 the expected read depth are extracted and scaffolded on Mo17 T2T and B73v5 as an alternate. Assemblies annotated with Liftoff. Repeat content assessed.
      https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/HiFiScaffolding.sh
    - 	CG119 correction: The incorrectly placed contig is moved with bedtools and remerged into a new assembly.
      https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/CG119_correction.sh

    - HiFi Reads Repeats, Length, Coverage , Quality : CCS reads are compared to repeat consensus sequences using blast. Total repeat content is then calculated.Then, contigs are generated and the same process is repeated. Finally, the contig v satellite repeat output is used to check for anchored ends on contigs, defined as contigs that do not have a satellite hit within 100bp of start or end bp.
      https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/RepContentAndContigAnchor.sh

    - biases
   
   		flow cell code: https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/flowcell_bias.sh
   
   		length bias : https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/length_bias.sh
      

4. Arrays
   **_Conserved Repeat Positions._**
   - Finding Arrays : Repeat hits within 10kb were merged using bedtools. Arrays with >=10 monomers and >= 10% repeat content were defined as repeat arrays. Array positions were then compared to 100N gaps using bedtools. Then, array positions were compared to core genes using bedtools closest.
     https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/ArrayPositions.sh	
     
   - Conserved Array Positions : Array positions were clustered together based on repeat type and shared up and down core genes.
     https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/Conserved_Arrays.R


5. Structure
   
   **_Repeat Structure._**
   
   - Model : Classification model creation from sample data. LDA model chosen to further analysis.
     https://github.com/dawelab/MaizeSatelliteEvolution/tree/main/Model_Classification
   - Raw reads : Raw repeat structure predicted from monomer similarity network toplogy using LDA model.
      https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/Raw_Structure.sh
   - Assembly structure : Assembly reoeat structure predicted using non-overlapping 10kb bins from monomer similarity network toplogy using LDA model.
     https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/ArrayStructure.sh
   - Assessing HOR structure: labeling monomers by cluster, checking character string for repeating kmers
     https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/HOR_re-eval.sh
     
6. Shared HOR patterns
   - Generating consensus monomers for all identified HOR's
     https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/HOR_reclustering.sh
   - clustering consensus monomers and looking for patterns
     https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/HOR_recluster_Pattern.sh
   - Purity of shared patterns and pulling HORs that occur multiple times within one array
     https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/functions/Shared_HOR_purity.R

7. Repeat Array Comparisons
   -   Grouping homologous arrays-- dot plots
	 https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/Array_Homolog_Grouping.sh
   -   Grouping homologous arrays-- pariwise comparisons
    	https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/Homolog_Similiarity.sh
    - Array dot plot generation and summary in R
      	https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/HomologousArrayComps.R

8. Chip-seq
  	-  Chipseq alignemnt to Mo17, using unique k-mers
    	https://github.com/dawelab/MaizeSatelliteEvolution/blob/main/src/chip_Mo17.sh
	-  plotting in R with patterns
	 <<<<< insert code  >>>>
