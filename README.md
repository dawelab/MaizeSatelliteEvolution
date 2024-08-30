# MaizeSatelliteEvolution


Repeat Consensus Sequences : https://github.com/dawelab/ETR_in_HiFi/tree/main/Consensus_Repeats

## Project outline:
1. Issues in OLD Assemblies
**_Repeat Content Estimation._**
	- Repeats in Raw Data, mapping : Long reads from Gapless Ab10 assemly are aligned to the gapless Ab10 assembly using minimap2. Read hits are filtered to a single positions (2308) and averaged in 10kb windows.
https://github.com/dawelab/ETR_in_HiFi/blob/main/analysis/issues_old_assemblies.sh


2. Issues in NEW Assemblies
   **_Repeat Content Estimation._**
    - HiFi Reads Prep : Raw HiFi are converted to CCS in 10 chunks, then bam files are merged together, filtered to a minimum quality of 99%, and converted to a fastq.
      https://github.com/dawelab/ETR_in_HiFi/blob/main/Additional_assemblies/HiFi_reads.sh
    - HiFi Reads Repeats, Length, Coverage , Quality : CCS reads are compared to repeat consensus sequences using blast. Total repeat content is then calculated.Then, contigs are generated and the same process is repeated. Finally, the contig v satellite repeat output is used to check for anchored ends on contigs, defined as contigs that do not have a satellite hit within 100bp of start or end bp.
      https://github.com/dawelab/ETR_in_HiFi/blob/main/analysis/issues_new_asemblies.sh
      
   **_HiFi Assembly and Annotation._**
    - HiFi Contig Generation, summary, scaffolding, annotation : Contigs with > 1/2 the expected read depth are extracted and scaffolded on Mo17 T2T and B73v5 as an alternate. Assemblies annotated with Liftoff. Repeat content assessed.
      https://github.com/dawelab/ETR_in_HiFi/blob/main/analysis/scaffolding.sh
    - 	CG119 correction: The incorrectly placed contig is moved with bedtools and remerged into a new assembly.
      https://github.com/dawelab/ETR_in_HiFi/blob/main/analysis/CG119_correction.sh

3. Arrays
   **_Conserved Repeat Positions._**
   - Finding Arrays : Repeat hits within 10kb were merged using bedtools. Arrays with >=10 monomers and >= 10% repeat content were defined as repeat arrays. Array positions were then compared to 100N gaps using bedtools. Then, array positions were compared to core genes using bedtools closest.
     https://github.com/dawelab/ETR_in_HiFi/blob/main/analysis/array_positions.sh
   - Conserved Array Positions : Array positions were clustered together based on repeat type and shared up and down core genes.
    https://github.com/dawelab/ETR_in_HiFi/blob/main/analysis/Conserved_Arrays.R

4. Structure 
   **_Repeat Structure._**
   - Model : Classification model creation from sample data. LDA model chosen to further analysis.
     https://github.com/dawelab/ETR_in_HiFi/tree/main/Classification_Model
   - Raw reads : Raw repeat structure predicted from monomer similarity network toplogy using LDA model.
      https://github.com/dawelab/ETR_in_HiFi/blob/main/analysis/raw_structure.sh
   - Assembly structure : Assembly reoeat structure predicted using non-overlapping 10kb bins from monomer similarity network toplogy using LDA model.
     https://github.com/dawelab/ETR_in_HiFi/blob/main/analysis/array_structure.sh
