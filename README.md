# Circadian Transcription Factors (CircTF)

The goal of this project is to evaluate computational approaches to identify conserved binding sites for transcription factors regulating circadian genes in Brassicaceae.

This file describes the scripts used as of summer of 2023. The scripts are organized by the order in which they were run.


============================================================
Script: 
    MetaCycle_v01.R 
------------------------------------------------------------
Input:
    - ../input/tpm_counts.txt
    - ../input/LIBRARIES_KEY.txt
    - ../input/JGI_syntenicHits/*.synHits
------------------------------------------------------------
Steps:
    - Preprocessing
        - Match counts to libraries keys
        - Keep 7 genotypes (L58|WO83|A03|VT123|Pcglu|O302V|R500)
        - Keep 7 timepoints (WTP1 - WTP7)
        - Add missing timepoints and reps as NA
    - MetaCycle analysis
        - Specify timepoints (17,21,25,29,33,37,41)
        - Run meta2d using JTK
        - Save each genotype as a separate txt file
    - Identifying circadian genes
        - List files from both meta2d and syntenic hits
        - Rename genotypes (WO83<-WO_83, Pcglu<-PCGglu, O302V<-O_302V)
        - Keep genes with FDR < 0.001
        - Match Metaclycle gene ID with syntenic hits' gene ID (except for R500)
        - Finish looping through all genotypes  
        - Remove duplicates from A03 (and maybe other genotypes?)
        - Save output to ../output/MetaCycle/CircadianGenes.csv
------------------------------------------------------------
Output:
    - ../input/MetaCycle/input.txt (per genotype)
    - ../output/MetaCycle/CircadianGenes.csv
============================================================




============================================================
Script: 
    get_2kb_fasta.sh
------------------------------------------------------------
Input:
    - ../input/Genomes/*.fa
------------------------------------------------------------
Steps:
    - Rename fasta files (for ease of use)
        - mv BrapaO_302V_711_v1.0.fa BrapaO_302V.fa
        - mv Brapassp_chinensisvar_communisPCGlu_712_v2.0.fa BrapaPCGlu.fa
        - mv Brapassp_chinensisvar_parachinensisL58_709_v1.0.fa BrapaL58.fa
        - mv Brapassp_oleiferavar_oleiferaWO_83_715_v1.0.fa BrapaWO_83.fa
        - mv Brapassp_pekinensisvar_pekinensisA03_708_v1.0.fa BrapaA03.fa
        - mv Brapassp_rapavar_rapaVT123_714_v1.0.fa BrapaVT123.fa
        - mv Brapassp_trilocularisR500_795_v2.0.fa BrapaR500.fa
    - Create index file (samtools faidx)
    - Get chromosome size
    - Rename gff3 file
    - Switch columns Name and Attributes
    - Flank the 2kb upstream sequence (bedtools flank)
    - Get 2kb sequence (bedtools getfasta)
    - Finish looping through all genotypes
------------------------------------------------------------
Output:
    - ./output/*.chrom.sizes
    - ./output/*.genes.gff3
    - ./output/*.promoters.bed
    - ./output/Regions/*.genes.2kb.promoters.bed
============================================================




============================================================
Script: 
    CircadianRegions.py
------------------------------------------------------------
Input:
    - ../output/GFF3_gene_only//*.gff3
    - ../output/Regions/*.2kb.promoters.bed
    - ../output/MetaCycle/CircadianGenes.csv
------------------------------------------------------------
Steps:
    - List gff3 and bd Files
    - Create dataframe with gene ID and promoter sequence
    - Match dataframe with gene ID from MetaCycle
    - Finish looping through all genotypes
    - Read in Metacycle output
    - Rename genes and merge with promoter sequence
    - Save combined output as ../output/CircadianRegions_2kb.csv
------------------------------------------------------------
Output:
    - ../output/CircadianRegions_2kb.csv
============================================================



