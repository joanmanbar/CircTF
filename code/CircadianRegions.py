#!/usr/bin/env python


# Import libraries
import pandas as pd
import glob



# List gff3 and bed files
gff3_files = glob.glob('../output/GFF3_gene_only//*.gff3', recursive=False)
bed_files = glob.glob('../output/Regions/*.2kb.promoters.bed', recursive=False)

# CReate empty dataframe
All_regions = pd.DataFrame()
# List the seven genotypes
genotypes = ['L58','R500','WO_83','PCGlu','VT123','A03','O_302V']

for g in genotypes:
    
    # Current genotype (cg)
    cg = g
    
    # Match genotype with gff3 file
    gff3_genes = next(x for x in gff3_files if cg in x)
    df_gff3 = pd.read_csv(gff3_genes, sep="\t", header=None, comment="#") # read
    # Define colnames
    gff3_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    df_gff3.columns = gff3_names
    # Get gene name
    df_gff3['gene'] = df_gff3['attributes'].apply(lambda x: x.split('Name=')[1])

    # Match genotype with bed file 
    region_genes = next(x for x in bed_files if cg in x)
    df_region = pd.read_csv(region_genes, sep="\t", header=None, comment="#") # read
    region_names = ["attributes", "region2kb"] # colnames
    df_region.columns = region_names
    # Get position
    df_region['position'] = df_region['attributes'].apply(lambda x: x.split('::')[1])
    # Get gene name
    df_region['gene'] = df_region['attributes'].apply(lambda x: x.split('::')[0])
    df_region['gene'] = df_region['gene'].apply(lambda x: x.split('Name=')[1])
    df_region = df_region.drop('attributes', axis=1) # remove unnecessary column

    # Merge gff3 with bed file
    merged_df = pd.merge(df_gff3, df_region, on='gene')
    
    # Append to df
    All_regions = pd.concat([All_regions,merged_df])
    


# Circadian genes

# Read output form MetaCycle
CircadianGenes = '../output/MetaCycle/CircadianGenes.csv'
CircadianGenes = pd.read_csv(CircadianGenes)
# Get gene name to match with promoter
CircadianGenes = CircadianGenes.rename(columns={'NewCycID': 'gene'})
# Merge files
CircRegions_2kb = pd.merge(CircadianGenes,All_regions, on='gene')
# Write as csv (probably not the ideal format)
CircRegions_2kb.to_csv('../output/CircadianRegions_2kb.csv')