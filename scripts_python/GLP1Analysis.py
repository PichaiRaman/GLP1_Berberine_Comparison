#!/usr/bin/env python3
"""
GLP-1 Agonist Differential Expression Analysis
Author: Pichai Raman
Date: 7/7/2025
Description: Analysis of GLP-1 agonist effects on gene expression
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pydeseq2 import preprocessing, dds, ds
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# SETUP OUTPUT FOLDER
# =============================================================================

# Create output folder
output_folder = "output"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    print(f"✓ Created output folder: {output_folder}")
else:
    print(f"✓ Output folder already exists: {output_folder}")

# =============================================================================
# DATA IMPORT AND PREPROCESSING
# =============================================================================

# Read in Sample Meta-Data
# Source: https://jhubiostatistics.shinyapps.io/recount3-study-explorer/
# Dataset: ERP070186 gene counts for gencode_v23
print("Reading sample metadata...")
sample_data = pd.read_csv("data/sra.sra.ERP070186.MD", sep='\t')
sample_data = sample_data[["external_id", "experiment_attributes"]]

# Assign treatment groups based on experiment attributes
sample_data["treatment"] = sample_data["experiment_attributes"].apply(
    lambda x: "Control" if "vehicle" in str(x).lower() else "GLP1"
)

# Read in gene expression data
print("Reading gene expression data...")
data = pd.read_csv("data/sra.gene_sums.ERP070186.M023", sep='\t', skiprows=2)
data = data[["gene_id"] + list(sample_data["external_id"])]

# =============================================================================
# QUALITY CONTROL AND FILTERING
# =============================================================================

# Filter genes with low expression (counts < 10 across all samples)
print("Filtering low-expression genes...")
count_cols = [col for col in data.columns if col != "gene_id"]
count_filter = data[count_cols].max(axis=1) > 10
data_filt = data[count_filter].copy()
print(f"Removed {sum(~count_filter)} low-expression genes")

# =============================================================================
# PRINCIPAL COMPONENT ANALYSIS (PCA)
# =============================================================================

print("Performing initial PCA analysis...")

# Extract count data for PCA and log-transform
count_data = data_filt[count_cols].values
log_count_data = np.log2(count_data + 1)  # Add 1 to avoid log(0)

# Perform PCA
scaler = StandardScaler()
log_count_data_scaled = scaler.fit_transform(log_count_data.T)
pca = PCA()
pca_result = pca.fit_transform(log_count_data_scaled)

# Calculate variance explained by each principal component
var_explained = pca.explained_variance_ratio_ * 100

# Create PCA plot using matplotlib
print("Creating initial PCA plot...")

# Prepare data for plotting
pca_data = pd.DataFrame({
    'PC1': pca_result[:, 0],
    'PC2': pca_result[:, 1],
    'Sample': sample_data["external_id"].values,
    'Treatment': sample_data["treatment"].values
})

# Create PCA plot
plt.figure(figsize=(12, 8))
colors = {'GLP1': '#E74C3C', 'Control': '#3498DB'}

for treatment in ['Control', 'GLP1']:
    mask = pca_data['Treatment'] == treatment
    plt.scatter(pca_data.loc[mask, 'PC1'], pca_data.loc[mask, 'PC2'], 
               c=colors[treatment], label=treatment, s=100, alpha=0.8)

# Add sample labels
for idx, row in pca_data.iterrows():
    plt.annotate(row['Sample'], (row['PC1'], row['PC2']), 
                xytext=(5, 5), textcoords='offset points', fontsize=8)

plt.xlabel(f'PC1 ({var_explained[0]:.1f}% variance)')
plt.ylabel(f'PC2 ({var_explained[1]:.1f}% variance)')
plt.title(f'PCA of GLP-1 Gene Expression Data (Log2 Counts)\n'
          f'PC1: {var_explained[0]:.1f}% variance, PC2: {var_explained[1]:.1f}% variance')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save the initial PCA plot
plt.savefig(os.path.join(output_folder, "PCA_GLP1_Initial.pdf"), dpi=300, bbox_inches='tight')
plt.close()

# =============================================================================
# SAMPLE QUALITY CONTROL - REMOVE OUTLIERS
# =============================================================================

print("Identifying and removing outlier samples...")

# Remove samples identified as outliers
bad_samples = ["ERR2104427", "ERR2104420"]
keep_samples = [s for s in sample_data["external_id"] if s not in bad_samples]
sample_data = sample_data[sample_data["external_id"].isin(keep_samples)].copy()

# Update data to include only good samples
data = data[["gene_id"] + keep_samples].copy()

# Re-filter data after removing outliers
count_cols = [col for col in data.columns if col != "gene_id"]
count_filter = data[count_cols].max(axis=1) > 10
data_filt = data[count_filter].copy()

# =============================================================================
# UPDATED PCA AFTER OUTLIER REMOVAL
# =============================================================================

print("Performing updated PCA after outlier removal...")

# Extract count data for updated PCA
count_data = data_filt[count_cols].values
log_count_data = np.log2(count_data + 1)

# Perform updated PCA
scaler = StandardScaler()
log_count_data_scaled = scaler.fit_transform(log_count_data.T)
pca = PCA()
pca_result = pca.fit_transform(log_count_data_scaled)
var_explained = pca.explained_variance_ratio_ * 100

# Prepare data for updated plotting
pca_data = pd.DataFrame({
    'PC1': pca_result[:, 0],
    'PC2': pca_result[:, 1],
    'Sample': sample_data["external_id"].values,
    'Treatment': sample_data["treatment"].values
})

# Create updated PCA plot
plt.figure(figsize=(12, 8))

for treatment in ['Control', 'GLP1']:
    mask = pca_data['Treatment'] == treatment
    plt.scatter(pca_data.loc[mask, 'PC1'], pca_data.loc[mask, 'PC2'], 
               c=colors[treatment], label=treatment, s=100, alpha=0.8)

# Add sample labels
for idx, row in pca_data.iterrows():
    plt.annotate(row['Sample'], (row['PC1'], row['PC2']), 
                xytext=(5, 5), textcoords='offset points', fontsize=8)

plt.xlabel(f'PC1 ({var_explained[0]:.1f}% variance)')
plt.ylabel(f'PC2 ({var_explained[1]:.1f}% variance)')
plt.title(f'PCA of GLP-1 Gene Expression Data (Log2 Counts) - After QC\n'
          f'PC1: {var_explained[0]:.1f}% variance, PC2: {var_explained[1]:.1f}% variance')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save the updated PCA plot
plt.savefig(os.path.join(output_folder, "PCA_GLP1_AfterQC.pdf"), dpi=300, bbox_inches='tight')
plt.close()

# =============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

print("Setting up differential expression analysis...")

# Prepare count matrix for DESeq2
count_matrix = data_filt[count_cols].values
gene_ids = data_filt["gene_id"].values

# Create sample information
sample_info = sample_data.set_index("external_id")["treatment"]

# Create sample metadata DataFrame for DESeq2
sample_metadata = pd.DataFrame({
    'condition': [sample_info[col] for col in count_cols]
}, index=count_cols)

# Create count matrix DataFrame - transpose so samples are rows and genes are columns
count_df = pd.DataFrame(count_matrix.T, index=count_cols, columns=gene_ids)

print("Running DESeq2 differential expression analysis...")

# Initialize DESeqDataSet
dds_obj = dds.DeseqDataSet(
    counts=count_df,
    metadata=sample_metadata,
    design_factors="condition"
)

# Run DESeq2 analysis
dds_obj.deseq2()

# Get results for GLP1 vs Control comparison
stat_res = ds.DeseqStats(dds_obj, contrast=["condition", "GLP1", "Control"])
stat_res.summary()
res_glp1_vs_control = stat_res.results_df

# Sort by adjusted p-value
res_glp1_vs_control = res_glp1_vs_control.sort_values('padj')

# Filter significant genes (padj < 0.05 and |log2FoldChange| > 0)
sig_glp1_vs_control = res_glp1_vs_control[
    (res_glp1_vs_control['padj'] < 0.05) & 
    (np.abs(res_glp1_vs_control['log2FoldChange']) > 0)
]

# Print summary statistics
print("\n=== GLP-1 DIFFERENTIAL EXPRESSION ANALYSIS SUMMARY ===")
print(f"Total genes analyzed: {len(res_glp1_vs_control)}")
print(f"Significant genes (padj < 0.05 & |log2FC| > 1): {len(sig_glp1_vs_control)}")
print(f"Up-regulated genes: {sum(sig_glp1_vs_control['log2FoldChange'] > 0)}")
print(f"Down-regulated genes: {sum(sig_glp1_vs_control['log2FoldChange'] < 0)}")

# =============================================================================
# GENE ANNOTATION FROM GTF FILE
# =============================================================================

print("Reading GTF file for gene annotation...")

# Read GTF file and extract Ensembl ID to gene name mapping
gtf_data = []
with open("data/mouse.gene_sums.M023.gtf", 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 9 and parts[2] == 'gene':
            # Extract gene_id and gene_name from attributes
            attrs = parts[8]
            gene_id = None
            gene_name = None
            
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1]
                elif attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1]
            
            if gene_id:
                gtf_data.append({'ensembl_id': gene_id, 'gene_name': gene_name})

gene_map = pd.DataFrame(gtf_data)
print("Gene annotation mapping complete.")

# =============================================================================
# GENE ANNOTATION AND OUTPUT
# =============================================================================

print("Adding gene annotations to results...")

# Reset index to gene_id for merging
res_glp1_vs_control = res_glp1_vs_control.reset_index().rename(columns={'index': 'gene_id'})

# Merge with gene_map on gene_id and ensembl_id
res_glp1_vs_control = res_glp1_vs_control.merge(
    gene_map, left_on='gene_id', right_on='ensembl_id', how='left'
)

# Fallback to Ensembl ID if no gene symbol found
res_glp1_vs_control['Mouse_Symbol'] = res_glp1_vs_control['gene_name'].fillna(res_glp1_vs_control['gene_id'])

# Set gene_id as index for consistency with R output
res_glp1_vs_control = res_glp1_vs_control.set_index('gene_id')

# Save annotated results
output_file = os.path.join(output_folder, "GLP1_vs_Control_DEGs_annotated.csv")
res_glp1_vs_control.to_csv(output_file)
print(f"Results saved to: {output_file}")

# =============================================================================
# VISUALIZATION: VOLCANO PLOT
# =============================================================================

print("Creating volcano plot...")

# Prepare labels for volcano plot
volcano_labels = res_glp1_vs_control['Mouse_Symbol'].fillna(pd.Series(res_glp1_vs_control.index, index=res_glp1_vs_control.index))

# Create volcano plot
plt.figure(figsize=(12, 8))

# Define significance threshold
sig_threshold = (res_glp1_vs_control['padj'] < 0.05) & (np.abs(res_glp1_vs_control['log2FoldChange']) > 1)

# Plot points
plt.scatter(res_glp1_vs_control.loc[~sig_threshold, 'log2FoldChange'], 
           -np.log10(res_glp1_vs_control.loc[~sig_threshold, 'padj']), 
           c='gray', alpha=0.6, s=20, label='Not significant')

plt.scatter(res_glp1_vs_control.loc[sig_threshold, 'log2FoldChange'], 
           -np.log10(res_glp1_vs_control.loc[sig_threshold, 'padj']), 
           c='red', alpha=0.8, s=30, label='Significant')

# Add threshold lines
plt.axhline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.7)
plt.axvline(-1, color='blue', linestyle='--', alpha=0.7)
plt.axvline(1, color='blue', linestyle='--', alpha=0.7)

# Add labels for top genes
top_genes = sig_threshold.sum()
up_genes = (sig_threshold & (res_glp1_vs_control['log2FoldChange'] > 0)).sum()
down_genes = (sig_threshold & (res_glp1_vs_control['log2FoldChange'] < 0)).sum()

plt.xlabel('Log2 Fold Change')
plt.ylabel('-log10(adjusted p-value)')
plt.title(f'GLP-1 vs Control Differential Expression\n'
          f'Significant genes: {top_genes} ({up_genes} up, {down_genes} down)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save volcano plot
plt.savefig(os.path.join(output_folder, "Volcano_GLP1_vs_Control.pdf"), dpi=300, bbox_inches='tight')
plt.close()

print("Volcano plot saved to: Volcano_GLP1_vs_Control.pdf")

# =============================================================================
# VISUALIZATION: HEATMAP OF TOP DIFFERENTIALLY EXPRESSED GENES
# =============================================================================

print("Creating heatmap of top differentially expressed genes...")

# Get top 50 differentially expressed genes (after setting index)
sig_glp1_vs_control_indexed = res_glp1_vs_control[
    (res_glp1_vs_control['padj'] < 0.05) & 
    (np.abs(res_glp1_vs_control['log2FoldChange']) > 1)
]
top_genes = sig_glp1_vs_control_indexed.head(50).index.tolist()

# Prepare data for heatmap
plot_data = data_filt[data_filt['gene_id'].isin(top_genes)].set_index('gene_id')[count_cols]

# Log2 transform and z-score normalization
plot_data_log = np.log2(plot_data + 1)
plot_data_z = stats.zscore(plot_data_log, axis=1)

# Set row names to gene symbols for better readability
gene_symbols = res_glp1_vs_control.loc[top_genes, 'Mouse_Symbol'].values
plot_data_z_df = pd.DataFrame(plot_data_z, index=gene_symbols, columns=count_cols)

# Create annotation for samples
annotation_col = pd.DataFrame({
    'Treatment': [sample_info[col] for col in count_cols]
}, index=count_cols)

# Create heatmap
plt.figure(figsize=(10, 12))

# Create color map for annotations
treatment_colors = {'Control': '#3498DB', 'GLP1': '#E74C3C'}
annotation_colors = annotation_col['Treatment'].map(treatment_colors)

# Create heatmap
sns.clustermap(plot_data_z_df, 
               cmap='RdBu_r', 
               center=0,
               col_colors=annotation_colors,
               xticklabels=True,
               yticklabels=True,
               figsize=(10, 12),
               cbar_kws={'label': 'Z-score'})

plt.title('Top 50 Differentially Expressed Genes: GLP-1 vs Control')
plt.tight_layout()

# Save heatmap
plt.savefig(os.path.join(output_folder, "Heatmap_TopDEGs_GLP1_vs_Control.pdf"), 
            dpi=300, bbox_inches='tight')
plt.close()

print("Heatmap saved to: Heatmap_TopDEGs_GLP1_vs_Control.pdf")

# =============================================================================
# ANALYSIS COMPLETION
# =============================================================================

print("\n=== GLP-1 ANALYSIS COMPLETE ===")
print("All results saved in 'output' folder:")
print("1. GLP1_vs_Control_DEGs_annotated.csv - Annotated differential expression results")
print("2. Volcano_GLP1_vs_Control.pdf - Volcano plot of differential expression")
print("3. Heatmap_TopDEGs_GLP1_vs_Control.pdf - Heatmap of top differentially expressed genes")
print("4. PCA_GLP1_Initial.pdf - Initial PCA plot")
print("5. PCA_GLP1_AfterQC.pdf - PCA plot after quality control")
print("\nAnalysis completed successfully!") 
