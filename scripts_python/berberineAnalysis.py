#!/usr/bin/env python3
"""
Berberine Differential Expression Analysis
Author: Pichai Raman
Date: 7/7/2025
Description: Analysis of Berberine effects on gene expression
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

# Read in Berberine gene expression data
print("Reading Berberine gene expression data...")
data = pd.read_excel("data/GSE264072_1_genes_fpkm_expression.xlsx")

# Filter genes with low expression (counts < 10) and low FPKM (< 1)
print("Filtering low-expression genes...")
count_cols = [col for col in data.columns if 'count' in col.lower()]
fpkm_cols = [col for col in data.columns if 'fpkm' in col.lower()]

count_filter = data[count_cols].max(axis=1) > 10
fpkm_filter = data[fpkm_cols].max(axis=1) > 1
data_filt = data[count_filter & fpkm_filter].copy()
print(f"Filtered to {len(data_filt)} genes.")

# =============================================================================
# PRINCIPAL COMPONENT ANALYSIS (PCA)
# =============================================================================

print("Performing PCA analysis...")

# Extract count data for PCA and log-transform
count_cols_pca = [col for col in data_filt.columns if 'count' in col.lower()]
count_data = data_filt[count_cols_pca].values

# Set row names to gene identifiers (keep Ensembl IDs as row names)
gene_ids = data_filt.iloc[:, 0].values

# Log transform the count data (add 1 to avoid log(0))
log_count_data = np.log2(count_data + 1)

# Create sample information for PCA
sample_names_pca = count_cols_pca
sample_info_pca = pd.DataFrame({
    'sample': sample_names_pca,
    'treatment': ['Berberine' if 'thb' in name.lower() else 'Control' for name in sample_names_pca]
}, index=sample_names_pca)

# Perform PCA
scaler = StandardScaler()
log_count_data_scaled = scaler.fit_transform(log_count_data.T)
pca = PCA()
pca_result = pca.fit_transform(log_count_data_scaled)
var_explained = pca.explained_variance_ratio_ * 100

# Prepare data for plotting
pca_data = pd.DataFrame({
    'PC1': pca_result[:, 0],
    'PC2': pca_result[:, 1],
    'Sample': sample_names_pca,
    'Treatment': sample_info_pca['treatment'].values
})

# Create PCA plot
plt.figure(figsize=(12, 8))
colors = {'Berberine': '#E74C3C', 'Control': '#3498DB'}

for treatment in ['Control', 'Berberine']:
    mask = pca_data['Treatment'] == treatment
    plt.scatter(pca_data.loc[mask, 'PC1'], pca_data.loc[mask, 'PC2'], 
               c=colors[treatment], label=treatment, s=100, alpha=0.8)

# Add sample labels
for idx, row in pca_data.iterrows():
    plt.annotate(row['Sample'], (row['PC1'], row['PC2']), 
                xytext=(5, 5), textcoords='offset points', fontsize=8)

plt.xlabel(f'PC1 ({var_explained[0]:.1f}% variance)')
plt.ylabel(f'PC2 ({var_explained[1]:.1f}% variance)')
plt.title(f'PCA of Berberine Gene Expression Data (Log2 Counts)\n'
          f'PC1: {var_explained[0]:.1f}% variance, PC2: {var_explained[1]:.1f}% variance')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save the PCA plot
plt.savefig(os.path.join(output_folder, "PCA_Berberine.pdf"), dpi=300, bbox_inches='tight')
plt.close()

# Print PCA summary
print("PCA Analysis Summary (Log2 Counts):")
print("==================================")
print(f"Total samples: {len(pca_result)}")
print(f"Total genes: {log_count_data.shape[1]}")
print(f"PC1 variance explained: {var_explained[0]:.1f}%")
print(f"PC2 variance explained: {var_explained[1]:.1f}%")

# =============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

print("Setting up differential expression analysis...")

# Prepare count matrix for DESeq2
count_matrix = data_filt[count_cols].values
gene_ids = data_filt.iloc[:, 0].values

# Create sample information data frame
sample_names = count_cols
sample_info = pd.DataFrame({
    'sample': sample_names,
    'treatment': ['Berberine' if 'thb' in name.lower() else 'Control' for name in sample_names]
}, index=sample_names)

# Create sample metadata DataFrame for DESeq2
sample_metadata = pd.DataFrame({
    'condition': sample_info['treatment']
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

# Get results for Berberine vs Control comparison
stat_res = ds.DeseqStats(dds_obj, contrast=["condition", "Berberine", "Control"])
stat_res.summary()
res_berb_vs_control = stat_res.results_df

# Sort by adjusted p-value
res_berb_vs_control = res_berb_vs_control.sort_values('padj')

# Filter significant genes (padj < 0.05 and |log2FoldChange| > 0)
sig_berb_vs_control = res_berb_vs_control[
    (res_berb_vs_control['padj'] < 0.05) & 
    (np.abs(res_berb_vs_control['log2FoldChange']) > 0)
]

# Print summary statistics
print("\n=== BERBERINE DIFFERENTIAL EXPRESSION ANALYSIS SUMMARY ===")
print(f"Total genes analyzed: {len(res_berb_vs_control)}")
print(f"Significant genes (padj < 0.05 & |log2FC| > 1): {len(sig_berb_vs_control)}")
print(f"Up-regulated genes: {sum(sig_berb_vs_control['log2FoldChange'] > 0)}")
print(f"Down-regulated genes: {sum(sig_berb_vs_control['log2FoldChange'] < 0)}")

# =============================================================================
# GENE SYMBOL ANNOTATION
# =============================================================================

print("Adding mouse gene symbols to results...")

# Check what columns are available in the results
print(f"Available columns in results: {list(res_berb_vs_control.columns)}")

# Add mouse gene symbols for labeling (if available)
if 'gene_name' in data_filt.columns:
    # Match Ensembl IDs (row names) to gene names using the first column
    gene_name_map = data_filt.set_index(data_filt.columns[0])['gene_name']
    res_berb_vs_control['Mouse_Symbol'] = res_berb_vs_control.index.map(gene_name_map)
elif 'Gene' in data_filt.columns:
    gene_name_map = data_filt.set_index(data_filt.columns[0])['Gene']
    res_berb_vs_control['Mouse_Symbol'] = res_berb_vs_control.index.map(gene_name_map)
elif 'Symbol' in data_filt.columns:
    gene_name_map = data_filt.set_index(data_filt.columns[0])['Symbol']
    res_berb_vs_control['Mouse_Symbol'] = res_berb_vs_control.index.map(gene_name_map)
else:
    res_berb_vs_control['Mouse_Symbol'] = res_berb_vs_control.index

# Set gene_id as index for consistency with R output
res_berb_vs_control = res_berb_vs_control.reset_index().rename(columns={'index': 'gene_id'}).set_index('gene_id')

# Save annotated results
output_file = os.path.join(output_folder, "Berberine_vs_Control_DEGs_annotated.csv")
res_berb_vs_control.to_csv(output_file)
print(f"Results saved to: {output_file}")

# =============================================================================
# VISUALIZATION: VOLCANO PLOT
# =============================================================================

print("Creating volcano plot...")

# Prepare labels for volcano plot
volcano_labels = res_berb_vs_control['Mouse_Symbol'].fillna(pd.Series(res_berb_vs_control.index, index=res_berb_vs_control.index))

# Create volcano plot
plt.figure(figsize=(12, 8))

# Define significance threshold
sig_threshold = (res_berb_vs_control['padj'] < 0.05) & (np.abs(res_berb_vs_control['log2FoldChange']) > 1)

# Plot points
plt.scatter(res_berb_vs_control.loc[~sig_threshold, 'log2FoldChange'], 
           -np.log10(res_berb_vs_control.loc[~sig_threshold, 'padj']), 
           c='gray', alpha=0.6, s=20, label='Not significant')

plt.scatter(res_berb_vs_control.loc[sig_threshold, 'log2FoldChange'], 
           -np.log10(res_berb_vs_control.loc[sig_threshold, 'padj']), 
           c='red', alpha=0.8, s=30, label='Significant')

# Add threshold lines
plt.axhline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.7)
plt.axvline(-1, color='blue', linestyle='--', alpha=0.7)
plt.axvline(1, color='blue', linestyle='--', alpha=0.7)

# Add labels for top genes
top_genes = sig_threshold.sum()
up_genes = (sig_threshold & (res_berb_vs_control['log2FoldChange'] > 0)).sum()
down_genes = (sig_threshold & (res_berb_vs_control['log2FoldChange'] < 0)).sum()

plt.xlabel('Log2 Fold Change')
plt.ylabel('-log10(adjusted p-value)')
plt.title(f'Berberine vs Control Differential Expression\n'
          f'Significant genes: {top_genes} ({up_genes} up, {down_genes} down)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save volcano plot
plt.savefig(os.path.join(output_folder, "Volcano_Berberine_vs_Control.pdf"), dpi=300, bbox_inches='tight')
plt.close()

print("Volcano plot saved to: Volcano_Berberine_vs_Control.pdf")

# =============================================================================
# VISUALIZATION: HEATMAP OF TOP DIFFERENTIALLY EXPRESSED GENES
# =============================================================================

print("Creating heatmap of top differentially expressed genes...")

# Get top 50 differentially expressed genes (after setting index)
sig_berb_vs_control_indexed = res_berb_vs_control[
    (res_berb_vs_control['padj'] < 0.05) & 
    (np.abs(res_berb_vs_control['log2FoldChange']) > 1)
]
top_genes = sig_berb_vs_control_indexed.head(50).index.tolist()

# Prepare data for heatmap
plot_data = data_filt[data_filt.iloc[:, 0].isin(top_genes)].set_index(data_filt.columns[0])[count_cols]

# Log2 transform and z-score normalization
plot_data_log = np.log2(plot_data + 1)
plot_data_z = stats.zscore(plot_data_log, axis=1)

# Set row names to gene symbols for better readability
gene_symbols = res_berb_vs_control.loc[top_genes, 'Mouse_Symbol'].values
plot_data_z_df = pd.DataFrame(plot_data_z, index=gene_symbols, columns=count_cols)

# Create annotation for samples
annotation_col = pd.DataFrame({
    'Treatment': [sample_info.loc[col, 'treatment'] for col in count_cols]
}, index=count_cols)

# Create heatmap
plt.figure(figsize=(10, 12))

# Create color map for annotations
treatment_colors = {'Control': '#3498DB', 'Berberine': '#E74C3C'}
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

plt.title('Top 50 Differentially Expressed Genes: Berberine vs Control')
plt.tight_layout()

# Save heatmap
plt.savefig(os.path.join(output_folder, "Heatmap_TopDEGs_Berberine_vs_Control.pdf"), 
            dpi=300, bbox_inches='tight')
plt.close()

print("Heatmap saved to: Heatmap_TopDEGs_Berberine_vs_Control.pdf")

# =============================================================================
# ANALYSIS COMPLETION
# =============================================================================

print("\n=== BERBERINE ANALYSIS COMPLETE ===")
print("All results saved in 'output' folder:")
print("1. Berberine_vs_Control_DEGs_annotated.csv - Annotated differential expression results")
print("2. Volcano_Berberine_vs_Control.pdf - Volcano plot of differential expression")
print("3. Heatmap_TopDEGs_Berberine_vs_Control.pdf - Heatmap of top differentially expressed genes")
print("4. PCA_Berberine.pdf - PCA plot")
print("\nAnalysis completed successfully!") 
