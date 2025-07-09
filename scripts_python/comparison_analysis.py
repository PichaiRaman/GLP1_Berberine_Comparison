#!/usr/bin/env python3
"""
GLP-1 vs Berberine Differential Expression Comparison
Author: Pichai Raman
Date: 2025
Description: Compare results from GLP-1 and Berberine differential expression analyses
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import hypergeom
import matplotlib_venn as venn
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
# DATA LOADING
# =============================================================================

print("Loading differential expression results...")

# Read GLP-1 results
glp1_data = None
try:
    glp1_data = pd.read_csv("output/GLP1_vs_Control_DEGs_annotated.csv", index_col=0)
    print(f"✓ GLP-1 data loaded: {len(glp1_data)} genes, {len(glp1_data.columns)} columns")
except FileNotFoundError:
    print("✗ output/GLP1_vs_Control_DEGs_annotated.csv not found")

# Read Berberine results
berberine_data = None
try:
    berberine_data = pd.read_csv("output/Berberine_vs_Control_DEGs_annotated.csv", index_col=0)
    print(f"✓ Berberine data loaded: {len(berberine_data)} genes, {len(berberine_data.columns)} columns")
except FileNotFoundError:
    print("✗ output/Berberine_vs_Control_DEGs_annotated.csv not found")

# =============================================================================
# DATA EXPLORATION
# =============================================================================

print("\n" + "="*60)
print("DATA EXPLORATION SUMMARY")
print("="*60)

if glp1_data is not None:
    print("\nGLP-1 vs Control Analysis:")
    print(f"  Total genes: {len(glp1_data)}")
    print(f"  Significant genes (padj < 0.05): {sum(glp1_data['padj'] < 0.05)}")
    print(f"  Up-regulated: {sum((glp1_data['padj'] < 0.05) & (glp1_data['log2FoldChange'] > 0))}")
    print(f"  Down-regulated: {sum((glp1_data['padj'] < 0.05) & (glp1_data['log2FoldChange'] < 0))}")
    print(f"  Columns: {', '.join(glp1_data.columns)}")
    
    print("\n  First 5 rows:")
    print(glp1_data.head())

if berberine_data is not None:
    print("\nBerberine vs Control Analysis:")
    print(f"  Total genes: {len(berberine_data)}")
    print(f"  Significant genes (padj < 0.05): {sum(berberine_data['padj'] < 0.05)}")
    print(f"  Up-regulated: {sum((berberine_data['padj'] < 0.05) & (berberine_data['log2FoldChange'] > 0))}")
    print(f"  Down-regulated: {sum((berberine_data['padj'] < 0.05) & (berberine_data['log2FoldChange'] < 0))}")
    print(f"  Columns: {', '.join(berberine_data.columns)}")
    
    print("\n  First 5 rows:")
    print(berberine_data.head())

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_significant_genes(data, padj_threshold=0.05, log2fc_threshold=0):
    """Get significant genes based on thresholds"""
    if data is None:
        return pd.DataFrame()
    
    significant = data[
        (data['padj'] < padj_threshold) & 
        (np.abs(data['log2FoldChange']) > log2fc_threshold)
    ].sort_values('padj')
    
    return significant

def get_upregulated_genes(data, padj_threshold=0.05, log2fc_threshold=0):
    """Get up-regulated genes"""
    if data is None:
        return pd.DataFrame()
    
    upregulated = data[
        (data['padj'] < padj_threshold) & 
        (data['log2FoldChange'] > log2fc_threshold)
    ].sort_values('padj')
    
    return upregulated

def get_downregulated_genes(data, padj_threshold=0.05, log2fc_threshold=0):
    """Get down-regulated genes"""
    if data is None:
        return pd.DataFrame()
    
    downregulated = data[
        (data['padj'] < padj_threshold) & 
        (data['log2FoldChange'] < -log2fc_threshold)
    ].sort_values('padj')
    
    return downregulated

def perform_hypergeometric_test(overlap_count, set1_count, set2_count, total_genes):
    """Perform hypergeometric test"""
    # Hypergeometric test parameters:
    # k = overlap_count (number of successes in sample)
    # M = total_genes (population size)
    # n = set1_count (number of successes in population)
    # N = set2_count (sample size)
    
    p_value = hypergeom.sf(overlap_count - 1, total_genes, set1_count, set2_count)
    return p_value

# =============================================================================
# GENE COMPARISON ANALYSIS
# =============================================================================

print("\n" + "="*60)
print("SIGNIFICANT GENES COMPARISON")
print("="*60)

# Get significant genes for each treatment
glp1_sig = get_significant_genes(glp1_data)
berberine_sig = get_significant_genes(berberine_data)

print(f"\nGLP-1 significant genes: {len(glp1_sig)}")
print(f"Berberine significant genes: {len(berberine_sig)}")

# Get up-regulated and down-regulated genes
glp1_up = get_upregulated_genes(glp1_data)
glp1_down = get_downregulated_genes(glp1_data)
berberine_up = get_upregulated_genes(berberine_data)
berberine_down = get_downregulated_genes(berberine_data)

print(f"\nGLP-1 up-regulated: {len(glp1_up)}")
print(f"GLP-1 down-regulated: {len(glp1_down)}")
print(f"Berberine up-regulated: {len(berberine_up)}")
print(f"Berberine down-regulated: {len(berberine_down)}")

# Calculate total number of genes for hypergeometric test
total_genes = max(len(glp1_data) if glp1_data is not None else 0, 
                 len(berberine_data) if berberine_data is not None else 0)

if len(glp1_sig) > 0 and len(berberine_sig) > 0:
    # Get gene symbols for overlap analysis
    glp1_symbols = glp1_sig['Mouse_Symbol'].dropna()
    glp1_symbols = glp1_symbols[glp1_symbols != '']
    berberine_symbols = berberine_sig['Mouse_Symbol'].dropna()
    berberine_symbols = berberine_symbols[berberine_symbols != '']
    
    # Get up-regulated and down-regulated gene symbols
    glp1_up_symbols = glp1_up['Mouse_Symbol'].dropna()
    glp1_up_symbols = glp1_up_symbols[glp1_up_symbols != '']
    glp1_down_symbols = glp1_down['Mouse_Symbol'].dropna()
    glp1_down_symbols = glp1_down_symbols[glp1_down_symbols != '']
    berberine_up_symbols = berberine_up['Mouse_Symbol'].dropna()
    berberine_up_symbols = berberine_up_symbols[berberine_up_symbols != '']
    berberine_down_symbols = berberine_down['Mouse_Symbol'].dropna()
    berberine_down_symbols = berberine_down_symbols[berberine_down_symbols != '']
    
    # Find overlapping genes
    overlap_genes = set(glp1_symbols) & set(berberine_symbols)
    glp1_only = set(glp1_symbols) - set(berberine_symbols)
    berberine_only = set(berberine_symbols) - set(glp1_symbols)
    
    # Find overlapping up-regulated genes
    overlap_up = set(glp1_up_symbols) & set(berberine_up_symbols)
    
    # Find overlapping down-regulated genes
    overlap_down = set(glp1_down_symbols) & set(berberine_down_symbols)
    
    print("\n=== OVERALL OVERLAP ANALYSIS ===")
    print(f"Overlapping significant genes: {len(overlap_genes)}")
    print(f"GLP-1 only: {len(glp1_only)}")
    print(f"Berberine only: {len(berberine_only)}")
    
    print("\n=== UP-REGULATED GENES OVERLAP ===")
    print(f"Overlapping up-regulated genes: {len(overlap_up)}")
    print(f"GLP-1 up only: {len(set(glp1_up_symbols) - set(berberine_up_symbols))}")
    print(f"Berberine up only: {len(set(berberine_up_symbols) - set(glp1_up_symbols))}")
    
    print("\n=== DOWN-REGULATED GENES OVERLAP ===")
    print(f"Overlapping down-regulated genes: {len(overlap_down)}")
    print(f"GLP-1 down only: {len(set(glp1_down_symbols) - set(berberine_down_symbols))}")
    print(f"Berberine down only: {len(set(berberine_down_symbols) - set(glp1_down_symbols))}")
    
    # Perform hypergeometric tests
    print("\n=== HYPERGEOMETRIC TESTS ===")
    
    # Overall overlap test
    if len(overlap_genes) > 0:
        overall_p = perform_hypergeometric_test(
            len(overlap_genes), 
            len(glp1_symbols), 
            len(berberine_symbols), 
            total_genes
        )
        print(f"Overall overlap p-value: {overall_p:.2e}")
    
    # Up-regulated overlap test
    if len(overlap_up) > 0:
        up_p = perform_hypergeometric_test(
            len(overlap_up), 
            len(glp1_up_symbols), 
            len(berberine_up_symbols), 
            total_genes
        )
        print(f"Up-regulated overlap p-value: {up_p:.2e}")
    
    # Down-regulated overlap test
    if len(overlap_down) > 0:
        down_p = perform_hypergeometric_test(
            len(overlap_down), 
            len(glp1_down_symbols), 
            len(berberine_down_symbols), 
            total_genes
        )
        print(f"Down-regulated overlap p-value: {down_p:.2e}")
    
    # Show some overlapping genes
    if len(overlap_genes) > 0:
        print("\nSample of overlapping genes (first 10):")
        overlap_sample = list(overlap_genes)[:10]
        
        for gene in overlap_sample:
            # Find the gene in both datasets
            glp1_idx = glp1_sig[glp1_sig['Mouse_Symbol'] == gene].index
            berberine_idx = berberine_sig[berberine_sig['Mouse_Symbol'] == gene].index
            
            if len(glp1_idx) > 0 and len(berberine_idx) > 0:
                glp1_fc = glp1_sig.loc[glp1_idx[0], 'log2FoldChange']
                berberine_fc = berberine_sig.loc[berberine_idx[0], 'log2FoldChange']
                
                print(f"  {gene}: GLP-1 FC={glp1_fc:.2f}, Berberine FC={berberine_fc:.2f}")

# =============================================================================
# VISUALIZATION
# =============================================================================

print("\nCreating comparison plots...")

# Create comparison plots
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# Plot 1: Log2FoldChange distributions
if glp1_data is not None and berberine_data is not None:
    axes[0, 0].hist(glp1_data['log2FoldChange'], bins=50, alpha=0.7, 
                    color='#E74C3C', label='GLP-1', density=True)
    axes[0, 0].hist(berberine_data['log2FoldChange'], bins=50, alpha=0.7, 
                    color='#3498DB', label='Berberine', density=True)
    axes[0, 0].set_title('Distribution of Log2 Fold Changes')
    axes[0, 0].set_xlabel('Log2 Fold Change')
    axes[0, 0].set_ylabel('Density')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

# Plot 2: P-value distributions
if glp1_data is not None and berberine_data is not None:
    axes[0, 1].hist(-np.log10(glp1_data['padj']), bins=50, alpha=0.7, 
                    color='#E74C3C', label='GLP-1', density=True)
    axes[0, 1].hist(-np.log10(berberine_data['padj']), bins=50, alpha=0.7, 
                    color='#3498DB', label='Berberine', density=True)
    axes[0, 1].set_title('Distribution of P-values')
    axes[0, 1].set_xlabel('-log10(adjusted p-value)')
    axes[0, 1].set_ylabel('Density')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

# Plot 3: GLP-1 volcano plot
if glp1_data is not None:
    sig_glp1 = (glp1_data['padj'] < 0.05) & (np.abs(glp1_data['log2FoldChange']) > 1)
    axes[1, 0].scatter(glp1_data.loc[~sig_glp1, 'log2FoldChange'], 
                      -np.log10(glp1_data.loc[~sig_glp1, 'padj']), 
                      c='gray', alpha=0.6, s=20)
    axes[1, 0].scatter(glp1_data.loc[sig_glp1, 'log2FoldChange'], 
                      -np.log10(glp1_data.loc[sig_glp1, 'padj']), 
                      c='red', alpha=0.8, s=30)
    axes[1, 0].axhline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.7)
    axes[1, 0].axvline(-1, color='blue', linestyle='--', alpha=0.7)
    axes[1, 0].axvline(1, color='blue', linestyle='--', alpha=0.7)
    axes[1, 0].set_title('GLP-1 vs Control')
    axes[1, 0].set_xlabel('Log2 Fold Change')
    axes[1, 0].set_ylabel('-log10(adjusted p-value)')
    axes[1, 0].grid(True, alpha=0.3)

# Plot 4: Berberine volcano plot
if berberine_data is not None:
    sig_berb = (berberine_data['padj'] < 0.05) & (np.abs(berberine_data['log2FoldChange']) > 1)
    axes[1, 1].scatter(berberine_data.loc[~sig_berb, 'log2FoldChange'], 
                      -np.log10(berberine_data.loc[~sig_berb, 'padj']), 
                      c='gray', alpha=0.6, s=20)
    axes[1, 1].scatter(berberine_data.loc[sig_berb, 'log2FoldChange'], 
                      -np.log10(berberine_data.loc[sig_berb, 'padj']), 
                      c='red', alpha=0.8, s=30)
    axes[1, 1].axhline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.7)
    axes[1, 1].axvline(-1, color='blue', linestyle='--', alpha=0.7)
    axes[1, 1].axvline(1, color='blue', linestyle='--', alpha=0.7)
    axes[1, 1].set_title('Berberine vs Control')
    axes[1, 1].set_xlabel('Log2 Fold Change')
    axes[1, 1].set_ylabel('-log10(adjusted p-value)')
    axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(output_folder, "GLP1_Berberine_Comparison.pdf"), 
            dpi=300, bbox_inches='tight')
plt.close()

print("✓ Comparison plots saved as 'GLP1_Berberine_Comparison.pdf'")

# =============================================================================
# VENN DIAGRAM
# =============================================================================

if len(glp1_sig) > 0 and len(berberine_sig) > 0:
    print("\nCreating Venn diagrams...")
    
    # Get gene symbols for Venn diagrams
    glp1_symbols = glp1_sig['Mouse_Symbol'].dropna()
    glp1_symbols = glp1_symbols[glp1_symbols != '']
    berberine_symbols = berberine_sig['Mouse_Symbol'].dropna()
    berberine_symbols = berberine_symbols[berberine_symbols != '']
    
    # Overall overlap
    overlap_genes = set(glp1_symbols) & set(berberine_symbols)
    
    # Create Venn diagrams
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Overall overlap Venn diagram
    venn.venn2([set(glp1_symbols), set(berberine_symbols)], 
               set_labels=('GLP-1', 'Berberine'), ax=axes[0])
    axes[0].set_title('Overall Significant Genes')
    
    # Up-regulated overlap Venn diagram
    glp1_up_symbols = glp1_up['Mouse_Symbol'].dropna()
    glp1_up_symbols = glp1_up_symbols[glp1_up_symbols != '']
    berberine_up_symbols = berberine_up['Mouse_Symbol'].dropna()
    berberine_up_symbols = berberine_up_symbols[berberine_up_symbols != '']
    overlap_up = set(glp1_up_symbols) & set(berberine_up_symbols)
    
    venn.venn2([set(glp1_up_symbols), set(berberine_up_symbols)], 
               set_labels=('GLP-1', 'Berberine'), ax=axes[1])
    axes[1].set_title('Up-regulated Genes')
    
    # Down-regulated overlap Venn diagram
    glp1_down_symbols = glp1_down['Mouse_Symbol'].dropna()
    glp1_down_symbols = glp1_down_symbols[glp1_down_symbols != '']
    berberine_down_symbols = berberine_down['Mouse_Symbol'].dropna()
    berberine_down_symbols = berberine_down_symbols[berberine_down_symbols != '']
    overlap_down = set(glp1_down_symbols) & set(berberine_down_symbols)
    
    venn.venn2([set(glp1_down_symbols), set(berberine_down_symbols)], 
               set_labels=('GLP-1', 'Berberine'), ax=axes[2])
    axes[2].set_title('Down-regulated Genes')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, "Venn_Diagram_Overlapping_Genes.pdf"), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✓ Venn diagrams saved as 'Venn_Diagram_Overlapping_Genes.pdf'")

# =============================================================================
# SAVE OVERLAPPING GENES
# =============================================================================

if len(glp1_sig) > 0 and len(berberine_sig) > 0:
    print("\nSaving overlapping genes...")
    
    # Get gene symbols
    glp1_symbols = glp1_sig['Mouse_Symbol'].dropna()
    glp1_symbols = glp1_symbols[glp1_symbols != '']
    berberine_symbols = berberine_sig['Mouse_Symbol'].dropna()
    berberine_symbols = berberine_symbols[berberine_symbols != '']
    overlap_genes = set(glp1_symbols) & set(berberine_symbols)
    
    if len(overlap_genes) > 0:
        # Create comparison dataframe for overall overlap
        comparison_data = []
        
        for gene in overlap_genes:
            # Find in GLP-1 data
            glp1_idx = glp1_sig[glp1_sig['Mouse_Symbol'] == gene].index
            berberine_idx = berberine_sig[berberine_sig['Mouse_Symbol'] == gene].index
            
            glp1_fc = glp1_sig.loc[glp1_idx[0], 'log2FoldChange'] if len(glp1_idx) > 0 else np.nan
            glp1_padj = glp1_sig.loc[glp1_idx[0], 'padj'] if len(glp1_idx) > 0 else np.nan
            berberine_fc = berberine_sig.loc[berberine_idx[0], 'log2FoldChange'] if len(berberine_idx) > 0 else np.nan
            berberine_padj = berberine_sig.loc[berberine_idx[0], 'padj'] if len(berberine_idx) > 0 else np.nan
            
            comparison_data.append({
                'Gene_Symbol': gene,
                'GLP1_log2FC': glp1_fc,
                'GLP1_padj': glp1_padj,
                'Berberine_log2FC': berberine_fc,
                'Berberine_padj': berberine_padj
            })
        
        comparison_df = pd.DataFrame(comparison_data)
        comparison_df = comparison_df.sort_values('GLP1_padj')
        
        comparison_df.to_csv(os.path.join(output_folder, "Overlapping_Significant_Genes.csv"), 
                            index=False)
        print(f"✓ Saved {len(comparison_df)} overlapping genes to 'Overlapping_Significant_Genes.csv'")
    
    # Save up-regulated overlapping genes
    glp1_up_symbols = glp1_up['Mouse_Symbol'].dropna()
    glp1_up_symbols = glp1_up_symbols[glp1_up_symbols != '']
    berberine_up_symbols = berberine_up['Mouse_Symbol'].dropna()
    berberine_up_symbols = berberine_up_symbols[berberine_up_symbols != '']
    overlap_up = set(glp1_up_symbols) & set(berberine_up_symbols)
    
    if len(overlap_up) > 0:
        up_comparison_data = []
        
        for gene in overlap_up:
            glp1_idx = glp1_up[glp1_up['Mouse_Symbol'] == gene].index
            berberine_idx = berberine_up[berberine_up['Mouse_Symbol'] == gene].index
            
            glp1_fc = glp1_up.loc[glp1_idx[0], 'log2FoldChange'] if len(glp1_idx) > 0 else np.nan
            glp1_padj = glp1_up.loc[glp1_idx[0], 'padj'] if len(glp1_idx) > 0 else np.nan
            berberine_fc = berberine_up.loc[berberine_idx[0], 'log2FoldChange'] if len(berberine_idx) > 0 else np.nan
            berberine_padj = berberine_up.loc[berberine_idx[0], 'padj'] if len(berberine_idx) > 0 else np.nan
            
            up_comparison_data.append({
                'Gene_Symbol': gene,
                'GLP1_log2FC': glp1_fc,
                'GLP1_padj': glp1_padj,
                'Berberine_log2FC': berberine_fc,
                'Berberine_padj': berberine_padj
            })
        
        up_comparison_df = pd.DataFrame(up_comparison_data)
        up_comparison_df = up_comparison_df.sort_values('GLP1_padj')
        
        up_comparison_df.to_csv(os.path.join(output_folder, "Overlapping_Upregulated_Genes.csv"), 
                               index=False)
        print(f"✓ Saved {len(up_comparison_df)} overlapping up-regulated genes to 'Overlapping_Upregulated_Genes.csv'")
    
    # Save down-regulated overlapping genes
    glp1_down_symbols = glp1_down['Mouse_Symbol'].dropna()
    glp1_down_symbols = glp1_down_symbols[glp1_down_symbols != '']
    berberine_down_symbols = berberine_down['Mouse_Symbol'].dropna()
    berberine_down_symbols = berberine_down_symbols[berberine_down_symbols != '']
    overlap_down = set(glp1_down_symbols) & set(berberine_down_symbols)
    
    if len(overlap_down) > 0:
        down_comparison_data = []
        
        for gene in overlap_down:
            glp1_idx = glp1_down[glp1_down['Mouse_Symbol'] == gene].index
            berberine_idx = berberine_down[berberine_down['Mouse_Symbol'] == gene].index
            
            glp1_fc = glp1_down.loc[glp1_idx[0], 'log2FoldChange'] if len(glp1_idx) > 0 else np.nan
            glp1_padj = glp1_down.loc[glp1_idx[0], 'padj'] if len(glp1_idx) > 0 else np.nan
            berberine_fc = berberine_down.loc[berberine_idx[0], 'log2FoldChange'] if len(berberine_idx) > 0 else np.nan
            berberine_padj = berberine_down.loc[berberine_idx[0], 'padj'] if len(berberine_idx) > 0 else np.nan
            
            down_comparison_data.append({
                'Gene_Symbol': gene,
                'GLP1_log2FC': glp1_fc,
                'GLP1_padj': glp1_padj,
                'Berberine_log2FC': berberine_fc,
                'Berberine_padj': berberine_padj
            })
        
        down_comparison_df = pd.DataFrame(down_comparison_data)
        down_comparison_df = down_comparison_df.sort_values('GLP1_padj')
        
        down_comparison_df.to_csv(os.path.join(output_folder, "Overlapping_Downregulated_Genes.csv"), 
                                 index=False)
        print(f"✓ Saved {len(down_comparison_df)} overlapping down-regulated genes to 'Overlapping_Downregulated_Genes.csv'")

# =============================================================================
# CORRELATION ANALYSIS
# =============================================================================

if len(glp1_sig) > 0 and len(berberine_sig) > 0:
    print("\nPerforming correlation analysis...")
    
    # Get gene symbols for correlation analysis
    glp1_symbols = glp1_sig['Mouse_Symbol'].dropna()
    glp1_symbols = glp1_symbols[glp1_symbols != '']
    berberine_symbols = berberine_sig['Mouse_Symbol'].dropna()
    berberine_symbols = berberine_symbols[berberine_symbols != '']
    overlap_genes = set(glp1_symbols) & set(berberine_symbols)
    
    if len(overlap_genes) > 0:
        # Get fold changes for overlapping genes
        glp1_fc = []
        berberine_fc = []
        
        for gene in overlap_genes:
            glp1_idx = glp1_sig[glp1_sig['Mouse_Symbol'] == gene].index
            berberine_idx = berberine_sig[berberine_sig['Mouse_Symbol'] == gene].index
            
            if len(glp1_idx) > 0 and len(berberine_idx) > 0:
                glp1_fc.append(glp1_sig.loc[glp1_idx[0], 'log2FoldChange'])
                berberine_fc.append(berberine_sig.loc[berberine_idx[0], 'log2FoldChange'])
        
        # Calculate correlation
        correlation = np.corrcoef(glp1_fc, berberine_fc)[0, 1]
        print(f"Pearson correlation between GLP-1 and Berberine fold changes: {correlation:.3f}")
        
        # Create correlation plot
        plt.figure(figsize=(8, 8))
        
        plt.scatter(glp1_fc, berberine_fc, alpha=0.8, c='#2E86AB', s=50)
        
        # Add regression line
        z = np.polyfit(glp1_fc, berberine_fc, 1)
        p = np.poly1d(z)
        plt.plot(glp1_fc, p(glp1_fc), "r--", linewidth=2)
        
        # Add diagonal line
        min_val = min(min(glp1_fc), min(berberine_fc))
        max_val = max(max(glp1_fc), max(berberine_fc))
        plt.plot([min_val, max_val], [min_val, max_val], 'gray', linestyle='--', alpha=0.7)
        
        plt.xlabel('GLP-1 Log2 Fold Change')
        plt.ylabel('Berberine Log2 Fold Change')
        plt.title(f'Fold Change Correlation\nr = {correlation:.3f}')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        plt.savefig(os.path.join(output_folder, "Fold_Change_Correlation.pdf"), 
                    dpi=300, bbox_inches='tight')
        plt.close()
        
        print("✓ Correlation plot saved as 'Fold_Change_Correlation.pdf'")

# =============================================================================
# ANALYSIS COMPLETION
# =============================================================================

print("\n" + "="*60)
print("ANALYSIS COMPLETE")
print("="*60)
print("All results saved in 'output' folder:")
print("- GLP1_Berberine_Comparison.pdf: Comparison plots")
print("- Venn_Diagram_Overlapping_Genes.pdf: Venn diagrams (overall, up-regulated, down-regulated)")
print("- Overlapping_Significant_Genes.csv: Overall overlapping genes")
print("- Overlapping_Upregulated_Genes.csv: Overlapping up-regulated genes")
print("- Overlapping_Downregulated_Genes.csv: Overlapping down-regulated genes")
print("- Fold_Change_Correlation.pdf: Correlation plot of fold changes")
print("\nAnalysis completed successfully!") 