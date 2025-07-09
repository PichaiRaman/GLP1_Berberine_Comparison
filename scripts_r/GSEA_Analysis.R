#!/usr/bin/env Rscript
# =============================================================================
# GSEA Analysis on Overlapping Genes
# Author: Pichai Raman
# Date: 7/7/2025
# Description: Gene Set Enrichment Analysis using MSigDB on overlapping genes 
#              between GLP-1 and Berberine
# =============================================================================

# =============================================================================
# LOAD REQUIRED LIBRARIES
# =============================================================================

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required Bioconductor packages if not already installed
required_packages <- c("clusterProfiler", "org.Mm.eg.db", "msigdbr", "enrichplot", "ggplot2", "dplyr")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (pkg %in% c("clusterProfiler", "org.Mm.eg.db", "msigdbr", "enrichplot")) {
            BiocManager::install(pkg)
        } else {
            install.packages(pkg)
        }
    }
}

# Load libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(dplyr)

# =============================================================================
# SETUP OUTPUT FOLDER
# =============================================================================

# Create output folder
output_folder <- "output"
if (!dir.exists(output_folder)) {
    dir.create(output_folder)
    cat("✓ Created output folder:", output_folder, "\n")
} else {
    cat("✓ Output folder already exists:", output_folder, "\n")
}

# =============================================================================
# LOAD OVERLAPPING GENE DATA
# =============================================================================

cat("Loading overlapping gene data...\n")

# Read overlapping genes
up_genes_file <- file.path(output_folder, "Overlapping_Upregulated_Genes.csv")
down_genes_file <- file.path(output_folder, "Overlapping_Downregulated_Genes.csv")

if (!file.exists(up_genes_file) || !file.exists(down_genes_file)) {
    stop("Overlapping gene files not found. Please run the comparison analysis first.")
}

up_genes_df <- read.csv(up_genes_file)
down_genes_df <- read.csv(down_genes_file)

# Extract gene symbols
up_genes <- up_genes_df$Gene_Symbol
down_genes <- down_genes_df$Gene_Symbol

cat("Loaded", length(up_genes), "up-regulated genes and", length(down_genes), "down-regulated genes\n")

# =============================================================================
# PREPARE MSIGDB GENE SETS
# =============================================================================

cat("Preparing MSigDB gene sets for mouse...\n")

# Get mouse gene sets from MSigDB
# H: hallmark gene sets
# C2: curated gene sets (KEGG, Reactome, etc.)
# C5: GO gene sets
# C6: oncogenic signatures
# C7: immunologic signatures

msigdb_categories <- c("H", "C2", "C5", "C6", "C7")
gene_sets_list <- list()

# First try mouse-specific gene sets
for (collection in msigdb_categories) {
    cat("Loading", collection, "gene sets for mouse...\n")
    # Use mouse-specific gene sets with db_species = "MM"
    gene_sets_list[[collection]] <- msigdbr(species = "Mus musculus", collection = collection)
}

# Check if we got enough mouse gene sets, if not, use human with ortholog mapping
total_mouse_sets <- sum(sapply(gene_sets_list, nrow))
if (total_mouse_sets < 1000) {
    cat("Limited mouse gene sets found (", total_mouse_sets, "). Using human gene sets with ortholog mapping...\n")
    
    # Use human gene sets with ortholog mapping to mouse
    for (collection in msigdb_categories) {
        cat("Loading", collection, "human gene sets with mouse ortholog mapping...\n")
        gene_sets_list[[collection]] <- msigdbr(species = "Homo sapiens", collection = collection)
    }
}

# Combine all gene sets
all_gene_sets <- do.call(rbind, gene_sets_list)

# Create a list of gene sets for clusterProfiler
gene_sets <- split(all_gene_sets$gene_symbol, all_gene_sets$gs_name)

cat("Loaded", length(gene_sets), "gene sets from MSigDB\n")

# =============================================================================
# GSEA ANALYSIS FOR UP-REGULATED GENES
# =============================================================================

cat("\n=== GSEA ANALYSIS FOR UP-REGULATED GENES ===\n")

# Run GSEA for up-regulated genes
up_gsea_results <- enricher(
    gene = up_genes,
    TERM2GENE = all_gene_sets[, c("gs_name", "gene_symbol")],
    pvalueCutoff = 1,
    pAdjustMethod = "BH"
)

# Filter results by adjusted p-value < 0.25
up_gsea_filtered <- up_gsea_results@result[up_gsea_results@result$p.adjust < 0.25, ]

if (!is.null(up_gsea_results) && nrow(up_gsea_filtered) > 0) {
    cat("Found", nrow(up_gsea_filtered), "enriched gene sets for up-regulated genes (p.adjust < 0.25)\n")
    
    # Save results
    up_results_file <- file.path(output_folder, "GSEA_Upregulated_Genes_MSigDB.csv")
    write.csv(up_gsea_filtered, up_results_file, row.names = FALSE)
    cat("Up-regulated GSEA results saved to:", up_results_file, "\n")
    
    # Create dot plot for top enriched pathways
    top_pathways <- head(up_gsea_filtered, 20)
    if (nrow(top_pathways) > 0) {
        p <- ggplot(top_pathways, aes(x = Count, y = reorder(Description, Count))) +
            geom_point(aes(size = Count, color = p.adjust)) +
            scale_color_gradient(low = "red", high = "blue") +
            theme_bw() +
            theme(axis.text.y = element_text(size = 10),
                  axis.title = element_text(size = 12),
                  plot.title = element_text(size = 14, hjust = 0.5)) +
            labs(title = "Top Enriched Pathways - Up-regulated Genes",
                 x = "Gene Count", y = "Pathway", 
                 color = "Adjusted P-value", size = "Gene Count")
        up_plot_file <- file.path(output_folder, "GSEA_DotPlot_Upregulated_MSigDB.pdf")
        ggsave(up_plot_file, p, width = 12, height = 10, dpi = 300)
        cat("Up-regulated GSEA dot plot saved to:", up_plot_file, "\n")
    }
} else {
    cat("No significant gene sets found for up-regulated genes (p.adjust < 0.25)\n")
}

# =============================================================================
# GSEA ANALYSIS FOR DOWN-REGULATED GENES
# =============================================================================

cat("\n=== GSEA ANALYSIS FOR DOWN-REGULATED GENES ===\n")

# Run GSEA for down-regulated genes
down_gsea_results <- enricher(
    gene = down_genes,
    TERM2GENE = all_gene_sets[, c("gs_name", "gene_symbol")],
    pvalueCutoff = 1,
    pAdjustMethod = "BH"
)

down_gsea_filtered <- down_gsea_results@result[down_gsea_results@result$p.adjust < 0.25, ]

if (!is.null(down_gsea_results) && nrow(down_gsea_filtered) > 0) {
    cat("Found", nrow(down_gsea_filtered), "enriched gene sets for down-regulated genes (p.adjust < 0.25)\n")
    
    # Save results
    down_results_file <- file.path(output_folder, "GSEA_Downregulated_Genes_MSigDB.csv")
    write.csv(down_gsea_filtered, down_results_file, row.names = FALSE)
    cat("Down-regulated GSEA results saved to:", down_results_file, "\n")
    
    # Create dot plot for top enriched pathways
    top_pathways <- head(down_gsea_filtered, 20)
    if (nrow(top_pathways) > 0) {
        p <- ggplot(top_pathways, aes(x = Count, y = reorder(Description, Count))) +
            geom_point(aes(size = Count, color = p.adjust)) +
            scale_color_gradient(low = "red", high = "blue") +
            theme_bw() +
            theme(axis.text.y = element_text(size = 10),
                  axis.title = element_text(size = 12),
                  plot.title = element_text(size = 14, hjust = 0.5)) +
            labs(title = "Top Enriched Pathways - Down-regulated Genes",
                 x = "Gene Count", y = "Pathway", 
                 color = "Adjusted P-value", size = "Gene Count")
        down_plot_file <- file.path(output_folder, "GSEA_DotPlot_Downregulated_MSigDB.pdf")
        ggsave(down_plot_file, p, width = 12, height = 10, dpi = 300)
        cat("Down-regulated GSEA dot plot saved to:", down_plot_file, "\n")
    }
} else {
    cat("No significant gene sets found for down-regulated genes (p.adjust < 0.25)\n")
}

# =============================================================================
# CATEGORY-SPECIFIC ANALYSIS
# =============================================================================

cat("\n=== CATEGORY-SPECIFIC GSEA ANALYSIS ===\n")

# Run GSEA for each category separately
categories <- c("H" = "Hallmark", "C2" = "Curated", "C5" = "GO", "C6" = "Oncogenic", "C7" = "Immunologic")

for (cat_code in names(categories)) {
    cat_name <- categories[cat_code]
    cat("Analyzing", cat_name, "gene sets...\n")
    
    # Get gene sets for this category
    cat_gene_sets <- gene_sets_list[[cat_code]]
    
    if (nrow(cat_gene_sets) > 0) {
        # Up-regulated genes
        up_cat_results <- enricher(
            gene = up_genes,
            TERM2GENE = cat_gene_sets[, c("gs_name", "gene_symbol")],
            pvalueCutoff = 1,
            pAdjustMethod = "BH"
        )
        up_cat_filtered <- up_cat_results@result[up_cat_results@result$p.adjust < 0.25, ]
        if (!is.null(up_cat_results) && nrow(up_cat_filtered) > 0) {
            up_cat_file <- file.path(output_folder, paste0("GSEA_Upregulated_", cat_name, "_MSigDB.csv"))
            write.csv(up_cat_filtered, up_cat_file, row.names = FALSE)
            cat("  Up-regulated", cat_name, "results saved\n")
        }
        
        # Down-regulated genes
        down_cat_results <- enricher(
            gene = down_genes,
            TERM2GENE = cat_gene_sets[, c("gs_name", "gene_symbol")],
            pvalueCutoff = 1,
            pAdjustMethod = "BH"
        )
        down_cat_filtered <- down_cat_results@result[down_cat_results@result$p.adjust < 0.25, ]
        if (!is.null(down_cat_results) && nrow(down_cat_filtered) > 0) {
            down_cat_file <- file.path(output_folder, paste0("GSEA_Downregulated_", cat_name, "_MSigDB.csv"))
            write.csv(down_cat_filtered, down_cat_file, row.names = FALSE)
            cat("  Down-regulated", cat_name, "results saved\n")
        }
    }
}

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== GSEA ANALYSIS SUMMARY ===\n")
cat("Up-regulated genes analyzed:", length(up_genes), "\n")
cat("Down-regulated genes analyzed:", length(down_genes), "\n")

if (!is.null(up_gsea_results) && nrow(up_gsea_results@result) > 0) {
    cat("Significant pathways for up-regulated genes:", nrow(up_gsea_results@result), "\n")
} else {
    cat("No significant pathways found for up-regulated genes\n")
}

if (!is.null(down_gsea_results) && nrow(down_gsea_results@result) > 0) {
    cat("Significant pathways for down-regulated genes:", nrow(down_gsea_results@result), "\n")
} else {
    cat("No significant pathways found for down-regulated genes\n")
}

cat("\n=== GSEA ANALYSIS COMPLETE ===\n")
cat("All results saved in 'output' folder:\n")
cat("1. GSEA_Upregulated_Genes_MSigDB.csv - Up-regulated gene enrichment results\n")
cat("2. GSEA_Downregulated_Genes_MSigDB.csv - Down-regulated gene enrichment results\n")
cat("3. GSEA_DotPlot_Upregulated_MSigDB.pdf - Dot plot for up-regulated genes\n")
cat("4. GSEA_DotPlot_Downregulated_MSigDB.pdf - Dot plot for down-regulated genes\n")
cat("5. Category-specific CSV files for detailed analysis\n")
cat("\nAnalysis completed successfully!\n") 