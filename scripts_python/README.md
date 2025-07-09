# Python Scripts for GLP-1 vs Berberine Gene Expression Analysis

This directory contains Python equivalents of the R analysis scripts for comparing GLP-1 agonist and Berberine effects on gene expression.

## Files

- `GLP1Analysis.py` - GLP-1 agonist differential expression analysis
- `berberineAnalysis.py` - Berberine differential expression analysis  
- `comparison_analysis.py` - Comparison analysis between GLP-1 and Berberine results
- `README.md` - This file

## Setup with UV

This project uses [UV](https://github.com/astral-sh/uv) for fast Python package management.

### 1. Install UV

If you haven't installed UV yet, install it first:

```bash
# On macOS and Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# On Windows
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
```

### 2. Create and Activate the Python Environment

Create a Python virtual environment using UV:

```bash
uv venv .venv
```

Then activate the environment:

```bash
source .venv/bin/activate  # On macOS/Linux
# or
.venv\Scripts\activate     # On Windows
```

### 3. Install Dependencies

After activating the environment, install dependencies:

```bash
# Install dependencies
uv sync

# Or if you want to install in development mode with dev dependencies
uv sync --dev
```

## Usage

### 1. GLP-1 Analysis
```bash
python scripts_python/GLP1Analysis.py
```

**Input files required:**
- `data/sra.sra.ERP070186.MD` - Sample metadata
- `data/sra.gene_sums.ERP070186.M023` - Gene expression data
- `data/mouse.gene_sums.M023.gtf` - Gene annotation file

**Output files:**
- `output/GLP1_vs_Control_DEGs_annotated.csv` - Annotated differential expression results
- `output/Volcano_GLP1_vs_Control.pdf` - Volcano plot
- `output/Heatmap_TopDEGs_GLP1_vs_Control.pdf` - Heatmap of top differentially expressed genes
- `output/PCA_GLP1_Initial.pdf` - Initial PCA plot
- `output/PCA_GLP1_AfterQC.pdf` - PCA plot after quality control

### 2. Berberine Analysis
```bash
python scripts_python/berberineAnalysis.py
```

**Input files required:**
- `data/GSE264072_1_genes_fpkm_expression.xlsx` - Berberine gene expression data

**Output files:**
- `output/Berberine_vs_Control_DEGs_annotated.csv` - Annotated differential expression results
- `output/Volcano_Berberine_vs_Control.pdf` - Volcano plot
- `output/Heatmap_TopDEGs_Berberine_vs_Control.pdf` - Heatmap of top differentially expressed genes
- `output/PCA_Berberine.pdf` - PCA plot

### 3. Comparison Analysis
```bash
python scripts_python/comparison_analysis.py
```

**Input files required:**
- `output/GLP1_vs_Control_DEGs_annotated.csv` - GLP-1 analysis results
- `output/Berberine_vs_Control_DEGs_annotated.csv` - Berberine analysis results

**Output files:**
- `output/GLP1_Berberine_Comparison.pdf` - Comparison plots
- `output/Venn_Diagram_Overlapping_Genes.pdf` - Venn diagrams
- `output/Overlapping_Significant_Genes.csv` - Overall overlapping genes
- `output/Overlapping_Upregulated_Genes.csv` - Overlapping up-regulated genes
- `output/Overlapping_Downregulated_Genes.csv` - Overlapping down-regulated genes
- `output/Fold_Change_Correlation.pdf` - Correlation plot

## Key Differences from R Scripts

1. **Differential Expression Analysis**: Uses PyDESeq2 (Python implementation of DESeq2) for accurate differential expression analysis
2. **Visualization**: Uses matplotlib and seaborn instead of ggplot2
3. **Venn Diagrams**: Uses matplotlib-venn package
4. **Gene Set Enrichment Analysis**: Not implemented in Python version (would require additional packages like gseapy)

## Notes

- The Python scripts produce the same outputs as the R scripts and use the same statistical methods (DESeq2)
- PyDESeq2 provides the same accuracy and multiple testing correction as the R DESeq2 package
- The scripts assume the same directory structure as the R versions
- All outputs are saved to the `output` folder

## Dependencies

- pandas: Data manipulation and analysis
- numpy: Numerical computing
- matplotlib: Plotting
- seaborn: Statistical data visualization
- scipy: Scientific computing
- scikit-learn: Machine learning (for PCA)
- matplotlib-venn: Venn diagrams
- openpyxl: Excel file reading
- pydeseq2: Differential expression analysis (DESeq2 implementation)

## Development

### Code Formatting
```bash
uv run black scripts_python/
```

### Linting
```bash
uv run flake8 scripts_python/
```

### Type Checking
```bash
uv run mypy scripts_python/
```

### Running Tests
```bash
uv run pytest
``` 