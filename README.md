# GLP1_Berberine_Comparison

## Project Idea

The core idea of this project is to investigate whether **berberine**, a natural compound often referred to as *“nature’s Ozempic,”* exerts molecular effects similar to those of **GLP-1 receptor agonists (GLP‑1RA)** — a class of FDA-approved therapeutics used in the treatment of obesity and type 2 diabetes.

Berberine has gained increasing attention for its metabolic benefits, including improved glycemic control, lipid regulation, and modest weight loss. Anecdotally and mechanistically, it has been compared to pharmacological GLP‑1 agonists like semaglutide, which powerfully regulate appetite, metabolism, and insulin sensitivity. However, it remains unclear whether **these two interventions share common transcriptomic signatures**, particularly at the level of gene expression programs in metabolically active tissues like adipose.

To explore this, I’m performing a **comparative transcriptomic analysis** using two publicly available RNA-seq datasets:

* **GSE264072**: A study profiling the effects of **tetrahydroberberrubine (THBru)**, a berberine derivative, on white adipose tissue (WAT) in mice with diet-induced obesity.
* **ERP070186**: In this study, **female ob/ob mice** were treated with either **liraglutide (a GLP-1 receptor agonist)** or vehicle via **subcutaneous osmotic mini-pumps** for 14 days to assess transcriptomic changes in peri-gonadal white adipose tissue (WAT). 

By analyzing differential gene expression, pathway enrichment, and key regulator activity, this project aims to identify:

* Which **genes** and **biological pathways** are similarly or uniquely affected by each treatment
* If berberine truly mirrors aspects of GLP‑1 action at the transcriptomic level — or if the similarity is more anecdotal than molecular

Ultimately, this comparison could reveal whether **berberine-based therapies** engage similar cellular programs as GLP‑1 pharmacotherapies, offering insight into their potential as **natural metabolic modulators** and informing future research in drug repurposing and combination strategies.


This project contains analysis scripts for comparing GLP-1 agonist and Berberine effects on gene expression. Both R and Python versions are available, along with Gene Set Enrichment Analysis (GSEA) using MSigDB.

## Project Structure

```
GLP1_Berberine_Analysis/
├── data/                          # Input data files
├── scripts_r/                     # R analysis scripts
│   └── GSEA_Analysis.R
├── scripts_python/                # Python analysis scripts
│   ├── GLP1Analysis.py
│   ├── berberineAnalysis.py
│   ├── comparison_analysis.py
│   └── README.md
├── output/                        # Analysis results (created by scripts)
├── pyproject.toml                 # Python project configuration
├── .python-version               # Python version specification
└── README.md                     # This file
```

## Quick Start

### Python Version (Recommended)

1. **Install UV** (if not already installed):
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   uv venv .venv
   source .venv/bin/activate
   ```

2. **Install dependencies**:
   ```bash
   uv sync
   ```

3. **Run analyses**:
   ```bash
   # GLP-1 analysis
   python scripts_python/GLP1Analysis.py
   
   # Berberine analysis
   python scripts_python/berberineAnalysis.py
   
   # Comparison analysis
   python scripts_python/comparison_analysis.py
   
   # GSEA analysis (R script, requires comparison analysis to be run first)
   Rscript scripts_r/GSEA_Analysis.R
   ```

## Python Environment Management with UV

This project uses [UV](https://github.com/astral-sh/uv) for fast, reproducible Python environment management. All dependencies are specified in `pyproject.toml` and locked in `uv.lock`.

- The environment is fully reproducible using the `uv.lock` file.
- To add or update dependencies, edit `pyproject.toml` and run `uv pip install <package>`.
- For development dependencies, use `uv sync --dev`.

## Analysis Workflow

1. **GLP-1 Analysis** (`GLP1Analysis.R` / `GLP1Analysis.py`)
   - Reads GLP-1 gene expression data
   - Performs quality control and outlier removal
   - Conducts differential expression analysis using DESeq2/PyDESeq2
   - Creates PCA plots, volcano plots, and heatmaps
   - Saves annotated results

2. **Berberine Analysis** (`berberineAnalysis.R` / `berberineAnalysis.py`)
   - Reads Berberine gene expression data
   - Performs differential expression analysis using DESeq2/PyDESeq2
   - Creates PCA plots, volcano plots, and heatmaps
   - Saves annotated results

3. **Comparison Analysis** (`comparison_analysis.R` / `comparison_analysis.py`)
   - Compares results from both treatments
   - Identifies overlapping genes (up and down-regulated separately)
   - Performs statistical tests (hypergeometric tests)
   - Creates Venn diagrams and correlation plots
   - Saves overlapping gene lists

4. **GSEA Analysis** (`GSEA_Analysis.R`)
   - Performs Gene Set Enrichment Analysis on overlapping genes
   - Uses MSigDB gene sets (mouse-specific with fallback to human orthologs)
   - Analyzes multiple categories: Hallmark, Curated, GO, Oncogenic, Immunologic
   - Creates enrichment dot plots and detailed CSV results
   - Requires comparison analysis to be completed first

## Input Data Requirements

### GLP-1 Analysis
- `data/sra.sra.ERP070186.MD` - Sample metadata
- `data/sra.gene_sums.ERP070186.M023` - Gene expression data
- `data/mouse.gene_sums.M023.gtf` - Gene annotation file

### Berberine Analysis
- `data/GSE264072_1_genes_fpkm_expression.xlsx` - Berberine gene expression data

## Output Files

All analyses save results to the `output/` folder:

### Differential Expression Results
- **GLP1_vs_Control_DEGs_annotated.csv** - GLP-1 differential expression results
- **Berberine_vs_Control_DEGs_annotated.csv** - Berberine differential expression results

### Visualizations
- **PCA plots** - Principal component analysis plots
- **Volcano plots** - Differential expression volcano plots
- **Heatmaps** - Top differentially expressed genes heatmaps
- **Venn diagrams** - Overlapping genes visualization
- **Correlation plots** - Fold change correlation between treatments

### Overlap Analysis
- **Overlapping_Upregulated_Genes.csv** - Genes up-regulated in both treatments
- **Overlapping_Downregulated_Genes.csv** - Genes down-regulated in both treatments
- **Overlapping_Significant_Genes.csv** - All overlapping significant genes

### GSEA Results
- **GSEA_Upregulated_Genes_MSigDB.csv** - Up-regulated gene enrichment results
- **GSEA_Downregulated_Genes_MSigDB.csv** - Down-regulated gene enrichment results
- **GSEA_DotPlot_Upregulated_MSigDB.pdf** - Dot plot for up-regulated genes
- **GSEA_DotPlot_Downregulated_MSigDB.pdf** - Dot plot for down-regulated genes
- **Category-specific CSV files** - Detailed results for each MSigDB category

## Key Features

- **Quality Control**: Outlier detection and removal
- **Differential Expression**: Statistical analysis using DESeq2/PyDESeq2
- **Gene Annotation**: Mouse gene symbol mapping
- **Overlap Analysis**: Identification of common genes between treatments
- **Statistical Testing**: Hypergeometric tests for overlap significance
- **Gene Set Enrichment**: MSigDB-based pathway analysis
- **Visualization**: Comprehensive plotting of results
- **Mouse-Specific Analysis**: Optimized for mouse gene expression data

## MSigDB Categories Used in GSEA

- **H (Hallmark)**: Well-defined biological states and processes
- **C2 (Curated)**: Gene sets from various sources including KEGG, Reactome
- **C5 (GO)**: Gene Ontology biological processes, molecular functions, cellular components
- **C6 (Oncogenic)**: Gene sets representing signatures of cellular perturbations
- **C7 (Immunologic)**: Gene sets representing cell states and perturbations within the immune system

## Development

### Python Development
```bash
# Install dev dependencies
uv sync --dev

# Code formatting
uv run black scripts_python/

# Linting
uv run flake8 scripts_python/

# Type checking
uv run mypy scripts_python/
```

### R Development
- Use RStudio or your preferred R IDE
- Install additional packages as needed from CRAN or Bioconductor
- For GSEA analysis, ensure Bioconductor packages are installed:
  ```r
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "msigdbr", "enrichplot"))
  ```

## Conclusions

Based on a **hypergeometric test**, the overlap between the **berberine** and **GLP‑1RA** transcriptomic responses is **statistically significant**, indicating that these treatments share more differentially expressed genes than would be expected by chance. Moreover, pathway enrichment analysis highlights several **common biological processes**, notably **GOMF\_MAGNESIUM\_ION\_TRANSMEMBRANE\_TRANSPORTER\_ACTIVITY** and **GOBP\_MAGNESIUM\_ION\_TRANSPORT**, which intriguingly intersect with known **GLP‑1 receptor signaling pathways**.


### Magnesium-Related Pathways & GLP‑1 Signaling

1. **Magnesium Ion Transport Activity (GOMF\_MAGNESIUM\_ION\_TRANSMEMBRANE\_TRANSPORTER\_ACTIVITY & GOBP\_MAGNESIUM\_ION\_TRANSPORT)**
   These pathways regulate cellular **Mg²⁺ homeostasis**, which is critical for enzymatic reactions, mitochondrial function, ion channel operations, and signal transduction ([shemed.co.uk][1]).

2. **GLP‑1/Electrolyte Interactions**
   GLP‑1 receptor agonists are known to promote **renal sodium excretion (natriuresis)** and can influence **electrolyte balance**, including magnesium, via altered fluid regulation and kidney function ([shemed.co.uk][1]).

3. **Magnesium & GLP‑1 Expression**
   In rodent models, **magnesium supplementation** alongside metformin was shown to **increase GLP‑1 levels and expression**, with concurrent improvements in glucose metabolism ([medcraveonline.com][2]).

---

### Why It Matters

The overlap of our **shared transcriptomic genes**, enriched in magnesium transport pathways, **suggests that both berberine and GLP‑1RA might influence magnesium handling** within adipose tissue. Given that magnesium:

* Acts as a **cofactor for ATP**, enzymes, and ion pumps
* Modulates **insulin sensitivity**, **mitochondrial efficiency**, and **hormone secretion** ([medcraveonline.com][2])
* May be involved in the **fluid and electrolyte dynamics** observed with GLP‑1RA ([shemed.co.uk][1])

This raises an intriguing mechanistic hypothesis: **Mg²⁺ transport regulation could be a shared mediator of metabolic effects triggered by both treatments.**


### Conclusion

The unexpected connection between magnesium transport and both transcriptional responses **merits further investigation**—it may reveal a novel axis linking natural compounds and pharmacologic GLP‑1 agonists in metabolic regulation.

[1]: https://www.shemed.co.uk/blog/electrolyte-balance-during-glp-1-assisted-weight-loss?utm_source=chatgpt.com "Electrolyte Balance During GLP-1-Assisted Weight Loss - SheMed"
[2]: https://medcraveonline.com/IJMBOA/oral-magnesium-supplementation-modulates-hepatic-and-intestinal-expression-of-some-carbohydrate-metabolizing-genes-in-type-2-diabetic-rats.html?utm_source=chatgpt.com "Oral magnesium supplementation modulates hepatic and intestinal ..."

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this analysis pipeline in your research, please cite:

```
Raman, P. (2025). GLP-1 vs Berberine Gene Expression Analysis. 
GitHub repository. https://github.com/username/GLP1_Berberine_Analysis
``` 
