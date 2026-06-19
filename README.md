# Analyze-Gene-Expression-Correlations
# Bioinformatics Sequence Analysis Pipeline

## Project Overview
This project provides an automated, modular framework for the **comparative genomic analysis** of gene sequences. By integrating bioinformatics tools with data science techniques, this pipeline converts raw FASTA-formatted sequences from the NCBI database into structured, statistically validated insights.

It is designed to assess how structural characteristics—specifically **GC content** and **sequence length**—correlate with the functional classification of genes within key biological pathways.

## Why This Pipeline Matters
In molecular biology, the sequence composition of a gene often reflects its functional constraints. This pipeline was designed to:
1. **Standardize Data Ingestion:** Automates the handling of non-uniform genomic data from NCBI datasets.
2. **Statistically Quantify Variation:** Uses **ANOVA (Analysis of Variance)** to determine if the differences in genomic composition between gene families are statistically significant ($p < 0.05$).
3. **Uncover Latent Patterns:** Employs **Principal Component Analysis (PCA)** to reduce the complexity of high-dimensional genomic features, revealing how genes cluster together in "feature space."
4. **Validate Biological Relevance:** Performs **Functional Enrichment Analysis** to verify if specific gene sets are statistically over-represented in biological processes (e.g., "Mitotic Cell Cycle").

## Technical Workflow
The pipeline follows a four-stage data lifecycle:

### 1. Data Processing & Feature Extraction
Utilizes `Biopython` to parse `.fna` files. The pipeline calculates:
* **Sequence Length:** Indicates the structural size of the gene.
* **GC Content:** $\text{GC\%} = \frac{\text{Count(G)} + \text{Count(C)}}{\text{Total Length}} \times 100$

### 2. Statistical Hypothesis Testing
To ensure findings are not due to noise, we perform an **ANOVA test** using `scipy.stats` to mathematically validate compositional differences between gene families.

### 3. Multivariate Dimensionality Reduction
Using `Scikit-Learn`, we apply a `StandardScaler` for feature normalization and project these features onto a 2D plane using **PCA**. This allows visual confirmation of how genes cluster by pathway.

### 4. Functional Annotation & Enrichment
Using `gseapy`, we query the **GO Biological Process database** to identify cellular roles where the input genes are most active.

## Getting Started

### Prerequisites
Ensure you have Python 3.x installed.

## Directory Layout
/data: Store your gene-specific folders here. Each must contain the standard NCBI ncbi_dataset/data/gene.fna structure.

/test: The output repository. The pipeline automatically creates this and stores generated enrichment_analysis.png plots and tabular reports.

## Execution
python Analyze_gene_expression_correlations.py

## Technologies Used
Bioinformatics: Biopython, GSEAPY, Bioservices
Data Science: Pandas, NumPy, Scikit-learn, SciPy
Visualization: Matplotlib, Seaborn
```bash

pip install -r requirements.txt
