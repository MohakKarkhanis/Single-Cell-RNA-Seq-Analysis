#Single Cell RNA sequence Analysis: Human Peripeharal Blood Mononuclear Cells

Single cell RNA sequence analysis is the process wherein gene expression profiles of each individual cell is analysed in depth which helps us to later classify them into clusters. These clusters represent identically functioning cells like for examples cells involved in sugar metabolism, water regulatory genes etc. are grouped in an individual cluster. 

This repository contains the code and methodology for processing, analysing, and clustering Single-Cell RNA Sequencing (scRNA-seq) data to identify distinct cellular populations based on their gene expression profiles.

The [Steps_Followed.docx](Steps_Followed.docx) file provides a comprehensive, step-by-step walkthrough of the entire pipeline, including code snippets, intermediate visual plots (Violin, Elbow Plots and UMAP), and the reason for each decision. Due to GitHub's file size limits, the raw HDF5 matrix files are not hosted here, but the data acquistion process is fully documented. 

**Tools**
1. Language: Python
2. Core scRNA-seq library: `Scanpy `
3. Data Processing: `NumPy`, `Pandas`
4. Clustering and Graphing: `igraph`, `leidenalg`
5. Visualization: `Matplotlib`, `Seaborn`

**Methodology**
1. Data Acquisition: Downloaded the cell matrix HDF5 (filtered) h5 [**10k Human PBMCs, 3' v3.1, Chromium Controller**](https://cf.10xgenomics.com/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_Controller/10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.h5) dataset from 10x Genomics.
2. Quality Control (QC) & Preprocessing: Calculated QC metrics to identify dead/dying cells via mitochondrial gene expression (`pct_counts_mt`).
**Filtering**: Removed low-quality cells (genes by counts < 200 and >2500) and cells with high mitochondrial content (> 5/0%). Filtered out genes expressed in fewer than 3 cells.
**Normalization**: Normalized total counts per cell to 10,000 and applied a log1p transformation.
**Feature Selection**: Identified and subsetted Highly Variable Genes (HVGs) to focus on the most biologically informative markers, followed by data scaling.
3. Dimensionality Reduction: Applied PCA to capture the most variance. An Elbow Plot was utilized to determine the optimal number of principal components. Also, applied UMAP on the PCA-reduced data to non-linearly embed the cells into a 2D space while preserving local and global data structures.
4. Graph-Based Clustering: Computed a K-nearest neighbour (KNN) graph and applied Leiden Clustering algorithm with a `resolution=0.5` to detect communities of transcriptionally similar cells, effectively segmenting the UMAP plot into distinct clusters.
5. Marker Gene Discovery & Annotation: Utilized the Wilcoxon rank-sum test (`rank_genes_groups`) to extract the top 10 defining marker genes for each Leiden Cluster. Conducted a robust literature survey using **PangloDB** and the **The Human Protein Atlas** to map these specific marker genes in `Blood` and `Immune system` part of the databases.

**Results**
The pipeline successfully cleaned a noisy raw scRNA-seq dataset, isolating high quality, viable cells. The combination of PCA, UMAP, and Leiden Clustering effectively grouped the Human PBMCs (Peripheral Blood Mononuclear Cells).

**Steps to run the model**
1. Download the Streamlit folder and set its directory location in `Terminal`. Ensure streamlit is installed or manually install it using `pip install streamlit`.
2. Run the Python script in Terminal using the command `streamlit run rna_analysis.py`.
