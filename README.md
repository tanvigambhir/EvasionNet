# EvasionNet: Spatial-Transcriptomic Microglial Evasion Networks

EvasionNet is an end-to-end single-nucleus RNA sequencing (snRNA-seq) pipeline designed to isolate and map the transcriptomic states of human microglia as they evade or respond to advancing Alzheimer's disease pathology. 

## Project Overview
Using spatial-transcriptomic data from the Allen Institute (Entorhinal Cortex), this project tracks the biological rewiring of immune cells across a continuous disease timeline (Braak Stages 0 -> 6). 

The pipeline identifies Disease-Associated Microglia (DAM) and mathematically isolates the primary gene drivers responsible for this state change, specifically highlighting the massive upregulation of non-coding RNAs like *MALAT1* and *ADAM7-AS1* alongside systemic mitochondrial gene collapse.

## Computational Architecture & Cross-Language Bridging
A major bottleneck in spatial-transcriptomics is format lock-in (e.g., Python-exclusive `.h5ad` architectures failing in R environments due to `SeuratDisk` dependency issues on ARM processors). 

To ensure a robust, reproducible pipeline, EvasionNet utilizes a cross-language extraction strategy:
1. **Python (`scanpy` / `pandas`):** Natively parses the raw `.h5ad` object, extracting and transposing the sparse counts matrix and metadata into universal CSVs.
2. **R (`Seurat`):** Rebuilds the object cleanly for high-performance Quality Control, PCA, UMAP dimensionality reduction, and Differential Gene Expression (DGE) analysis.

## Key Visualizations

### 1. Disease State Clustering (UMAP)
Microglial transcriptomes were reduced to 2D space. The visualization reveals a distinct, physically isolated "island" of Alzheimer's microglia (left) completely detached from the homeostatic core (center), proving a fundamental transcriptomic state change.
*(See `/plots/UMAP_DiseaseState.png`)*

### 2. Temporal Pathology Mapping (Ridge Plot)
By utilizing Braak Staging as a continuous pseudo-time variable, the pipeline maps the activation of the primary evasion marker, *MALAT1*. Expression remains dormant in healthy brains (Stage 0) and violently spikes as tau pathology reaches terminal severity (Stage 6).
*(See `/plots/RidgePlot_MALAT1.png`)*

## How to Run
1. Clone the repository.
2. Run `scripts/extract_h5ad.py` in your standard ML/Python environment to bypass dependency bottlenecks and generate the matrices.
3. Run `scripts/pipeline.R` in RStudio to execute the Seurat workflow, generate visualizations, and compile the final `.rds` object.

## Tech Stack
* **Languages:** Python, R
* **Libraries:** Scanpy, Pandas, Seurat, ggplot2
* **Biology:** snRNA-seq, Differential Gene Expression, Neuroimmunology, Alzheimer's Pathology (Braak Staging)
