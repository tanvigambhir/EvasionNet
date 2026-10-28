# EvasionNet: Microglial States Across Alzheimer's Disease Progression

EvasionNet is an end-to-end single-nucleus RNA sequencing (snRNA-seq) analysis pipeline that tracks how human microglia, the brain's resident immune cells, change as Alzheimer's disease pathology advances across Braak stages 0 to 6.

## Project Overview

- **Question:** How do microglial transcriptional states shift as Alzheimer's pathology progresses?
- **Data:** Publicly available Allen Institute human brain data from the entorhinal cortex, one of the earliest regions affected in Alzheimer's disease. The source dataset is part of a spatial transcriptomics resource; this analysis uses the single-nucleus expression data.
- **Approach:** Quality control, dimensionality reduction (PCA, UMAP), clustering, and differential gene expression, with Braak stage treated as a continuous pseudo-time variable for disease progression.

## Key Observations

- Alzheimer's-associated microglia form a distinct cluster on the UMAP, separate from the homeostatic microglial core (`UMAP_DiseaseState.png`).
- Expression of the long non-coding RNAs **MALAT1** and **ADAM7-AS1** increases at later Braak stages (`RidgePlot_MALAT1.png`).
- Mitochondrial gene expression decreases with disease stage.

## Workflow

The pipeline is split across Python and R to avoid SeuratDisk `.h5ad` conversion issues on ARM (Apple Silicon) machines.

1. **`extract_h5ad.ipynb` (Python):** Loads the `.h5ad` object with Scanpy and exports the sparse counts matrix and cell metadata to CSV with Pandas.
2. **`ProxiMap-AD.R` (R):** Rebuilds the Seurat object from the CSVs and runs quality control, normalization, PCA, UMAP, clustering, and differential expression. Generates the figures and saves the final `.rds` object.

## Figures

| File | Description |
|---|---|
| `QC_ViolinPlot.png` | Quality-control metrics per nucleus |
| `UMAP_DiseaseState.png` | UMAP of microglia colored by disease state |
| `FeaturePlot_Markers.png` | Expression of microglial marker genes |
| `RidgePlot_MALAT1.png` | MALAT1 expression across Braak stages |

## How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/tanvigambhir/EvasionNet.git
   cd EvasionNet
   ```
2. Download the Allen Institute `.h5ad` file and update the file path in `extract_h5ad.ipynb`.
3. Run `extract_h5ad.ipynb` to export the counts matrix and metadata to CSV.
4. Open `ProxiMap-AD.R` in RStudio and run it to generate the figures and the final `.rds` object.

## Limitations and Next Steps

- **MALAT1** is highly abundant in nuclear RNA, so its signal in snRNA-seq should be treated as a hypothesis requiring orthogonal validation.
- **Mitochondrial gene expression** is also used as a quality-control metric in single-nucleus data; the observed decrease should be checked against QC filtering before biological interpretation.
- **Next steps:** incorporate spatial coordinates to test whether disease-associated microglia localize near pathology, and validate candidate genes in an independent cohort.

## Tech Stack

Python (Scanpy, Pandas) · R (Seurat, ggplot2) · snRNA-seq · Differential Expression · Neuroimmunology · Braak Staging
