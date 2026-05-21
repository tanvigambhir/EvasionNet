# EvasionNet: Single-Cell RNA-Seq Pipeline for Microglial Alzheimer's Pathology
# Author: Tanvi Gambhir


# Loading Libraries
library(Seurat)

# Load the CSVs that were just made
counts <- read.csv("/Users/tanvigambhir/Documents/neuro-ml/brain-spatial-data/ec_counts.csv", row.names = 1)
metadata <- read.csv("/Users/tanvigambhir/Documents/neuro-ml/brain-spatial-data/ec_metadata.csv", row.names = 1)

# Initialize the Seurat Object
ec_microglia <- CreateSeuratObject(counts = counts, meta.data = metadata, project = "ProxiMap-AD")

# View the metadata
head(ec_microglia@meta.data)

# 1. Check the balance of Disease States (exact number of Normal vs. Alzheimer's cells)
print("Disease Breakdown:")
table(ec_microglia$disease)

# 2. Check the breakdown of Pathology progression (how many cells we have at each stage of the disease)
print("Braak Stage Breakdown:")
table(ec_microglia$BraakStage)

# 3. Check the Microglial Subclusters
print("Microglial Substates:")
table(ec_microglia$seurat.subclusters)

# 4. Generate the Standard Quality Control (QC) Plot (creates a Violin Plot comparing the number of genes found in normal vs AD microglia)
VlnPlot(ec_microglia, features = c("nFeature_RNA", "nCount_RNA"), group.by = "disease", pt.size = 0.1)

# 5. The QC Filter (filter doublets (>1500) and debris (<200))
ec_microglia <- subset(ec_microglia, subset = nFeature_RNA > 200 & nFeature_RNA < 1500)

# 6. Normalization and Feature Selection
ec_microglia <- NormalizeData(ec_microglia)
ec_microglia <- FindVariableFeatures(ec_microglia, selection.method = "vst", nfeatures = 2000)

# 7. Dimensionality Reduction (PCA & UMAP)
dim(ec_microglia)
ec_microglia <- ScaleData(ec_microglia)
ec_microglia <- RunPCA(ec_microglia, npcs = 30)
ec_microglia <- RunUMAP(ec_microglia, dims = 1:20)

# 8. Visualize global disease state clustering
DimPlot(ec_microglia, reduction = "umap", group.by = "disease", pt.size = 0.5) + 
  ggtitle("UMAP of EC Microglia: Disease State")

# 9. Differential Gene Expression (DGE)
Idents(ec_microglia) <- "disease"
evasion_markers <- FindMarkers(ec_microglia, 
                               ident.1 = "Alzheimer disease", 
                               ident.2 = "normal", 
                               min.pct = 0.25)


#Sort by adjusted p-value for top biological drivers
evasion_markers <- evasion_markers[order(evasion_markers$p_val_adj), ]

# 10. Disease Progression Visualization
head(evasion_markers)
#                  p_val avg_log2FC pct.1 pct.2     p_val_adj
#ENSG00000251562 1.163196e-174   1.096540 0.998 0.994 3.721880e-170 -> MALAT1
#ENSG00000198938 1.460845e-113  -2.553291 0.289 0.654 4.674267e-109 -> MT-CO3
#ENSG00000253535  3.054858e-83   3.404956 0.531 0.062  9.774629e-79 -> ADAM7-AS1
#ENSG00000198804  1.564842e-67  -2.383239 0.217 0.485  5.007026e-63 -> MT-CO1
#ENSG00000198712  5.677552e-59  -2.170634 0.275 0.522  1.816646e-54 -> MT-CO2
#ENSG00000198763  4.451549e-53  -2.257483 0.205 0.439  1.424362e-48 -> MT-ND2

# 11. Set the exact Ensembl IDs for MALAT1, APOE, and TREM2
# ENSG00000251562 = MALAT1
# ENSG00000130203 = APOE
# ENSG00000095970 = TREM2
# Included APOE and TREM2 alongside top hit (MALAT1) to act as positive controls. 
# Because APOE and TREM2 are the most heavily validated genetic risk factors for 
# Alzheimer's Disease-Associated Microglia (DAM), mapping them visually proves 
# that the isolated evasion cluster aligns with established AD pathology.

library(patchwork)

p1 <- FeaturePlot(ec_microglia, features = "ENSG00000251562", pt.size = 0.5, cols = c("lightgrey", "darkred")) + ggtitle("MALAT1")
p2 <- FeaturePlot(ec_microglia, features = "ENSG00000130203", pt.size = 0.5, cols = c("lightgrey", "darkred")) + ggtitle("APOE")
p3 <- FeaturePlot(ec_microglia, features = "ENSG00000095970", pt.size = 0.5, cols = c("lightgrey", "darkred")) + ggtitle("TREM2")

clean_feature_plot <- p1 | p2 | p3
clean_feature_plot

# Chronologically order Braak Stages (0 -> 6)
ec_microglia$BraakStage <- factor(ec_microglia$BraakStage, levels = sort(unique(ec_microglia$BraakStage)))

# Visualize MALAT1 expression across pathology timeline
RidgePlot(ec_microglia, features = "ENSG00000251562", group.by = "BraakStage", fill.by = "ident") +
  ggtitle("MALAT1 Evasion Activation Across Alzheimer's Progression")

# 12. Save the finalized, fully computed Seurat Object
saveRDS(ec_microglia, file = "/Users/tanvigambhir/Documents/neuro-ml/brain-spatial-data/EvasionNet_Final.rds")
