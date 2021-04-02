library(Seurat)
library(CellID)
library(data.table)
library(tidyverse)
library(SingleCellExperiment)

set.seed(1)
# Cite Seq ----------------------------------------------------------------
cbmc.rna <- as.sparse(read.csv(file = "data/cite_seq/cite_seq.csv.gz", sep = ",", header = TRUE, row.names = 1))
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)
cbmc.adt <- as.sparse(read.csv(file = "data/cite_seq/cite_seq_protein.csv.gz", sep = ",", header = TRUE, row.names = 1))
cbmc.adt <- cbmc.adt[setdiff(rownames(x = cbmc.adt), c("CCR5", "CCR7", "CD10")), ]
cbmc <- CreateSeuratObject(counts = cbmc.rna)
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)
ElbowPlot(cbmc, ndims = 50)
cbmc <- FindNeighbors(cbmc, dims = 1:25)
cbmc <- FindClusters(cbmc, resolution = 0.8)
new.cluster.ids <- c("CD4 T", "CD14 Mono", "CD4 T", "NK", "CD14 Mono", "Mouse", "B", 
                     "CD8 T", "CD16 Mono", "T/Mono doublets", "T/Mono doublets", "CD34", "Multiplets", "Mouse", "Eryth", "Mk", 
                     "Mouse", "DC", "pDC")
names(new.cluster.ids) <- levels(cbmc)
cbmc <- RenameIdents(cbmc, new.cluster.ids)
cbmc[["ADT"]] <- CreateAssayObject(counts = cbmc.adt)
cbmc <- NormalizeData(cbmc, assay = "ADT", normalization.method = "CLR")
cbmc <- ScaleData(cbmc, assay = "ADT")
cbmc <- SubsetData(cbmc, ident.remove = c("Mouse", "Multiplets", "T/Mono doublets"))
genes <- rownames(cbmc@assays$RNA@counts)
genes <- genes[genes %in% HgProteinCodingGenes ]
genes <- names(which(scater::nexprs(cbmc@assays$RNA@counts[genes,],byrow = T) >= 5))
seurat_cbmc_filtered <- CreateSeuratObject(counts = cbmc@assays$RNA@counts[genes,])
seurat_cbmc_filtered@active.ident <- cbmc@active.ident
seurat_cbmc_filtered <- NormalizeData(seurat_cbmc_filtered)
seurat_cbmc_filtered <- FindVariableFeatures(seurat_cbmc_filtered)
seurat_cbmc_filtered <- ScaleData(seurat_cbmc_filtered, features = rownames(seurat_cbmc_filtered))
seurat_cbmc_filtered <- RunPCA(seurat_cbmc_filtered, verbose = FALSE)
seurat_cbmc_filtered <- RunUMAP(seurat_cbmc_filtered, dims = 1:50)
seurat_cbmc_filtered <- RunTSNE(seurat_cbmc_filtered, dims = 1:50, num_threads = 16)

seurat_cbmc_filtered <- RunMCA(seurat_cbmc_filtered, nmcs = 50)
write_rds(seurat_cbmc_filtered, "data/cite_seq/cite_seurat.rds")

# Reap Seq ----------------------------------------------------------------



download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100501&format=file",destfile = "data/reap_seq/raw.tar")
untar("data/reap_seq/raw.tar",exdir = "data/reap_seq/")
file1 <- as.matrix(fread("data/reap_seq/GSM2685238_mRNA_2_PBMCs_matrix.txt.gz"), rownames = "V1")
file2 <- as.matrix(fread("data/reap_seq/GSM2685239_mRNA_3_PBMCs_matrix.txt.gz"), rownames = "V1")
raw_mat <- cbind(file1,file2)
raw_mat <- raw_mat[rownames(raw_mat) %in% HgProteinCodingGenes,]  
reap_seurat <- CreateSeuratObject(raw_mat, min.cells = 5)
reap_seurat <- NormalizeData(reap_seurat)
reap_seurat <- FindVariableFeatures(reap_seurat)
reap_seurat <- ScaleData(reap_seurat)
reap_seurat <- RunPCA(reap_seurat)
reap_seurat <- RunUMAP(reap_seurat, dims = 1:50)
reap_seurat <- FindNeighbors(reap_seurat,dims = 1:25)
reap_seurat <- FindClusters(reap_seurat, resolution = 1.2)
new.cluster.ids <- c("CD4 T", "CD14 Mono", "CD14 Mono", "CD8 T", "CD16 Mono", "B", 
                     "NK", "B", "unknown", "DC", "pDC", "Mk", "unknown")
names(new.cluster.ids) <- levels(reap_seurat)
reap_seurat <- RenameIdents(reap_seurat, new.cluster.ids)
file1ADT <- as.matrix(fread("data/reap_seq/GSM2685243_protein_2_PBMCs_matrix.txt.gz"), rownames = "V1")
file2ADT <- as.matrix(fread("data/reap_seq/GSM2685244_protein_3_PBMCs_matrix.txt.gz"), rownames = "V1")
raw_mat_adt <- cbind(file1ADT,file2ADT)
reap_seurat[["ADT"]] <- CreateAssayObject(raw_mat_adt)
reap_seurat <- NormalizeData(reap_seurat, assay = "ADT", normalization.method = "CLR")
reap_seurat <- ScaleData(reap_seurat, assay = "ADT")
reap_seurat <- RunMCA(reap_seurat, nmcs = 50)
reap_seurat$cell_type1 <- reap_seurat@active.ident
reap_seurat$cell_type1 <- factor(reap_seurat$cell_type1, levels = c("Mk", "B", "CD4 T", "CD8 T", "NK",  "CD14 Mono", "CD16 Mono", "DC", "pDC", "unknown"))
write_rds(reap_seurat,"data/reap_seq/reap_seurat.rds")