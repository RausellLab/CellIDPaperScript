source("R/utilitary_function.R")
# RNA 10X -----------------------------------------------------------------
list_seurat_10X <- pbmclapply(str_subset(list.files("data/10X/",full.names = T), "droplet"), mc.cores = 16, function(x) {load(x); return(UpdateSeuratObject(tiss))})
seurat_10X <- merge(x = list_seurat_10X[[1]], y = list_seurat_10X[-1])
seurat_10X <- subset(seurat_10X, features = rownames(seurat_10X)[rownames(seurat_10X) %in% MgProteinCodingGenes])
seurat_10X <- subset(seurat_10X, features = rownames(seurat_10X)[nexprs(seurat_10X@assays$RNA@counts, byrow = T)>5])
seurat_10X$cell_type1 <- seurat_10X$cell_ontology_class
seurat_10X$cell_type1[is.na(seurat_10X$cell_type1)] <- "unknown"
write_rds(seurat_10X,"data/10X/seurat_10X.rds")

# RNA SmartSeq ------------------------------------------------------------
list_seurat_SS <- pbmclapply(str_subset(list.files("data/SmartSeq/",full.names = T), "facs"), mc.cores = 16, function(x) {load(x); return(UpdateSeuratObject(tiss))})
seurat_SS <- merge(x = list_seurat_SS[[1]], y = list_seurat_SS[-1])
seurat_SS <- subset(seurat_SS, features = rownames(seurat_SS)[rownames(seurat_SS) %in% MgProteinCodingGenes])
seurat_SS <- subset(seurat_SS, features = rownames(seurat_SS)[nexprs(seurat_SS@assays$RNA@counts, byrow = T)>5])
seurat_SS$cell_type1 <- seurat_SS$cell_ontology_class
seurat_SS$cell_type1[is.na(seurat_SS$cell_type1)] <- "unknown"
write_rds(seurat_SS,"data/SmartSeq/seurat_SS.rds")

# ATAC Seq ----------------------------------------------------------------
ATAC <- read_rds("data/ATAC/activity_scores.quantitative.rds")
ATAC <- ATAC[nexprs(ATAC, byrow = T)>=5,]
rownames(ATAC) <- str_to_title(rownames(ATAC))
ATAC <- ATAC[rownames(ATAC) %in% MgProteinCodingGenes,]
metadata_ATAC <- fread("data/ATAC/cell_metadata.txt")
seurat_ATAC <- CreateSeuratObject(ATAC, meta.data = as.data.frame(metadata_ATAC) %>% set_rownames(metadata_ATAC$cell))
seurat_ATAC$cell_type1 <- seurat_ATAC$cell_label
write_rds(seurat_ATAC,"data/ATAC/seurat_ATAC.rds")

# -------------------------------------------------------------------------
seurat_SS$tissue[seurat_SS$tissue == "Large_Intestine"] <- "LargeIntestine"
seurat_ATAC$tissue[seurat_ATAC$tissue == "BoneMarrow"] <- "Marrow"
seurat_10X$tissue[seurat_10X$tissue == "Heart_and_Aorta"] <- "Heart"
int_ATAC_tissue_SS <- intersect(unique(seurat_SS$tissue),unique(seurat_ATAC$tissue))
int_ATAC_tissue_10X <- intersect(unique(seurat_10X$tissue),unique(seurat_ATAC$tissue))

seurat_10X_sub <- subset(seurat_10X, cells = names(seurat_10X$tissue)[seurat_10X$tissue %in% int_ATAC_tissue_10X])
seurat_10X_sub <- subset(seurat_10X_sub, features =names(which(nexprs(seurat_10X_sub@assays$RNA@counts, byrow = T)>=5)))
seurat_10X_sub <- RunSeuratBasics(seurat_10X_sub)
seurat_10X_sub <- RunMCA(seurat_10X_sub, nmcs = 50)
seurat_10X_sub$cell_type1[is.na(seurat_10X_sub$cell_type1)] <- "unknown"

seurat_ATAC_sub_10X <- subset(seurat_ATAC, cells = names(seurat_ATAC$tissue)[seurat_ATAC$tissue %in% int_ATAC_tissue_10X])
seurat_ATAC_sub_10X <- subset(seurat_ATAC_sub_10X, features =names(which(nexprs(seurat_ATAC_sub2@assays$RNA@counts, byrow = T)>=5)))
seurat_ATAC_sub_10X <- RunSeuratBasics(seurat_ATAC_sub_10X)
seurat_ATAC_sub_10X <- RunMCA(seurat_ATAC_sub_10X, nmcs =50)

seurat_SS_sub <- subset(seurat_SS, cells = names(seurat_SS$tissue)[seurat_SS$tissue %in% int_ATAC_tissue_SS])
seurat_SS_sub <- subset(seurat_SS_sub, features =names(which(nexprs(seurat_SS_sub@assays$RNA@counts, byrow = T)>=5)))
seurat_SS_sub <- RunSeuratBasics(seurat_SS_sub)
seurat_SS_sub <- RunMCA(seurat_SS_sub, nmcs = 50)
seurat_SS_sub$cell_type1[is.na(seurat_SS_sub$cell_type1)] <- "unknown"

seurat_ATAC_sub_SS <- subset(seurat_ATAC, cells = names(seurat_ATAC$tissue)[seurat_ATAC$tissue %in% int_ATAC_tissue_SS])
seurat_ATAC_sub_SS <- subset(seurat_ATAC_sub, features =names(which(nexprs(seurat_ATAC_sub@assays$RNA@counts, byrow = T)>=5)))
seurat_ATAC_sub_SS <- RunSeuratBasics(seurat_ATAC_sub)
seurat_ATAC_sub_SS <- RunMCA(seurat_ATAC_sub, nmcs =50)
seurat_ATAC_sub_SS@reductions$umap@cell.embeddings[,2] <- seurat_ATAC_sub_SS@reductions$umap@cell.embeddings[,2]*-1

write_rds(seurat_ATAC_sub_SS, path = "data/ATAC/seurat_ATAC_sub_SS.rds")
write_rds(seurat_ATAC_sub_10X, path = "data/ATAC/seurat_ATAC_sub_10X.rds")
write_rds(seurat_SS_sub, path = "data/SmartSeq/seurat_SS_sub.rds")
write_rds(seurat_10X_sub, path = "data/10X/seurat_10X_sub.rds")