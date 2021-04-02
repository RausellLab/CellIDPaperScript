source("R/utilitary_function.R")
seurat_ATAC_sub_SS <- read_rds("data/ATAC/seurat_ATAC_sub_SS.rds")
seurat_ATAC_sub_10X <- read_rds("data/ATAC/seurat_ATAC_sub_10X.rds")
seurat_10X_sub <- read_rds("data/10X/seurat_10X_sub.rds")
seurat_SS_sub <- read_rds("data/SmartSeq/seurat_SS_sub.rds")

# SS Ref ----------------------------------------------------------------------

GSG_SS <- GetGroupGeneSet(seurat_SS_sub, group.by = "cell_type1", n.features = 50)
GSC_SS <- GetCellGeneSet(seurat_SS_sub, n.features = 50)

HGTg <- RunCellHGT(seurat_ATAC_sub_SS, pathways = GSG_SS, log.trans = T) 
pred <- rownames(HGTg)[apply(HGTg,2, which.max)]
pred <- ifelse(apply(HGTg,2, max) >2, pred, "unassigned")
seurat_ATAC_sub_SS$pred_CellID_G <- pred

HGTc <- RunCellHGT(seurat_ATAC_sub_SS, pathways = GSC_SS)
pred <- seurat_SS_sub$cell_type1[rownames(HGTc)[apply(HGTc,2, which.max)]]
pred <- ifelse(apply(HGTc,2, max)>2, pred, "unassigned")
seurat_ATAC_sub_SS$pred_CellID_C <- pred

seurat_ATAC_sub_SS$pred_Seurat        <- pred_Seurat(seurat_ATAC_sub_SS,ref_seurat = seurat_SS_sub)
seurat_ATAC_sub_SS$pred_scmap_cluster <- pred_scmap_cluster(seurat_ATAC_sub_SS,ref_seurat = seurat_SS_sub)
seurat_ATAC_sub_SS$pred_scmap_cell    <- pred_scmap_cell(seurat_ATAC_sub_SS,ref_seurat = seurat_SS_sub)
seurat_ATAC_sub_SS$pred_SCN           <- pred_SCN(x = seurat_ATAC_sub_SS,ref_seurat = seurat_SS_sub)
seurat_ATAC_sub_SS$pred_CHETAH        <- pred_CHETAH(x = seurat_ATAC_sub_SS,ref_seurat = seurat_SS_sub)
seurat_ATAC_sub_SS$pred_CaSTLe        <- pred_CaSTLe(x = seurat_ATAC_sub_SS,ref_seurat = seurat_SS_sub)
seurat_ATAC_sub_SS$pred_MNN           <- pred_MNN(x = seurat_ATAC_sub_SS,ref_seurat = seurat_SS_sub)
seurat_ATAC_sub_SS$pred_scPred        <- pred_scPred(x = seurat_ATAC_sub_SS,ref_seurat = seurat_SS_sub)
seurat_ATAC_sub_SS$pred_SingleR       <- pred_SingleR(x = seurat_ATAC_sub_SS,ref_seurat = seurat_SS_sub)
seurat_ATAC_sub_SS$pred_scID          <- pred_scID(x = seurat_ATAC_sub_SS,ref_seurat = seurat_SS_sub)


# 10X Ref -----------------------------------------------------------------


GSG_10X <- head(GetGroupGeneSet(seurat_10X_sub, group.by = "cell_type1",  dims = 1:50,  n.features = 50), -1)
GSC_10X <- GetCellGeneSet(seurat_10X_sub,  n.features = 50, dims = 1:50)

HGTg <- RunCellHGT(seurat_ATAC_sub_10X, pathways = GSG_10X, log.trans = T) 
pred <- rownames(HGTg)[apply(HGTg,2, which.max)]
pred <- ifelse(apply(HGTg,2, max) >2, pred, "unassigned")
seurat_ATAC_sub_10X$pred_CellID_G <- pred

HGTc <- RunCellHGT(seurat_ATAC_sub_10X, pathways = GSC_10X)
pred <- seurat_10X_sub$cell_type1[rownames(HGTc)[apply(HGTc,2, which.max)]]
pred <- ifelse(apply(HGTc,2, max)>2, pred, "unassigned")
seurat_ATAC_sub_10X$pred_CellID_C <- pred

seurat_ATAC_sub_10X$pred_Seurat        <- pred_Seurat(seurat_ATAC_sub_10X,ref_seurat = seurat_10X_sub)
seurat_ATAC_sub_10X$pred_scmap_cluster <- pred_scmap_cluster(seurat_ATAC_sub_10X,ref_seurat = seurat_10X_sub)
seurat_ATAC_sub_10X$pred_scmap_cell    <- pred_scmap_cell(seurat_ATAC_sub_10X,ref_seurat = seurat_10X_sub)
seurat_ATAC_sub_10X$pred_SCN           <- pred_SCN(x = seurat_ATAC_sub_10X,ref_seurat = seurat_10X_sub)
seurat_ATAC_sub_10X$pred_CHETAH        <- pred_CHETAH(x = seurat_ATAC_sub_10X,ref_seurat = seurat_10X_sub)
seurat_ATAC_sub_10X$pred_CaSTLe        <- pred_CaSTLe(x = seurat_ATAC_sub_10X,ref_seurat = seurat_10X_sub)
seurat_ATAC_sub_10X$pred_MNN           <- pred_MNN(x = seurat_ATAC_sub_10X,ref_seurat = seurat_10X_sub)
seurat_ATAC_sub_10X$pred_scPred        <- pred_scPred(x = seurat_ATAC_sub_10X,ref_seurat = seurat_10X_sub)
seurat_ATAC_sub_10X$pred_SingleR       <- pred_SingleR(x = seurat_ATAC_sub_10X,ref_seurat = seurat_10X_sub)
seurat_ATAC_sub_10X$pred_scID          <- pred_scID(x = seurat_ATAC_sub_10X,ref_seurat = seurat_10X_sub)

write_rds(seurat_ATAC_sub_SS, path = "data/ATAC/seurat_ATAC_sub_SS.rds")
write_rds(seurat_ATAC_sub_10X, path = "data/ATAC/seurat_ATAC_sub_10X.rds")