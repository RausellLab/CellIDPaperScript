library(CellID)
library(SCINA)
library(AUCell)
library(foreach)
library(tidyverse)
library(data.table)
library(gtable)
library(grid)
library(ggpubr)
library(ggforce)
library(patchwork)
library(RColorBrewer)

#  Load Data --------------------------------------------------------------
immune_levels <-  c("CD34", "Eryth", "Mk", "B", "CD4 T", "CD8 T", "NK",  "CD14 Mono", "CD16 Mono", "DC", "pDC")
cite_seurat <- readRDS(file = "data/cite_seq/cite_seurat.rds")
cite_seurat$cell_type1 <- cite_seurat@active.ident
cite_seurat$cell_type1 <- factor(cite_seurat$cell_type1, immune_levels)
reap_seurat <- readRDS(file = "data/reap_seurat.rds")
immune_signature <- read_rds("data/immune_signature/immune_signature.rds")
immune_cell <- c("HSC","MPP","CMP","GMP", "MEP", "Erythrocytes", "Megakaryocytes","Platelets", "B-cells","CD4 T-cells","CD8 T-cells","NK cells",  "Basophils","Eosinophils","Neutrophils","CD14 Monocytes","CD16 Monocytes","Macrophages","DC","cDC","pDC")
immune_signature <- immune_signature[immune_cell]

getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
palCol <- getPalette(11)
palCol <- setNames(palCol, immune_levels)

DimPlot(cite_seurat, reduction = "mca")

#   ____________________________________________________________________________
#   CiteSeq                                                                 ####

# HyperGeometric Test ------------------------------------------------------
{
    set.seed(1)
    cite_HGT <-
        RunCellHGT(
            X = cite_seurat,
            reduction = "mca",
            pathways = immune_signature,
            log.trans = T, minSize = 5
        )
    cite_seurat@assays[["ImSig"]] <- CreateAssayObject(data = cite_HGT)
    rows <- sort(factor(rownames(cite_HGT), levels = immune_cell))
    pred <- rownames(cite_HGT)[apply(cite_HGT, 2, which.max)]
    cite_seurat@meta.data$pred_CellID <-
        ifelse(2 < apply(cite_HGT, 2, max), pred, "unassigned")
    cite_HGT <- cite_HGT[immune_cell,names(sort(cite_seurat@active.ident))]
}

# Aucell ------------------------------------------------------------------

{
    set.seed(1)
    library(AUCell)
    AURanking <- AUCell_buildRankings(cite_seurat@assays$RNA@counts, nCores = 16)
    cells_AUC <- AUCell_calcAUC(immune_signature, AURanking, nCores =  16)
    pred_aucell <- rownames(cells_AUC)[apply(AUCell::getAUC(cells_AUC), 2,which.max)]
    pred_aucell <- ifelse(apply(AUCell::getAUC(cells_AUC), 2,max) >= 0.1, pred_aucell, "unassigned")
    cite_seurat@meta.data$pred_AUCell <- pred_aucell
}
# SCINA -------------------------------------------------------------------

{
    set.seed(1)
    results_scina <- SCINA(cite_seurat@assays$RNA@data, immune_signature, allow_unknown = T, max_iter = 200, convergence_n = 20, convergence_rate = 0.999)
    cite_seurat$pred_SCINA <- results_scina$cell_labels
    cite_seurat$pred_SCINA[cite_seurat$pred_SCINA == "unknown"] <- "unassigned"
}

write_rds(cite_seurat, path = "data/cite_seq/cite_seurat.rds")

cite_seurat <- RunTSNE(cite_seurat, dims = 1:50, num_threads = 16, perplexity = 20)
cite_seurat <- RunUMAP(cite_seurat, dims = 1:50)
DimPlot(cite_seurat, reduction = "umap")
FeaturePlot(cite_seurat, "MEP", reduction = "tsne", min.cutoff = 2)
#   ____________________________________________________________________________
#   Reap Seq####
# HyperGeometric Test ------------------------------------------------------
{
    set.seed(1)
    reap_HGT <-
        RunCellHGT(
            X = reap_seurat,
            reduction = "mca",
            pathways = immune_signature,
            dims = 1:50,
            n.features = 200,
            minSize = 5
        )
    reap_seurat@assays[["ImSig"]] <- CreateAssayObject(data = reap_HGT)
    rows <- sort(factor(rownames(reap_HGT), levels = immune_cell))
    pred <- rownames(reap_HGT)[apply(reap_HGT, 2, which.max)]
    reap_seurat@meta.data$pred_CellID <-
        ifelse(2 < apply(reap_HGT, 2, max), pred, "unassigned")
}

# Aucell ------------------------------------------------------------------

{
    set.seed(1)
    library(AUCell)
    AURanking <- AUCell_buildRankings(reap_seurat@assays$RNA@counts, nCores = 16)
    cells_AUC <- AUCell_calcAUC(immune_signature, AURanking, nCores =  16)
    pred_aucell <- rownames(cells_AUC)[apply(AUCell::getAUC(cells_AUC), 2,which.max)]
    pred_aucell <- ifelse(apply(AUCell::getAUC(cells_AUC), 2,max) >= 0.1, pred_aucell, "unassigned")
    reap_seurat@meta.data$pred_AUCell <- pred_aucell
}
# SCINA -------------------------------------------------------------------

{
    set.seed(1)
    results_scina <- SCINA(reap_seurat@assays$RNA@data, immune_signature, allow_unknown = T,max_iter = 200, convergence_n = 20,convergence_rate = 0.999)
    reap_seurat$pred_SCINA <- results_scina$cell_labels
    reap_seurat$pred_SCINA[reap_seurat$pred_SCINA == "unknown"] <- "unassigned"
}

# -------------------------------------------------------------------------
# Mapping
# -------------------------------------------------------------------------


map_cite <- list(c("B-cells"), c("CD14 Monocytes"), c("CD16 Monocytes"), c("HSC","MPP","GMP","MEP","CMP"),c("CD4 T-cells"), c("CD8 T-cells"),  c("cDC","DC"), c("Erythrocytes","MEP"),c("Megakaryocytes", "Platelets"),c("NK cells"), c("pDC")) %>% set_names(unique(sort(as.vector(cite_seurat@active.ident))))
map_reap <- list(c("B-cells"), c("CD14 Monocytes"), c("CD16 Monocytes"), c("CD4 T-cells"), c("CD8 T-cells"),  c("cDC","DC"), c("Megakaryocytes", "Platelets"),c("NK cells"), c("pDC")) %>% set_names(head(unique(sort(as.vector(reap_seurat@active.ident))),-1))

ImmuneBenchCell <-
    foreach(
        Seurat = list(cite_seurat, subset(
            reap_seurat, idents = "unknown", invert = T
        )),
        mapping = list(map_cite, map_reap),
        data = c("CITE", "REAP"),
        .combine = rbind
    ) %:% foreach(
        predictions = c("pred_CellID", "pred_SCINA", "pred_AUCell"),
        .combine = rbind
    ) %do% {
        TruthDF <-
            tibble(
                Author = as.vector(Seurat$cell_type1),
                Predictions = as.vector(unlist(Seurat[[predictions]])),
                map_label = as.vector(sapply(mapping[as.vector(Seurat$cell_type1)], function(x)
                    paste0(x, collapse = ", ")))
            )
        TruthDF <- TruthDF %>% filter(map_label != "removed")
        TruthDF <-
            TruthDF %>% mutate(positive = mapply(
                x = .$map_label,
                y = .$Predictions,
                FUN = function(x, y) {
                    x %like% y
                }
            ))
        TruthDF <-
            TruthDF %>% mutate(final_map = ifelse(positive, Author, Predictions))
        PredDF <- sapply(unique(TruthDF$Author), function(x, truth) {
            TruthDFpos <- TruthDF %>% filter(Author == x)
            TruthDFneg <- TruthDF %>% filter(Author != x)
            TP <- sum(TruthDFpos$Author == TruthDFpos$final_map)
            FN <- sum(TruthDFpos$Author != TruthDFpos$final_map)
            FP <- sum(TruthDFneg$final_map %in% mapping[[x]])
            Recall <- TP / (TP + FN)
            Precision <- TP / (FP + TP)
            F1 <- 2 * (Recall * Precision) / (Recall + Precision)
            return(c(
                Precision = Precision,
                Recall = Recall,
                F1 = F1
            ))
        },
        truth = TruthDF) %>% as.data.frame() %>%  rownames_to_column(var = "metrics") %>%  gather("cell_type", "value",-1) %>%  mutate(methods = predictions) %>%  mutate(data = data) %>%  mutate(value = ifelse(is.na(value), 0, value))
    } %>% dplyr::select(data, methods, cell_type, everything()) %>%  mutate(methods = str_remove(methods, "pred_")) %>% mutate(methods = str_replace(methods, "CellID_G", "CellID(G)")) %>% mutate(methods = str_replace(methods, "CellID_C", "CellID(C)")) %>%  mutate(cell_type = factor(cell_type, immune_levels))


ImmuneBenchOverall <- ImmuneBenchCell  %>%  group_by(methods, data, metrics) %>%  summarise(macro = mean(value), sd = sd(value))

ImmuneBenchCITE <-ImmuneBenchCell %>% filter(data == "CITE") %>%  dplyr::select(data, everything())
ImmuneBenchREAP <-ImmuneBenchCell %>% filter(data == "REAP") %>%  dplyr::select(data, everything())
ImmuneBenchOverallCITE <- ImmuneBenchOverall %>% filter(data == "CITE") %>%  dplyr::select(data, everything())
ImmuneBenchOverallREAP <- ImmuneBenchOverall %>% filter(data == "REAP") %>%  dplyr::select(data, everything())

xlsx::write.xlsx(as.data.frame(ImmuneBenchCITE), sheetName = "Cell_Population_CITE", "../FinalTable/SupTable3.xlsx", append = F)
xlsx::write.xlsx(as.data.frame(ImmuneBenchOverallCITE), sheetName = "Overall_CITE", "../FinalTable/SupTable3.xlsx", append = T)
xlsx::write.xlsx(as.data.frame(ImmuneBenchREAP), sheetName = "Cell_Population_REAP", "../FinalTable/SupTable3.xlsx", append = T)
xlsx::write.xlsx(as.data.frame(ImmuneBenchOverallREAP), sheetName = "Overall_REAP", "../FinalTable/SupTable3.xlsx", append = T)