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
immune_levels <-  c("CD34", "Eryth", "Mk", "B", "CD4 T", "CD8 T", "NK", "CD14 Mono", "CD16 Mono", "DC", "pDC")
immune_cell <- c("HSC","MPP","CMP","GMP", "MEP", "Erythrocytes", "Megakaryocytes", "Platelets", "B-cells", "CD4 T-cells", "CD8 T-cells", "NK cells", "Basophils", "Eosinophils", "Neutrophils", "CD14 Monocytes", "CD16 Monocytes", "Macrophages", "DC", "cDC", "pDC")
cite_seurat2 <- readRDS(file = "data/cite_seq/cite_seurat.rds")
reap_seurat <- readRDS(file = "data/reap_seq/reap_seurat.rds")
cite_seurat$cell_type1 <- factor(cite_seurat@active.ident, immune_levels)
immune_signature <- read_rds("data/immune_signature/immune_signature.rds")
immune_signature <- immune_signature[immune_cell]

getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
palCol <- setNames(getPalette(11), immune_levels)

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

cite_seurat <- RunUMAP(cite_seurat, dims = 1:50)

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


# Hybrid Methods ----------------------------------------------------------

#AUCell ranking hypergeometric test
AUCellID <- function (X, pathways, reduction = "mca", n.features = 200, features = NULL, dims = seq(50), minSize = 10, log.trans = TRUE, p.adjust = TRUE) 
{
    DT <- getRanking(AUCell_buildRankings(X@assays$RNA@counts, nCores = 16))
    message("ranking genes")
    features <- rownames(DT)
    cells <- colnames(DT)
    i <- pbapply::pbapply(DT, 2, order)[seq(n.features), ]
    j <- rep(seq(ncol(DT)), each = n.features)
    TargetMatrix <- sparseMatrix(i, j, x = 1, dims = c(length(features), 
                                                       length(cells)), dimnames = list(features, cells))
    pathways <- lapply(pathways, function(x) x[x %fin% features])
    pathways <- pathways[sapply(pathways, function(x) length(x) >= 
                                    minSize)]
    message("calculating number of success\n")
    PathwayMat <- pbapply::pbsapply(pathways, function(x) which(features %fin% 
                                                                    x), simplify = F)
    PathwayLen <- unlist(lapply(PathwayMat, length))
    j <- rep(seq(length(PathwayMat)), times = PathwayLen)
    PathwayMatrix <- sparseMatrix(unlist(PathwayMat), j, x = 1, 
                                  dims = c(length(features), length(PathwayMat)), dimnames = list(features, 
                                                                                                  names(PathwayMat)))
    q <- as.data.frame((t(TargetMatrix) %*% PathwayMatrix) - 
                           1)
    m <- sapply(pathways, function(x) sum(x %fin% features))
    n <- sapply(m, function(x) length(features) - x)
    k <- n.features
    message("performing hypergeometric test\n")
    A <- pbapply::pbmapply(FUN = function(q, m, n, k) {
        listhyper <- phyper(seq(-1, max(q)), m, n, k, lower.tail = F)[q + 
                                                                          2]
        return(listhyper)
    }, q = q, m = m, n = n, k = k)
    rownames(A) <- rownames(q)
    A <- t(A)
    if (p.adjust) {
        A <- apply(A, 2, function(x) p.adjust(x, "BH"))
    }
    if (log.trans) {
        A <- as.sparse(-log10(A))
    }
    return(A)
}

#Naive logFC ranking & hypergeometric test
NaiveCellID <- function (X, pathways, reduction = "mca", n.features = 200, features = NULL, dims = seq(50), minSize = 10, log.trans = TRUE, p.adjust = TRUE) 
{
    RM <- (rowMeans(X@assays$RNA@counts)+1)
    DT <- apply(X@assays$RNA@counts, 2, function(x) frank(log((RM+1)/(x+1)), ties.method = "random"))
    message("ranking genes")
    features <- rownames(X@assays$RNA@counts)
    cells <- colnames(X@assays$RNA@counts)
    i <- pbapply::pbapply(DT, 2, order)[seq(n.features), ]
    j <- rep(seq(ncol(DT)), each = n.features)
    TargetMatrix <- sparseMatrix(i, j, x = 1, dims = c(length(features), 
                                                       length(cells)), dimnames = list(features, cells))
    pathways <- lapply(pathways, function(x) x[x %fin% features])
    pathways <- pathways[sapply(pathways, function(x) length(x) >= 
                                    minSize)]
    message("calculating number of success\n")
    PathwayMat <- pbapply::pbsapply(pathways, function(x) which(features %fin% 
                                                                    x), simplify = F)
    PathwayLen <- unlist(lapply(PathwayMat, length))
    j <- rep(seq(length(PathwayMat)), times = PathwayLen)
    PathwayMatrix <- sparseMatrix(unlist(PathwayMat), j, x = 1, 
                                  dims = c(length(features), length(PathwayMat)), dimnames = list(features, 
                                                                                                  names(PathwayMat)))
    q <- as.data.frame((t(TargetMatrix) %*% PathwayMatrix) - 
                           1)
    m <- sapply(pathways, function(x) sum(x %fin% features))
    n <- sapply(m, function(x) length(features) - x)
    k <- n.features
    message("performing hypergeometric test\n")
    A <- pbapply::pbmapply(FUN = function(q, m, n, k) {
        listhyper <- phyper(seq(-1, max(q)), m, n, k, lower.tail = F)[q + 
                                                                          2]
        return(listhyper)
    }, q = q, m = m, n = n, k = k)
    rownames(A) <- rownames(q)
    A <- t(A)
    if (p.adjust) {
        A <- apply(A, 2, function(x) p.adjust(x, "BH"))
    }
    if (log.trans) {
        A <- as.sparse(-log10(A))
    }
    return(A)
}


# Cite Seq ----------------------------------------------------------------

AUCellID_Cite <- AUCellID(X = cite_seurat, reduction = "mca",n.features = 200, pathways = immune_signature, log.trans = T, minSize = 5)
AUCellID_Cite_pred <- rownames(AUCellID_Cite)[apply(AUCellID_Cite, 2, which.max)]
cite_seurat@meta.data$pred_CellID_AUC <- ifelse(2 < apply(AUCellID_Cite, 2, max), AUCellID_Cite_pred, "unassigned")

NaiveCite <- NaiveCellID(X = cite_seurat, reduction = "mca", n.features = 200, pathways = immune_signature, log.trans = T, minSize = 5)
AUcite_HGT_naive <- rownames(NaiveCite)[apply(NaiveCite, 2, which.max)]
cite_seurat@meta.data$pred_Naive <-ifelse(2 < apply(NaiveCite, 2, max), AUcite_HGT_naive, "unassigned")


# Reap Seq ----------------------------------------------------------------

AUCellID_Reap <- AUCellID(X = reap_seurat, reduction = "mca", n.features = 200, pathways = immune_signature, log.trans = T, minSize = 5)
AUCellID_Reap_pred <- rownames(AUCellID_Reap)[apply(AUCellID_Reap, 2, which.max)]
reap_seurat@meta.data$pred_CellID_AUC <- ifelse(2 < apply(AUCellID_Reap, 2, max), AUCellID_Reap_pred, "unassigned")

NaiveReap <- NaiveCellID(X = reap_seurat, reduction = "mca",n.features = 200, pathways = immune_signature, log.trans = T, minSize = 5)
AUreap_HGT_naive <- rownames(NaiveReap)[apply(NaiveReap, 2, which.max)]
reap_seurat@meta.data$pred_Naive <- ifelse(2 < apply(NaiveReap, 2, max), AUreap_HGT_naive, "unassigned")


cite_seurat@assays$ImSigAUC <- CreateAssayObject(AUCellID_Cite)
cite_seurat@assays$ImSigNaive <- CreateAssayObject(NaiveCite)


# Calculate Metrics -------------------------------------------------------

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
        predictions = c("pred_CellID", "pred_SCINA", "pred_AUCell", "pred_Naive", "pred_CellID_AUC"),
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

# Overall assessment
ImmuneBenchOverall <- (ImmuneBenchCell  %>%  group_by(methods, data, metrics) %>%  summarise(`standard deviation` = sd(value), value = mean(value)))[,c(1,2,3,5,4)] %>%  mutate(value = round(value, digits = 2), `standard deviation` = round(`standard deviation`, digits = 2))


# Split Table -------------------------------------------------------------

ImmuneBenchCITE <-ImmuneBenchCell %>% filter(data == "CITE") %>%  dplyr::select(data, everything()) %>%  mutate(value = round(value, digits = 2))
ImmuneBenchREAP <-ImmuneBenchCell %>% filter(data == "REAP") %>%  dplyr::select(data, everything()) %>%  mutate(value = round(value, digits = 2))
ImmuneBenchOverallCITE <- ImmuneBenchOverall %>% filter(data == "CITE") %>%  dplyr::select(data, everything())
ImmuneBenchOverallREAP <- ImmuneBenchOverall %>% filter(data == "REAP") %>%  dplyr::select(data, everything())

# write supplementary table 7 ---------------------------------------------

xlsx::write.xlsx(x = as.data.frame(ImmuneBenchCITE), sheetName = "ST7.CellPopulation_CITE", file = "../FinalTable/SupTable7.xlsx", append = F, row.names = F)
xlsx::write.xlsx(x = as.data.frame(ImmuneBenchOverallCITE), sheetName = "ST7.Overall_CITE", file = "../FinalTable/SupTable7.xlsx", append = T, row.names = F)
xlsx::write.xlsx(x = as.data.frame(ImmuneBenchREAP), sheetName = "ST7.CellPopulation_REAP", file = "../FinalTable/SupTable7.xlsx", append = T, row.names = F)
xlsx::write.xlsx(x = as.data.frame(ImmuneBenchOverallREAP), sheetName = "ST7.Overall_REAP", file = "../FinalTable/SupTable7.xlsx", append = T, row.names = F)


# save seurat -------------------------------------------------------------

write_rds(reap_seurat,"data/reap_seq/reap_seurat.rds")
write_rds(cite_seurat, "data/cite_seq/cite_seurat.rds")
