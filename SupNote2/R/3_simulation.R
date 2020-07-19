library(pheatmap)
library(data.table)
library(foreach)
library(tidyverse)
library(cowplot)
library(gtable)
library(grid)
library(ggpubr)
library(SCINA)
library(AUCell)
library(Seurat)
library(ggforce)
library(CellID)
library(doMC)
registerDoMC(cores = 50)

setDTthreads(1)
cbmc_seurat <- readRDS(file = "data/cite_seq/seurat_cbmc_filtered.rds")
reap_seurat <- readRDS("data/reap_seurat.rds")
median_UMI_rep <- median(reap_seurat$nCount_RNA)
median_UMI_cite <- median(cbmc_seurat$nCount_RNA)
wilcox.test(cbmc_seurat$nCount_RNA,reap_seurat$nCount_RNA)
immune_signature <- read_rds("data/immune_signature/immune_signature.rds")
immune_cell <- c("HSC","MPP","CMP","GMP","B-cells","CD4 T-cells","CD8 T-cells","NK cells","pDC","Megakaryocytes", "Platelets","Erythrocytes","MEP","Basophils","Eosinophils","Neutrophils","CD14 Monocytes","CD16 Monocytes","Macrophages","DC","cDC")
immune_signature <- immune_signature[immune_cell]

# -------------------------------------------------------------------------

#   ____________________________________________________________________________
#   Cite PDC                                                                ####

# Generate Simulation Data ------------------------------------------------
simulation_data_pdc <- foreach(y= 1:100, .final = function(x) setNames(x,paste0("data",1:100))) %dopar% {
  set.seed(y)
  pDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident == "pDC"])
  nonpDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident != "pDC"])
  A <- sample(pDCs,1)
  B <- sample(nonpDCs, 999)
  cells <- c(A,B)
  sub <- SubsetData(cbmc_seurat,cells = cells)
  sub <- RunMCA(sub)
  sub@meta.data$rare <- rownames(sub@meta.data) %in% A
  HGT <- RunCellHGT(sub, pathways = immune_signature, log.trans = T, p.adjust = F)
}

simulation_data_pdc[[1]]["pDC",1000]
TP <- sapply(simulation_data_pdc, function(x) {x["pDC",1000]>2}) %>% sum
FP <- sapply(simulation_data_pdc, function(x){sum(x["pDC",-1000]>2)})
FN <- sapply(simulation_data_pdc, function(x){sum((!x$predictions[x$rare] %in% c("pDC")))})
Precision <- (sum(TP)/(sum(TP)+sum(FP)))
Recall <- sum(TP)/sum(FN+TP)
F1 <- (2 * Precision * Recall)/(Precision + Recall)



simulation_data_cd34 <- foreach(y= 1:100, .final = function(x) setNames(x,paste0("data",1:100))) %dopar% {
  set.seed(y)
  pDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident == "CD34"])
  nonpDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident != "CD34"])
  A <- sample(pDCs,1)
  B <- sample(nonpDCs, 999)
  cells <- c(A,B)
  sub <- SubsetData(cbmc_seurat,cells = cells)
  HGT <- RunCellHGT(sub, pathways = immune_signature, log.trans = T, p.adjust = F, minSize =5)
}

TP <- sapply(simulation_data_cd34, function(x) {any(x[c("HSC", "CMP", "GMP", "MPP", "MEP"),cbmc_seurat@active.ident[colnames(x)] == "CD34"]>2)}) %>% sum
FN <- sapply(simulation_data_cd34, function(x) {(colSums(x[c("HSC", "CMP", "GMP", "MPP", "MEP"),cbmc_seurat@active.ident[colnames(x)] != "CD34"]>2))!=0}) %>% colSums()
FN <- sapply(simulation_data_cd34, function(x){sum((!x$predictions[x$rare] %in% c("pDC")))})
Precision <- (sum(TP)/(sum(TP)+sum(FP)))
Recall <- sum(TP)/sum(FN+TP)
F1 <- (2 * Precision * Recall)/(Precision + Recall)


simulation_data_ery <- foreach(y= 1:100, .final = function(x) setNames(x,paste0("data",1:100))) %dopar% {
  set.seed(y)
  pDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident == "Eryth"])
  nonpDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident != "Eryth"])
  A <- sample(pDCs,1)
  B <- sample(nonpDCs, 999)
  cells <- c(A,B)
  sub <- SubsetData(cbmc_seurat,cells = cells)
  sub <- RunMCA(sub)
  HGT <- (RunCellHGT(sub, pathways = immune_signature, log.trans = T, p.adjust = F))[,cells]
}

sapply(simulation_data_ery, function(x) head(x$pred_AUCell,1)) %>% table


cite_seurat$cell_type1[simulation_data_ery[[1]][,1:2] %>% colnames]

TP <- sapply(simulation_data_ery, function(x) {any(x[c("MEP", "Erythrocytes"),1]>2)})
FP <- sapply(simulation_data_ery, function(x) {colSums(x[c("MEP", "Erythrocytes"),-1]>2)}) %>%  colSums()
FP <- sapply(simulation_data_ery, function(x){sum(x$predictions[!x$rare] %in% c("MEP", "Erythrocytes"))})
FN <- sapply(simulation_data_ery, function(x){sum(!x$predictions[x$rare] %in% c("MEP", "Erythrocytes"))})
Precision <- (sum(TP)/(sum(TP)+sum(FP)))
Recall <- sum(TP)/sum(FN+TP)
F1 <- (2 * Precision * Recall)/(Precision + Recall)



# -------------------------------------------------------------------------

simulation_data_pdc <- foreach(y= 1:100, .final = function(x) setNames(x,paste0("data",1:100))) %dopar% {
  set.seed(y)
  
  
  A <- sample(pDCs,1)
  subset()
  cells <- c(A,B)
  sub <- SubsetData(cbmc_seurat,cells = cells)
  sub <- RunMCA(sub)
  sub@meta.data$rare <- rownames(sub@meta.data) %in% A
  HGT <- RunCellHGT(sub, pathways = immune_signature, log.trans = T, p.adjust = T)
}


pDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident == "pDC"])
nonpDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident != "pDC"])

test <- simulation_data_pdc <- foreach(y= pDCs) %dopar% {
  keep <- c(nonpDCs, y)
  sub <- subset(x = cbmc_seurat, cells = keep)
  sub <- RunMCA(sub)
  sub@meta.data$rare <- rownames(sub@meta.data) %in% y
  HGT <- RunCellHGT(sub, pathways = immune_signature, log.trans = T, p.adjust = T)
}

sum(unlist(lapply(seq(49), function(x)(test[[x]]["pDC",colnames(test[[x]]) %in% pDCs])))>2)

rownames(test[[1]])[apply(test[[1]], 2 , which.max)]
