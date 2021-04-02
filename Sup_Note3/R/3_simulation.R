library(pheatmap)
library(data.table)
library(foreach)
library(tidyverse)
library(magrittr)
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
cbmc_seurat <- readRDS("data/cite_seq/cite_seurat.rds")
reap_seurat <- readRDS("data/reap_seq/reap_seurat.rds")
median_UMI_reap <- median(reap_seurat$nCount_RNA)
median_UMI_cite <- median(cbmc_seurat$nCount_RNA)
wilcox.test(cbmc_seurat$nCount_RNA,reap_seurat$nCount_RNA)
immune_signature <- read_rds("data/immune_signature/immune_signature.rds")
immune_cell <- c("HSC","MPP","CMP","GMP","B-cells","CD4 T-cells","CD8 T-cells","NK cells","pDC","Megakaryocytes", "Platelets","Erythrocytes","MEP","Basophils","Eosinophils","Neutrophils","CD14 Monocytes","CD16 Monocytes","Macrophages","DC","cDC")
immune_signature <- immune_signature[immune_cell]

# -------------------------------------------------------------------------

#   ____________________________________________________________________________
#   Cite PDC                                                                ####

# Generate Simulation Data ------------------------------------------------
simulation_data_pdc <- foreach(y = 1:100, .final = function(x) setNames(x,paste0("data",1:100))) %dopar% {
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


pDCTP <- sapply(simulation_data_pdc, function(x) {x["pDC",1000]>2})
pDCFP <- sapply(simulation_data_pdc, function(x){sum(x["pDC",-1000]>2)})
pDCFN <- 1-pDCTP
pDCPrecision <- mean(na.omit(pDCTP/(pDCTP+pDCFP)))
pDCRecall <- mean(pDCTP/(pDCFN+pDCTP))
pDCF1 <- (2 * pDCPrecision * pDCRecall)/(pDCPrecision + pDCRecall)



simulation_data_cd34 <- foreach(y = 1:100, .final = function(x) setNames(x,paste0("data",1:100))) %dopar% {
  set.seed(y)
  pDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident == "CD34"])
  nonpDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident != "CD34"])
  A <- sample(pDCs,1)
  B <- sample(nonpDCs, 999)
  cells <- c(A,B)
  sub <- SubsetData(cbmc_seurat,cells = cells)
  HGT <- RunCellHGT(sub, pathways = immune_signature[c("HSC", "CMP", "GMP", "MPP", "MEP")], log.trans = T, p.adjust = T, minSize = 5)
}

CD34TP <- sapply(simulation_data_cd34, function(x) {any(x[c("HSC", "CMP", "GMP", "MPP", "MEP"),cbmc_seurat@active.ident[colnames(x)] == "CD34"]>2)}) %>% as.numeric
CD34FP <- sapply(simulation_data_cd34, function(x) {(colSums(x[c("HSC", "CMP", "GMP", "MPP", "MEP"),cbmc_seurat@active.ident[colnames(x)] != "CD34"]>2))!=0}) %>% colSums()
CD34FN <- 1 - CD34TP
CD34Precision <- mean(na.omit(CD34TP/(CD34TP+CD34FP)))
CD34Recall <- mean(CD34TP/(CD34FN+CD34TP))
CD34F1 <- (2 * CD34Precision * CD34Recall)/(CD34Precision + CD34Recall)


simulation_data_ery <- foreach(y = 1:100, .final = function(x) setNames(x,paste0("data",1:100))) %dopar% {
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



EryTP <- sapply(simulation_data_ery, function(x) {any(x[c("Erythrocytes"),1]>2)})
EryFP <- sapply(simulation_data_ery, function(x) {(x[c("Erythrocytes"),-1]>2)}) %>%  colSums()
EryFN <- 1-EryTP
EryPrecision <- mean(na.omit(EryTP/(EryTP+EryFP)))
EryRecall <- mean(EryTP/(EryFN+EryTP))
EryF1 <- (2 * EryPrecision * EryRecall)/(EryPrecision + EryRecall)

SupTable8 <- do.call(rbind,list(c(pDCPrecision,pDCRecall,pDCF1),
                                c(CD34Precision,CD34Recall,CD34F1),
                                c(EryPrecision,EryRecall,EryF1))) %>%  
  set_rownames(c("pDC", "CD34+", "Erythrocytes")) %>%  
  set_colnames(c("Precision", "Recall", "F1")) %>%  
  as.data.frame() %>%  
  rownames_to_column("cell_type") %>%  
  gather("metric", 'values', 2:4) %>%  
  arrange(cell_type) %>% 
  mutate(values = round(values, digits = 3)) %>%  
  tibble()

xlsx::write.xlsx(SupTable8, file = "../../../Definitive/CellIDPaperScript/FinalTable/SupTable8_original.xlsx")



# Competitive -------------------------------------------------------------

predictions_ery <- lapply(simulation_data_ery, function(x) ifelse(apply(x, 2, max) >2, rownames(x)[apply(x, 2, which.max)], "unassigned"))

EryTP <- sapply(predictions_ery, function(x) {x[1]  == "Erythrocytes"})
EryFP <- sapply(predictions_ery, function(x) {x[-1] == "Erythrocytes"}) %>%  colSums()
EryFN <- 1-EryTP
EryPrecision <- mean(na.omit(EryTP/(EryTP+EryFP)))
EryRecall <- mean(EryTP/(EryFN+EryTP))
EryF1 <- (2 * EryPrecision * EryRecall)/(EryPrecision + EryRecall)


simulation_data_cd34 <- foreach(y= 1:100, .final = function(x) setNames(x,paste0("data",1:100))) %dopar% {
  set.seed(y)
  pDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident == "CD34"])
  nonpDCs <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident != "CD34"])
  A <- sample(pDCs,1)
  B <- sample(nonpDCs, 999)
  cells <- c(A,B)
  sub <- SubsetData(cbmc_seurat,cells = cells)
  HGT <- RunCellHGT(sub, pathways = immune_signature, log.trans = T, p.adjust = T, minSize = 5)
}

predictions_CD34 <- lapply(simulation_data_cd34, function(x) ifelse(apply(x, 2, max) >2, rownames(x)[apply(x, 2, which.max)], "unassigned"))
CD34cells <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident %in% "CD34"])
Erycells <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident %in% "Eryth"])
CD34TP <- sapply(predictions_CD34, function(x) {x[names(x) %in% CD34cells] %in% c("HSC", "CMP", "GMP", "MPP", "MEP")}) 
CD34FP <- sapply(predictions_CD34, function(x) {(x[!(names(x) %in% CD34cells)] %in% c("HSC", "CMP", "GMP", "MPP", "MEP")) & (!x[(names(x) %in% Erycells)] %in% c("MEP"))}) %>% colSums()
CD34FN <- 1 - CD34TP
CD34Precision <- mean(na.omit(CD34TP/(CD34TP+CD34FP)))
CD34Recall <- mean(CD34TP/(CD34FN+CD34TP))
CD34F1 <- (2 * CD34Precision * CD34Recall)/(CD34Precision + CD34Recall)

predictions_pdc <- lapply(simulation_data_pdc, function(x) ifelse(apply(x, 2, max) >2, rownames(x)[apply(x, 2, which.max)], "unassigned"))
pDCcells <- names(cbmc_seurat@active.ident[cbmc_seurat@active.ident %in% "pDC"])
pDCTP <- sapply(predictions_pdc, function(x) {x[names(x) %in% pDCcells] %in% c("pDC")})
pDCFP <- sapply(predictions_pdc, function(x){x[!names(x) %in% pDCcells] %in% c("pDC")}) %>%  colSums()
pDCFN <- 1-pDCTP
pDCPrecision <- mean(na.omit(pDCTP/(pDCTP+pDCFP)))
pDCRecall <- mean(pDCTP/(pDCFN+pDCTP))
pDCF1 <- (2 * pDCPrecision * pDCRecall)/(pDCPrecision + pDCRecall)


SupTable8 <- do.call(rbind,list(c(pDCPrecision,pDCRecall,pDCF1),
                                c(CD34Precision,CD34Recall,CD34F1),
                                c(EryPrecision,EryRecall,EryF1))) %>%  
  set_rownames(c("pDC", "CD34+", "Erythrocytes")) %>%  
  set_colnames(c("Precision", "Recall", "F1")) %>%  
  as.data.frame() %>%  
  rownames_to_column("cell_type") %>%  
  gather("metric", 'values', 2:4) %>%  
  arrange(cell_type) %>% mutate(values = round(values, digits = 3)) %>%  tibble()


xlsx::write.xlsx(SupTable8, file = "../../../Definitive/CellIDPaperScript/FinalTable/SupTable8.xlsx")
