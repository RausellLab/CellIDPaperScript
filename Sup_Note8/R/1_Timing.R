source("R/utilitary_function.R")
library(tidyverse)
library(CellID)
library(foreach)
library(profmem)
library(ggpubr)

# Load Data ---------------------------------------------------------------
Seurat10X <- readRDS("../SupNote5/data/10X/seurat_10X.rds")
SeuratSS <- readRDS("../SupNote5/data/SmartSeq/seurat_SS.rds")

# Run MCA -----------------------------------------------------------------

Seurat10X <- RunMCA(Seurat10X)
SeuratSS <- RunMCA(SeuratSS)

# Prepare Subset ----------------------------------------------------------


set.seed(1)
Seurat10XDS <- lapply(c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000), function(x) subset(Seurat10X, cells = sample(colnames(Seurat10X), x)))
GS.SS       <- GetCellGeneSet(SeuratSS, cells = sample(colnames(SeuratSS), 10000))
GS.SS.DS    <- lapply(c(10, 100, 200, 500, 1000, 2000, 5000, 10000), function(x) sample(GS.SS, x))
keep_celltype <- table(SeuratSS$cell_type1) %>%  sort %>%  tail(10) %>%  names
SeuratSS@active.ident <- factor(SeuratSS$cell_type1)
SeuratSSCommonCell <- subset(SeuratSS, ident = keep_celltype)
SeuratSSDS <- lapply(c(200, 500, 1000, 2000, 5000, 10000, 20000), function(x) subset(SeuratSSCommonCell, cells = sample(colnames(SeuratSSCommonCell), x)))


# -------------------------------------------------------------------------
#Query variable cells number
query_var <- foreach(func = pred_func, y = pred_func_name) %do% {
  print(y)
  pbmcapply::pbmclapply(Seurat10XDS, function(x){ 
    gc()
    tt <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14])
    time <- system.time(func(x, SeuratSSDS[[5]]))
    mem <- sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt
    return(c(time = unname(time[3]), memory = mem))
  }, mc.cores = 8
  )
}
write_rds(query_var, "data/query_var.rds")

QueryVarDF <- data.table::rbindlist(mapply(x = query_var, y = pred_func_name, FUN = function(x,y) {
  A <- as.data.frame(t(sapply(x, function(x) x)))
  A$ncells <- c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000)
  A$methods <- y
  return(A)
}, SIMPLIFY = F)) %>%  mutate(methods = str_remove(methods, "pred_"))

QueryVarDF$methods[QueryVarDF$methods == "CellID_G"] <- "CellID(G)"
QueryVarDF$methods[QueryVarDF$methods == "CellID_C"] <- "CellID(C)"
QueryVarDF$methods <- factor(QueryVarDF$methods, levels = c("CellID(G)","CellID(C)", "scmap_cluster","scmap_cell","Seurat","MNN","SCN","scPred","CHETAH","CaSTLe","scID","SingleR"))
write_rds(QueryVarDF, "data/QueryVarDF.rds")



# -------------------------------------------------------------------------
#Reference variable cells number
ref_var <- foreach(func = pred_func, y = pred_func_name) %do% {
  print(y)
  pbmcapply::pbmclapply(SeuratSSDS, function(x){ 
    gc()
    tt <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14])
    time <- system.time(func(Seurat10XDS[[5]], x))
    mem <- sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt
    return(c(time = unname(time[3]), memory = mem))
  }, mc.cores = 8
  )
}
write_rds(ref_var, "data/ref_var.rds")

RefVarDF <- data.table::rbindlist(mapply(x = ref_var, y = pred_func_name, FUN = function(x,y) {
  A <- as.data.frame(t(sapply(x, function(x) x)))
  A$ncells <- c(200, 500, 1000, 2000, 5000, 10000, 20000)
  A$methods <- y
  return(A)
  }, SIMPLIFY = F)) %>%  mutate(methods = str_remove(methods, "pred_"))

RefVarDF$methods[RefVarDF$methods == "CellID_G"] <- "CellID(G)"
RefVarDF$methods[RefVarDF$methods == "CellID_C"] <- "CellID(C)"
RefVarDF$methods <- factor(RefVarDF$methods, levels = c("CellID(G)","CellID(C)", "scmap_cluster","scmap_cell","Seurat","MNN","SCN","scPred","CHETAH","CaSTLe","scID","SingleR"))
write_rds(RefVarDF, "data/RefVarDF.rds")
