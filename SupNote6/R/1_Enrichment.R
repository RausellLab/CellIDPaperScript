source("R/utilitary_function.R")
library(tidyverse)
library(CellID)
library(foreach)
library(doMC)
library(profmem)
library(ggpubr)
ThemeSupFig15 <- theme(
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    face = "bold",
    size = 8
  ),
  axis.text.y = element_text(face = "bold", size =8),
  aspect.ratio = 1,
  axis.title = element_text(face = "bold", size =8),
  legend.text = element_text(face = "bold"),
  legend.title = element_text(face = "bold")
)
# Load Data ---------------------------------------------------------------
Seurat10X <- readRDS("../SupNote5/data/10X/seurat_10X.rds")
SeuratSS <- readRDS("../SupNote5/data/SmartSeq/seurat_SS.rds")

# Run MCA -----------------------------------------------------------------

Seurat10X <- RunMCA(Seurat10X)
SeuratSS <- RunMCA(SeuratSS)
registerDoMC(cores = 8)

# Prepare Subset ----------------------------------------------------------


set.seed(1)
Seurat10XDS <- lapply(c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000), function(x) subset(Seurat10X, cells = sample(colnames(Seurat10X), x)))
GS.SS       <- GetCellGeneSet(SeuratSS, cells = sample(colnames(SeuratSS), 10000))
GS.SS.DS    <- lapply(c(10, 100, 200, 500, 1000, 2000, 5000, 10000), function(x) sample(GS.SS, x))
keep_celltype <- table(SeuratSS$cell_type1) %>%  sort %>%  tail(10) %>%  names
SeuratSS@active.ident <- factor(SeuratSS$cell_type1)
SeuratSSCommonCell <- subset(SeuratSS, ident = keep_celltype)
SeuratSSDS <- lapply(c(200, 500, 1000, 2000, 5000, 10000, 20000), function(x) subset(SeuratSSCommonCell, cells = sample(colnames(SeuratSSCommonCell), x)))


# Bench --------------------------------------------------------------------
SupFig15Col <- setNames(RColorBrewer::brewer.pal(name = "Paired", n = 12), c("scmap_cluster","scmap_cell","Seurat","MNN", "CellID(G)","CellID(C)","SCN","scPred","CHETAH","CaSTLe","scID","SingleR"))
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

query_var <- read_rds("data/query_var.rds")
QueryVarDF <- data.table::rbindlist(mapply(x = query_var, y = pred_func_name, FUN = function(x,y) {
  A <- as.data.frame(t(sapply(x, function(x) x)))
  A$ncells <- c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000)
  A$methods <- y
  return(A)
}, SIMPLIFY = F)) %>%  mutate(methods = str_remove(methods, "pred_"))

QueryVarDF$methods[QueryVarDF$methods == "CellID_G"] <- "CellID(G)"
QueryVarDF$methods[QueryVarDF$methods == "CellID_C"] <- "CellID(C)"
QueryVarDF$methods <- factor(QueryVarDF$methods, levels = c("CellID(G)","CellID(C)", "scmap_cluster","scmap_cell","Seurat","MNN","SCN","scPred","CHETAH","CaSTLe","scID","SingleR"))

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

ref_var <- read_rds("data/ref_var.rds")
RefVarDF <- data.table::rbindlist(mapply(x = ref_var, y = pred_func_name, FUN = function(x,y) {
  A <- as.data.frame(t(sapply(x, function(x) x)))
  A$ncells <- c(200, 500, 1000, 2000, 5000, 10000, 20000)
  A$methods <- y
  return(A)
  }, SIMPLIFY = F)) %>%  mutate(methods = str_remove(methods, "pred_"))
RefVarDF$methods[RefVarDF$methods == "CellID_G"] <- "CellID(G)"
RefVarDF$methods[RefVarDF$methods == "CellID_C"] <- "CellID(C)"
RefVarDF$methods <- factor(RefVarDF$methods, levels = c("CellID(G)","CellID(C)", "scmap_cluster","scmap_cell","Seurat","MNN","SCN","scPred","CHETAH","CaSTLe","scID","SingleR"))



# Time --------------------------------------------------------------------

EnrichTime <- foreach(x = Seurat10XDS, .combine = rbind) %:% foreach(y = GS.SS.DS, .combine = c) %dopar% system.time(RunCellHGT(x, y, log.trans = F, p.adjust = F))[3]
write_rds(EnrichTime, "data/EnrichTime.rds")
EnrichTime <- readRDS("data/EnrichTime.rds")

colnames(EnrichTime) <- paste0(c(10, 100, 200, 500, 1000, 2000, 5000, 10000), " gene sets")
EnrichTimeDF <- as_tibble(EnrichTime)
EnrichTimeDF$ncells <- c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000)
EnrichTimeDF <- (gather(EnrichTimeDF, "ngenesets", "value", -ncells))
EnrichTimeDF$ngenesets <- str_remove(EnrichTimeDF$ngenesets, " gene sets")
EnrichTimeDF$ngenesets <- factor(EnrichTimeDF$ngenesets, c(10, 100, 200, 500, 1000, 2000, 5000, 10000))
colnames(EnrichTimeDF)[3] <- "Time"


# Memory -----------------------------------------------------------------
EnrichMemory <- foreach(x = Seurat10XDS, .combine = rbind) %:% foreach(y = GS.SS.DS, .combine = c) %do% {
  gc()
  tt <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14])
  mem <- RunCellHGT(x, y, log.trans = F, p.adjust = F)
  sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt
}
write_rds(EnrichMemory, "data/EnrichMemory.rds")
EnrichMemory <- read_rds("data/EnrichMemory.rds")

EnrichMemoryDF <- as_tibble(EnrichMemory)
colnames(EnrichMemoryDF) <- paste0(c(10, 100, 200, 500, 1000, 2000, 5000, 10000), " gene sets")
EnrichMemoryDF$ncells <- c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000)
EnrichMemoryDF <- (gather(EnrichMemoryDF, "ngenesets", "value", -ncells))
EnrichMemoryDF$ngenesets <- str_remove(EnrichMemoryDF$ngenesets, " gene sets")
EnrichMemoryDF$ngenesets <- factor(EnrichMemoryDF$ngenesets, c(10, 100, 200, 500, 1000, 2000, 5000, 10000))
SupFig15Col2 <- setNames(rev(RColorBrewer::brewer.pal(name = "Spectral", n = 8)), levels(EnrichMemoryDF$ngenesets))
colnames(EnrichMemoryDF)[3] <- "Memory"






# Figure ------------------------------------------------------------------


SupFig15A <- ggplot(QueryVarDF, aes(x = ncells, y = time, color = methods)) + geom_line() + geom_point() + scale_color_manual(values = SupFig15Col) + scale_x_continuous(breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000), trans = "log10") + ThemeSupFig15  + ylab("Time(s)") + xlab("Number of cells") + guides(color = guide_legend(ncol = 1))+ yscale("log10")
SupFig15B <- ggplot(QueryVarDF, aes(x = ncells, y = memory, color = methods)) + geom_line() + geom_point() + scale_color_manual(values = SupFig15Col) + scale_x_continuous(breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000), trans = "log10") + ThemeSupFig15  + ylab("Peak Memory (Mb)") + xlab("Number of cells") + guides(color = guide_legend(ncol = 1))+ yscale("log10")
SupFig15C <- ggplot(RefVarDF, aes(x = ncells, y = time, color = methods)) + geom_line() + geom_point() + scale_color_manual(values = SupFig15Col) + scale_x_continuous(breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000), trans = "log10") + ThemeSupFig15 + ylab("Time(s)") + xlab("Number of cells") + guides(color = guide_legend(ncol = 1))+ yscale("log10")
SupFig15D <- ggplot(RefVarDF, aes(x = ncells, y = memory, color = methods)) + geom_line() + geom_point() + scale_color_manual(values = SupFig15Col) + scale_x_continuous(breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000), trans = "log10") + ThemeSupFig15 + ylab("Peak Memory (Mb)") + xlab("Number of cells") + guides(color = guide_legend(ncol = 1))+ yscale("log10")
SupFig15E <- ggplot(EnrichTimeDF, aes(x = ncells, y = Time, color = ngenesets)) + geom_line() + geom_point() + xscale("log10") + yscale("log10") + theme_bw() + scale_x_continuous(breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000), trans = "log10") + ylab("Time (s)") + xlab("Number of cells") + ThemeSupFig15 + guides(color = guide_legend("n gene sets")) + scale_color_manual(values = SupFig15Col2) + theme(legend.margin = margin(r =30, l = 5))
SupFig15F <- ggplot(EnrichMemoryDF, aes(x = ncells, y = Memory, color = ngenesets)) + geom_line() + geom_point() + xscale("log10") + yscale("log10") + theme_bw() + scale_x_continuous(breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000), trans = "log10") + ylab("Peak Memory (Mb)") + xlab("Number of cells") + ThemeSupFig15 + guides(color = guide_legend("n gene sets")) + scale_color_manual(values = SupFig15Col2) + theme(legend.margin = margin(r =30, l = 5))



AD <- ggarrange(SupFig15A,
          SupFig15B, 
          SupFig15C, 
          SupFig15D, common.legend = T, 
          labels = "AUTO", ncol = 2, nrow = 2, legend = "right")
EF <- ggarrange(SupFig15E,
          SupFig15F, common.legend = T,labels = c("E","F"), ncol = 2, nrow = 1, legend = "right")


# Gather ------------------------------------------------------------------
SupFig15 <- ggarrange(AD,EF, ncol = 1, heights = c(2,1)) 
ggsave(filename = "../FinalFigure/SupFig15.pdf", SupFig15, dpi = 600, units = "mm", width = 180, height = 210)
ggsave(filename = "../FinalFigure/SupFig15.png", SupFig15, dpi = 600, units = "mm", width = 180, height = 210)
