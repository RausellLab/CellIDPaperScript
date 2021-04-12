source("R/utilitary_function.R")
library(CellID)

seurat_10X <- read_rds("data/10X/seurat_10X.rds")
seurat_SS <- read_rds("data/SmartSeq/seurat_SS.rds")

seurat_10X$tissue[seurat_10X$tissue == "Heart_and_Aorta"] <- "Heart"
seurat_SS <- subset(seurat_SS, cells = colnames(seurat_SS)[(seurat_SS$tissue) %in% unique(seurat_10X$tissue)])
seurat_10X <- subset(seurat_10X, cells = colnames(seurat_10X)[unique(seurat_10X$tissue) %in% (seurat_SS$tissue)])

int_cell_type <- intersect(seurat_SS$cell_type1,seurat_10X$cell_type1) 

seurat_10X <- subset(seurat_10X, cells = colnames(seurat_10X)[(seurat_10X$cell_type1) %in% int_cell_type])
seurat_SS <- subset(seurat_SS, cells = colnames(seurat_SS)[(seurat_SS$cell_type1) %in% int_cell_type])

cell_type_sup <- names(table(seurat_SS$cell_type1)[table(seurat_SS$cell_type1) > 200])
seurat_10X <- subset(seurat_10X, cells = colnames(seurat_10X)[(seurat_10X$cell_type1) %in% cell_type_sup])
seurat_SS <- subset(seurat_SS, cells = colnames(seurat_SS)[(seurat_SS$cell_type1) %in% cell_type_sup])
tst <- table(seurat_10X$tissue, seurat_10X$cell_type1)>0
apply(tst, 1, FUN = function(x){colnames(tst)[x]})

seurat_10X <- RunMCA(seurat_10X) 
seurat_SS <- RunMCA(seurat_SS) 


seurat_SS$cell_type1 %>% table

# Geneset SS --------------------------------------------------------------


GS_All_SS   <- GetCellGeneSet(seurat_SS)

Tissue_Split_SS <- SplitObject(seurat_SS, split.by = "tissue")
GS_Tissue_SS <- pbmclapply(Tissue_Split_SS, function(x) GetCellGeneSet(RunMCA(x, nmcs = 50)), mc.cores = 20)
GS_Tissue_SS <- GS_Tissue_SS %>% unname() %>%  unlist(recursive = F)

CT_Split_SS <- SplitObject(seurat_SS, split.by = "cell_type1")
GS_CT_SS <- pbmclapply(CT_Split_SS, function(x) GetCellGeneSet(RunMCA(x, nmcs = 50)), mc.cores = 40)
GS_CT_SS <- GS_CT_SS %>% unname() %>%  unlist(recursive = F)


# Geneset 10 --------------------------------------------------------------

GS_All_10X   <- GetCellGeneSet(seurat_10X)

Tissue_Split_10X <- SplitObject(seurat_10X, split.by = "tissue")
GS_Tissue_10X <- pbmclapply(Tissue_Split_10X, function(x) GetCellGeneSet(RunMCA(x, nmcs = 50)), mc.cores = 20)
GS_Tissue_10X <- GS_Tissue_10X %>% unname() %>%  unlist(recursive = F)

CT_Split_10X <- SplitObject(seurat_10X, split.by = "cell_type1")
GS_CT_10X <- pbmclapply(CT_Split_10X, function(x) GetCellGeneSet(RunMCA(x, nmcs = 10), dims =1:10), mc.cores = 40)
GS_CT_10X <- GS_CT_10X %>% unname() %>%  unlist(recursive = F)



# intersection ------------------------------------------------------------

ncell_type <- apply(table(seurat_10X$tissue,seurat_10X$cell_type1),1, function(x) x !=0) %>% colSums()

median_cell_type <- function(x,y){
  as.data.frame(mapply(x = x[names(y)], y = y, function(x,y){sum(x %in% y)})) %>%  
    rownames_to_column() %>%  
    set_colnames(c("cell", "overlap")) %>%  
    mutate(Tissue = seurat_10X$tissue[cell]) %>% 
    group_by(Tissue) %>%  
    summarise(median_overlap = median(overlap)) %>%  
    mutate(ncell_type = ncell_type[Tissue])
}

median_overall <- function(x,y){median(mapply(x = x[names(y)], y = y, function(x,y){sum(x %in% y)}))}
sd_overall <- function(x,y){sd(mapply(x = x[names(y)], y = y, function(x,y){sum(x %in% y)}))}


Int_GS_All_Tissue_SS <- median_overall(GS_All_SS, GS_Tissue_SS)
Int_GS_All_CT_SS <- median_overall(GS_All_SS, GS_CT_SS)
Int_GS_Tissue_CT_SS <- median_overall(GS_Tissue_SS, GS_CT_SS)

SD_GS_All_Tissue_SS <- round(sd_overall(GS_All_SS, GS_Tissue_SS), digits = 2)
SD_GS_All_CT_SS <- round(sd_overall(GS_All_SS, GS_CT_SS), digits = 2)
SD_GS_Tissue_CT_SS <- round(sd_overall(GS_Tissue_SS, GS_CT_SS), digits = 2)


Int_GS_All_Tissue_10X <- median_overall(GS_All_10X, GS_Tissue_10X)
Int_GS_All_CT_10X <- median_overall(GS_All_10X, GS_CT_10X)
Int_GS_Tissue_CT_10X <- median_overall(GS_Tissue_10X, GS_CT_10X)

SD_GS_All_Tissue_10X <- round(sd_overall(GS_All_10X, GS_Tissue_10X),digits = 2)
SD_GS_All_CT_10X <- round(sd_overall(GS_All_10X, GS_CT_10X),digits = 2)
SD_GS_Tissue_CT_10X <- round(sd_overall(GS_Tissue_10X, GS_CT_10X),digits = 2)



overlap_table_SS <- data.frame(All = c(NA, Int_GS_All_Tissue_SS, Int_GS_All_CT_SS),
                               Tissue = c(Int_GS_All_Tissue_SS, NA, Int_GS_Tissue_CT_SS),
                               Cell_Type = c(Int_GS_All_CT_SS, Int_GS_Tissue_CT_SS, NA), 
                               row.names = c("All","Tissue","Cell_Type")
)

overlap_table_10X <- data.frame(All = c(NA, Int_GS_All_Tissue_10X, Int_GS_All_CT_10X),
                                Tissue = c(Int_GS_All_Tissue_10X, NA, Int_GS_Tissue_CT_10X),
                                Cell_Type = c(Int_GS_All_CT_10X, Int_GS_Tissue_CT_10X, NA), 
                                row.names = c("All","Tissue","Cell_Type")
)



All_HGT <- RunCellHGT(seurat_10X, pathways = GS_All_SS, dims = 1:50, n.features = 200)
All_pred <- as.vector(seurat_SS$cell_type1[rownames(All_HGT)[(apply(All_HGT, 2, which.max))]])
All_prediction <- ifelse((apply(All_HGT, 2, max))>2, All_pred, "unassigned")


Tissue_HGT <- RunCellHGT(seurat_10X, pathways = GS_Tissue_SS, dims = 1:50, n.features = 200)
Tissue_pred <- as.vector(seurat_SS$cell_type1[rownames(Tissue_HGT)[(apply(Tissue_HGT, 2, which.max))]])
Tissue_prediction <- ifelse((apply(Tissue_HGT, 2, max))>2, Tissue_pred, "unassigned")


CT_HGT <- RunCellHGT(seurat_10X, pathways = GS_CT_SS, dims = 1:50, n.features = 200)
CT_pred <- as.vector(seurat_SS$cell_type1[rownames(CT_HGT[,names(reference_label)])[(apply(CT_HGT[,names(reference_label)], 2, which.max))]])
CT_prediction <- ifelse((apply(CT_HGT, 2, max))>2, CT_pred, "unassigned")



reference_label <- factor(seurat_10X$cell_type1)
levels(reference_label) <- c(levels(reference_label), "unassigned")

CT_label <- factor(CT_prediction, levels = levels(reference_label))
Tissue_label <- factor(Tissue_prediction, levels = levels(reference_label))
All_Label <- factor(All_prediction, levels = levels(reference_label))
CT_bench <- caret::confusionMatrix(CT_label,reference = reference_label)
Tissue_bench <- caret::confusionMatrix(Tissue_label,reference = reference_label)
All_bench <- caret::confusionMatrix(All_Label,reference = reference_label)



CT_bench <- CT_bench$byClass[-nrow(CT_bench$byClass), c("Precision", "Recall", "F1")] %>%  round(digits = 3)
Tissue_bench <- Tissue_bench$byClass[-nrow(Tissue_bench$byClass), c("Precision", "Recall", "F1")] %>%  round(digits = 3)
All_bench <- All_bench$byClass[-nrow(All_bench$byClass), c("Precision", "Recall", "F1")] %>%  round(digits = 3)



CT_bench[is.na(CT_bench)] <- 0
Tissue_bench[is.na(Tissue_bench)] <- 0
All_bench[is.na(All_bench)] <- 0

write_rds(CT_bench, "data/CT_bench.rds")
write_rds(Tissue_bench, "data/Tissue_bench.rds")
write_rds(All_bench, "data/All_bench.rds")


CT_bench <- read_rds("data/CT_bench.rds")
Tissue_bench <- read_rds("data/Tissue_bench.rds")
All_bench <- read_rds("data/All_bench.rds")

rownames(CT_bench) <- str_remove(rownames(CT_bench), "Class: ")
rownames(Tissue_bench) <- str_remove(rownames(Tissue_bench), "Class: ")
rownames(All_bench) <- str_remove(rownames(All_bench), "Class: ")

df <- data.frame(
  All = colMeans(All_bench),
  Tissue = colMeans(Tissue_bench),
  CellType = colMeans(CT_bench)) %>%  round(digits = 3)


xlsx::write.xlsx(All_bench, file = "../FinalTable/SupTable12.xlsx", sheetName = "All", append = F)
xlsx::write.xlsx(Tissue_bench, file = "../FinalTable/SupTable12.xlsx", sheetName = "Tissue", append = T)
xlsx::write.xlsx(CT_bench, file = "../FinalTable/SupTable12.xlsx", sheetName = "CellType", append = T)
xlsx::write.xlsx(df, file = "../FinalTable/SupTable12.xlsx", sheetName = "global", append = T)

features <- rownames(seurat_10X)
PathwayMat <- pbapply::pbsapply(GS_Tissue_SS, function(x) which(features %fin%  x), simplify = F)
PathwayLen <- unlist(lapply(PathwayMat, length))
j <- rep(seq(length(PathwayMat)), times = PathwayLen)
PathwayMatrix <- sparseMatrix(unlist(PathwayMat), j, x = 1, 
                              dims = c(length(features), length(PathwayMat)), dimnames = list(features, 
                                                                                              names(PathwayMat)))

TargetMat <- pbapply::pbsapply(GS_Tissue_10X, function(x) which(features %fin%  x), simplify = F)
TargetLen <- unlist(lapply(TargetMat, length))
j <- rep(seq(length(TargetMat)), times = TargetLen)
TargetMatrix <- sparseMatrix(unlist(TargetMat), j, x = 1, 
                             dims = c(length(features), length(TargetMat)), dimnames = list(features, 
                                                                                            names(TargetMat)))

q <- as.data.frame((t(TargetMatrix) %*% PathwayMatrix) - 1)
m <- sapply(GS_Tissue_SS, function(x) sum(x %fin% features))
n <- sapply(m, function(x) length(features) - x)
k <- 200
message("performing hypergeometric test\n")
A <- pbapply::pbmapply(FUN = function(q, m, n, k) {
  listhyper <- phyper(seq(-1, max(q)), m, n, k, lower.tail = F)[q + 
                                                                  2]
  return(listhyper)
}, q = q, m = m, n = n, k = k)
rownames(A) <- rownames(q)
A <- t(A)
A <- apply(A, 2, function(x) p.adjust(x, "BH"))
B <- as.sparse(-log10(A))

pml <- rownames(B)[apply(B, 2, which.max)]
pred <- factor(as.vector(seurat_SS$cell_type1[pml]), levels = levels(reference_label))
prout <- apply(B, 2, which.max)>2
pred[!prout] <- "unassigned"
pred <- factor(pred, levels = levels(reference_label))
gobiden <- confusionMatrix(pred, reference = reference_label)
Tissue_Tissue_df <- gobiden$byClass[-28,c("Precision", "Recall", "F1")] %>%  round(digits=3)
gobiden$byClass[-28,c("Precision", "Recall", "F1")] %>%  colMeans()

features <- rownames(seurat_10X)
PathwayMat <- pbapply::pbsapply(GS_CT_SS, function(x) which(features %fin%  x), simplify = F)
PathwayLen <- unlist(lapply(PathwayMat, length))
j <- rep(seq(length(PathwayMat)), times = PathwayLen)
PathwayMatrix <- sparseMatrix(unlist(PathwayMat), j, x = 1, 
                              dims = c(length(features), length(PathwayMat)), dimnames = list(features, 
                                                                                              names(PathwayMat)))

TargetMat <- pbapply::pbsapply(GS_CT_10X, function(x) which(features %fin%  x), simplify = F)
TargetLen <- unlist(lapply(TargetMat, length))
j <- rep(seq(length(TargetMat)), times = TargetLen)
TargetMatrix <- sparseMatrix(unlist(TargetMat), j, x = 1, 
                             dims = c(length(features), length(TargetMat)), dimnames = list(features, 
                                                                                            names(TargetMat)))

q <- as.data.frame((t(TargetMatrix) %*% PathwayMatrix) - 1)
m <- sapply(GS_CT_SS, function(x) sum(x %fin% features))
n <- sapply(m, function(x) length(features) - x)
k <- 200
message("performing hypergeometric test\n")
A1 <- pbapply::pbmapply(FUN = function(q, m, n, k) {
  listhyper <- phyper(seq(-1, max(q)), m, n, k, lower.tail = F)[q + 
                                                                  2]
  return(listhyper)
}, q = q, m = m, n = n, k = k)
rownames(A1) <- rownames(q)
A1 <- t(A1)
A1 <- apply(A1, 2, function(x) p.adjust(x, "BH"))
C <- as.sparse(-log10(A1))


CT_CT <- rownames(C[,reference_label%>%  names])[apply(C[,reference_label%>%  names], 2, which.max)]
lel <- apply(C[,reference_label%>%  names], 2, max)>2
CT_CT[!lel] <- "unassigned"
table(CT_CT)
pred_CT_CT <- factor(as.vector(seurat_SS$cell_type1[CT_CT]), levels = levels(reference_label))
pred_CT_CT[!lel] <- "unassigned"
table(pred_CT_CT)
gotrump <- confusionMatrix(pred_CT_CT, reference = reference_label)
CT_CT_DF <- gotrump$byClass[-28,c("Precision", "Recall", "F1")] %>%  round(digits=3)
CT_CT_global <- colMeans(CT_CT_DF[,3:5])
Tissue_Tissue_global <- colMeans(Tissue_Tissue_df[,3:5])
rownames(CT_CT_DF) <- str_remove(rownames(CT_CT_DF), "Class: ")
rownames(Tissue_Tissue_df) <- str_remove(rownames(CT_CT_DF), "Class: ")

df <- data.frame(
  All = colMeans(All_bench),
  Tissue = colMeans(Tissue_bench),
  CellType = colMeans(CT_bench)) %>%  round(digits = 3)


xlsx::write.xlsx(All_bench, file = "../FinalTable/SupTable12.xlsx", sheetName = "All", append = F)
xlsx::write.xlsx(Tissue_bench, file = "../FinalTable/SupTable12.xlsx", sheetName = "Tissue", append = T)
xlsx::write.xlsx(CT_bench, file = "../FinalTable/SupTable12.xlsx", sheetName = "CellType", append = T)
xlsx::write.xlsx(df, file = "../FinalTable/SupTable12.xlsx", sheetName = "global", append = T)

tab_tis <- table(seurat_10X$tissue, seurat_10X$cell_type1)
Tissue <- apply(tab_tis,2, function(x){paste0(rownames(tab_tis)[x !=0], collapse = ", ")})

All_bench <- All_bench %>% as.data.frame(row.names = rownames(All_bench))  %>%  rownames_to_column(var = "cell_type") %>%  mutate(Tissue = Tissue)
Tissue_bench <- Tissue_bench %>% as.data.frame(row.names = rownames(Tissue_bench)) %>%  rownames_to_column(var = "cell_type")  %>%  mutate(Tissue = Tissue)
CT_bench <- CT_bench %>% as.data.frame(row.names = rownames(CT_bench)) %>%  rownames_to_column(var = "cell_type")  %>%  mutate(Tissue = Tissue)
CT_CT_DF <- CT_CT_DF %>% as.data.frame(row.names = rownames(CT_CT_DF)) %>%  rownames_to_column(var = "cell_type")  %>%  mutate(Tissue = Tissue)
Tissue_Tissue_df <- Tissue_Tissue_df %>% as.data.frame(row.names = rownames(Tissue_Tissue_df)) %>%  rownames_to_column(var = "cell_type")  %>%  mutate(Tissue = Tissue)


All_bench <- All_bench %>% select(cell_type, Tissue, everything())
Tissue_bench <- Tissue_bench %>% select(cell_type, Tissue, everything())
CT_bench <- CT_bench %>% select(cell_type, Tissue, everything())
CT_CT_DF <- CT_CT_DF %>% select(cell_type, Tissue, everything())
Tissue_Tissue_df <- Tissue_Tissue_df %>% select(cell_type, Tissue, everything())


df <- data.frame(
  All_All = colMeans(All_bench[,3:5]),
  Tissue_All = colMeans(Tissue_bench[,3:5]),
  CellType_All = colMeans(CT_bench[,3:5]),
  Tissue_Tissue = colMeans(Tissue_Tissue_df[,3:5]),
  CellType_CellType = colMeans(CT_CT_DF[,3:5])
) %>%  round(digits = 3)


xlsx::write.xlsx(All_bench, file = "../FinalTable/SupTable12.xlsx", sheetName = "All_All", append = F)
xlsx::write.xlsx(Tissue_bench, file = "../FinalTable/SupTable12.xlsx", sheetName = "Tissue_All", append = T)
xlsx::write.xlsx(CT_bench, file = "../FinalTable/SupTable12.xlsx", sheetName = "CellType_All", append = T)
xlsx::write.xlsx(CT_CT_DF, file = "../FinalTable/SupTable12.xlsx", sheetName = "Tissue_Tissue", append = T)
xlsx::write.xlsx(Tissue_Tissue_df, file = "../FinalTable/SupTable12.xlsx", sheetName = "CellType_CellType", append = T)
xlsx::write.xlsx(df, file = "../FinalTable/SupTable12.xlsx", sheetName = "global", append = T)
xlsx::write.xlsx(overlap_table_10X, file = "../FinalTable/SupTable12.xlsx", sheetName = "overlap_table_10X", append = T)
xlsx::write.xlsx(overlap_table_SS, file = "../FinalTable/SupTable12.xlsx", sheetName = "overlap_table_SS", append = T)
