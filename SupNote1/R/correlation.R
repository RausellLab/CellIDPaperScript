set.seed(1)
library(CellID)
library(Seurat)
library(splatter)
library(scater)
library(tidyverse)
library(foreach)
library(irlba)
library(magrittr)
library(ComplexHeatmap)
library(FastKNN)
library(cowplot)
library(ggpubr)
library(grid)
library(gtable)
library(pbmcapply)
library(doMC)
registerDoMC(cores =50)

# Simulate Data Fig1 ------------------------------------------------------
simulation <- function(x){
    print(x)
    params <- splatter::newSplatParams()
    sim <- splatter::splatSimulate(params, nGenes = 5000, batchCells = 1000, group.prob = c(0.3, 0.15, 0.25, 0.2, 0.1), method = "groups", dropout.type = "experiment", dropout.mid = 3, de.facLoc = 2,seed=x)
    sim <- sim[nexprs(sim,byrow = T)>=5,]
    SingleCellExperiment::logcounts(sim) <- log1p(SingleCellExperiment::counts(sim))
    sim <- scater::runPCA(sim,ncomponents = 10,feature_set = rownames(sim), scale_features = T)
    sim <- CellID::RunMCA(sim,nmcs = 10)
    correlation <- cor(reducedDim(sim,type = "PCA"),reducedDim(sim,"MCA"),method = "pearson")
return(correlation)
}


sim100 <- pbmclapply(X = 1:100, FUN = simulation, mc.cores = 50)


# Heatmap Correlation -----------------------------------------------------

Reduce(cbind, lapply(sim100, function(x) diag(abs(x)))) %>% apply(1,sd)
Reduce(cbind,lapply(sim100, function(x) diag(abs(x)))) %>% apply(1,median)
Mediancor_Matrix <- apply(simplify2array(lapply(sim100,abs)), c(1,2), median)

# Cells -------------------------------------------------------------------
params <- newSplatParams()
sim <- splatSimulate(params, nGenes = 5000, batchCells= 1000, group.prob = c(0.3, 0.15, 0.25, 0.2, 0.1), method = "groups", dropout.type = "experiment" ,dropout.mid= 3, de.facLoc= 2,seed = 1)
logcounts(sim) <- log1p(counts(sim))
sim <- sim[nexprs(sim, byrow = T)>=5]
print(sum((logcounts(sim))==0)/(dim(logcounts(sim))[1]*dim(logcounts(sim))[2]))
sim <- runPCA(sim,ncomponents = 10, subset_row = rownames(sim),scale = T)
sim <- RunMCA(sim,nmcs = 10)
sim_seurat <- as.Seurat(sim)
sim_seurat@reductions$MCA@feature.loadings <- attributes(x = reducedDim(sim,"MCA"))$genesCoordinates
sim_seurat@reductions$MCA@misc$mca.flag <- TRUE
sim_seurat@reductions$PCA@cell.embeddings <- - sim_seurat@reductions$PCA@cell.embeddings

ggPCA <- DimPlot(sim_seurat, reduction = "PCA",  pt.size = 0.8, cols ="gray") + ggpubr::theme_classic2()  + theme(legend.position="none", aspect.ratio = 1)
ggMCA <- DimPlot(sim_seurat, reduction = "MCA",  pt.size = 0.8, cols ="gray") + ggpubr::theme_classic2() +theme(legend.position="none", aspect.ratio = 1) 
ggMCA_gene_space <- Loadings(sim_seurat,"MCA") %>%  as.data.frame() %>% ggplot(aes(x= MCA_1, y=MCA_2)) + geom_point(size=0.1, shape = 4) + ggpubr::theme_classic2() + theme(legend.position = "none", aspect.ratio = 1)
SupFig1A <- ggarrange(ggPCA, legend = "none", labels = "A", vjust = 12, font.label = list(size =18))
SupFig1BC <- ggarrange(plotlist = list(ggMCA, ggMCA_gene_space), ncol = 1,legend = "none", labels = c("B", "C"), font.label = list(size =18))
SupFig1ABC <- ggarrange(SupFig1A,SupFig1BC, font.label = list(size =18),legend = "none")
SupFig1Dgb <- ComplexHeatmap::Heatmap(abs(Mediancor_Matrix),cluster_rows = F, cluster_columns = F, name = "Median\nSpearman\ncorrelation", 
                                col = circlize::colorRamp2(c(0, 0.5, 1), c("#4575B4", "white", "#D73027")), 
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                    grid.text(sprintf("%.2f", abs(Mediancor_Matrix)[i, j]), x, y, gp = gpar(fontsize = 11, fontface = "bold"))
                                },border =T,rect_gp = gpar(col = "white", lwd = 1),column_names_rot = 45,row_names_gp = gpar(fontface = "bold"),column_names_gp = gpar(fontface = "bold"),heatmap_legend_param = list(legend_height = unit(6, "cm"),border = "black",grid_width = unit(1, "cm"),labels_gp = gpar(fontface = "bold"),title_position = "topcenter"))
SupFig1D <- ggarrange(ggdraw(grid.grabExpr(draw(SupFig1Dgb)))+ theme(plot.margin = unit(c(1,0,0,1), "cm")), labels = "D",font.label = list(size =18))
SupFig1 <- ggarrange(SupFig1ABC,SupFig1D)

ggsave(filename = "../FinalFigure/SupFig1.pdf", plot = SupFig1, width = 12.5 ,height = 6)
ggsave(filename = "../FinalFigure/SupFig1.png", plot = SupFig1, width = 12.5 ,height = 6, dpi = 600)

# NEarest Neighbour fold change -----------------------
drop <- function(sim,Dim,k,x){
    print(x)
    nameCell <- colnames(sim)[x]
    k.near <- colnames(sim)[k.nearest.neighbors(x,distance_matrix = distance, k=k)]
    #Plot for one cell
    non.k.near <- colnames((sim))[!colnames((sim)) %in% c(k.near,nameCell)]
    gene_rank <- names(GetCellGeneRanking(sim, dims = 1:Dim)[[x]])
    norm_matrix <- exp(logcounts(sim))
    fold_change_knear <- log(rowMeans(norm_matrix[,k.near])/(rowMeans(norm_matrix[,non.k.near])))
    #Plot for one cell
    bin <- rep(seq(20),each=round(length(gene_rank)/20))
    bin <- bin[(seq(length(rank(fold_change_knear))))]
    bin[is.na(bin)] <- 20
    DF <- tibble(genes =gene_rank,
                 logFC = fold_change_knear[gene_rank],
                 zero = counts(sim)[gene_rank,paste0("Cell",x)]==0,
                 bin = factor(bin)) 
    DF$zero[DF$zero==FALSE] <- "detected gene"
    DF$zero[DF$zero==TRUE] <- "non detected gene"
    return(DF)
}


dummy <- function(sim,x,k){
    nameCell <- colnames(sim)[x]
    norm <- exp(logcounts(sim))
    logfold_cell <- log(sort(norm[,x]/rowMeans(norm[,-x]), decreasing = T))
    rank <- names(logfold_cell)
    distance <- fields::rdist(t(logcounts(sim)),compact=F)
    colnames(distance) <- colnames(sim)
    rownames(distance) <- colnames(sim)
    k.near <- colnames(sim)[k.nearest.neighbors(x,distance_matrix = distance, k=k)]
    non.k.near <- colnames((sim))[!colnames((sim)) %in% c(k.near,nameCell)]
    bin <- rep(seq(20),each=round(length(rank)/20))
    bin <- bin[(seq(length(rank)))]
    bin[is.na(bin)] <- 20
    logFC <- log(rowMeans(norm[,k.near] + 1)/rowMeans(norm[,non.k.near] + 1))[rank]
    return(data.frame(rank,logFC,bin, zero =ifelse(counts(sim)[rank,x] == 0, yes = "non detected gene", "detected gene")))
}


aucell <- function(sim,Dim,k,x){
    print(x)
    nameCell <- colnames(sim)[x]
    colnames(distance) <- colnames(sim)
    rownames(distance) <- colnames(sim)
    k.near <- colnames(sim)[k.nearest.neighbors(x,distance_matrix = distance, k=k)]
    #Plot for one cell
    non.k.near <- colnames((sim))[!colnames((sim)) %in% c(k.near,nameCell)]
    gene_rank <- names(sort(AUCell::getRanking(AUranking)[,x]))
    norm_matrix <- exp(logcounts(sim))
    fold_change_knear <- log(rowMeans(norm_matrix[,k.near])/(rowMeans(norm_matrix[,non.k.near])))
    #Plot for one cell
    bin <- rep(seq(20),each=round(length(gene_rank)/20))
    bin <- bin[(seq(length(rank(fold_change_knear))))]
    bin[is.na(bin)] <- 20
    DF <- tibble(genes =gene_rank,
                 logFC = fold_change_knear[gene_rank],
                 zero = counts(sim)[gene_rank,paste0("Cell",x)]==0,
                 bin = factor(bin)) 
    DF$zero[DF$zero==FALSE] <- "detected gene"
    DF$zero[DF$zero==TRUE] <- "non detected gene"
    return(DF)
}

params <- splatter::newSplatParams()
sim <- splatSimulate(params, nGenes = 5000, batchCells= 1000, group.prob = c(0.3, 0.15, 0.25, 0.2, 0.1), method = "groups", dropout.type="experiment", dropout.mid= 5, de.facLoc = 2, seed=1)
sim <- splatSimulate(params, nGenes = 5000, batchCells= 1000, group.prob = c(0.3, 0.15, 0.25, 0.2, 0.1), method = "groups", dropout.type = "experiment",  dropout.mid = 3, de.facLoc = 2, seed=1)
logcounts(sim) <- log1p(counts(sim))
sim <- as.Seurat(sim)
sim <- NormalizeData(sim)
sim <- as.SingleCellExperiment(sim)
sim <- sim[nexprs(sim,byrow = T)>= 5,]
print(sum((counts(sim))==0)/(dim(counts(sim))[1]*dim(counts(sim))[2]))
sim <- RunMCA(sim, nmcs = 50,log.transformed = F)


distance <- fields::rdist(t(logcounts(sim)),compact=F)
AUranking <- AUCell::AUCell_buildRankings(exprMat = logcounts(sim))


DF <- drop(sim = sim, Dim = 10, 50, x = 1)
DUM <- dummy(sim,x = 1,50)
AUCELL <- aucell(sim = sim, Dim = 10, 50, x = 1)
cellID_one_cell <- ggboxplot(DF, x = "bin", y = "logFC", color = "black",fill ="zero",outlier.size = 0, add = "jitter",add.params = list(size=0.1)) + geom_hline(yintercept = 0) + labs(fill="") +theme(legend.key.size =  unit(0.5, "in"))
dummy_one_cell <- ggboxplot(DUM, x = "bin", y = "logFC", color = "black",fill ="zero",outlier.size = 0, add = "jitter",add.params = list(size=0.1)) + geom_hline(yintercept = 0) + labs(fill="") +theme(legend.key.size =  unit(0.5, "in"))
aucell_one_cell <- ggboxplot(AUCELL, x = "bin", y = "logFC", color = "black",fill ="zero",outlier.size = 0, add = "jitter",add.params = list(size=0.1)) + geom_hline(yintercept = 0) + labs(fill="") +theme(legend.key.size =  unit(0.5, "in"))


cor_cellid <- foreach(x=1:1000,.combine = rbind) %dopar% {
    A <- drop(sim = sim, Dim = 10, k = 50, x = x)
    A$rank <- seq(1:nrow(A))
    corel <- cor.test(-A$rank, A$logFC, method = "spearman",exact = T)
    c(pval=corel$p.value, rho = corel$estimate)
}

cor_aucell <- foreach(x=1:1000,.combine = rbind) %dopar% {
    A <- aucell(sim = sim, Dim = 10, k = 50, x = x)
    A$rank <- seq(1:nrow(A))
    corel <- cor.test(-A$rank, A$logFC, method = "spearman",exact = T)
    c(pval=corel$p.value, rho = corel$estimate)
}

cor_dummy <- foreach(x=1:1000,.combine = rbind) %dopar% {
    A <- dummy(sim = sim, k = 50, x = x)
    A$rank <- seq(1:nrow(A))
    corel <- cor.test(-A$rank, A$logFC, method = "spearman",exact = T)
    c(pval=corel$p.value, rho = corel$estimate)
}

median_cor_cellid <- apply(cor_cellid, 2, median)
median_cor_aucell <- apply(cor_aucell, 2, median)
median_cor_dummy <- apply(cor_dummy, 2, median)
corr_df <- do.call(rbind,list(median_cor_cellid,
                   median_cor_aucell,
                   median_cor_dummy))
colnames(corr_df) <- c("pval", "rho")
corr_df <- data.frame(corr_df)
corr_df$method <- c("cellID", "AUCell", "Naive")
corr_df <- corr_df[,c(3,1,2)]
write_csv((corr_df), path = "figure/Fig2Correlation.csv")

# AUcell Ranking ----------------------------------------------------------
median_DF_cellID <- foreach(x=1:1000,.combine = rbind) %dopar% {drop(sim, Dim = 50, k = 50,x = x) %>% mutate(data=x)}
median_DF_aucell <- foreach(x=1:1000,.combine = rbind) %dopar% {aucell(sim, Dim = 10, k = 50,x = x) %>% mutate(data=x)}
median_DF_dummy  <- foreach(x=1:1000,.combine = rbind) %dopar% {dummy(sim, k = 50,x = x) %>% mutate(data=x)}

# Other -------------------------------------------------------------------

median_drop <- function(x){DF1 <- x %>% group_by(bin,data) %>% summarise(median_logFC=median(logFC))
DF1 <- DF1 %>% dplyr::mutate(zero="none") %>% select(bin,data,zero,median_logFC) %>% ungroup
DF2 <- x %>% group_by(bin,data,zero) %>% summarise(median_logFC=median(logFC)) %>% ungroup
DFall <- rbind(DF1,DF2)
return(DFall)
}

DFcellID <- median_drop(median_DF_cellID)
DFcellID$zero[DFcellID$zero == "zero"] <- "non detected gene"
DFcellID$zero[DFcellID$zero == "non zero"] <- "detected gene"
gg_median_cellid <- ggboxplot(dplyr::filter(DFcellID,zero!="none"), x = "bin", y = "median_logFC", color = "black",fill ="zero",outlier.size = 0) + geom_hline(yintercept = 0)+ coord_cartesian(ylim = c(-4,2))+ labs(fill="") +theme(legend.key.size =  unit(0.5, "in"))

DFaucell <- median_drop(median_DF_aucell)
DFaucell$zero[DFaucell$zero == "zero"] <- "non detected gene"
DFaucell$zero[DFaucell$zero == "non zero"] <- "detected gene"
gg_median_aucell <- ggboxplot(dplyr::filter(DFaucell,zero!="none"), x = "bin", y = "median_logFC", color = "black",fill ="zero",outlier.size = 0) + geom_hline(yintercept = 0)+ coord_cartesian(ylim = c(-4,2))+ labs(fill="") +theme(legend.key.size =  unit(0.5, "in"))

DFdummy <- median_drop(median_DF_dummy)
DFaucell$zero[DFaucell$zero == "zero"] <- "non detected gene"
DFaucell$zero[DFaucell$zero == "non zero"] <- "detected gene"
gg_median_dummy <- ggboxplot(dplyr::filter(DFdummy,zero!="none"), x = "bin", y = "median_logFC", color = "black",fill ="zero",outlier.size = 0) + geom_hline(yintercept = 0)+ coord_cartesian(ylim = c(-4,2))+ labs(fill="") +theme(legend.key.size =  unit(0.5, "in"))

supFig2 <- ggarrange(plotlist = list(cellID_one_cell,gg_median_cellid,dummy_one_cell, gg_median_dummy, aucell_one_cell, gg_median_aucell), ncol = 2, nrow = 3, common.legend = T, labels = c("A","B","C","D","E","F"))
ggsave(plot = supFig2,filename = "../FinalFigure/supFig2.pdf", width = 12, height = 10)
ggsave(plot = supFig2,filename = "../FinalFigure/supFig2.png", width = 12, height = 10, dpi = 600)
