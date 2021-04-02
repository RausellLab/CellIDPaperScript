source("R/0_library.R")
# Simulate Data Fig1 ------------------------------------------------------
simulation <- function(x) {
  print(x)
  params <- splatter::newSplatParams()
  sim <- splatter::splatSimulate(params, nGenes = 5000, batchCells = 1000, group.prob = c(0.3, 0.15, 0.25, 0.2, 0.1), method = "groups", dropout.type = "experiment", dropout.mid = 3, de.facLoc = 2, seed = x)
  sim <- sim[nexprs(sim, byrow = T) >= 25, ]
  SingleCellExperiment::logcounts(sim) <- log1p(SingleCellExperiment::counts(sim))
  sim <- scater::runPCA(sim, ncomponents = 10, feature_set = rownames(sim), scale_features = T)
  sim <- CellID::RunMCA(sim, nmcs = 10)
  correlation <- cor(reducedDim(sim, type = "PCA"), reducedDim(sim, "MCA"), method = "spearman")
  return(correlation)
}


sim100 <- pbmclapply(X = 1:100, FUN = simulation, mc.cores = 50)


# Heatmap Correlation -----------------------------------------------------

Reduce(cbind, lapply(sim100, function(x) diag(abs(x)))) %>% apply(1, sd)
Reduce(cbind, lapply(sim100, function(x) diag(abs(x)))) %>% apply(1, median)
Mediancor_Matrix <- apply(simplify2array(lapply(sim100, abs)), c(1, 2), median)

# Cells -------------------------------------------------------------------

# Create simulated data
params <- newSplatParams()
sim <- splatSimulate(params,
  nGenes = 5000,
  batchCells = 1000, group.prob = c(0.3, 0.15, 0.25, 0.2, 0.1), method = "groups", dropout.type = "experiment", dropout.mid = 3, de.facLoc = 2, seed = 1
)
logcounts(sim) <- log1p(counts(sim))
sim <- sim[nexprs(sim, byrow = T) >= 5]
print(sum((logcounts(sim)) == 0) / (dim(logcounts(sim))[1] * dim(logcounts(sim))[2]))
sim <- runPCA(sim, ncomponents = 10, subset_row = rownames(sim), scale = T)
sim <- RunMCA(sim, nmcs = 10)
sim_seurat <- as.Seurat(sim)
sim_seurat@reductions$MCA@feature.loadings <- attributes(x = reducedDim(sim, "MCA"))$genesCoordinates
sim_seurat@reductions$MCA@misc$mca.flag <- TRUE
sim_seurat@reductions$PCA@cell.embeddings <- -sim_seurat@reductions$PCA@cell.embeddings

ggPCA <- DimPlot(sim_seurat, reduction = "PCA", pt.size = 0.8, cols = "gray") + ggpubr::theme_classic2() + theme(legend.position = "none", aspect.ratio = 1)
ggMCA <- DimPlot(sim_seurat, reduction = "MCA", pt.size = 0.8, cols = "gray") + ggpubr::theme_classic2() + theme(legend.position = "none", aspect.ratio = 1)
ggMCA_gene_space <- Loadings(sim_seurat, "MCA") %>%
  as.data.frame() %>%
  ggplot(aes(x = MCA_1, y = MCA_2)) +
  geom_point(size = 0.1, shape = 4) +
  ggpubr::theme_classic2() +
  theme(legend.position = "none", aspect.ratio = 1)
SupFig1A <- ggarrange(ggPCA, legend = "none", labels = "A", vjust = 12, font.label = list(size = 18))
SupFig1BC <- ggarrange(plotlist = list(ggMCA, ggMCA_gene_space), ncol = 1, legend = "none", labels = c("B", "C"), font.label = list(size = 18))
SupFig1ABC <- ggarrange(SupFig1A, SupFig1BC, font.label = list(size = 18), legend = "none")
SupFig1Dgb <- ComplexHeatmap::Heatmap(abs(Mediancor_Matrix),
  cluster_rows = F, cluster_columns = F, name = "Median\nSpearman\ncorrelation",
  col = circlize::colorRamp2(c(0, 0.5, 1), c("#4575B4", "white", "#D73027")),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", abs(Mediancor_Matrix)[i, j]), x, y, gp = gpar(fontsize = 11, fontface = "bold"))
  }, border = T, rect_gp = gpar(col = "white", lwd = 1), column_names_rot = 45, row_names_gp = gpar(fontface = "bold"), column_names_gp = gpar(fontface = "bold"), heatmap_legend_param = list(legend_height = unit(6, "cm"), border = "black", grid_width = unit(1, "cm"), labels_gp = gpar(fontface = "bold"), title_position = "topcenter")
)
SupFig1D <- ggarrange(ggdraw(grid.grabExpr(draw(SupFig1Dgb))) + theme(plot.margin = unit(c(1, 0, 0, 1), "cm")), labels = "D", font.label = list(size = 18))
SupFig1 <- ggarrange(SupFig1ABC, SupFig1D)

ggsave(filename = "../FinalFigure/SupFig1.pdf", plot = SupFig1, width = 12.5, height = 6)
ggsave(filename = "../FinalFigure/SupFig1.png", plot = SupFig1, width = 12.5, height = 6, dpi = 600)
