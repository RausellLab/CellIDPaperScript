source("R/0_library.R")
# PCA Biplots-----------------------------------------------------------------

simulation2 <- function(x) {
  params <- newSplatParams()
  sim <- splatSimulate(params, nGenes = 5000, batchCells = 1000, group.prob = c(0.3, 0.15, 0.25, 0.2, 0.1), method = "groups", dropout.type = "experiment", dropout.mid = 3, de.facLoc = 2, seed = x)
  sim <- sim[nexprs(sim, byrow = T) >= 5]
  sim <- normalize(sim)
  sim <- as.Seurat(sim)
  sim <- NormalizeData(sim)
  sim <- RunMCA(sim, nmcs = 10)
  test_matrix <- t(sim@assays$RNA@data)
  test_matrix_centered <- scale(test_matrix, center = TRUE, scale = TRUE)
  X <- test_matrix_centered
  I <- nrow(X)
  J <- ncol(X)
  Q <- diag(1 / I, I)
  Q_sqrt <- diag(sqrt(1 / I), I)
  Q_sqrt_inv <- diag(1 / sqrt(1 / I), I)
  T <- (Q_sqrt) %*% X
  T_SVD <- svd(T, nv = 10, nu = 10) # small matrix
  D_sing <- diag(T_SVD$d[1:10])
  D_sing_inv <- diag(1 / T_SVD$d[1:10])
  D_eign <- D_sing^2
  D_eign_inv <- diag((1 / ((T_SVD$d^2)[1:10])))
  D_eign_inv_sqrt <- diag(sqrt(1 / ((T_SVD$d^2)[1:10])))
  U <- T_SVD$u
  V <- T_SVD$v
  Std_Row_coordiantes <- diag(1 / sqrt(1 / I), I) %*% U
  Ppal_Col_coordiantes <- (V %*% D_sing) %>% set_colnames(paste0("pca_", 1:10))
  cor(sim@reductions$mca@feature.loadings, Ppal_Col_coordiantes)
}

gene_cor <- pbmclapply(X = 1:100, FUN = simulation2, mc.cores = 50)
Reduce(cbind, lapply(gene_cor, function(x) diag(abs(x)))) %>% apply(1, sd)
Reduce(cbind, lapply(gene_cor, function(x) diag(abs(x)))) %>% apply(1, median)
Mediangenecor_Matrix <- apply(simplify2array(lapply(gene_cor, abs)), c(1, 2), median)
diag(Mediangenecor_Matrix)
