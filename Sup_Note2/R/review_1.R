Seuratsim <- runPCA(sim)
PCA <- SingleCellExperiment::reducedDim(sim, type = "PCA")
pca_distance <- fields::rdist(PCA)
rownames(pca_distance) <- colnames(sim)
colnames(pca_distance) <- colnames(sim)

cosineDist <- function(x){
    as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}

cosine_distance <- as.matrix(cosineDist(t(logcounts(sim))))
cosine_distance

params <- newSplatParams()
sim <- splatSimulate(params, nGenes = 5000, batchCells= 1000, group.prob = c(0.3, 0.15, 0.25, 0.2, 0.1), method = "groups", dropout.type = "experiment" ,dropout.mid= 3, de.facLoc= 2,seed = 1)
sim <- sim[nexprs(sim, byrow = T)>=5]
sim <- normalize(sim)
sim <- as.Seurat(sim)
sim <- NormalizeData(sim)
sim <- RunMCA(sim)



library(irlba)
library(magrittr)
library(multicon)

test_matrix <- t(sim@assays$RNA@data)
test_matrix_centered <- scale(test_matrix, center = TRUE, scale = FALSE) 
# test_matrix_centered <- scale(test_matrix, center = TRUE, scale = TRUE) 

X <- test_matrix_centered
I <- nrow(X)
J <- ncol(X)
Q <- diag(1/I,I)
Q_sqrt <- diag(sqrt(1/I),I)
Q_sqrt_inv <- diag(1/sqrt(1/I),I)
T <- (Q_sqrt) %*% X
T_SVD <- svd(T, nv = 50, nu = 50) # small matrix
# T = U D_sing V'
D_sing <-diag(T_SVD$d[1:50])
D_sing_inv <-diag(1/T_SVD$d[1:50])
D_eign <-D_sing^2
D_eign_inv <-diag((1/((T_SVD$d^2)[1:50])))
D_eign_inv_sqrt <-diag(sqrt(1/((T_SVD$d^2)[1:50])))



U<-T_SVD$u 
V<-T_SVD$v 



# The row-metric preserving biplot is given by:
Ppal_Row_coordiantes = diag(1/sqrt(1/I),I) %*% U %*% D_sing 
Std_Col_coordiantes = V


# The column-metric preserving biplot (aka Covariation biplot) is given by:
Std_Row_coordiantes = diag(1/sqrt(1/I),I) %*% U
Ppal_Col_coordiantes = V %*% D_sing
#====================




plot(sim@reductions$mca@cell.embeddings[,c(1,3)], col = "red")
points(sim@reductions$mca@feature.loadings[,c(1,3)], col = "black")
plot(U[,c(1,3)], col = "black")
points(V[,c(1,3)], col = "red")
plot(PCA$x[,c(1,3)])
points(PCA$rotation[,c(1,3)], col ="red")



MCA <- ggplot() + geom_point(as.data.frame(sim@reductions$mca@cell.embeddings[,c(1,2)]) %>% set_colnames(c("dim1","dim2")), mapping = aes(x = dim1, y= dim2)) + 
    geom_point(as.data.frame(sim@reductions$mca@feature.loadings[,c(1,2)]) %>% set_colnames(c("dim1","dim2")), mapping = aes(x = dim1, y= dim2), color = "red", shape = 4, size =0.1) + theme_classic() +theme(aspect.ratio = 1) + ggtitle("MCA")


Biplots <- ggplot() + geom_point(as.data.frame(-Std_Row_coordiantes[,c(1,2)]) %>% set_colnames(c("dim1","dim2")), mapping = aes(x = dim1, y= dim2)) + 
    geom_point(as.data.frame(-Ppal_Col_coordiantes[,c(1,2)]), mapping = aes(x = V1, y= V2), color = "red", shape = 4, size =0.1) + theme_classic() +theme(aspect.ratio = 1)  + ggtitle("Biplots") 


library(patchwork)

MCA +Biplots 

MCA_dist <- as.vector(fields::rdist(sim@reductions$mca@cell.embeddings[,1:10], sim@reductions$mca@feature.loadings[,1:10]))
biplots_dist <- as.vector(fields::rdist(Std_Row_coordiantes[,1:10], Ppal_Col_coordiantes[,1:10]))

correlation <- cor.test(MCA_dist,biplots_dist, method = "spearman")




# Heatmap -----------------------------------------------------------------

simulation2 <- function(x){
    params <- newSplatParams()
    sim <- splatSimulate(params, nGenes = 5000, batchCells= 1000, group.prob = c(0.3, 0.15, 0.25, 0.2, 0.1), method = "groups", dropout.type = "experiment" ,dropout.mid= 3, de.facLoc= 2,seed = x)
    sim <- sim[nexprs(sim, byrow = T)>=5]
    sim <- normalize(sim)
    sim <- as.Seurat(sim)
    sim <- NormalizeData(sim)
    sim <- RunMCA(sim, nmcs = 10)
    test_matrix <- t(sim@assays$RNA@data)
    test_matrix_centered <- scale(test_matrix, center = TRUE, scale = TRUE) 
    X <- test_matrix_centered
    I <- nrow(X)
    J <- ncol(X)
    Q <- diag(1/I,I)
    Q_sqrt <- diag(sqrt(1/I),I)
    Q_sqrt_inv <- diag(1/sqrt(1/I),I)
    T <- (Q_sqrt) %*% X
    T_SVD <- svd(T, nv = 10, nu = 10) # small matrix
    D_sing <-diag(T_SVD$d[1:10])
    D_sing_inv <-diag(1/T_SVD$d[1:10])
    D_eign <-D_sing^2
    D_eign_inv <-diag((1/((T_SVD$d^2)[1:10])))
    D_eign_inv_sqrt <-diag(sqrt(1/((T_SVD$d^2)[1:10])))
    U<-T_SVD$u 
    V<-T_SVD$v 
    Std_Row_coordiantes = diag(1/sqrt(1/I),I) %*% U
    Ppal_Col_coordiantes = (V %*% D_sing) %>% set_colnames(paste0("pca_",1:10))
    cor(sim@reductions$mca@feature.loadings, Ppal_Col_coordiantes)
}

gene_cor <- pbmclapply(X = 1:100, FUN = simulation2, mc.cores = 50)
Reduce(cbind, lapply(gene_cor, function(x) diag(abs(x)))) %>% apply(1,sd)
Reduce(cbind,lapply(gene_cor, function(x) diag(abs(x)))) %>% apply(1,median)
Mediangenecor_Matrix <- apply(simplify2array(lapply(gene_cor,abs)), c(1,2), median)
SupFig1Dgene_cor <- ComplexHeatmap::Heatmap(abs(Mediangenecor_Matrix),cluster_rows = F, cluster_columns = F, name = "Median\nSpearman\ncorrelation", 
                                      col = circlize::colorRamp2(c(0, 0.5, 1), c("#4575B4", "white", "#D73027")), 
                                      cell_fun = function(j, i, x, y, width, height, fill) {
                                          grid.text(sprintf("%.2f", abs(Mediangenecor_Matrix)[i, j]), x, y, gp = gpar(fontsize = 11, fontface = "bold"))
                                      },border =T,rect_gp = gpar(col = "white", lwd = 1),column_names_rot = 45,row_names_gp = gpar(fontface = "bold"),column_names_gp = gpar(fontface = "bold"),heatmap_legend_param = list(legend_height = unit(6, "cm"),border = "black",grid_width = unit(1, "cm"),labels_gp = gpar(fontface = "bold"),title_position = "topcenter"))

ggsave(filename = "../FinalFigure/SupFig1Dgene_cor.pdf", plot = ggdraw(grid.grabExpr(draw(SupFig1Dgene_cor))), width = 7 ,height = 6)
ggsave(filename = "../FinalFigure/SupFig1Dgene_cor.png", plot = ggdraw(grid.grabExpr(draw(SupFig1Dgene_cor))), width = 7 ,height = 6, dpi = 600)
