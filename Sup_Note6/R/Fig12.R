#   ____________________________________________________________________________
#   DimPlot 10X                                                             ####
formating <- function(x) x +scale_colour_manual(values = SupFig12Col, drop =F) + theme(plot.margin = unit(c(0.1,0,0,0.1), "lines"))+
  scale_y_continuous(breaks = seq(from = -40 , to = 40, by =20)) +
  scale_x_continuous(breaks = seq(from = -40 , to = 40, by =20)) +
  coord_fixed(ylim=c(-50, 50),xlim=c(-50, 50)) + guides(col = guide_legend(ncol = 1, override.aes = list(size = 5)))

cell_type_plot_10X <- c(head(rowterms,-8), "Enterocytes")
cells_to_plot <- colnames(seurat_ATAC_sub_10X)[seurat_ATAC_sub_10X$cell_type1 %in% cell_type_plot_10X]
SupFig12Col <- colorRampPalette(colors = RColorBrewer::brewer.pal(8,'Set1'))(length(cell_type_plot_10X))
SupFig12Col <- c(setNames(SupFig12Col,cell_type_plot_10X))
SupFig12Col <- c(SupFig12Col, other = "grey30", unassigned = "lightgrey")
Dim10X <- DimPlot(seurat_ATAC_sub_10X, group.by = "cell_type1",cells = cells_to_plot, reduction = "tsne")
Dim10X$data$cell_type1 <- factor(Dim10X$data$cell_type1, levels = names(SupFig12Col))
Dim10X <- Dim10X %>%  formating

map10XATAClabel <- map10X[seurat_ATAC_sub_10X$cell_type1]

ok <- mapply(x = map10XATAClabel, y = seurat_ATAC_sub_10X$pred_CellID_G, z = seurat_ATAC_sub_10X$cell_type1, function(x,y,z){ifelse(y%in%x, z , last(names(map10X)[sapply(map10X, function(x,y){any(x %in% y)},y = y,simplify = T)]))})
ok[seurat_ATAC_sub_10X$pred_CellID_G == "unassigned"] <- "unassigned" 
ok <- ok[seurat_ATAC_sub_10X$cell_type1 %in% cell_type_plot_10X]
Dim10Xbis <- DimPlot(seurat_ATAC_sub_10X, group.by = "cell_type1", cells = cells_to_plot, reduction = "tsne")
Dim10Xbis$data$cell_type1 <- factor(ok, levels = c(cell_type_plot_10X,"unassigned","other"))
Dim10Xbis$data$cell_type1[is.na(Dim10Xbis$data$cell_type1)] <- "other"
Dim10Xbis <- Dim10Xbis %>% formating

ok <- mapply(x = map10XATAClabel, y = seurat_ATAC_sub_10X$pred_CellID_C, z = seurat_ATAC_sub_10X$cell_type1, function(x,y,z){ifelse(y%in%x, z , last(names(map10X)[sapply(map10X, function(x,y){any(x %in% y)},y = y,simplify = T)]))})
ok[seurat_ATAC_sub_10X$pred_CellID_C == "unassigned"] <- "unassigned" 
ok <- ok[seurat_ATAC_sub_10X$cell_type1 %in% cell_type_plot_10X]
Dim10Xbis2 <- DimPlot(seurat_ATAC_sub_10X, group.by = "cell_type1", cells = cells_to_plot, reduction = "tsne")
Dim10Xbis2$data$cell_type1 <- factor(ok, levels = c(cell_type_plot_10X,"unassigned","other"))
Dim10Xbis2$data$cell_type1[is.na(Dim10Xbis2$data$cell_type1)] <- "other"
Dim10Xbis2 <- Dim10Xbis2 %>% formating
SupFig12A <- ggarrange(Dim10X,Dim10Xbis,Dim10Xbis2, labels = c("Original","CellID(G)","CellID(C)"),common.legend = T,legend = "right", ncol = 3, font.label = list(size =21))



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### DimPlot SS                                                             ####


cell_type_plot_SS <- c(head(rowterms,-8), "Enterocytes")
cells_to_plot <- colnames(seurat_ATAC_sub_SS)[seurat_ATAC_sub_SS$cell_type1 %in% cell_type_plot_SS]
DimSS <- DimPlot(seurat_ATAC_sub_SS, group.by = "cell_type1",cells = cells_to_plot,order = cells_to_plot, reduction = "tsne")
DimSS$data$cell_type1 <- factor(DimSS$data$cell_type1, levels = c(cell_type_plot_SS,"unassigned","other"))
DimSS <- DimSS %>%  formating

mapSSATAClabel <- map10X[seurat_ATAC_sub_SS$cell_type1]
ok <- mapply(x = mapSSATAClabel, y = seurat_ATAC_sub_SS$pred_CellID_G, z = seurat_ATAC_sub_SS$cell_type1, function(x,y,z){ifelse(y%in%x, z , last(names(mapSS)[sapply(mapSS, function(x,y){any(x %in% y)},y = y,simplify = T)]))})
ok[seurat_ATAC_sub_SS$pred_CellID_G == "unassigned"] <- "unassigned" 
ok <- ok[seurat_ATAC_sub_SS$cell_type1 %in% cell_type_plot_SS]
DimSSbis <- DimPlot(seurat_ATAC_sub_SS, group.by = "cell_type1", cells = cells_to_plot, reduction = "tsne")
DimSSbis$data$cell_type1 <- factor(ok, levels = c(cell_type_plot_SS,"unassigned","other"))
DimSSbis$data$cell_type1[is.na(DimSSbis$data$cell_type1)] <- "other"
DimSSbis <- DimSSbis %>% formating

ok <- mapply(x = mapSSATAClabel, y = seurat_ATAC_sub_SS$pred_CellID_C, z = seurat_ATAC_sub_SS$cell_type1, function(x,y,z){ifelse(y%in%x, z , last(names(mapSS)[sapply(mapSS, function(x,y){any(x %in% y)},y = y,simplify = T)]))})
ok[seurat_ATAC_sub_SS$pred_CellID_C == "unassigned"] <- "unassigned" 
ok <- ok[seurat_ATAC_sub_SS$cell_type1 %in% cell_type_plot_SS]
DimSSbis2 <- DimPlot(seurat_ATAC_sub_SS, group.by = "cell_type1", cells = cells_to_plot, reduction = "tsne")
DimSSbis2$data$cell_type1 <- factor(ok, levels = c(cell_type_plot_SS,"unassigned","other"))
DimSSbis2$data$cell_type1[is.na(DimSSbis2$data$cell_type1)] <- "other"
DimSSbis2 <- DimSSbis2 %>% formating

SupFig12 <- ggarrange(DimSS,DimSSbis,DimSSbis2,Dim10X,Dim10Xbis,Dim10Xbis2,nrow = 2,  ncol = 3,common.legend = T,legend = "right", labels = "AUTO", font.label = list(size =24)) + theme(legend.text = element_text(face = "bold"), axis.title = element_blank())
ggsave("../FinalFigure/SupFig12.pdf", SupFig12,width = 183,height = 130, units = "mm", scale = 2)
ggsave("../FinalFigure/SupFig12.png", SupFig12,width = 183,height = 130, units = "mm", scale = 2)



