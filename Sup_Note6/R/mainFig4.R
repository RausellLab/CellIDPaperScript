library(ComplexHeatmap)
library(circlize)
library(patchwork)

# Main A ------------------------------------------------------------------
formating <- function(x) x +scale_colour_manual(values = SupFig12Col, drop =F) + theme(plot.margin = unit(c(0.1,0,0,0.1), "lines"))+ guides(col = guide_legend(ncol = 1, override.aes = list(size = 4)))
cell_type_plot_10X <- c(head(rowterms,-8), "Enterocytes")
cells_to_plot <- colnames(seurat_ATAC_sub_10X)[seurat_ATAC_sub_10X$cell_type1 %in% cell_type_plot_10X]
SupFig12Col <- colorRampPalette(colors = RColorBrewer::brewer.pal(8,'Set1'))(length(cell_type_plot_10X))
SupFig12Col <- c(setNames(SupFig12Col,cell_type_plot_10X))
SupFig12Col <- c(SupFig12Col, other = "grey30", unassigned = "lightgrey")
Main4A1 <- DimPlot(seurat_ATAC_sub_10X, group.by = "cell_type1",cells = cells_to_plot, reduction = "tsne")
Main4A1$data$cell_type1 <- factor(Main4A1$data$cell_type1, levels = names(SupFig12Col))
Main4A1 <- Main4A1 %>%  formating


Main4A1 <- DimPlot(seurat_10X_sub, group.by = "cell_type1", reduction = "tsne")+ guides(color = F)
Main4A1$data$cell_type1 <- factor(Main4A1$data$cell_type1, levels = names(SupFig12Col)) 
Main4A1 <- Main4A1 %>%  formating


ok <- pbapply::pblapply(X = unique(seurat_10X_sub$cell_type1), FUN = function(y) names(map10X)[sapply(map10X, function(x) y %in% x)]) %>% set_names(unique(seurat_10X_sub$cell_type1))
ok[sapply(ok, function(x) is_empty(x))] <- "other"
ok$`immature T cell` <- "T cells"
ok$`T cell` <- "T cells"
ok$`lung endothelial cell` <- "Type I pneumocytes"
ok$`DN1 thymic pro-T cell` <- "T cells"
ok$`immature B cell` <- "Immature B cells"
ok$`B cell` <- "B cells"
ok$`early pro-B cell` <- "Immature B cells"
ok$`late pro-B cell` <- "Immature B cells"
ok$`duct epithelial cell` <- "Collecting duct"
ok$`kidney capillary endothelial cell` <- "Endothelial I (glomerular)"
ok$`kidney proximal straight tubule epithelial cell` <- "Proximal tubule"
ok$`kidney collecting duct epithelial cell` <- "Collecting duct"
ok$`Fraction A pre-pro B cell` <- "Immature B cells"
ok$promonocyte <- "Monocytes"
ok$`endothelial cell` <- "Endothelial I cells"
ok$`kidney cell` <- "Endothelial I (glomerular)"
ok$`type II pneumocyte` <- "Type II pneumocytes"
mp <- ok[seurat_10X_sub$cell_type1]
conv <- sapply(mp, function(x) sample(x,1))

seurat_10X_sub$cell_type2 <- conv

Main4A1 <- DimPlot(seurat_10X_sub, group.by = "cell_type2", reduction = "tsne")+ guides(color = F)
Main4A1$data$cell_type2 <- factor(conv, levels = c(cell_type_plot_10X,"unassigned","other"))
Main4A1 <- Main4A1 %>% formating



ok <- mapply(x = map10XATAClabel, y = seurat_ATAC_sub_10X$pred_CellID_G, z = seurat_ATAC_sub_10X$cell_type1, function(x,y,z){ifelse(y%in%x, z , last(names(map10X)[sapply(map10X, function(x,y){any(x %in% y)},y = y,simplify = T)]))})
ok[seurat_ATAC_sub_10X$pred_CellID_G == "unassigned"] <- "unassigned" 
ok <- ok[seurat_ATAC_sub_10X$cell_type1 %in% cell_type_plot_10X]
Main4A2 <- DimPlot(seurat_ATAC_sub_10X, group.by = "cell_type1", cells = cells_to_plot, reduction = "tsne")
Main4A2$data$cell_type1 <- factor(ok, levels = c(cell_type_plot_10X,"unassigned","other"))
Main4A2$data$cell_type1[is.na(Main4A2$data$cell_type)] <- "other"
Main4A2 <- Main4A2 %>% formating
Main4ALeg <- ggdraw(get_legend(Main4A2+ theme(legend.text = element_text(size =7, face = "bold"), legend.margin = margin(0,0,0,0), legend.key.size = unit(0,"mm"))+ guides(color = guide_legend(override.aes = list(size =2), ncol =1))))
Main4A1 <- Main4A1 +theme_bw() + theme(panel.grid = element_blank()) + theme(legend.position = "none",aspect.ratio = 1,axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(face = "bold", size =10))
Main4A2 <- Main4A2 +theme_bw() + theme(panel.grid = element_blank()) + theme(legend.position = "none",aspect.ratio = 1,axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(face = "bold", size =10))
void <- ggplot() + geom_blank() + theme_void()
Main4A <- ggarrange(Main4A1, void, Main4A2, nrow = 1, widths = c(2.5,1,2.5)) 
# Main B ------------------------------------------------------------------

Main4BDF <- ATACBenchOverall %>%  filter(data == "10X", metrics == "F1")
Main4B <- ggbarplot(Main4BDF, x = "methods", y = "value", fill = "methods", label = "value", lab.nb.digits = 2, lab.size = 2.5) + theme(
  axis.text.y            = element_text(size = 8, face = "bold", hjust = 5),
  axis.text.x            = element_text(size =8 , face = "bold", angle = 45, hjust = 1),
  axis.title             = element_blank(),
  legend.text            = element_text(size = 8, face = "bold"),
  legend.title           = element_blank(),
  legend.position = "none"
) + coord_cartesian(ylim = c(0, 1)) + scale_fill_manual(values = SupFig14Col) + guides(fill = guide_legend(nrow=4, override.aes = list(size = 0.5)), color = guide_legend(override.aes = list(size =0.5))) + geom_vline(xintercept =2.5)
Main4B$layers[[2]]$aes_params$fontface <- "bold"

Main4BLeg <- ggdraw(get_legend(ggscatter(Main4BDF, x = "methods", y = "value", fill = "methods", shape = 22)+ scale_fill_manual(values = SupFig14Col) + guides(fill = guide_legend(nrow=4, override.aes = list(size = 2))) + theme(legend.key.height = unit(0,"mm"), legend.key.width = unit(0,"mm"), legend.text = element_text(face = "bold", size = 8), legend.title = element_blank()))) 
# Main C ------------------------------------------------------------------

{
  contingency <- table(seurat_ATAC_sub_10X$cell_type1, factor(seurat_ATAC_sub_10X$pred_CellID_G, levels = unique(seurat_ATAC_sub_10X$pred_CellID_G)))
  contingency <- as.matrix.data.frame(contingency,rownames.force = T) %>% set_colnames(colnames(contingency))
  contingency <- (contingency/rowSums(contingency))
  contingency <- contingency[rowterms, colterms]
  colnames(contingency) <- str_remove(colnames(contingency), "cell")
  colnames(contingency) <- str_remove(colnames(contingency), " $")
  colnames(contingency) <- str_replace(colnames(contingency), "precursor", "pre")
  colnames(contingency) <- str_replace(colnames(contingency), "progenitors", "pro")
  rownames(contingency) <- str_replace(rownames(contingency), "progenitors", "pro")
  colnames(contingency) <- str_replace(colnames(contingency), "kidney loop of Henle ascending limb epithelial", "kidney loop of Henle")
  colnames(contingency) <- str_replace(colnames(contingency), "kidney proximal straight tubule epithelial", "proximal straight tubule")
  colnames(contingency) <- str_replace(colnames(contingency), "kidney collecting duct epithelial", "collecting duct")
  colnames(contingency) <- str_replace(colnames(contingency), "endothelial", "endo.")
  colnames(contingency) <- str_replace(colnames(contingency), "monocyte", "mono")
  rownames(contingency) <- str_replace(rownames(contingency), "Endothelial I \\(glomerular\\)", "Endo glomerular")
  rownames(contingency) <- str_replace(rownames(contingency), "Distal convoluted tubule", "DCT")
  rownames(contingency) <- str_replace(rownames(contingency), "Alveolar macrophages", "Alveolar Macro")
  rownames(contingency) <- str_replace(rownames(contingency), "Type I pneumocytes", "Pneumocytes I")
  rownames(contingency) <- str_replace(rownames(contingency), "Type II pneumocytes", "Pneumocytes II")
  rownames(contingency) <- str_remove(rownames(contingency), "cells")
  Main4C <-
    ComplexHeatmap::Heatmap(
      name = "%",
      contingency[1:26, ],
      cluster_rows = F,
      cluster_columns = F, column_names_rot = 45,
      col = rev(heat.colors(101)),
      row_title = "Original Curated Labels",
      row_title_gp = gpar(fontsize = 10, fontface ="bold"),
      column_title = "CellID predictions",
      column_title_gp = gpar(fontsize = 10, fontface ="bold"),
      row_names_gp = gpar(fontsize = 7, fontface ="bold"),
      column_names_gp = gpar(fontsize = 7,fontface ="bold"),
      row_names_side = "left",
      border = T,
      rect_gp = gpar(col = "grey", lwd = 1), show_heatmap_legend = F,
      column_title_side = "bottom", heatmap_legend_param = list()
    )
}

lgd = Legend(col_fun = colorRamp2(seq(0,1, length.out = 101), rev(heat.colors(101))), labels_gp = gpar(fontsize = 7,fontface ="bold"), at = c(0,  0.2,0.4,0.6,0.8, 1), border = "black", legend_height = unit(4, "cm"), title = "%")
Leg4C <- ggdraw(grid.grabExpr(draw(lgd)))

dev.off()
# Main all ----------------------------------------------------------------
annotA <- ggplot() + geom_text(aes(x = 0,y = 0, label = "A", fontface = "bold"), size = 6) + theme_void()
annotB <- ggplot() + geom_text(aes(x = 0,y = 0, label = "B", fontface = "bold"), size = 6) + theme_void()
annotC <- ggplot() + geom_text(aes(x = 0,y = 0, label = "C", fontface = "bold"), size = 6) + theme_void()
annotD <- ggplot() + geom_text(aes(x = 0,y = 0, label = "D", fontface = "bold"), size = 6) + theme_void()
annotF1 <- ggplot() + geom_text(aes(x = 0,y = 0, label = "F1 Score", fontface = "bold"), size = 3.7) + theme_void()
annottSNE1 <- ggplot() + geom_text(aes(x = 0,y = 0, label = "tSNE1", fontface = "bold"), size = 3.7) + theme_void()
annottSNE2 <- ggplot() + geom_text(aes(x = 0,y = 0, label = "tSNE2", fontface = "bold"), size = 3.7, angle = 90) + theme_void()
annotLine <- ggplot() + geom_hline(yintercept =1) + theme_void()

void <- ggplot(Main4BDF) + geom_blank() + theme_void()
design  <-
  c(
    Main4A = area(t = 1, l = 2, b = 10, r = 16),
    Main4ALeg = area(t = 1,l = 18,b = 10,r = 19),
    Main4B = area(t = 14, l = 15, b = 19, r = 20),
    Main4BLeg = area(t = 13, l = 15, b = 13, r = 20),
    Leg4C = area(t = 13, l = 14, b = 18, r = 14),
    Main4C = area(t = 13, l = 0, b = 20, r = 13),
    annotA = area(t = 1, l = 1, b = 1,r = 1),
    annotB = area(t = 1, l = 10, b = 1,r = 11),
    annotC = area(t = 11, l = 1, b = 12,r = 1),
    annotD = area(t = 11, l = 14, b = 12,r = 15),
    annotF1 = area(t = 20, l = 17, b = 20,r = 18),
    annottSNE1 = area(t = 9, l = 9, b = 11,r = 10),
    annottSNE2 = area(t = 1, l = 1, b = 10,r = 2),
    annotLine = area(t = 9, l = 4, b = 10,r = 14)
  )
Main4 <- (  
            Main4A + 
            Main4ALeg + 
            Main4B + 
            Main4BLeg + 
            Leg4C + 
            grid.grabExpr(draw(Main4C))  + 
            annotA + 
            annotB + 
            annotC +
            annotD +   
            annotF1 + 
            annottSNE1 + 
            annottSNE2 + 
            annotLine) + plot_layout(design = design)
ggsave("../FinalFigure/Main/Main4.pdf", plot = Main4, width =230, height = 200, units= "mm")
ggsave("../FinalFigure/Main/Main4.png", plot = Main4, width =230, height = 200, units= "mm", dpi = 600)
ggsave("../FinalFigure/Main/Main4.svg", plot = Main4, width =230, height = 200, units= "mm", dpi = 600)
