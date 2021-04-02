# Contingency table -------------------------------------------------------
rowterms <- c(
  "Hematopoietic progenitors",
  "Alveolar macrophages",
  "Dendritic cells",
  "Macrophages",
  "Monocytes",
  "Immature B cells", 
  "Activated B cells",
  "B cells",
  "T cells",
  "Regulatory T cells",
  "NK cells",
  "Erythroblasts",
  "Endothelial I cells",
  "Endothelial II cells",
  "Hepatocytes",
  "Endothelial I (glomerular)",
  "Collecting duct",
  "DCT/CD",
  "Distal convoluted tubule",
  "Proximal tubule",
  "Proximal tubule S3",
  "Loop of henle",
  "Cardiomyocytes",
  "Type I pneumocytes",
  "Type II pneumocytes",
  "Podocytes",
  "Sperm",
  "Purkinje cells",
  "Microglia",
  "Ex. neurons SCPN",
  "Cerebellar granule cells",
  "Oligodendrocytes",
  "Inhibitory neurons",
  "Enterocytes"
)
colterms <-
  c(  "hematopoietic precursor cell",
      "leukocyte",
      "granulocytopoietic cell",
      "granulocyte",
      "myeloid cell",
      "alveolar macrophage",
      "basophil",
      "mast cell",
      "non-classical monocyte",
      "classical monocyte",
      "monocyte",
      "macrophage",
      "early pro-B cell",
      "late pro-B cell",
      "immature B cell",
      "Fraction A pre-pro B cell",
      "B cell",
      "immature T cell",
      "T cell",
      "DN1 thymic pro-T cell",
      "natural killer cell",
      "erythroblast",
      'proerythroblast',
      "dendritic cell",
      "endothelial cell",
      "lung endothelial cell",
      "hepatocyte",
      "endothelial cell of hepatic sinusoid" ,
      "kidney cell",
      "kidney capillary endothelial cell",
      "duct epithelial cell",
      "kidney collecting duct epithelial cell",
      "kidney proximal straight tubule epithelial cell",
      "kidney loop of Henle ascending limb epithelial cell",
      "cardiac muscle cell",
      "endocardial cell",
      "type II pneumocyte",
      "stromal cell",
      "fibroblast",
      "mesangial cell",
      "unassigned"
  )

{
  contingency <- table(seurat_ATAC_sub_10X$cell_type1, factor(seurat_ATAC_sub_10X$pred_CellID_G, levels = unique(seurat_ATAC_sub_10X$pred_CellID_G)))
  contingency <- as.matrix.data.frame(contingency,rownames.force = T) %>% set_colnames(colnames(contingency))
  contingency <- (contingency/rowSums(contingency))
  contingency <- contingency[rowterms, colterms]
  colnames(contingency) <- str_remove(colnames(contingency), "cell")
  SupFig13A <-  ggdraw(pheatmap(contingency,cluster_rows = F,cluster_cols = F,color = rev(heat.colors(101)),angle_col = 45, cellheight = 7, cellwidth = 7, silent = T,fontsize_col = 6.5,fontsize_row = 6.5,legend = F)$gtable)
  }

{ 
  contingency <- table(seurat_ATAC_sub_10X$cell_type1, factor(seurat_ATAC_sub_10X$pred_CellID_C,levels = unique(seurat_ATAC_sub_10X$pred_CellID_G)))
  contingency <- as.matrix.data.frame(contingency,rownames.force = T) %>% set_colnames(colnames(contingency))
  contingency <- (contingency/rowSums(contingency))
  contingency <- contingency[c(head(rowterms,-8), "Enterocytes"), colterms]
  colnames(contingency) <- str_remove(colnames(contingency), "cell")
  SupFig13B <-  ggdraw(pheatmap(contingency,cluster_rows = F,cluster_cols = F,color = rev(heat.colors(101)),angle_col = 45,cellheight = 7, cellwidth = 7, silent = T,fontsize_col = 6.5,fontsize_row = 6.5,legend = F)$gtable)
}


lgd = Legend(col_fun = colorRamp2(seq(0,1, length.out = 101), rev(heat.colors(101))), labels_gp = gpar(fontsize = 7,fontface ="bold"), at = c(0,  0.2,0.4,0.6,0.8, 1), border = "black", legend_height = unit(4, "cm"), title = "%")
Leg4C <- ggdraw(grid.grabExpr(draw(lgd)))

design <- c(area(1,2,7,10), area(1,11,5,11))
SupFig13 <- SupFig13B + Leg4C + plot_layout(design = design)

ggsave(filename = "../FinalFigure/SupFig13.pdf", SupFig13, width = 7, height = 5)
ggsave(filename = "../FinalFigure/SupFig13.png", SupFig13, width = 7, height = 5)
