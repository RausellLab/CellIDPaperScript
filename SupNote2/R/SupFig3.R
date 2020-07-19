theme_supFig3 <- theme(
  title = element_text(size = 10, face = "bold.italic"),
  legend.text = element_text(size = 12, face = "bold"),
  legend.title = element_text(size = 12, face = "bold"),
  axis.title.x = element_text(size = 10, face = "bold"),
  axis.title.y = element_text(size = 10, face = "bold"),
  strip.text = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(size = 10, face = "bold"),
  axis.text.y = element_text(size = 10, face = "bold"),
  legend.title.align = 0.5,
  legend.position = "top",
  plot.margin = unit(x = c(0.5,0,0.3,1), units = "cm")
)


citeCol <- palCol[unique(sort(ImmuneBenchCITE$cell_type))]
names(citeCol) <- unique(sort(ImmuneBenchCITE$cell_type))
SupFig3A  <- DimPlot(cite_seurat, label = T,group.by = "cell_type1",label.size = 3) + scale_color_manual(values = c(citeCol, unknown = "darkgray"))  + theme(aspect.ratio=1) + theme_supFig3
SupFig3A$layers[[1]]$show.legend <- F
SupFig3A$layers[[2]]$aes_params$fontface <-"bold"
cite_match_gs_df <- map_cite %>% unlist %>% as.data.frame() %>% rownames_to_column()
colnames(cite_match_gs_df) <- c("type", "gs")
cite_match_gs_df$type <- str_remove(cite_match_gs_df$type, "\\d$")
cite_converter <- c(setNames(cite_match_gs_df$type,cite_match_gs_df$gs), unassigned = "unassigned")
cite_seurat$SupFig3BCol <- cite_converter[cite_seurat$pred_CellID]
cite_seurat$SupFig3BCol[is.na(cite_seurat$SupFig3BCol)] <- "other"
cite_seurat$SupFig3BCol <- factor(cite_seurat$SupFig3BCol, levels = c(names(citeCol),"other","unassigned"))
SupFig3B <- DimPlot(cite_seurat, group.by = "SupFig3BCol", label = T, label.size = 3) + theme(aspect.ratio=1) + scale_color_manual(values = c(citeCol, other = "black", unassigned = "lightgrey"))+ theme_supFig3
SupFig3B$layers[[2]]$data <- SupFig3B$layers[[2]]$data %>% filter(!SupFig3BCol %in% c("unassigned", "other")) 
SupFig3B$layers[[2]]$aes_params$fontface <- "bold"
SupFig3Legend <- ggarrange(get_legend(SupFig3B))
SupFig3B$layers[[1]]$show.legend <- F

SupFig3C <-
  Seurat::DoHeatmap(
    cite_seurat,
    features = rownames(cite_HGT),
    slot = "data",
    assay = "ImSig",
    raster = T,
    group.by = "cell_type1" ,
    group.colors = citeCol,
    disp.max = 10,
    size = 3.5,lines.width = 110
  ) + scale_fill_gradientn(values=c(1, .2, .1, 0), colours=c("#D73027",  "#FEDF8F", "#4E7DB8", "#4575B4"), breaks = seq(0,10,2), na.value = 'lightgrey') + theme(plot.margin = margin(0.8,0.5,0,0.5, "cm")) + theme(
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    legend.position = "bottom"
  ) 

SupFig3C$labels$fill <- "-Log10 pval"
SupFig3C$layers[[2]]$show.legend <- F
SupFig3C$layers[[4]]$aes_params$fontface <-  "bold"
SupFig3C$layers[[3]]$geom_params$raster <- ifelse(SupFig3C$layers[[3]]$geom_params$raster == "#FFFFFF", "lightgrey",SupFig3C$layers[[3]]$geom_params$raster)
SupFig3C$layers[[3]]$geom_params$ymin <- 21.45
SupFig3C$layers[[4]]$data$label.x.pos <- SupFig3C$layers[[4]]$data$label.x.pos - 110
SupFig3C <- SupFig3C + geom_rect(mapping=aes(xmin=length(SupFig3C$layers[[3]]$geom_params$raster)-110, xmax=length(SupFig3C$layers[[3]]$geom_params$raster), ymin=0.5, ymax=22.5), fill = "white")
SupFig3C$layers <- SupFig3C$layers[c(1,2,3,5,4)]

SupFig3 <- ggarrange(ggarrange(SupFig3A, SupFig3B), SupFig3C, ncol = 1) + geom_text(data = data.frame(label = c("A","B","C"), x = c(0.05,0.55,0.05), y = c(0.95,0.95,0.5)), aes(x=x, y=y, label = label), size = 7, fontface = "bold")

ggsave("../FinalFigure/supFig3.pdf",SupFig3,width = 10, height = 10, dpi = 320)
ggsave("../FinalFigure/supFig3.png",SupFig3,width = 10, height = 10, dpi = 320)
