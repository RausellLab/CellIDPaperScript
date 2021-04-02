theme_supFig4 <- theme(
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


reapCol <- palCol[unique(sort(ImmuneBenchREAP$cell_type))]
names(reapCol) <- unique(sort(ImmuneBenchREAP$cell_type))
SupFig4A  <- DimPlot(reap_seurat, label = T,group.by = "cell_type1",label.size = 3) + scale_color_manual(values = c(reapCol, unknown = "black"))  + theme(aspect.ratio=1) + theme_supFig4
SupFig4A$layers[[1]]$show.legend <- F
SupFig4A$layers[[2]]$aes_params$fontface <-"bold"
SupFig4A$layers[[2]]$data <- SupFig4A$layers[[2]]$data %>% filter(cell_type1 != "unknown")
reap_match_gs_df <- rownames_to_column(as.data.frame(unlist(map_reap)))
colnames(reap_match_gs_df) <- c("type", "gs")
reap_match_gs_df$type <- str_remove(reap_match_gs_df$type, "\\d$")
reap_converter <- c(setNames(reap_match_gs_df$type,reap_match_gs_df$gs), unassigned = "unassigned")
reap_seurat$SupFig4BCol <- reap_converter[reap_seurat$pred_CellID]
reap_seurat$SupFig4BCol[is.na(reap_seurat$SupFig4BCol)] <- "other"
reap_seurat$SupFig4BCol <- factor(reap_seurat$SupFig4BCol, levels = c(names(reapCol),"other","unassigned"))
SupFig4B <- DimPlot(reap_seurat, group.by = "SupFig4BCol", label = T, label.size = 3) + theme(aspect.ratio=1) + scale_color_manual(values = c(reapCol, other = "black", unknown = "black", unassigned = "lightgrey"))+ theme_supFig4
SupFig4B$layers[[2]]$data <- SupFig4B$layers[[2]]$data %>% filter(!SupFig4BCol %in% c("unassigned", "other")) 
SupFig4B$layers[[2]]$aes_params$fontface <- "bold"
SupFig4Legend <- ggarrange(get_legend(SupFig4B))
SupFig4B$layers[[1]]$show.legend <- F

SupFig4C <-
  Seurat::DoHeatmap(
    reap_seurat,
    features = rownames(reap_HGT),
    slot = "data",
    assay = "ImSig",
    raster = T,
    group.by = "cell_type1" ,
    group.colors = c(reapCol, "black"),
    disp.max = 10,
    size = 3.5,lines.width = 110
  ) + scale_fill_gradientn(values=c(1, .2, .1, 0), colours=c("#D73027",  "#FEDF8F", "#4E7DB8", "#4575B4"), breaks = seq(0,10,2), na.value = 'lightgrey') + theme(plot.margin = margin(0.8,0.5,0,0.5, "cm")) + theme(
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    legend.position = "bottom"
  ) 

SupFig4C$labels$fill <- "-Log10 pval"
SupFig4C$layers[[2]]$show.legend <- F
SupFig4C$layers[[4]]$aes_params$fontface <-  "bold"
SupFig4C$layers[[3]]$geom_params$raster <- ifelse(SupFig4C$layers[[3]]$geom_params$raster == "#FFFFFF", "lightgrey",SupFig4C$layers[[3]]$geom_params$raster)
SupFig4C$layers[[3]]$geom_params$ymin <- 19.5
SupFig4C$layers[[4]]$data$label.x.pos <- SupFig4C$layers[[4]]$data$label.x.pos - 80
SupFig4C <- SupFig4C + geom_rect(mapping=aes(xmin=length(SupFig4C$layers[[3]]$geom_params$raster)-110, xmax=length(SupFig4C$layers[[3]]$geom_params$raster), ymin=0.5, ymax=20.5), fill = "white")
SupFig4C$layers <- SupFig4C$layers[c(1,2,3,5,4)]

SupFig4 <- ggarrange(ggarrange(SupFig4A, SupFig4B), SupFig4C, ncol = 1) + geom_text(data = data.frame(label = c("A","B","C"), x = c(0.05,0.55,0.05), y = c(0.95,0.95,0.5)), aes(x=x, y=y, label = label), size = 7, fontface = "bold")

ggsave("../FinalFigure/supFig4.pdf",SupFig4,width = 10, height = 10, dpi = 320)
ggsave("../FinalFigure/supFig4.png",SupFig4,width = 10, height = 10, dpi = 320)
