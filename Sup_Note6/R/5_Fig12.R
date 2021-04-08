SupFig12Col <- RColorBrewer::brewer.pal(name = "Paired", n = 12)
SupFig12Col <- setNames(SupFig12Col, c("scmap_cluster","scmap_cell","Seurat","MNN", "CellID(G)","CellID(C)","SCN","scPred","CHETAH","CaSTLe","scID","SingleR"))

SupFig12 <- ggbarplot(ATACBenchOverall, x = "methods", y = "value", facet.by = c("data", "metrics"), fill = "methods",label   = T,
          lab.size = 2.5, lab.nb.digits = 2) + theme(
    aspect.ratio           = 1,
    axis.text.y            = element_text(size = 10, face = "bold"),
    axis.text.x            = element_blank(),
    axis.ticks.x           = element_blank(),
    axis.title             = element_blank(),
    strip.text             = element_text(size = 12, face = "bold"),
    legend.text            = element_text(size = 12, face = "bold"),
    legend.title           = element_blank()
  ) + coord_cartesian(ylim = c(0, 1)) + scale_fill_manual(values = SupFig12Col) + guides(fill = guide_legend(nrow=2)) + geom_vline(xintercept =2.5)
SupFig12$layers[[2]]$aes_params$fontface <- "bold"


ggsave("../FinalFigure/SupFig12.pdf", SupFig12,  dpi = 600, units = "mm", width = 180, height =180)
ggsave("../FinalFigure/SupFig12.png", SupFig12,  dpi = 600, units = "mm", width = 180, height =180)

