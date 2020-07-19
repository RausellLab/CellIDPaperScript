ThemeSupFig15 <- theme(
  axis.text.x = element_text(angle = 45,hjust = 1,face = "bold",size = 8),axis.text.y = element_text(face = "bold", size =8),aspect.ratio = 1,axis.title = element_text(face = "bold", size =8),legend.text = element_text(face = "bold"),legend.title = element_text(face = "bold"))
SupFig15Col <- setNames(RColorBrewer::brewer.pal(name = "Paired", n = 12), c("scmap_cluster","scmap_cell","Seurat","MNN", "CellID(G)","CellID(C)","SCN","scPred","CHETAH","CaSTLe","scID","SingleR"))

read_rds("data/QueryVarDF.rds")
read_rds("data/RefVarDF.rds")
SupFig15A <- ggplot(QueryVarDF, aes(x = ncells, y = time, color = methods)) + geom_line() + geom_point() + scale_color_manual(values = SupFig15Col) + scale_x_continuous(breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000), trans = "log10") + ThemeSupFig15  + ylab("Time(s)") + xlab("Number of cells") + guides(color = guide_legend(ncol = 1))+ yscale("log10")
SupFig15B <- ggplot(QueryVarDF, aes(x = ncells, y = memory, color = methods)) + geom_line() + geom_point() + scale_color_manual(values = SupFig15Col) + scale_x_continuous(breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000), trans = "log10") + ThemeSupFig15  + ylab("Peak Memory (Mb)") + xlab("Number of cells") + guides(color = guide_legend(ncol = 1))+ yscale("log10")
SupFig15C <- ggplot(RefVarDF, aes(x = ncells, y = time, color = methods)) + geom_line() + geom_point() + scale_color_manual(values = SupFig15Col) + scale_x_continuous(breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000), trans = "log10") + ThemeSupFig15 + ylab("Time(s)") + xlab("Number of cells") + guides(color = guide_legend(ncol = 1))+ yscale("log10")
SupFig15D <- ggplot(RefVarDF, aes(x = ncells, y = memory, color = methods)) + geom_line() + geom_point() + scale_color_manual(values = SupFig15Col) + scale_x_continuous(breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000), trans = "log10") + ThemeSupFig15 + ylab("Peak Memory (Mb)") + xlab("Number of cells") + guides(color = guide_legend(ncol = 1))+ yscale("log10")
SupFig15 <- ggarrange(SupFig15A,SupFig15B, SupFig15C, SupFig15D, common.legend = T, 
                      labels = "AUTO", ncol = 2, nrow = 2, legend = "right")


# Gather ------------------------------------------------------------------

ggsave(filename = "../FinalFigure/SupFig15.pdf", SupFig15, dpi = 600, units = "mm", width = 180, height = 210)
ggsave(filename = "../FinalFigure/SupFig15.png", SupFig15, dpi = 600, units = "mm", width = 180, height = 210)
