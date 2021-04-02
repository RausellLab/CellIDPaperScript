library(patchwork)
WuDF <- FetchData(SeuratWu, vars = c("UMAP_1", "UMAP_2", "Airway.Brush.Tuft", "Intestinal.Brush.Tuft", "Gnat3", "Il25", "leukotriene biosynthetic process (GO:0019370)---GO Biological Process", "interleukin-17-mediated signaling pathway (GO:0097400)---GO Biological Process"))
FletcherDF <- FetchData(SeuratFletcher, vars = c("UMAP_1", "UMAP_2", "Airway.Brush.Tuft", "Intestinal.Brush.Tuft", "Gnat3", "Il25", "leukotriene biosynthetic process (GO:0019370)---GO Biological Process", "interleukin-17-mediated signaling pathway (GO:0097400)---GO Biological Process"))
colnames(WuDF)[7:8] <- c("Leukotriene_biosynthesis", "IL17_signaling")
colnames(FletcherDF)[7:8] <- c("Leukotriene_biosynthesis", "IL17_signaling")


# Legend ------------------------------------------------------------------

leg <-
  ggscatter(
    FletcherDF,
    x = "UMAP_1",
    y = "UMAP_2",
    size = 0.5,
    color = "Gnat3"
  ) + theme(
    aspect.ratio = 1,
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.size = unit(1.5, "line"),
    legend.text = element_text(size = 12, face = "bold")
  ) + scale_color_gradient(
    low = "grey",
    high = "red",
    labels = c("Low", "High"),
    breaks = c(0, 3)
  )
ggleg <- as_ggplot(get_legend(leg))


# Zoom Function -----------------------------------------------------------

ggzoomFig11 <- function(df, x, xmin, xmax, ymin, ymax, rect, line1, line2) {
  main <- ggscatter(df[order(df[[x]]), ], x = "UMAP_1", y = "UMAP_2", size = 0.5, color = x) + theme(aspect.ratio = 1, legend.position = "none") + scale_color_gradient(low = "grey", high = "red", limits = c(1, NA), na.value = "grey") 
  zoom <- main + coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) + ggplot2::theme_void() + theme(legend.position = "none", aspect.ratio = 1, plot.background = element_rect(colour = "black", size = 2))
  main <- main + annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = 0.2, color = "black", fill = "white")
  zoom$layers[[1]]$inherit.aes <- F
  gg <- main + annotation_custom(ggplotGrob(zoom), xmin = rect[1], xmax = rect[2], ymin = rect[3], ymax = rect[4])
  gg <- gg + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), title = element_text(size = 8, face = "bold"))
  gg <- gg + geom_segment(aes(x = line1[1], y = line1[2], xend = line1[3], yend = line1[4]), size = 0.2, alpha = 0.5)
  gg <- gg + geom_segment(aes(x = line2[1], y = line2[2], xend = line2[3], yend = line2[4]), size = 0.2, alpha = 0.5) + ggtitle(x)
  return(gg)
}

ggFletcher <- lapply(c("Airway.Brush.Tuft", "Gnat3", "Leukotriene_biosynthesis", "Intestinal.Brush.Tuft",  "Il25",  "IL17_signaling"), function(x) ggzoomFig11(FletcherDF, x, xmin = -2.6, xmax = -2.1, ymin = 2, ymax = 2.3, rect = c(4, 10, -2.5, -5), line1 = c(-2.1, 2, 5, -2.5), line2 = c(-2.1, 2, 5, -5)))
ggWu <- lapply(c("Airway.Brush.Tuft", "Gnat3", "Leukotriene_biosynthesis", "Intestinal.Brush.Tuft", "Il25", "IL17_signaling"), function(x) ggzoomFig11(WuDF, x, xmin = 3.4, xmax = 4.2, ymin = -10.6, ymax = -11.6, c(-16, -10, 14, 24), line1 = c(3.4, -10.6, -10, 24), line2 = c(3.4, -10.6, -10, 14)))

ggFletcher[[1]] <- ggFletcher[[1]] + geom_label(
  data = data.frame(label = paste0("n = ",c(37)), x = c(7), y = c(-5.4)),
  aes(x = x, y = y, label = label), label.size = 0.3, size = 2.5, label.padding = unit(0.1, "lines")
)

ggWu[[1]] <- ggWu[[1]] + geom_label(
  data = data.frame(label = paste0("n = ",c(5)), x = c(-13), y = c(13.5)),
  aes(x = x, y = y, label = label), label.size = 0.3, size = 2.5, label.padding = unit(0.1, "lines")
)

A <- ggpubr::ggarrange(plotlist = ggFletcher, nrow = 2, ncol = 3, labels = "A", font.label = list(size = 16))
B <- ggpubr::ggarrange(plotlist = ggWu, nrow = 2, ncol = 3, labels = "B", font.label = list(size = 16))
void <- ggplot() + geom_blank() + theme_void()

SupFig11 <- (A/void/B/ggleg) +plot_layout(ncol = 1, heights = c(1,0.1,1,0.2))

ggsave("../FinalFigure/SupFig11.pdf", SupFig11, dpi = 320, units = "mm", width = 180, height = 210)
ggsave("../FinalFigure/SupFig11.png", SupFig11, dpi = 320, units = "mm", width = 180, height = 210)
