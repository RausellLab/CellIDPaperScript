
IntestinalBenchCell <- readRDS("data/IntestinalBenchCell.rds")
IntestinalBenchOverall <- readRDS("data/IntestinalBenchOverall.rds")


ggRareCell <- IntestinalBenchOverall %>%  filter(metrics == "F1", data == "Plasschaert Mouse") %>%  ggbarplot(x = "methods", y = "value", fill = "lightgrey", label = T, lab.size = 2, lab.nb.digits = 2) + theme_bw() + theme(
  aspect.ratio = 1,
  axis.title = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1),
  strip.text.x = element_text(size = 8.2, face = "bold"),
  strip.text.y = element_text(size = 8.2, face = "bold"),
  axis.text = element_text(size = 7, face = "bold", colour = "black"),
  legend.text = element_text(size = 8.2, face = "bold"),
  legend.title = element_blank()
)  + scale_y_continuous(
  breaks =
    c(0, 0.25, 0.5, 0.75, 1),
  limits = c(0, 1.01)
) + geom_vline(xintercept = 2.5) + ggtitle("Prediction F1 score")

ggRejection <- IntestinalBenchCell %>% filter(cell_type == "Rejection") %>% filter(metrics == "Recall", data == "Plasschaert Mouse") %>%  ggbarplot(x = "methods", y = "value", fill = "lightgrey", label = T, lab.size = 2, lab.nb.digits = 2) + theme_bw() + theme(
  aspect.ratio = 1,
  axis.title = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1),
  strip.text.x = element_text(size = 8.2, face = "bold"),
  strip.text.y = element_text(size = 8.2, face = "bold"),
  axis.text = element_text(size = 7, face = "bold", colour = "black"),
  legend.text = element_text(size = 8.2, face = "bold"),
  legend.title = element_blank()
)  + scale_y_continuous(
  breaks =
    c(0, 0.25, 0.5, 0.75, 1),
  limits = c(0, 1.01)
) + geom_vline(xintercept = 2.5) + ggtitle("Succesful rejection rate")

ggsave(plot = ggRareCell, filename = "../FinalFigure/Main3D_review.pdf", width = 5, height = 5, dpi = 300)
ggsave(plot = ggRareCell, filename = "../FinalFigure/Main3D_review.png", width = 5, height = 5, dpi = 300)

# Update 10C ---------------------------------------------------------------
SupFig10BarCol <- RColorBrewer::brewer.pal(name = "Paired", n = 12)
SupFig10BarCol <- setNames(SupFig10BarCol, c("scmap_cluster", "scmap_cell", "Seurat", "MNN", "CellID(G)", "CellID(C)", "SCN", "scPred", "CHETAH", "CaSTLe", "scID", "SingleR"))

SupFig10Bar <-
  ggbarplot(
    IntestinalBenchOverall,
    x = "methods",
    y = "value",
    fill = "methods",
    facet.by = c("metrics", "data"),
    label = T, lab.nb.digits = 2, 
    lab.size = 2.5,
  ) + theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    strip.text.x = element_text(size = 8.2, face = "bold"),
    strip.text.y = element_text(size = 8.2, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 8.2, face = "bold"),
    legend.text = element_text(size = 8.2, face = "bold"),
    legend.title = element_blank()
  ) + scale_fill_manual(values = SupFig10BarCol) + scale_y_continuous(
    breaks =
      c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1.05)
  ) + geom_vline(xintercept = 2.5) + guides(fill = guide_legend(nrow = 2))


df_Rejection <- IntestinalBenchCell%>%
  mutate(methods = factor(methods, c("CellID(G)", "CellID(C)", "scmap_cluster", "scmap_cell", "Seurat", "MNN", "SCN", "scPred", "CHETAH", "CaSTLe", "scID", "SingleR"))) %>% filter(cell_type == "Rejection") %>% filter(metrics == "Recall") 
df_Rejection$cell_type <- NULL
df_Rejection$metrics <- "Rejection Rate"

SupFig10C_review <- ggbarplot(
  rbind(IntestinalBenchOverall, df_Rejection),
  x = "methods",
  y = "value",
  fill = "methods",
  facet.by = c("metrics", "data"),
  label = T, lab.nb.digits = 2, 
  lab.size = 2.5,
) + theme(
  axis.title = element_blank(),
  axis.text.x = element_blank(),
  strip.text.x = element_text(size = 8.2, face = "bold"),
  strip.text.y = element_text(size = 8.2, face = "bold"),
  axis.ticks.x = element_blank(),
  axis.text = element_text(size = 8.2, face = "bold"),
  legend.text = element_text(size = 8.2, face = "bold"),
  legend.title = element_blank()
) + scale_fill_manual(values = SupFig10BarCol) + scale_y_continuous(
  breaks =
    c(0, 0.25, 0.5, 0.75, 1),
  limits = c(0, 1.05)
) + geom_vline(xintercept = 2.5) + guides(fill = guide_legend(nrow = 2))


ggsave(plot = SupFig10C_review, filename = "../FinalFigure/SupFig10C_review.pdf", width = 10, height = 8, dpi = 300)
ggsave(plot = SupFig10C_review, filename = "../FinalFigure/SupFig10C_review.png", width = 10, height = 8, dpi = 300)
