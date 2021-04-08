source("R/utilitary_function.R")
SeuratEpithelial <- readRDS("../SupNote3/data/SeuratEpithelial.rds")
SeuratIntestinal <- readRDS("data/SeuratIntestinal.rds")
IntestinalBenchOverall <- readRDS("data/IntestinalBenchOverall.rds")
haber_umap <- SeuratEpithelial[[4]]@reductions$tsne@cell.embeddings %>%  as.data.frame() %>%
  mutate(
    cell_type1 = SeuratEpithelial[[4]]$cell_type1,
    `Montoro_ref_CellID(G)` = SeuratIntestinal$pred_CellID_G__mon,
    `Montoro_ref_CellID(C)` = SeuratIntestinal$pred_CellID_C__mon,
    `Plasschaert_Mouse_ref_CellID(G)` = SeuratIntestinal$pred_CellID_G__pla,
    `Plasschaert_Mouse_ref_CellID(C)` = SeuratIntestinal$pred_CellID_C__pla
  ) %>%
  gather("methods", "predictions", -(1:3)) %>%
  separate("methods",
    into = c("reference", "methods"),
    sep = "_ref_"
  ) %>%
  as_tibble()

# UMAP --------------------------------------------------------------------
Epithelial_celltype <- c("Basal", "Krt4.13", "Secretory", "Goblet", "Ciliated", "Brush.Tuft", "Brush+PNEC", "PNEC", "Ionocytes", "unassigned", "other")
Epithelial_color <- c("firebrick1", "purple", "dodgerblue", "cyan4", "darkgreen", "chocolate4", "chocolate2", "orange", "cyan", "lightgrey", "black")
EpithelialPalette <- setNames(Epithelial_color, Epithelial_celltype)
SupFig8UMAPCol <- c(Enterocyte = "grey10", TA = "grey40", Stem = "grey70", Paneth = "grey90", Endocrine = "orange", EpithelialPalette)

gg_haber_predictions <- ggscatter(
  haber_umap,
  x = "tSNE_1",
  y =  "tSNE_2",
  color = "predictions",
  size = 0.01,
  facet.by = c("methods", "reference")
) +
  theme(
    aspect.ratio = 1,
    strip.text.x = element_text(size = 8.2, face = "bold"),
    strip.text.y = element_text(size = 8.2, face = "bold"),
    axis.title = element_text(size = 8.2, face = "bold"),
    axis.text = element_text(size = 8.2, face = "bold"),
    legend.text = element_text(size = 8.2, face = "bold"),
    legend.title = element_blank(),
    legend.key.height  =  unit(1, "mm")
  ) + scale_color_manual(values = SupFig8UMAPCol, drop = F) + guides(color = guide_legend(override.aes = list(size = 3), nrow = 3))

haber_umap$cell_type1 <- factor(haber_umap$cell_type1, c("Brush.Tuft","Enterocyte", "TA", "Endocrine", "Stem", "Paneth", "Goblet"))
gg_haber_author <- ggscatter(
  haber_umap,
  x = "tSNE_1",
  y =  "tSNE_2",
  color = "cell_type1",
  size = 0.01
) +
  theme(
    aspect.ratio = 1,
    strip.text.x = element_text(size = 8.2, face = "bold"),
    strip.text.y = element_text(size = 8.2, face = "bold"),
    axis.title = element_text(size = 8.2, face = "bold"),
    axis.text = element_text(size = 8.2, face = "bold"),
    legend.text = element_text(size = 8.2, face = "bold"),
    legend.title = element_blank(),
    legend.key.height  =  unit(1, "mm")
  ) + scale_color_manual(values = SupFig8UMAPCol, drop = F) + guides(color = guide_legend(override.aes = list(size = 3), nrow = 3)) 


# Barplot -----------------------------------------------------------------

SupFig8BarCol <- RColorBrewer::brewer.pal(name = "Paired", n = 12)
SupFig8BarCol <- setNames(SupFig8BarCol, c("scmap_cluster", "scmap_cell", "Seurat", "MNN", "CellID(G)", "CellID(C)", "SCN", "scPred", "CHETAH", "CaSTLe", "scID", "SingleR"))

SupFig8Bar <-
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
  ) + scale_fill_manual(values = SupFig8BarCol) + scale_y_continuous(
    breaks =
      c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1.05)
  ) + geom_vline(xintercept = 2.5) + guides(fill = guide_legend(nrow = 2))


df_Rejection <- IntestinalBenchCell%>%
  mutate(methods = factor(methods, c("CellID(G)", "CellID(C)", "scmap_cluster", "scmap_cell", "Seurat", "MNN", "SCN", "scPred", "CHETAH", "CaSTLe", "scID", "SingleR"))) %>% filter(cell_type == "Rejection") %>% filter(metrics == "Recall") 
df_Rejection$cell_type <- NULL
df_Rejection$metrics <- "Rejection Rate"

SupFig8C_review <- ggbarplot(
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
) + scale_fill_manual(values = SupFig8BarCol) + scale_y_continuous(
  breaks =
    c(0, 0.25, 0.5, 0.75, 1),
  limits = c(0, 1.05)
) + geom_vline(xintercept = 2.5) + guides(fill = guide_legend(nrow = 2))


ggsave(
  "../FinalFigure/SupFig8.pdf",
  SupFig8,
  width = 183,
  height = 240,
  units = "mm",
  dpi = 320
)

ggsave(
  "../FinalFigure/SupFig8.png",
  SupFig8,
  width = 183,
  height = 240,
  units = "mm",
  dpi = 900
)

