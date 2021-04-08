# -------------------------------------------------------------------------
# SupFig7 -----------------------------------------------------------------
# -------------------------------------------------------------------------

source("R/utilitary_function.R")
SeuratPancreas <- readRDS("data/SeuratPancreas.rds")
SeuratEpithelial <- readRDS("data/SeuratEpithelial.rds")

Pancreas_celltype <- c("alpha","beta", "delta", "gamma", "epsilon", "acinar", "ductal",  "PSC", "endothelial", "macrophage", "mast", "t_cell", "schwann", "unassigned", "other")
Pancreas_color <- c("firebrick1", "coral", "dodgerblue", "cyan", "cyan4","chocolate", "chocolate4",  "olivedrab","lightgreen","purple", "violet", "pink", "gold", "lightgrey", "black")
PancreasPalette <- setNames(Pancreas_color, Pancreas_celltype)

Fig7_muraro <- FetchData(SeuratPancreas[[2]], vars = c("tSNE_1", "tSNE_2"))
Fig7_muraro$Author <- factor(SeuratPancreas[[2]]$cell_type1, Pancreas_celltype)
Fig7_muraro$CellID_C <- factor(SeuratPancreas[[2]]$pred_CellID_G, Pancreas_celltype)
Fig7_muraro$CellID_G <- factor(SeuratPancreas[[2]]$pred_CellID_C, Pancreas_celltype)
Fig7_muraro$data <- "Muraro"
Fig7_muraro <- Fig7_muraro %>%  gather("num","cell_type",3:5)

Fig7_seger <- FetchData(SeuratPancreas[[3]], vars = c("tSNE_1", "tSNE_2"))
Fig7_seger$Author <- factor(SeuratPancreas[[3]]$cell_type1, Pancreas_celltype)
Fig7_seger$CellID_C <- factor(SeuratPancreas[[3]]$pred_CellID_G, Pancreas_celltype)
Fig7_seger$CellID_G <- factor(SeuratPancreas[[3]]$pred_CellID_C, Pancreas_celltype)

Fig7_seger$data <- "Segerstolpe"
Fig7_seger <- Fig7_seger %>%  gather("num","cell_type",3:5)

Fig7_pancreas <- rbind(Fig7_muraro,Fig7_seger)
Fig7_pancreas$cell_type[is.na(Fig7_pancreas$cell_type)] <- "other"
Fig7_pancreas$cell_type <- factor(Fig7_pancreas$cell_type, levels = Pancreas_celltype)
Fig7_pancreas$num <- str_replace(Fig7_pancreas$num, "CellID_C", "CellID(C)")
Fig7_pancreas$num <- str_replace(Fig7_pancreas$num, "CellID_G", "CellID(G)")

Fig7_pancreas$num <- factor(Fig7_pancreas$num, levels = c("Author", "CellID(G)", "CellID(C)"))
SupFig7_pancreas <- ggscatter(
  Fig7_pancreas,
  x = "tSNE_1" ,
  y = "tSNE_2",
  size = 0.01,
  color = "cell_type",
  facet.by = c("data", "num")
) + theme(aspect.ratio = 1) + scale_color_manual(values = PancreasPalette,
                                                 drop = F,
                                                 na.value = "black") + guides(color = guide_legend(override.aes = list(size = 2))) + theme(
                                                   strip.text.x    = element_text(face = "bold", size = 12),
                                                   strip.text.y    = element_text(face = "bold", size = 12),
                                                   legend.title    = element_blank(),
                                                   axis.title      = element_text(size = 10 , face = "bold"),
                                                   axis.text       = element_text(size = 10 , face = "bold"),
                                                   legend.text     = element_text(size = 10 , face = "bold"),
                                                   legend.key.size = unit(0.1, "cm"),
                                                   panel.spacing = unit(1.5, "lines")
                                                 )


Epithelial_celltype <- c("Basal", "Krt4.13", "Secretory", "Goblet", "Ciliated", "Brush.Tuft", "Brush+PNEC","PNEC", "Ionocytes", "unassigned", "other")
Epithelial_color <- c("firebrick1", "purple", "dodgerblue", "cyan4", "darkgreen", "chocolate4", "chocolate2","orange", "cyan","lightgrey", "black")
EpithelialPalette <- setNames(Epithelial_color, Epithelial_celltype)


Fig7_human          <- FetchData(SeuratEpithelial[[2]], vars = c("tSNE_1", "tSNE_2"))
Fig7_human$Author   <- factor(SeuratEpithelial[[2]]$cell_type1, Epithelial_celltype)
Fig7_human$CellID_C <- factor(SeuratEpithelial[[2]]$pred_CellID_C, Epithelial_celltype)
Fig7_human$CellID_G <- factor(SeuratEpithelial[[2]]$pred_CellID_G, Epithelial_celltype)
Fig7_human$data     <- "Plasschaert Human"
Fig7_human          <- Fig7_human %>%  gather("num","cell_type",3:5)

Fig7_montoro          <- FetchData(SeuratEpithelial[[3]], vars = c("tSNE_1", "tSNE_2"))
Fig7_montoro$Author   <- factor(SeuratEpithelial[[3]]$cell_type1, Epithelial_celltype)
Fig7_montoro$CellID_C <- factor(SeuratEpithelial[[3]]$pred_CellID_G, Epithelial_celltype)
Fig7_montoro$CellID_G <- factor(SeuratEpithelial[[3]]$pred_CellID_C, Epithelial_celltype)
Fig7_montoro$data     <- "Montoro"
Fig7_montoro          <- Fig7_montoro %>%  gather("num","cell_type",3:5)

Fig7_Epithelial <- rbind(Fig7_human,Fig7_montoro)
Fig7_Epithelial$cell_type[is.na(Fig7_Epithelial$cell_type)] <- "other"
Fig7_Epithelial$cell_type <- factor(Fig7_Epithelial$cell_type, levels = Epithelial_celltype)
Fig7_Epithelial$num <- str_replace(Fig7_Epithelial$num, "CellID_C", "CellID(C)")
Fig7_Epithelial$num <- str_replace(Fig7_Epithelial$num, "CellID_G", "CellID(G)")
Fig7_Epithelial$num <- factor(Fig7_Epithelial$num, c("Author","CellID(G)","CellID(C)"))
Fig7_Epithelial$importance <- Fig7_Epithelial$cell_type %in% c("other", "unassigned")
Fig7_Epithelial <- arrange(Fig7_Epithelial, importance)
Fig7_Epithelial <- Fig7_Epithelial[rev(seq(nrow(Fig7_Epithelial))),]
SupFig7_Epithelial <- ggscatter(
  Fig7_Epithelial,
  x = "tSNE_1" ,
  y = "tSNE_2",
  size = 0.01,
  color = "cell_type",
  facet.by = c("data", "num")
) + theme(aspect.ratio = 1) + scale_color_manual(values = EpithelialPalette,
                                                 drop = F,
                                                 na.value = "black") + guides(color = guide_legend(override.aes = list(size = 2)))  + theme(
                                                   strip.text      = element_text(face = "bold", size = 10),
                                                   legend.title    = element_blank(),
                                                   axis.title      = element_text(size = 10 , face = "bold"),
                                                   axis.text       = element_text(size = 10 , face = "bold"),
                                                   legend.text     = element_text(size = 10 , face = "bold"),
                                                   legend.key.size = unit(0.1, "cm"),
                                                   panel.spacing = unit(1.5, "lines")
                                                 )


leg_pancreas   <- ggdraw(get_legend(SupFig7_pancreas + theme(legend.position = "right")))
leg_epithelial <- ggdraw(get_legend(SupFig7_Epithelial + theme(legend.position = "right")))
SupFig7_A <- SupFig7_pancreas + theme(legend.position = "none")
SupFig7_B <- SupFig7_Epithelial + theme(legend.position = "none")
design <- c(
  area(t = 1, l = 1, b = 20, r = 32),
  area(t = 21, l = 1, b = 40, r = 32),
  area(t = 1, l = 34, b = 18, r = 40),
  area(t = 21, l = 33, b = 40, r = 40)
)

SupFig7 <- SupFig7_A  + SupFig7_B + leg_pancreas + leg_epithelial + plot_layout(design = design) + geom_text(data = data.frame(label = c("A","B"), x = c(-4.5,-4.5), y= c(2.35,1.12)), aes(x = x, y = y, label = label, fontface = "bold"), size= 6)
SupFig7 <- SupFig7 + theme(plot.margin = unit(c(0,0,0,0),"cm"))
ggsave("../FinalFigure/SupFig7.pdf", SupFig7, width = 183, height = 210, units = "mm")
ggsave("../FinalFigure/SupFig7.png", SupFig7, width = 183, height = 210, units = "mm")
