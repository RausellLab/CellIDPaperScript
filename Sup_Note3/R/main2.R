ImmuneBenchOverall <- ImmuneBenchOverall %>%  group_by(data, metrics) %>%  mutate(colorMacro = max(macro) == macro, colorSD = min(sd) == sd)
ImmuneBenchCell$metrics <- factor(ImmuneBenchCell$metrics, levels = c("Precision", "Recall", "F1"))
ImmuneBenchOverall$metrics <- factor(ImmuneBenchOverall$metrics, levels = c("Precision", "Recall", "F1"))
ImmuneBenchCell$methods <- factor(ImmuneBenchCell$methods, levels = c("CellID", "SCINA", "AUCell"))

Main2BDF <- as_tibble(ImmuneBenchCell %>%  filter(data == "CITE", metrics  == "F1")) %>% arrange(methods, cell_type) %>%  mutate(x = c(seq(0.8,1.2, length.out = 11), 1 + seq(0.8,1.2, length.out = 11), 2 + seq(0.8,1.2, length.out = 11))) %>%  mutate(methods)

citeCol <- palCol[unique(sort(ImmuneBenchCITE$cell_type))]
names(citeCol) <- unique(sort(ImmuneBenchCITE$cell_type))
Main2B <- ggboxplot(
  Main2BDF,
  x = "methods",
  y = "value",
  width = 0.5
) + geom_point(aes(fill = cell_type, x = x, y = value), shape = 21, size = 3) + 
  geom_text(data = ImmuneBenchOverall %>%  filter(data == "CITE", metrics  == "F1"), aes(x = methods, y= 1.12, label = round(macro, digits = 2), color = colorMacro), fontface = "bold", size = 2.5) + 
  geom_text(data = ImmuneBenchOverall %>%  filter(data == "CITE", metrics  == "F1"), aes(x = methods, y= 1.06, label = round(sd, digits = 2), color = colorSD),  fontface = "bold", size = 2.5) + 
  theme(
    aspect.ratio = 1,
    legend.title = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold")
  ) + grids(linetype = "dashed") + scale_fill_manual(values = citeCol) + scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1.2)) + guides(color =F) + scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"))

Main2D <- SupFig3C + theme(legend.position = "none")

Main2A <- SupFig3B + theme(axis.text = element_blank(),
                           axis.ticks = element_blank(), 
                           axis.title.x = element_text(size = 8),
                           axis.title.y = element_text(size = 8))

labelDF <- data.frame(UMAP_1 = c(-5, -2.5, 2, -12, -8, -7, -6, 3, 7, 5, -11),
                      UMAP_2 = c(-20.4744696, 0, -18.1897492, -4, 8, -5, -8, 7, 0.5, -13, -19),
                      SupFig3BCol = c("Eryth","Mk","CD34","DC","CD14 Mono","CD16 Mono","B","CD4 T","CD8 T","NK","pDC"))

segmentDF <- labelDF %>% rename(UMAP_1start=UMAP_1, UMAP_2start=UMAP_2)
segmentDF <- inner_join(Main2A$layers[[2]]$data, segmentDF)

Main2A$layers[[1]]$show.legend <- F
Main2A$layers[[2]] <- NULL
Main2A <- Main2A + geom_segment(data = segmentDF, aes(x = UMAP_1,y= UMAP_2, xend = UMAP_1start, yend = UMAP_2start)) + geom_label(data = labelDF, aes(x = UMAP_1, y = UMAP_2, label = SupFig3BCol), size =2, fontface ="bold") 




set.seed(12345)
Main2C <- ggscatter(
  Fig6B_df[sample(seq(nrow(Fig6B_df))),],
  x = "UMAP_1",
  y = "UMAP_2",
  size = 0.2,
  color = "pred", 
  alpha = 0.9,
) + coord_cartesian(xlim =  c(-2.0, -1.1), ylim = c(-22,-17)) + scale_color_gradientn(
  values = c(1, .2, .1, 0),
  limits = c(0, 10),
  colours = c("#D73027", "#FEDF8F", "#4E7DB8", "#4575B4"),
  breaks = seq(0, 10, 2),
  na.value = 'lightgrey',
  oob = scales::squish,
  guide = guide_colorbar(label = TRUE,
                         frame.colour = "black")
) + facet_wrap(~ ident, ncol = 6) + theme(
  legend.position = "none",
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  aspect.ratio = 1,
  strip.text = element_text(size = 8, face =  "bold"),
  axis.title.x = element_text(size = 8, face  = "bold"),
  axis.title.y = element_blank()
)


Leg <- DimPlot(cite_seurat, group.by = "SupFig3BCol", label = T, label.size = 3)

Main2Leg1 <- ggarrange(get_legend(ggplot(data.frame(ident = factor(c(immune_levels,"unassigned"),c(immune_levels,"unassigned")))) + geom_point(aes(x = ident, y = ident, fill = ident), shape =21, color = "black", size = 3) + theme_classic() +  scale_fill_manual(values = c(citeCol,  unassigned = "lightgrey"))+ theme(legend.text =element_text(size =8, face ="bold"), legend.key.size = unit(0,"mm"), title = element_blank())))
Main2Leg2 <- ggarrange(get_legend(ggscatter(
  Fig6B_df[sample(seq(nrow(Fig6B_df))),],
  x = "UMAP_1",
  y = "UMAP_2",
  size = 0.2,
  color = "pred", 
  alpha = 0.9,
) + coord_cartesian(xlim =  c(-2.0, -1.1), ylim = c(-22,-17)) + scale_color_gradientn(
  values = c(1, .2, .1, 0),
  limits = c(0, 10),
  colours = c("#D73027", "#FEDF8F", "#4E7DB8", "#4575B4"),
  breaks = seq(0, 10, 2),
  na.value = 'lightgrey',
  oob = scales::squish,
  guide = guide_colorbar(direction = "vertical",label = TRUE,
                         frame.colour = "black", title = "-Log10 pval", title.theme = element_text(size = 8, face = "bold"), title.position = "top", title.hjust = 0, label.theme = element_text(size = 8, face = "bold"))
)))

Main2D <- SupFig3C + theme(legend.position = "none")
Main2D$layers[[5]]$aes_params$size <- 2
Main2D <- Main2D + theme(axis.text.y = element_text(size = 6))


library(patchwork)
Void <- ggplot() + geom_blank() + theme_void()
Main2LabelDF_A <- ggplot(tibble(label = "A", x = -1, y = 0, )) + geom_text(aes(x = x , y = y, label = label,fontface = "bold"), size = 5) + theme_void()
Main2LabelDF_B <- ggplot(tibble(label = "B", x = -1, y = 0, )) + geom_text(aes(x = x , y = y, label = label,fontface = "bold"), size = 5) + theme_void()
Main2LabelDF_C <- ggplot(tibble(label = "C", x = -1, y = 0, )) + geom_text(aes(x = x , y = y, label = label,fontface = "bold"), size = 5) + theme_void()
Main2LabelDF_D <- ggplot(tibble(label = "D", x = -1, y = 0, )) + geom_text(aes(x = x , y = y, label = label,fontface = "bold"), size = 5) + theme_void()
design <-
  c(Main2A = area(2, 2, 19, 9),
    Main2B = area(2, 11, 19, 18),
    Main2D = area(27, 2, 42, 19),
    Main2C = area(20, 2, 26, 19),
    Main2Leg1 = area(9, 19, 14, 21),
    Main2Leg2 = area(24, 20, 30, 21),
    Main2LabelDF_A = area(1, 1, 2, 2),
    Main2LabelDF_B = area(1, 10, 2, 11),
    Main2LabelDF_C = area(19, 1, 20, 2),
    Main2LabelDF_D = area(26, 1, 27, 2)
  )

Main2 <- Main2A + Main2B + Main2D + Main2C + Main2Leg1 + Main2Leg2 + Main2LabelDF_A + Main2LabelDF_B + Main2LabelDF_C + Main2LabelDF_D + plot_layout(design = design)
ggsave("../FinalFigure/Main/Main2.pdf", Main2, units = "mm", dpi = 320, width = 210, height = 180)
ggsave("../FinalFigure/Main/Main2.svg", Main2, units = "mm", dpi = 320, width = 210, height = 180)
ggsave("../FinalFigure/Main/Main2.png", Main2, units = "mm", dpi = 320, width = 210, height = 180)
