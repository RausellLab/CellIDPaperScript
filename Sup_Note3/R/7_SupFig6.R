

theme_supFig6 <- theme(
  title = element_text(size = 10, face = "bold.italic"),
  legend.text = element_text(size = 12, face = "bold"),
  legend.title = element_text(size = 12, face = "bold"),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  strip.text = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(size = 9, face = "bold"),
  axis.text.y = element_text(size = 9, face = "bold"),
  legend.title.align = 0.5,
  legend.position = "right",
  plot.margin = unit(x = c(0.5, 0, 0.3, 1), units = "cm")
)

show_cell <- c("CD34", "Eryth")
cite_seurat$Fig6A_group <-
  ifelse(cite_seurat$cell_type1 %in% show_cell,
         as.vector(cite_seurat$cell_type1),
         "other")
SupFig6A <-
  DimPlot(cite_seurat, group.by = "Fig6A_group")  + scale_color_manual(values = c(citeCol[show_cell], other = "grey")) + theme_supFig6 + theme(legend.position = "top")
SupFig6A$data$zoom <- F
alt_df <- SupFig6A$data
alt_df$zoom <- T
SupFig6A <-
  SupFig6A + geom_point(
    data = alt_df,
    aes(x = UMAP_1, y = UMAP_2, fill = Fig6A_group),
    shape = 21,
    size = 3
  ) + scale_fill_manual(values = c(citeCol[show_cell], other = "grey")) + guides(fill =
                                                                                   F)
SupFig6A <-
  SupFig6A + facet_zoom(
    xlim =  c(-2.1, -1),
    ylim = c(-22, -17),
    horizontal = F,
    zoom.size = 1,
    zoom.data = zoom
  )

SupFig6A <-
  ggarrange(
    SupFig6A,
    labels = " ",
    vjust = 0,
    font.label = list(size = 20)
  ) + theme(aspect.ratio = 2)

Fig6B_df <-
  FetchData(
    cite_seurat,
    vars = c(
      "UMAP_1",
      "UMAP_2",
      "HSC",
      "MPP",
      "CMP",
      "MEP",
      "GMP",
      "Erythrocytes"
    )
  ) %>%  gather("ident", "pred",-seq(2)) %>%  mutate(ident = factor(
    ident,
    levels = c("HSC", "MPP", "CMP",  "GMP", "MEP", "Erythrocytes")
  )) %>% group_by(ident) %>%  arrange(pred)


SupFig6B <- ggscatter(
  Fig6B_df,
  x = "UMAP_1",
  y = "UMAP_2",
  fill = "pred",
  shape = 21,
  size = 1,
  colour = "black"
) + coord_cartesian(xlim =  c(-2.1,-1), ylim = c(-22, -17)) + scale_fill_gradientn(
  values = c(1, .2, .1, 0),
  limits = c(0, 10),
  colours = c("#D73027",  "#FEDF8F", "#4E7DB8", "#4575B4"),
  breaks = seq(0, 10, 2),
  na.value = 'lightgrey',
  oob = scales::squish,
  guide = guide_colorbar(label = TRUE,
                         frame.colour = "black")) + theme_supFig6  + facet_wrap( ~ ident, ncol = 2) + theme(aspect.ratio = 1, legend.position = "top")
  SupFig6B <- ggarrange(SupFig6B,  font.label = list(size = 20))
  

SupFig6Design <- c(area(
  t = 1,
  l = 1,
  b = 19,
  r = 9
),
area(
  t = 1,
  l = 10,
  b = 18,
  r = 20
))

SupFig6 <-
  SupFig6A + SupFig6B + plot_layout(design = SupFig6Design) + geom_text(mapping = aes(
    x = -0.70,
    y = 0.95,
    label = "A",
    fontface = "bold"
  ),
  size = 7) + geom_text(mapping = aes(
    x = 0.1,
    y = 0.95,
    label = "B",
    fontface = "bold"
  ),
  size = 7)
ggsave(
  SupFig6,
  filename = "../FinalFigure/SupFig6.pdf",
  width = 180,
  height = 180 ,
  units = "mm"
)
ggsave(
  SupFig6,
  filename = "../FinalFigure/SupFig6.png",
  width = 180,
  height = 180 ,
  dpi = 320,
  units = "mm"
)
