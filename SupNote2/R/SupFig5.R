# SupFig5A ----------------------------------------------------------------

ImmuneBenchOverall <- ImmuneBenchOverall %>%  group_by(data, metrics) %>%  mutate(colorMacro = max(macro) == macro, colorSD = min(sd) == sd)
ImmuneBenchCell$metrics <- factor(ImmuneBenchCell$metrics, levels = c("Precision", "Recall", "F1"))
ImmuneBenchOverall$metrics <- factor(ImmuneBenchOverall$metrics, levels = c("Precision", "Recall", "F1"))

SupFig5A <- ggboxplot(
  ImmuneBenchCell %>%  filter(data == "CITE"),
  x = "methods",
  y = "value",
  width = 0.5, 
  facet.by = "metrics",
  outlier.shape = NA,
) + geom_jitter(aes(fill = cell_type), shape = 21, size = 4, width = 0.2) + 
  geom_text(data = ImmuneBenchOverall %>%  filter(data == "CITE"), aes(x = methods, y= 1.12, label = round(macro, digits = 2), color = colorMacro), fontface = "bold", size = 3) + 
  geom_text(data = ImmuneBenchOverall %>%  filter(data == "CITE"), aes(x = methods, y= 1.06, label = round(sd, digits = 2), color = colorSD),  fontface = "bold", size = 3) + 
  theme(
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    strip.text.x = element_text(size = 10, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold")
  ) + grids(linetype = "dashed") + scale_fill_manual(values = citeCol) + scale_y_continuous(breaks = seq(0, 1, by = 0.25)) + guides(color =F)+ scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"))
SupFig5A$layers[[2]]$position <- position_dodge(0.7)

SupFig5B <- ggboxplot(
  ImmuneBenchCell %>%  filter(data == "REAP"),
  x = "methods",
  y = "value",
  width = 0.5, 
  facet.by = "metrics",
  outlier.shape = NA
) + geom_jitter(aes(fill = cell_type), shape = 21, size = 4, width = 0.2) + 
  geom_text(data = ImmuneBenchOverall %>%  filter(data == "REAP"), aes(x = methods, y= 1.12, label = round(macro, digits = 2), color = colorMacro), fontface = "bold", size = 3) + 
  geom_text(data = ImmuneBenchOverall %>%  filter(data == "REAP"), aes(x = methods, y= 1.06, label = round(sd, digits = 2), color = colorSD),  fontface = "bold", size = 3) + 
  theme(
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    strip.text.x = element_text(size = 10, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold")
  ) + grids(linetype = "dashed") + scale_fill_manual(values = citeCol) + scale_y_continuous(breaks = seq(0, 1, by = 0.25)) + guides(color =F) + scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"))
SupFig5B$layers[[2]]$position <- position_dodge(0.7)


supFig5 <- ggarrange(SupFig5A, SupFig5B, labels = c("A", "B"), ncol = 1, font.label = list(size = 16), vjust = 0, common.legend = T)
ggsave("../FinalFigure/supFig5.pdf", plot = supFig5, height = 180, width = 180, units = "mm")
ggsave("../FinalFigure/supFig5.png", plot = supFig5, height = 180, width = 180, dpi = 320, units = "mm")
