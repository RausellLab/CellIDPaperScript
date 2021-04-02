SeuratPancreas[[2]]$predicted <- ifelse(SeuratPancreas[[2]]$pred_CellID_G == "schwann", yes = "schwann", no = "other")
SeuratPancreas[[3]]$predicted <- ifelse(SeuratPancreas[[3]]$pred_CellID_G == "schwann", yes = "schwann", no = "other")

schwann_genes <- c("SOX10", "S100B", "CRYAB", "NGFR", "CDH19", "PMP22")

schwann_muraro_df <- Seurat::FetchData(SeuratPancreas[[2]], vars = schwann_genes) %>%  rownames_to_column(var = "cell") %>%  gather("genes", "expression_level",-1) %>% inner_join(rownames_to_column(SeuratPancreas[[2]]@meta.data,"cell")[,c("cell","predicted")])
schwann_segers_df <- Seurat::FetchData(SeuratPancreas[[3]], vars = schwann_genes) %>%  rownames_to_column(var = "cell") %>%  gather("genes", "expression_level",-1) %>% inner_join(rownames_to_column(SeuratPancreas[[3]]@meta.data,"cell")[,c("cell","predicted")])


gg_schwann_muraro <-
    ggboxplot(
        schwann_muraro_df,
        x = "predicted",
        y = "expression_level",
        color = "predicted",
        add = "jitter",
        facet.by = "genes"
    ) +
    stat_compare_means(
        label.y = 5.5,
        comparisons = list(c("other", "schwann")),
        paired = F,
        label = "p.format"
    ) + ylim(c(-0.1, 6)) + ylab("expression level") +
    theme(
        strip.text.x = element_text(face = "bold"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 12 , face = "bold"),
        axis.title.y =  element_text(size = 12 , face = "bold"),
        legend.text = element_text(size = 12 , face = "bold"),
        legend.key.size = unit(1.2, "cm")
    ) + scale_color_manual(values = c(other = "grey", schwann = "gold"))

gg_schwann_muraro <-
    ggboxplot(
        schwann_muraro_df,
        x = "predicted",
        y = "expression_level",
        color = "predicted",
        add = "jitter",
        facet.by = "genes"
    ) +
    geom_signif(
        comparisons = list(c("other", "schwann")),
        map_signif_level = function(p) sprintf("p = %.2g", p),y_position = 5.5
    ) + ylim(c(-0.1, 6)) + ylab("expression level")  +
    theme(
        strip.text.x = element_text(face = "bold"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 12 , face = "bold"),
        axis.title.y =  element_text(size = 12 , face = "bold"),
        legend.text = element_text(size = 12 , face = "bold"),
        legend.key.size = unit(1.2, "cm")
    ) + scale_color_manual(values = c(other = "grey", schwann = "gold"))

gg_schwann_segers <-
    ggboxplot(
        schwann_segers_df,
        x = "predicted",
        y = "expression_level",
        color = "predicted",
        add = "jitter", 
        
        facet.by = "genes"
    ) +
    geom_signif(
        comparisons = list(c("other", "schwann")),
        map_signif_level = function(p) sprintf("p = %.2g", p),y_position = 5.5
    ) + ylim(c(-0.1, 6)) + ylab("expression level")  +
    theme(
        strip.text.x = element_text(face = "bold"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 12 , face = "bold"),
        axis.title.y =  element_text(size = 12 , face = "bold"),
        legend.text = element_text(size = 12 , face = "bold"),
        legend.key.size = unit(1.2, "cm")
    ) + scale_color_manual(values = c(other = "grey", schwann = "gold"))

gg_schwann_segers$layers[[1]]$data <- gg_schwann_segers$layers[[1]]$data %>%  filter(predicted == "other")

SupFig9 <- ggarrange(gg_schwann_muraro, gg_schwann_segers, nrow = 2, legend = "top",common.legend = T, labels = c("A", "B"), font.label = list(size = 16), vjust = 0)
ggsave(plot = SupFig9, filename = "../FinalFigure/SupFig9.pdf", width = 18, height = 20, units = "cm")
ggsave(plot = SupFig9, filename = "../FinalFigure/SupFig9.png", width = 18, height = 20, units = "cm", dpi = 600)
