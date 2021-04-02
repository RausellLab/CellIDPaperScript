source("R/utilitary_function.R")
library(foreach)

PancreasBenchOverall <- read_rds("data/PancreasBenchOverall.rds")
EpithelialBenchOverall <- read_rds("data/EpithelialBenchOverall.rds")

PancreasBenchCell <- read_rds("data/PancreasBenchCell.rds")
EpithelialBenchCell <- read_rds("data/EpithelialBenchCell.rds")

overall <- rbindlist(list(PancreasBenchOverall, EpithelialBenchOverall))
overall$value[is.na(overall$value)] <- 0
rare <- rbindlist(list(PancreasBenchCell,EpithelialBenchCell))


# Legend ------------------------------------------------------------------

rarePanCol <- c("lightcoral", "tomato1", "firebrick1", "red4") %>% set_names("epsilon", "macrophage", 'mast', "endothelial")
rareAirCol <- c("cyan", "deepskyblue", "deepskyblue4", "darkblue") %>% set_names("Ionocytes", "Brush.Tuft", 'Brush+PNEC', "PNEC")
rareCol <- c(rarePanCol, rareAirCol) 


ggLegShape <- 
    ggdraw(get_legend(
        ggplot(data.frame(
            x = seq(4),
            y = seq(4),
            data = c("22", "23", "24", "25")
        )) + geom_point(aes(
            x = x, y = y, shape = data
        ), size = 2.5) + scale_shape_manual(
            values = c(22, 23, 24, 25),
            name = "data",
            labels = c("Muraro", "Segerstolpe", "Plasschaert Human", "Montoro")
        )  + guides(shape = guide_legend(nrow = 1, override.aes = list(color = c("red", "red", "blue", "blue"), fill ="grey"))) + theme_bw() + theme(
            legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_text(size = 5, face = "bold", color = "black")
        ) 
    )) 



ggLegCol <-
    ggdraw(get_legend(
        ggplot(data.frame(
            x = seq(8),
            y = seq(8),
            cell_type = factor(
                c(
                    "PNEC",
                    "Brush.Tuft",
                    "Brush+PNEC",
                    "Ionocytes",
                    "epsilon",
                    "endothelial",
                    "mast",
                    "macrophage"
                ),
                levels = c(
                    "epsilon",
                    "macrophage",
                    "mast",
                    "endothelial",
                    "Ionocytes",
                    "Brush.Tuft",
                    "Brush+PNEC",
                    "PNEC"
                )
            )
        )) + geom_point(
            aes(x = x, y = y, fill = cell_type),
            shape = 21,
            size = 2.5
        ) + scale_fill_manual(values = rareCol, name = "cell type") + theme_bw() +
            theme(
                legend.position = "top",
                legend.title = element_blank(),
                legend.text = element_text(size = 5, face = "bold", color = "black"),legend.key.height = unit(0, "mm")
                
            )
    ))

# Global ------------------------------------------------------------------
overall$shape <-
    c(
        Muraro = 22,
        Segerstolpe = 23,
        `Plasschaert Human` = 24,
        Montoro = 25
    )[overall$data]
Globalrank <-
    foreach(
        x = unique(overall$metrics),
        .final = function(x)
            setNames(x, unique(overall$metrics))
    ) %do% {
        A <-
            overall %>% filter(metrics == x) %>% group_by(methods) %>% summarise(mean = mean(value)) %>% arrange(-mean) %>%  use_series(methods) %>% str_subset('CellID', negate = T)
        B <- c("CellID(G)", "CellID(C)")
        order.methods <- c(B, A)
        overall$methods <- factor(overall$methods, order.methods)
        mean_sd <-
            overall %>% filter(metrics == x)  %>% group_by(methods) %>% summarise(mean = round(mean(value), digits = 2),
                                                                                  sd = round(sd(value), digits = 2)) %>% mutate(minsd = sd == min(sd),
                                                                                                                                maxmean = mean == max(mean)) %>% arrange(-mean)
        p <-
            overall %>% filter(metrics == x)  %>%  arrange(methods) %>%  mutate(x = as.vector(sapply(seq(0,11,length.out = 12), function(x) x + c(seq(0.70, 0.90, length.out = 2), seq(1.1, 1.30, length.out = 2))))) %>% ggboxplot(
                x = "methods",
                y = "value",
                outlier.shape=  NA,
                ggtheme = theme_bw() 
            ) + scale_shape_identity() + scale_color_manual(values = c(Muraro = "red", Segerstolpe = "red", Montoro = "blue",  `Plasschaert Human` = "blue", `FALSE` = '#F8766D', `TRUE` = '#00BFC4'))  + theme(
                axis.title = element_blank(),
                axis.text.x =  element_text(
                    angle = 45,
                    hjust = 1, 
                    color = "black",
                    size = 7, 
                    face = "bold"
                ),
                axis.text.y = element_text(size = 5, face = "bold", color = "black")
            ) + geom_vline(xintercept = 2.5, size = 0.7)  + geom_point(aes(x=x, y= value, shape = shape, color = data), fill = "lightgrey") + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1)) + guides(color = FALSE, fill = FALSE) 
    }

SupFig8All <- ggarrange(
    plotlist = list(Globalrank$Precision, Globalrank$Recall),
    ncol = 1,
    nrow = 2,
    labels = "AUTO",
    vjust = 0
)



# Rare --------------------------------------------------------------------
rare$methods <- str_replace(rare$methods, "CellID_G", "CellID(G)")
rare$methods <- str_replace(rare$methods, "CellID_C", "CellID(C)")
rare$shape <- c(Muraro = 22, Segerstolpe=23, `Plasschaert Human` = 24, Montoro = 25)[rare$data]

Rarerank <-
    foreach(
        x = c("Precision", "Recall", "F1"),
        .final = function(x)
            setNames(x, c("Precision", "Recall", "F1"))
    ) %do% {
        rare_df <-
            rare %>% filter(
                cell_type %in% c(
                    "Ionocytes",
                    "PNEC",
                    "Brush.Tuft",
                    "Brush+PNEC",
                    "epsilon",
                    "endothelial",
                    "mast",
                    "macrophage"
                ),
                metrics == x
            )
        A <-
            rare_df %>% group_by(methods) %>% summarise(mean = mean(value)) %>% arrange(-mean) %>%  use_series(methods) %>% str_subset('CellID', negate = T)
        B <- c("CellID(G)", "CellID(C)")
        order.methods <- c(B, A)
        rare_df$methods <- factor(rare_df$methods, order.methods)
        rare_df_mean_sd <-
            rare_df  %>% group_by(methods) %>% summarise(mean = round(mean(value), digits = 2),
                                                         sd = round(sd(value), digits = 2)) %>% mutate(minsd = sd == min(sd),
                                                                                                       maxmean = mean == max(mean)) %>% arrange(-mean)
        
        p <-
            rare_df %>% ggboxplot(
                x = "methods",
                y = "value",
                add = "jitter", outlier.shape=  NA,
                add.params = list(
                    size = 1,
                    fill = "cell_type",
                    shape = "shape"
                ),
                ggtheme = theme_bw()
            ) + theme(
                axis.title = element_blank(),
                axis.text.x = element_text(
                    angle = 45,
                    hjust = 1, 
                    color = "black",
                    size = 7, 
                    face = "bold"
                ),
                axis.text.y = element_text(size = 7, face = "bold", color = "black")
            ) + geom_vline(xintercept = 2.5, size = 0.7) + guides(color = FALSE, fill = FALSE) +
            scale_shape_identity()
        p$layers[[2]]$position <- position_dodge(0.6)
        p + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1)) + scale_fill_manual(values = rareCol)
    }
SupFig8Rare <- ggarrange(plotlist = list(Rarerank$Precision,Rarerank$Recall), nrow = 2, labels = c("C", "D"), vjust = 0)



# Integrate ---------------------------------------------------------------


design <- c(area(1,1,1,5),
            area(1,6,1,10),
            area(2,1,10,5),
            area(2,6,10,10)
            )
SupFig8 <- ggLegShape + ggLegCol + SupFig8All + SupFig8Rare + plot_layout(design = design)
ggsave(filename = "../FinalFigure/SupFig8.pdf", SupFig8, units = "mm", width = 180, height = 200, dpi = 320)
ggsave(filename = "../FinalFigure/SupFig8.png", SupFig8, units = "mm", width = 180, height = 200, dpi = 320)
