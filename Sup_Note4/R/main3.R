source("R/utilitary_function.R")
PancreasBenchOverall <- read_rds("data/PancreasBenchOverall.rds")
EpithelialBenchOverall <- read_rds("data/EpithelialBenchOverall.rds")
PancreasBenchCell <- read_rds("data/PancreasBenchCell.rds")
EpithelialBenchCell <- read_rds("data/EpithelialBenchCell.rds")

overall <- rbindlist(list(PancreasBenchOverall, EpithelialBenchOverall))
overall$value[is.na(overall$value)] <- 0
rare <- rbindlist(list(PancreasBenchCell,EpithelialBenchCell))

rarePanCol <- c("lightcoral", "tomato1", "firebrick1", "red4") %>% set_names("epsilon", "macrophage", 'mast', "endothelial")
rareAirCol <- c("cyan", "deepskyblue", "deepskyblue4", "darkblue") %>% set_names("Ionocytes", "Brush.Tuft", 'Brush+PNEC', "PNEC")
rareCol <- c(rarePanCol, rareAirCol) 


# -------------------------------------------------------------------------

Pancreas_celltype <- c("alpha","beta", "delta", "gamma", "epsilon", "acinar", "ductal",  "PSC", "endothelial", "macrophage", "mast", "t_cell", "schwann", "unassigned", "other")
Pancreas_color <- c("firebrick1", "coral", "dodgerblue", "cyan", "cyan4","chocolate", "chocolate4",  "olivedrab","lightgreen","purple", "violet", "pink", "gold", "lightgrey", "black")
PancreasPalette <- setNames(Pancreas_color, Pancreas_celltype)

Epithelial_celltype <- c("Basal", "Krt4.13", "Secretory", "Goblet", "Ciliated", "Brush.Tuft", "Brush+PNEC","PNEC", "Ionocytes", "unassigned", "other")
Epithelial_color <- c("firebrick1", "purple", "dodgerblue", "cyan4", "darkgreen", "chocolate4", "chocolate2","orange", "cyan","lightgrey", "black")
EpithelialPalette <- setNames(Epithelial_color, Epithelial_celltype)


ggLegShape <- 
  ggdraw(get_legend(
    ggplot(data.frame(
      x = seq(4),
      y = seq(4),
      data = c("22", "23", "24", "25")
    )) + geom_point(aes(
      x = x, y = y, shape = data
    ), size = 2) + scale_shape_manual(
      values = c(22, 23, 24, 25),
      name = "data",
      labels = c("Muraro", "Segerstolpe", "Plasschaert Human", "Montoro")
    )  + guides(shape = guide_legend(nrow = 1, override.aes = list(color = c("red", "red", "blue", "blue"), fill ="grey"))) + theme_bw() + theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 7, face = "bold", color = "black")
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
      size = 2
    ) + scale_fill_manual(values = rareCol, name = "cell type") + theme_bw() +
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 7, face = "bold", color = "black"),legend.key.height = unit(0, "mm")
        
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

Main3B1 <-{A <-
      overall %>% filter(metrics == "F1") %>% group_by(methods) %>% summarise(mean = mean(value)) %>% arrange(-mean) %>%  use_series(methods) %>% str_subset('CellID', negate = T)
    B <- c("CellID(G)", "CellID(C)")
    order.methods <- c(B, A)
    overall$methods <- factor(overall$methods, order.methods)
    p <-
      overall %>% filter(metrics == x)  %>%  arrange(methods) %>%  mutate(x = as.vector(sapply(seq(0,11,length.out = 12), function(x) x + c(seq(0.70, 0.90, length.out = 2), seq(1.1, 1.30, length.out = 2))))) %>% ggboxplot(
        x = "methods",
        y = "value",
        outlier.shape = NA,
        ggtheme = theme_bw() 
      ) + scale_shape_identity() + scale_color_manual(values = c(Muraro = "red", Segerstolpe = "red", Montoro = "blue",  `Plasschaert Human` = "blue", `FALSE` = '#F8766D', `TRUE` = '#00BFC4'))  + theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7, face = "bold", color = "black")
      ) + geom_vline(xintercept = 2.5, size = 0.7)  + geom_point(aes(x=x, y= value, shape = shape, color = data), fill = "lightgrey") + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1), limits = c(0,1)) + guides(color = FALSE, fill = FALSE) 
}


# Rare --------------------------------------------------------------------

# Rare --------------------------------------------------------------------
rare$methods <- str_replace(rare$methods, "CellID_G", "CellID(G)")
rare$methods <- str_replace(rare$methods, "CellID_C", "CellID(C)")
rare$shape <- c(Muraro = 22, Segerstolpe = 23, `Plasschaert Human` = 24, Montoro = 25)[rare$data]

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
    p <-
      rare_df %>% ggboxplot(
        x = "methods",
        y = "value",
        add = "jitter",
        add.params = list(
          size = 2,
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
      ) + geom_vline(xintercept = 2.5, size = 0.7) + geom_text(data = rare_df_mean_sd, aes(
        x = methods,
        y = 1.14,
        label = mean,
        colour = maxmean,
        fontface = "bold"
      ), size = 2.5) + guides(color = FALSE, fill = FALSE) +
      scale_shape_identity()
    p$layers[[2]]$position <- position_dodge(0.6)
    p + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1)) + scale_fill_manual(values = rareCol)
  }

main_rare_df <- rare_df %>%  filter(metrics == "F1") %>%  mutate(cell_type = factor(cell_type, c("epsilon","endothelial","mast","macrophage", "Ionocytes","PNEC","Brush.Tuft","Brush+PNEC"))) %>%  arrange(methods, cell_type) %>%  mutate(x = as.vector(sapply(seq(0,11,length.out = 12), function(x) x + c(seq(0.70, 0.90, length.out = 6), seq(1.1, 1.30, length.out = 5)))))  
Main3B2 <- ggboxplot(data = main_rare_df, x = "methods", y = "value", ggtheme = theme_bw(), outlier.shape = NA) + geom_point(aes(x = x, y = value, shape = shape, fill = cell_type)) + scale_shape_identity() + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1), limits = c(0,1)) + scale_fill_manual(values = rareCol)  + scale_color_manual(values = c(Muraro = "red", Segerstolpe = "red", Montoro = "blue",  `Plasschaert Human` = "blue", `FALSE` = '#F8766D', `TRUE` = '#00BFC4')) + theme(
  axis.title = element_blank(),
  axis.text.x = element_text(
    angle = 45,
    hjust = 1, 
    color = "black",
    size = 8, 
    face = "bold"
  ),
  axis.text.y = element_text(size = 8, face = "bold", color = "black")
) + geom_vline(xintercept = 2.5, size = 0.7)+ guides(color = FALSE, fill = FALSE) + theme(
  axis.title = element_blank(),
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    size = 7,
    face = "bold",
    color = "black"
  ),
  axis.text.y = element_text(size = 7, face = "bold", color = "black")
)

# Intestinal --------------------------------------------------------------
SeuratIntestinal <- read_rds("../SupNote4/data/SeuratIntestinal.rds")
IntestinalBenchOverall <- read_rds("../SupNote4/data/IntestinalBenchOverall.rds")
SeuratIntestinal$Main3 <- as.vector(SeuratIntestinal$pred_CellID_G__pla)
SeuratIntestinal$Main3 <- ifelse(SeuratIntestinal$Main3 == "Secretory", yes = "Goblet", no = SeuratIntestinal$Main3)
SeuratIntestinal$Main3 <- ifelse(SeuratIntestinal$Main3 == "PNEC", yes = "Endocrine", no = SeuratIntestinal$Main3)
SeuratIntestinal$Main3 <- ifelse(SeuratIntestinal$Main3 %in% c("Ionocytes", "Basal", "Ciliated", "Krt4.13"), yes = "other", no = SeuratIntestinal$Main3)
Main3Intes1 <- DimPlot(SeuratIntestinal, group.by = 'Main3', reduction = "tsne", order = c("Brush.Tuft", "Endocrine", "Goblet", "other")) + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size =8 , face = "bold")) + scale_color_manual(values = c(other = "black", unassigned = "lightgrey", Endocrine = unname(rareCol["PNEC"]), Brush.Tuft = unname(rareCol["Brush.Tuft"]),  Goblet = "cyan4")) + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size =8 , face = "bold"),  legend.position = "none", aspect.ratio = 1)

Main3Intes1Leg <- ggdraw(get_legend(ggplot(data = data.frame(name = factor(c("Goblet", "unassigned","Endocrine" , "other", "Brush.Tuft"), levels = c("Brush.Tuft", "unassigned","Endocrine" , "other","Goblet" )), x= seq(5), y = seq(5))) + geom_point(aes(x=x, y=y, fill = name), shape=21) + theme_classic()+ guides(fill = guide_legend(nrow = 2, override.aes = list(size = 2))) + scale_fill_manual(values = c(other = "black", unassigned = "lightgrey", Endocrine = unname(rareCol["PNEC"]), Brush.Tuft = unname(rareCol["Brush.Tuft"]),  Goblet = "cyan4")) + theme(legend.title = element_blank())+ theme(
  legend.position = "top",
  legend.title = element_blank(),
  legend.text = element_text(size = 7, face = "bold", color = "black"),
  legend.key.size = unit(0, "mm")
) ))


Main3Intes2 <- IntestinalBenchOverall %>%  filter(metrics == "F1", data == "Plasschaert Mouse") %>%  ggbarplot(x = "methods", y = "value", fill = "lightgrey", label = T, lab.size = 2, lab.nb.digits = 2) + theme_bw() + theme(
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
) + geom_vline(xintercept = 2.5) 

# Olfactory ---------------------------------------------------------------
SeuratWu <- read_rds("../SupNote4/data/SeuratWu.rds")
DF <- FetchData(SeuratWu, vars = c("Airway.Brush.Tuft", "Intestinal.Brush.Tuft"))
SeuratWu$SCC <- ifelse((DF$Airway.Brush.Tuft >20) & (DF$Intestinal.Brush.Tuft >20), "putative SCC", "other")
Main3Olfa1 <- DimPlot(SeuratWu, group.by = "SCC", reduction = "tsne") + scale_color_manual(values = c(other= 'grey', `putative SCC` = unname(rareCol["Brush.Tuft"]))) + annotate("rect", xmin = 19, xmax = 25, ymin = -41, ymax = -48, alpha = 0.2, color = "Gold", fill = "white") + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size =8 , face = "bold"), legend.position = "none",  aspect.ratio = 1)

LegOlfa1 <- ggdraw(get_legend(ggplot(data = data.frame(name = factor(c("putative SCC", "unassigned"), levels = c("putative SCC", "unassigned")), x= seq(2), y = seq(2))) + geom_point(aes(x=x, y=y, fill = name), shape=21) + theme_classic()+ guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2))) + scale_fill_manual(values = c(unassigned = 'grey', `putative SCC` = unname(rareCol["Brush.Tuft"]))) + theme(legend.title = element_blank())+ theme(
  legend.position = "top",
  legend.title = element_blank(),
  legend.text = element_text(size = 7, face = "bold", color = "black")
) ))

DF2 <- FetchData(SeuratWu, vars = c("Airway.Brush.Tuft", "Intestinal.Brush.Tuft"))


Main3Olfa2DF <- FetchData(SeuratWu, vars = c("leukotriene biosynthetic process (GO:0019370)---GO Biological Process", 
                                             "interleukin-17-mediated signaling pathway (GO:0097400)---GO Biological Process",
                             "Eicosanoid Synthesis WP318---WikiPathways",
                             "Cholinergic synapse---KEGG",
                             "Transmission across Chemical Synapses Homo sapiens R-HSA-112315---Reactome"))[(DF$Airway.Brush.Tuft >20) & (DF$Intestinal.Brush.Tuft >20),] %>%  set_colnames(c("GO:0019370", "GO:0097400", "WP318", "KEGG-hsa04725", "R-HSA-112315")) %>%  gather("pathway","value")


Main3Olfa2 <- ggviolin(Main3Olfa2DF, x = "pathway", y = "value", size= 0.5, add = "jitter", fill = "lightgoldenrod1",add.params = list(size=0.1)) + ylab("-Log10 pvalue") + xlab("") + theme(aspect.ratio = 1, axis.text.x = element_text(size = 7 , face = "bold", angle = 45, hjust = 1), axis.text.y = element_text(size = 7 , face = "bold"), axis.title.y = element_text(size = 7 , face = "bold")) + scale_y_continuous(limits = c(0,8)) 


TitleMain2Intes <- ggplot() + theme_void() + geom_text(aes(x = 0,  y = 0, label ="F1 Score", fontface = "bold"),size = 3)
TitleMain2Olfa <- ggplot() + theme_void() + geom_text(aes(x = 0,  y = 0, label ="Pathway Enrichment", fontface = "bold"),size = 3)
AnnotA <- ggplot() + theme_void() + geom_text(aes(x = 0,  y = 0, label ="A", fontface = "bold"),size = 6)
AnnotB <- ggplot() + theme_void() + geom_text(aes(x = 0,  y = 0, label ="B", fontface = "bold"),size = 6)
AnnotC <- ggplot() + theme_void() + geom_text(aes(x = 0,  y = 0, label ="C", fontface = "bold"),size = 6)
AnnotD <- ggplot() + theme_void() + geom_text(aes(x = 0,  y = 0, label ="D", fontface = "bold"),size = 6)
AnnotE <- ggplot() + theme_void() + geom_text(aes(x = 0,  y = 0, label ="E", fontface = "bold"),size = 6)
AnnotF <- ggplot() + theme_void() + geom_text(aes(x = 0,  y = 0, label ="F", fontface = "bold"),size = 6)
void <- ggplot() + geom_blank() + theme_void()

# All ---------------------------------------------------------------------
library(patchwork)
design  <-
  c(
    ggLegShape = area(4, 2, 4, 19),
    Main3B1 = area(5, 2, 11, 19),
    ggLegCol = area(12, 2, 12, 19),
    Main3B2 = area(13, 2, 19, 19),
    Main3Intes1 = area(5, 21, 11, 30),
    Main3Intes1Leg = area(4, 21, 4, 30),
    Main3Olfa1 = area(5, 32, 11, 41),
    LegOlfa1 = area(4, 32, 4, 41),
    Main3Intes2 = area(13, 21, 19, 30),
    TitleMain2Intes = area(12, 21, 12, 30),
    Main3Olfa2 = area(13, 32, 19, 41),
    TitleMain2Olfa = area(12, 32, 12, 41),
    AnnotA = area(2, 1, 3, 1),
    AnnotB = area(12, 1, 12, 1),
    AnnotC = area(2, 20, 3, 21),
    AnnotD = area(12, 20, 12, 21),
    AnnotE = area(2, 31, 3, 32),
    AnnotF = area(12, 31, 12, 32)
  )

Main3 <-
  ggLegShape + Main3B1 + 
  ggLegCol + Main3B2 + 
  Main3Intes1 + Main3Intes1Leg + 
  Main3Olfa1 + LegOlfa1 + 
  Main3Intes2 + TitleMain2Intes + 
  Main3Olfa2 + TitleMain2Olfa + AnnotA + AnnotB + AnnotC + AnnotD + AnnotE + AnnotF + plot_layout(design = design)
ggsave("../FinalFigure/Main/Main3.pdf", Main3, units = "mm", dpi = 600, width = 245, height = 160)
ggsave("../FinalFigure/Main/Main3.svg", Main3, units = "mm", dpi = 600, width = 245, height = 160)
ggsave("../FinalFigure/Main/Main3.png", Main3, units = "mm", dpi = 1200, width = 245, height = 160)
