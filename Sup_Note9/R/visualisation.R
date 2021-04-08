library(patchwork)
source("../SupNote3/R/utilitary_function.R")
Pancreas <- readRDS("../SupNote3/data/SeuratPancreas.rds")[[1]]
Airway <- readRDS("../SupNote3/data/SeuratEpithelial.rds")[[1]]

# Pancreas ----------------------------------------------------------------

#color palette
Pancreas_celltype <- c("alpha","beta", "delta", "gamma", "epsilon", "acinar", "ductal",  "PSC", "endothelial", "macrophage", "mast", "t_cell", "schwann", "unassigned", "other")
Pancreas_color <- c("firebrick1", "coral", "dodgerblue", "cyan", "cyan4","chocolate", "chocolate4",  "olivedrab","lightgreen","purple", "violet", "pink", "gold", "lightgrey", "black")
PancreasPalette <- setNames(Pancreas_color, Pancreas_celltype)
ThemeSupFig16 <- theme(plot.margin = unit(c(0,0,0,0), "cm"),aspect.ratio = 1, title = element_text(size = 6), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.position = "none", plot.title = element_text(margin = margin(t = 10, b = -10)))

#axis label change
Pancreas@reductions$mca@key <- "MCA_"
colnames(Pancreas@reductions$mca@cell.embeddings) <- colnames(Pancreas@reductions$mca@cell.embeddings) %>% str_replace("mca", "MCA")

#order 
Pancreas$cell_type1 <-  factor(Pancreas$cell_type1, levels = Pancreas_celltype)

#Plot
AuthorPancreas <- DimPlotMC(Pancreas, group.by = "cell_type1",features = c("VWF", "GCG", "CPA1", "KRT19"), as.text = T) + theme(aspect.ratio = 1, title = element_text(size = 10),legend.text = element_text(size =8, face = "bold"),axis.text = element_text(size = 6, face = "bold"), axis.title = element_text(size =8, face = "bold"), legend.key.size = unit(0.1, "mm")) + guides(color= guide_legend(override.aes = list(size =1))) +scale_color_manual(values = PancreasPalette)
FPlotPancreas  <- FeaturePlot(Pancreas, features = c("VWF", "GCG", "CPA1", "KRT19"), reduction = "mca",combine = F, order =T) 
FPlotPancreas <- lapply(FPlotPancreas, function(x){x + ThemeSupFig16 + scale_color_gradient(low = "lightgrey", high = "red")})


# Airway ------------------------------------------------------------------

#color palette
Epithelial_celltype <- c("Basal", "Krt4.13", "Secretory", "Goblet", "Ciliated", "Brush.Tuft", "Brush+PNEC","PNEC", "Ionocytes", "unassigned", "other")
Epithelial_color <- c("firebrick1", "purple", "dodgerblue", "cyan4", "darkgreen", "chocolate4", "chocolate2","orange", "cyan","lightgrey", "black")
EpithelialPalette <- setNames(Epithelial_color, Epithelial_celltype)


#axis label change
Airway@reductions$mca@key <- "MCA_"
colnames(Airway@reductions$mca@cell.embeddings) <- colnames(Airway@reductions$mca@cell.embeddings) %>% str_replace("mca", "MCA")

#order 
Airway$cell_type1 <-  factor(Airway$cell_type1, levels = Epithelial_celltype)

AuthorAirway <- DimPlotMC(Airway, features = c("Chga", "Krt5", "Scgb1a1", "Foxj1"), as.text = T) + theme(aspect.ratio = 1, title = element_text(size = 10),legend.text = element_text(size =8, face = "bold"),axis.text = element_text(size = 6, face = "bold"), axis.title = element_text(size =8, face = "bold"), legend.key.size = unit(0.1, "mm"),) + guides(color= guide_legend(override.aes = list(size =1))) +scale_color_manual(values = EpithelialPalette)
FPlotAirway  <- FeaturePlot(Airway, features = c("Chga", "Krt5", "Scgb1a1", "Foxj1"), reduction = "mca",combine = F, order =T) 
FPlotAirway <-  lapply(FPlotAirway, function(x){x + ThemeSupFig16 + scale_color_gradient(low = "lightgrey", high = "red")})

LegSupFig14A <- get_legend(AuthorAirway)
LegSupFig14B <- get_legend(AuthorPancreas)
SupFig14A <- AuthorAirway + theme(legend.position = "none")
SupFig14B <- AuthorPancreas + theme(legend.position = "none")


SupFig14C <- ggarrange(plotlist = FPlotAirway)
SupFig14D <- ggarrange(plotlist = FPlotPancreas)


# Pathways ----------------------------------------------------------------

#Get Pathways
download.file("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2019_Mouse", "data/pathway/Wiki")
download.file("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse", "data/pathway/KEGG")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016", "data/pathway/Reactome")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018", "data/pathway/GOBP")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2018", "data/pathway/GOCC")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Molecular_Function_2018", "data/pathway/GOMF")

#Load Pathways
Wiki      <- fgsea::gmtPathways("data/pathway/Wiki")
KEGG      <- fgsea::gmtPathways("data/pathway/KEGG")
Reactome  <- fgsea::gmtPathways("data/pathway/Reactome")
GOBP      <- fgsea::gmtPathways("data/pathway/GOBP")
GOCC      <- fgsea::gmtPathways("data/pathway/GOCC")
GOMF      <- fgsea::gmtPathways("data/pathway/GOMF")

#Convert Pathways
PathHuman <- c(Wiki,KEGG,Reactome,GOBP,GOCC,GOMF)
PathMouse <- lapply(c(Wiki,KEGG,Reactome,GOBP,GOCC,GOMF), str_to_title)

#Run HGT
PancreasHGT <- RunCellHGT(Pancreas, PathHuman, log.trans = T, p.adjust = F)
AirwayHGT <- RunCellHGT(Airway, PathMouse, log.trans = T, p.adjust = F)

#Integrate
Pancreas@assays[["Pathway"]] <- CreateAssayObject(PancreasHGT)
Airway@assays[["Pathway"]] <- CreateAssayObject(AirwayHGT)


D1 <- FeaturePlot(Pancreas, order = T, "vasculogenesis (GO:0001570)", reduction = "mca", max.cutoff = 5, min.cutoff = 2) + ThemeSupFig16 + scale_color_gradient(low = "lightgrey", high = "red")
D2 <- FeaturePlot(Pancreas, order = T, "Pancreatic secretion", reduction = "mca", max.cutoff = 5, min.cutoff = 2) + ThemeSupFig16 + scale_color_gradient(low = "lightgrey", high = "red")
D3 <- FeaturePlot(Airway, order = T, "nervous system development (GO:0007399)", reduction = "mca", max.cutoff = 5, min.cutoff = 2) + ThemeSupFig16 + scale_color_gradient(low = "lightgrey", high = "red")
D4 <- FeaturePlot(Airway, order = T, "cilium movement (GO:0003341)", reduction = "mca", max.cutoff = 5, min.cutoff = 2) + ThemeSupFig16 + scale_color_gradient(low = "lightgrey", high = "red")

D1$labels$title <- "vasculogenesis"
D2$labels$title <- "pancreatic secretion"
D3$labels$title <- "nervous system development"
D4$labels$title <- "cilium movement"

SupFig14E <- ggarrange(plotlist = list(D1,D2))
SupFig14F <- ggarrange(plotlist = list(D3,D4))


#  patchwork --------------------------------------------------------------
AnnotA <- ggplot(mapping = aes(x = 0, y = 0, label = "A")) + geom_text(fontface = "bold", size = 6) + theme_void()
AnnotB <- ggplot(mapping = aes(x = 0, y = 0, label = "B")) + geom_text(fontface = "bold", size = 6) + theme_void()
AnnotC <- ggplot(mapping = aes(x = 0, y = 0, label = "C")) + geom_text(fontface = "bold", size = 6) + theme_void()
AnnotD <- ggplot(mapping = aes(x = 0, y = 0, label = "D")) + geom_text(fontface = "bold", size = 6) + theme_void()
AnnotE <- ggplot(mapping = aes(x = 0, y = 0, label = "E")) + geom_text(fontface = "bold", size = 6) + theme_void()
AnnotF <- ggplot(mapping = aes(x = 0, y = 0, label = "F")) + geom_text(fontface = "bold", size = 6) + theme_void()
design <- c(
  SupFig14A = area(2, 2, 20, 20),
  SupFig14B = area(22, 2, 40, 20),
  SupFig14C = area(2, 22, 20, 40),
  SupFig14D = area(22, 22, 40, 40),
  SupFig14E = area(42, 3, 55, 19),
  SupFig14F = area(42, 23, 55, 39),
  AnnotA = area(1,1,2,2),
  AnnotB = area(1,21,2,22),
  AnnotC = area(21,1,22,2),
  AnnotD = area(21,21,22,22),
  AnnotE = area(41,1,42,2),
  AnnotF = area(41,21,42,22)
)

SupFig16 <- SupFig14A + SupFig14B + SupFig14C + SupFig14D + SupFig14E + SupFig14F + AnnotA + AnnotB + AnnotC + AnnotD + AnnotE + AnnotF + plot_layout(design = design)

ggsave("../FinalFigure/SupFig14.pdf", plot = SupFig16, units = "mm", width = 180, height = 210, dpi = 600)
ggsave("../FinalFigure/SupFig14.png", plot = SupFig16, units = "mm", width = 180, height = 210, dpi = 600)
