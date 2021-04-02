source("R/utilitary_function.R")
SeuratPancreas <- readRDS("data/SeuratPancreas.rds")

dir.create("data/pathway")
download.file("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2019_Human", "data/pathway/Wiki.gmt")
download.file("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human", "data/pathway/KEGG.gmt")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016", "data/pathway/Reactome.gmt")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018", "data/pathway/GOBP.gmt")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2018", "data/pathway/GOCC.gmt")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Molecular_Function_2018", "data/pathway/GOMF.gmt")

Wiki     <- fgsea::gmtPathways("data/pathway/Wiki.gmt")
KEGG     <- fgsea::gmtPathways("data/pathway/KEGG.gmt")
Reactome <- fgsea::gmtPathways("data/pathway/Reactome.gmt")
GOBP     <- fgsea::gmtPathways("data/pathway/GOBP.gmt")
GOCC     <- fgsea::gmtPathways("data/pathway/GOCC.gmt")
GOMF     <- fgsea::gmtPathways("data/pathway/GOMF.gmt")

pathways_onto <- c(Wiki,KEGG,Reactome,GOBP,GOCC,GOMF)
label <- mapply(x = list(Wiki,KEGG,Reactome,GOBP,GOCC,GOMF), 
                y= c("WikiPathways", "KEGG", "Reactome","GO Biological Process","GO Cellular Component","GO Molecular Function"), 
                FUN = function(x,y){rep(y,length(x))}, 
                SIMPLIFY = T) %>% unlist
names(pathways_onto) <- paste0(names(pathways_onto),"___", label)

HGT_Muraro <- RunCellHGT(SeuratPancreas[[2]], log.trans = T, pathways = pathways_onto,  minSize = 5, p.adjust = F)
HGT_Segers <- RunCellHGT(SeuratPancreas[[3]], log.trans = T, pathways = pathways_onto,  minSize = 5, p.adjust = F)
MuraroCell <- names(SeuratPancreas[[2]]$pred_CellID_G[SeuratPancreas[[2]]$pred_CellID_G == "schwann"])
SegersCell <- names(SeuratPancreas[[3]]$pred_CellID_G[SeuratPancreas[[3]]$pred_CellID_G == "schwann"])
HGT_Muraro_summary <- HGT_Muraro[,MuraroCell] %>%  as.data.frame() %>%  rownames_to_column("geneset") %>%  gather("cell", "value", -1) %>%  group_by(geneset) %>%  summarise(median_log_10_pval = median(value), mean_log_10_pval = mean(value), sd_log_10_pval = sd(value))
HGT_Segers_summary <- HGT_Segers[,SegersCell] %>%  as.data.frame() %>%  rownames_to_column("geneset") %>%  gather("cell", "value", -1) %>%  group_by(geneset) %>%  summarise(median_log_10_pval = median(value), mean_log_10_pval = mean(value), sd_log_10_pval = sd(value)) %>%  arrange(-mean_log_10_pval)
HGT_Muraro_summary_filtered <- HGT_Muraro_summary %>%  arrange(-median_log_10_pval) %>%  filter(median_log_10_pval>2) %>%  separate(col = "geneset", into = c("geneset", "database"), sep = "___")
HGT_Segers_summary_filtered <- HGT_Segers_summary %>%  arrange(-median_log_10_pval) %>%  filter(median_log_10_pval>2) %>%  separate(col = "geneset", into = c("geneset", "database"), sep = "___")
GSMuraro <- GetCellGeneSet(SeuratPancreas[[2]], cells = MuraroCell)
GSSegers <- GetCellGeneSet(SeuratPancreas[[3]], cells =  SegersCell)
HGT_Muraro_summary_filtered$`50 pecent overlap` <- unlist(foreach(x = paste0(HGT_Muraro_summary_filtered$geneset, "___", HGT_Muraro_summary_filtered$database)) %:% foreach(y = GSMuraro, .combine = c, .final = function(x) paste(names(table(x))[table(x)/length(GSMuraro)>0.5], collapse = ", ")) %do% {pathways_onto[[x]][pathways_onto[[x]] %in% y]})
HGT_Segers_summary_filtered$`50 pecent overlap` <- unlist(foreach(x = paste0(HGT_Segers_summary_filtered$geneset, "___", HGT_Segers_summary_filtered$database)) %:% foreach(y = GSSegers, .combine = c, .final = function(x) paste(names(table(x))[table(x)/length(GSSegers)>0.5], collapse = ", ")) %do% {pathways_onto[[x]][pathways_onto[[x]] %in% y]})
PancreasPathwayTable <- list(Muraro_pathways = HGT_Muraro_summary_filtered, Segers_pathways = HGT_Segers_summary_filtered)

segersxlsx <- PancreasPathwayTable$Segers_pathways %>%  arrange(database, -median_log_10_pval) %>%  mutate(median_pvalue = 10^(-median_log_10_pval)) %>%  select(database, geneset, median_pvalue, everything())
muraroxlsx <- PancreasPathwayTable$Muraro_pathways %>%  arrange(database, -median_log_10_pval) %>%  mutate(median_pvalue = 10^(-median_log_10_pval)) %>%  select(database, geneset, median_pvalue, everything())
write_rds(segersxlsx, "data/segersxlsx.rds")
write_rds(muraroxlsx, "data/muraroxlsx.rds")

xlsx::write.xlsx(PancreasPathwayTable$Segers_pathways, file = "../FinalTable/SupTable6.xlsx", sheetName = "schwann_Muraro", append = F)
xlsx::write.xlsx(PancreasPathwayTable$Segers_pathways, file = "../FinalTable/SupTable6.xlsx", sheetName = "schwann_Segerstolpe", append =T)

SeuratPancreas[[2]]@assays[["Pathway"]] <- CreateAssayObject(data = HGT_Muraro) 
SeuratPancreas[[3]]@assays[["Pathway"]] <- CreateAssayObject(data = HGT_Segers) 

FeaturePlot(SeuratPancreas[[2]], features = "Neural Crest Differentiation WP2064---WikiPathways", min.cutoff = 2) + scale_color_gradient(low = "grey", high = "red")
FeaturePlot(SeuratPancreas[[3]], features = "Neural Crest Differentiation WP2064---WikiPathways", min.cutoff = 2) + scale_color_gradient(low = "grey", high = "red")

FeaturePlot(SeuratPancreas[[2]], features = "Collagen formation Homo sapiens R-HSA-1474290---Reactome", min.cutoff = 2) + scale_color_gradient(low = "grey", high = "red")
FeaturePlot(SeuratPancreas[[3]], features = "Collagen formation Homo sapiens R-HSA-1474290---Reactome", min.cutoff = 2) + scale_color_gradient(low = "grey", high = "red")

FeaturePlot(SeuratPancreas[[2]], features = "positive regulation of ossification (GO:0045778)---GO Biological Process", min.cutoff = 2, order= T) + scale_color_gradient(low = "grey", high = "red")
FeaturePlot(SeuratPancreas[[3]], features = "positive regulation of ossification (GO:0045778)---GO Biological Process", min.cutoff = 2, order= T) + scale_color_gradient(low = "grey", high = "red")

FeaturePlot(SeuratPancreas[[2]], features = "Neural Crest Differentiation WP2064---WikiPathways", min.cutoff = 2) + scale_color_gradient(low = "grey", high = "red")
FeaturePlot(SeuratPancreas[[3]], features = "Neural Crest Differentiation WP2064---WikiPathways", min.cutoff = 2) + scale_color_gradient(low = "grey", high = "red")
