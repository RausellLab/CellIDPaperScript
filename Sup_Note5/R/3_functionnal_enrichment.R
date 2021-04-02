
# Loading Data ------------------------------------------------------------

source("R/utilitary_function.R")
SeuratFletcher <- read_rds("data/SeuratFletcher.rds")
SeuratWu <- read_rds("data/SeuratWu.rds")

download.file("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2019_Mouse", "data/pathway/Wiki")
download.file("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse", "data/pathway/KEGG")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016", "data/pathway/Reactome")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018", "data/pathway/GOBP")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2018", "data/pathway/GOCC")
download.file("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Molecular_Function_2018", "data/pathway/GOMF")

Wiki <- fgsea::gmtPathways("data/pathway/Wiki")
BioPlanet <- fgsea::gmtPathways("data/pathway/BioPlanet")
KEGG <- fgsea::gmtPathways("data/pathway/KEGG")
BioCarta <- fgsea::gmtPathways("data/pathway/BioCarta")
HumanCyc <- fgsea::gmtPathways("data/pathway/HumanCyc")
Reactome <- fgsea::gmtPathways("data/pathway/Reactome")
GOBP <- fgsea::gmtPathways("data/pathway/GOBP")
GOCC <- fgsea::gmtPathways("data/pathway/GOCC")
GOMF <- fgsea::gmtPathways("data/pathway/GOMF")

pathways_onto <- c(Wiki, KEGG, Reactome, GOBP, GOCC, GOMF)
pathways_onto <- lapply(pathways_onto, str_to_title)
label <- mapply(
  x = list(Wiki, KEGG, Reactome, GOBP, GOCC, GOMF),
  y = c("WikiPathways", "KEGG", "Reactome", "GO Biological Process", "GO Cellular Component", "GO Molecular Function"),
  FUN = function(x, y) {
    rep(y, length(x))
  },
  SIMPLIFY = T
) %>% unlist()
names(pathways_onto) <- paste0(names(pathways_onto), "___", label)

# Enrichment --------------------------------------------------------------

HGT_Wu <- RunCellHGT(SeuratWu, pathways = pathways_onto, log.trans = T, n.features = 200, minSize = 5, p.adjust = F)
HGT_Fletcher <- RunCellHGT(SeuratFletcher, pathways = pathways_onto, log.trans = T, n.features = 200, minSize = 5, p.adjust = F)
SeuratWu@assays[["FUNC"]] <- CreateAssayObject(HGT_Wu)
SeuratFletcher@assays[["FUNC"]] <- CreateAssayObject(HGT_Fletcher)
write_rds(SeuratWu, "data/SeuratWu.rds")
write_rds(SeuratFletcher, "data/SeuratFletcher.rds")


# Summary DF --------------------------------------------------------------

#calculate log10 mean median

HGT_Wu_summary <- HGT_Wu[, WuCell] %>%
  as.data.frame() %>%
  rownames_to_column("geneset") %>%
  gather("cell", "value", -1) %>%
  group_by(geneset) %>%
  summarise(median_log_10_pval = median(value), mean_log_10_pval = mean(value), sd_log_10_pval = sd(value))
HGT_Fletcher_summary <- HGT_Fletcher[, FletcherCell] %>%
  as.data.frame() %>%
  rownames_to_column("geneset") %>%
  gather("cell", "value", -1) %>%
  group_by(geneset) %>%
  summarise(median_log_10_pval = median(value), mean_log_10_pval = mean(value), sd_log_10_pval = sd(value))

# geneset size
HGT_Wu_summary$size <- sapply(HGT_Wu_summary$geneset, function(x) {
  length(pathways_onto[[x]])
})
HGT_Fletcher_summary$size <- sapply(HGT_Fletcher_summary$geneset, function(x) {
  length(pathways_onto[[x]])
})

# filter high pvalue
HGT_Wu_summary_filtered <- HGT_Wu_summary %>%
  arrange(-median_log_10_pval) %>%
  filter(median_log_10_pval > 2) %>%
  separate(col = "geneset", into = c("geneset", "database"), sep = "___")
HGT_Fletcher_summary_filtered <- HGT_Fletcher_summary %>%
  arrange(-median_log_10_pval) %>%
  filter(median_log_10_pval > 2) %>%
  separate(col = "geneset", into = c("geneset", "database"), sep = "___")


# Calculate Overlap -------------------------------------------------------

GSWu <- GetCellGeneSet(SeuratWu, n.features = 200, dims = 1:50, cells = WuCell)
GSFletcher <- GetCellGeneSet(SeuratFletcher, n.features = 200, dims = 1:50, cells = FletcherCell)
HGT_Wu_summary_filtered$`50 pecent cells overlap` <- unlist(foreach(x = paste0(HGT_Wu_summary_filtered$geneset, "___", HGT_Wu_summary_filtered$database)) %:% foreach(y = GSWu, .combine = c, .final = function(x) paste(names(table(x))[table(x) / length(GSWu) > 0.5], collapse = ", ")) %do% {
  pathways_onto[[x]][pathways_onto[[x]] %in% y]
})
HGT_Fletcher_summary_filtered$`50 pecent cells overlap` <- unlist(foreach(x = paste0(HGT_Fletcher_summary_filtered$geneset, "___", HGT_Fletcher_summary_filtered$database)) %:% foreach(y = GSFletcher, .combine = c, .final = function(x) paste(names(table(x))[table(x) / length(GSFletcher) > 0.5], collapse = ", ")) %do% {
  pathways_onto[[x]][pathways_onto[[x]] %in% y]
})


# Organize Table ----------------------------------------------------------

OlfaPathwayTable <- list(Wu_pathways = HGT_Wu_summary_filtered, Fletcher_pathways = HGT_Fletcher_summary_filtered)

Wulsxlsx <- OlfaPathwayTable$Wu_pathways %>%
  arrange(database, -median_log_10_pval) %>%
  mutate(median_pvalue = 10^(-median_log_10_pval)) %>%
  dplyr::select(database, geneset, median_pvalue, everything())
Fletcherxlsx <- OlfaPathwayTable$Fletcher_pathways %>%
  arrange(database, -median_log_10_pval) %>%
  mutate(median_pvalue = 10^(-median_log_10_pval)) %>%
  dplyr::select(database, geneset, median_pvalue, everything())
write_rds(Wulsxlsx, "data/Wulsxlsx.rds")
write_rds(Fletcherxlsx, "data/Fletcherxlsx.rds")
xlsx::write.xlsx(Wulsxlsx, file = "../FinalTable/SupTable6.xlsx", sheetName = "SCC_Wu", append = T)
xlsx::write.xlsx(Fletcherxlsx, file = "../FinalTable/SupTable6.xlsx", sheetName = "SCC_Fletcher", append = T)

# Overlap -----------------------------------------------------------------

median(unlist(lapply(GSWu, function(x) {
  sum(x %in% GSAirway$Brush.Tuft) / 200
})))
sd(unlist(lapply(GSWu, function(x) {
  sum(x %in% GSAirway$Brush.Tuft) / 200
})))
median(unlist(lapply(GSFletcher, function(x) {
  sum(x %in% GSAirway$Brush.Tuft) / 200
})))
sd(unlist(lapply(GSFletcher, function(x) {
  sum(x %in% GSAirway$Brush.Tuft) / 200
})))

median(unlist(lapply(GSWu, function(x) {
  sum(x %in% GSIntestinal$Brush.Tuft) / 200
})))
sd(unlist(lapply(GSWu, function(x) {
  sum(x %in% GSIntestinal$Brush.Tuft) / 200
})))
median(unlist(lapply(GSFletcher, function(x) {
  sum(x %in% GSIntestinal$Brush.Tuft) / 200
})))
sd(unlist(lapply(GSFletcher, function(x) {
  sum(x %in% GSIntestinal$Brush.Tuft) / 200
})))

write_rds(SeuratWu, "data/SeuratWu.rds")
write_rds(SeuratFletcher, "data/SeuratFletcher.rds")
