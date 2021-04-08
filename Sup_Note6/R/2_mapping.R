source("R/utilitary_function.R")
seurat_ATAC_sub_SS      <- read_rds(path = "data/ATAC/seurat_ATAC_sub_SS.rds")
seurat_ATAC_sub_10X     <- read_rds(path = "data/ATAC/seurat_ATAC_sub_10X.rds")
seurat_SS_sub           <- read_rds(path = "data/SmartSeq/seurat_SS_sub.rds")
seurat_10X_sub          <- read_rds(path = "data/10X/seurat_10X_sub.rds")
n10X_sub_cell_type      <- seurat_10X_sub$cell_type1 %>% unique %>%  length()
nSS_sub_cell_type       <- seurat_SS_sub$cell_type1 %>% unique %>%  length()
nATAC_sub_SS_cell_type  <- seurat_ATAC_sub_SS$cell_type1 %>% unique %>%  length()
nATAC_sub_10X_cell_type <- seurat_ATAC_sub_10X$cell_type1 %>% unique %>%  length()

# mapping SS-----------------------------------------------------------------
ATAC_cell_type_plot_SS <- sort(unique(seurat_ATAC_sub_SS$cell_type1)) %>% str_subset("Collisions|Unknown",negate = T)

mapSS <- list(
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("B cell")),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("macro")),
  "removed",
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("B cell")),
  ("cardiac muscle cell"),
  ("removed"),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("duct|tubule")),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("duct|tubule")),
  (c("dendritic cells")),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("tubule")),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("endothelial")),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("endothelial")),
  c((sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("endothelial")), "fibroblast", "stromal cell"),
  sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("enterocyte|large intestine"),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("erythro")),
  ("removed"),
  (c("hematopoietic precursor cell", "granulocytopoietic cell", "early pro-B cell", "late pro-B cell", "granulocyte", "leukocyte", "promonocyte")),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("hepato")),
  c("immature B cell","late pro-B cell","precursor B cell"),
  ("removed"),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("kidney")),
  ("removed"),
  ("removed"),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("mono")),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("natural")),
  ("removed"),
  ("removed"),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("proximal")),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("proximal")),
  ("removed"),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("T cell")),
  ("removed"),
  (sort(unique(seurat_SS_sub$cell_type1)) %>% str_subset("T cell")),
  "epithelial cell of lung",
  "epithelial cell of lung"
) %>%  set_names(ATAC_cell_type_plot_SS)

mapSStable <- tibble(ATAC_cell_type_plot_SS, sapply(mapSS, function(x) paste0(x,  collapse = ", ")))
write_rds(mapSStable, "data/mapSStable.rds")
xlsx::write.xlsx(x = mapSStable, file = "../FinalTable/SupTable10.xlsx", append = T, sheetName = "ATAC_SS.RNA")

# mapping 10X-----------------------------------------------------------------

ATAC_10X_cell_type1 <- sort(unique(seurat_ATAC_sub_10X$cell_type1)) %>% str_subset("Collisions|Unknown",negate = T)
map10X <- list(
  (sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("B cell")),
  ("alveolar macrophage"),
  ("B cell"),
  ("cardiac muscle cell"),
  ("removed"),
  (sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("duct|tubule")),
  (sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("duct|tubule")),
  (c("dendritic cells")),
  (sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("duct|tubule")),
  (sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("endothelial")),
  (sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("endothelial")),
  c(sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("endothelial"), "fibroblast", "stromal cell"),
  ("removed"),
  (sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("erythro")),
  ("removed"),
  (c("hematopoietic precursor cell", "granulocytopoietic cell", "early pro-B cell", "late pro-B cell", "granulocyte", "leukocyte", "promonocyte")),
  c("hepatocyte"),
  (c("immature B cell","early pro-B cell","late pro-B cell", "Fraction A pre-pro B cell")),
  ("removed"),
  ("kidney loop of Henle ascending limb epithelial cell"),
  ("macrophage"),
  ("removed"),
  c(sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("mono"),"leukocyte"),
  ("natural killer cell"),
  ("removed"),
  ("removed"),
  (sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("proximal")),
  (sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("proximal")),
  ("removed"),
  (sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("T cell")),
  ("removed"),
  (sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("T cell")),
  c(sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("pneumocyte"), "lung endothelial cell"),
  (sort(unique(seurat_10X_sub$cell_type1)) %>% str_subset("pneumocyte"))
) %>%  set_names(ATAC_10X_cell_type1)


map10Xtable <- tibble(ATAC_10X_cell_type1, sapply(map10X, function(x) paste0(x,  collapse = ", "))) %>%  set_colnames(c("ATAC", "10X.RNA"))
write_rds(map10Xtable, "data/map10Xtable.rds")
xlsx::write.xlsx(x = map10Xtable, file = "../FinalTable/SupTable10.xlsx", append = T, sheetName = "ATAC_10X.RNA")



ATACBenchCell <- foreach(Seurat = list(seurat_ATAC_sub_10X, seurat_ATAC_sub_SS), mapping = list(map10X, mapSS), data = c("10X", "SS"), .combine = rbind) %:% foreach(predictions = pred_func_name, .combine = rbind) %do% {
  TruthDF <- tibble(Author = as.vector(Seurat$cell_type1), Predictions = as.vector(unlist(Seurat[[predictions]])), map_label = as.vector(sapply(mapping[Seurat$cell_type1], function(x) paste0(x, collapse =", " ))))
  TruthDF <- TruthDF %>%  filter(str_detect(Author,"Collisions|Unknown",negate = T))
  TruthDF <- TruthDF %>% filter(map_label != "removed")
  TruthDF <- TruthDF %>% mutate(positive = mapply(x = .$map_label, y = .$Predictions, FUN =function(x,y){x %like% y}))
  TruthDF <- TruthDF %>% mutate(final_map = ifelse(positive, Author, Predictions))
  PredDF <- sapply(unique(TruthDF$Author), function(x, truth){
    TruthDFpos <- TruthDF %>% filter(Author == x)
    TruthDFneg <- TruthDF %>% filter(Author != x)
    TP <- sum(TruthDFpos$Author == TruthDFpos$final_map)
    FN <- sum(TruthDFpos$Author != TruthDFpos$final_map)
    FP <- sum((TruthDFneg$final_map %in% mapping[[x]]) & (!TruthDFneg$positive))
    Recall <- TP/ (TP + FN)
    Precision <- TP/ (FP + TP)
    F1 <- 2*(Recall*Precision)/(Recall + Precision) 
    return(c(Precision = Precision, Recall = Recall, F1 = F1))
  },
  truth = TruthDF
  ) %>% as.data.frame() %>%  rownames_to_column(var = "metrics") %>%  gather("cell_type", "value", -1) %>%  mutate(methods = predictions) %>%  mutate(data = data) %>%  mutate(value = ifelse(is.na(value), 0, value))
} %>% select(data, methods, cell_type, everything()) %>%  mutate(methods = str_remove(methods, "pred_")) %>% mutate(methods = str_replace(methods, "CellID_G", "CellID(G)")) %>% mutate(methods = str_replace(methods, "CellID_C", "CellID(C)"))



df <- rbind(as.data.frame(table(seurat_ATAC_sub_10X$cell_type1)) %>%  mutate(data = "10X"), as.data.frame(table(seurat_ATAC_sub_SS$cell_type1))%>%  mutate(data = "SS")) %>% set_colnames(c("cell_type", "Freq", "data"))
ATACBenchOverall <- ATACBenchCell %>%  inner_join(df, by = c("data","cell_type")) %>% mutate(weighted = value*Freq)  %>%  group_by(methods, data, metrics) %>%  summarise(value = sum(weighted)/sum(Freq)) %>%  select(data, everything()) %>% ungroup %>%   mutate(metrics = factor(metrics, levels = c("Precision", "Recall", "F1"))) %>%  mutate(value = round(value, digits =3)) %>%
  mutate(methods = factor(methods, c("CellID(G)", "CellID(C)", "scmap_cluster", "scmap_cell", "Seurat", "MNN", "SCN", "scPred", "CHETAH", "CaSTLe", "scID", "SingleR")))

ATACBenchOverall <- ATACBenchCell %>%  inner_join(df, by = c("data","cell_type"))   %>%  group_by(methods, data, metrics) %>% summarise(value = mean(value)) %>%    select(data, everything()) %>% ungroup %>%   mutate(metrics = factor(metrics, levels = c("Precision", "Recall", "F1"))) %>%  mutate(value = round(value, digits =3)) %>%
  mutate(methods = factor(methods, c("CellID(G)", "CellID(C)", "scmap_cluster", "scmap_cell", "Seurat", "MNN", "SCN", "scPred", "CHETAH", "CaSTLe", "scID", "SingleR")))

Bench1 <- ATACBenchCell %>% mutate(value = round(value, digits =3)) %>%   spread(methods, value) %>%  rename(`Reference Data` = data)
Bench2 <- ATACBenchOverall %>% mutate(value = round(value, digits =3)) %>%   spread(methods, value) %>%  rename(`Reference Data` = data)
xlsx::write.xlsx(Bench1, file = "../FinalTable/SupTable5.xlsx", sheetName = "Cell_Population", append = T)
xlsx::write.xlsx(Bench2, file = "../FinalTable/SupTable5.xlsx", sheetName = "Overall", append = T)


write_rds(ATACBenchOverall, "data/ATACBenchOverall.rds")
write_rds(ATACBenchCell, "data/ATACBenchCell.rds")
