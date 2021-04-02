source("R/utilitary_function.R")
SeuratPancreas <- readRDS("data/SeuratPancreas.rds")
SeuratEpithelial <- readRDS("data/SeuratEpithelial.rds")

map_baron_muraro <-
    list(
        c("alpha"),
        c("beta"),
        c("delta"),
        c("gamma"),
        c("epsilon"),
        c("acinar"),
        c("ductal"),
        c("PSC"),
        c("endothelial"),
        c("removed")
    ) %>% set_names(c("alpha", "beta", "delta", "gamma", "epsilon", "acinar", "ductal", "PSC", "endothelial",  "unclassified"))

map_baron_segerstolpe <- list(
    c("alpha"),
    c("beta"),
    c("removed"),
    c("delta"),
    c("gamma"),
    c("epsilon"),
    c("acinar"),
    c("ductal"),
    c("PSC"),
    c("endothelial"),
    c("macrophage"),
    c("mast"),
    c("removed")
) %>% set_names(c("alpha", "beta", "co-expression", "delta", "gamma", "epsilon", "acinar", "ductal", "PSC", "endothelial", "macrophage", "mast",  "unclassified"))

write.xlsx(x = data.frame(Muraro = names(map_baron_muraro), Baron  = as.vector(sapply(map_baron_muraro, function(x) paste0(x, collapse =", " )))), file = "../FinalTable/SupTable10.xlsx", sheetName = "Baron_Muraro", append = F)
write.xlsx(x = data.frame(Segerstolpe = names(map_baron_segerstolpe), Baron  = as.vector(sapply(map_baron_segerstolpe, function(x) paste0(x, collapse =", " )))), file = "../FinalTable/SupTable10.xlsx", sheetName = "Baron_Segerstolpe", append = T)


# Cell Type --------------------------------------------------------------

PancreasBenchCell <- foreach(Seurat = SeuratPancreas[2:3], mapping = list(map_baron_muraro, map_baron_segerstolpe), data = c("Muraro", "Segerstolpe"), .combine = rbind) %:% foreach(predictions = pred_func_name, .combine = rbind) %do% {
    TruthDF <- tibble(Author = as.vector(Seurat$cell_type1), Predictions = as.vector(unlist(Seurat[[predictions]])), map_label = as.vector(sapply(mapping[Seurat$cell_type1], function(x) paste0(x, collapse =", " ))))
    TruthDF <- TruthDF %>% filter(map_label != "removed")
    Unassigned1 <- sum(TruthDF$Predictions == "unassigned")/length(TruthDF$Predictions)
    TruthDF <- TruthDF %>% mutate(positive = mapply(x = .$map_label, y = .$Predictions, FUN =function(x,y){x %like% y}))
    TruthDF <- TruthDF %>% mutate(final_map = ifelse(positive, Author, Predictions))
    PredDF <- sapply(unique(TruthDF$Author), function(x, truth){
        TruthDFpos <- TruthDF %>% filter(Author == x)
        Unassigned2 <- sum(TruthDFpos$Predictions == "unassigned")/length(TruthDFpos$Predictions)
        TruthDFneg <- TruthDF %>% filter(Author != x)
        TP <- sum(TruthDFpos$Author == TruthDFpos$final_map)
        FN <- sum(TruthDFpos$Author != TruthDFpos$final_map)
        FP <- sum(TruthDFneg$final_map %in% mapping[[x]])
        Recall <- TP/ (TP + FN)
        Precision <- TP/ (FP + TP)
        F1 <- 2*(Recall*Precision)/(Recall + Precision) 
        return(c(Precision = Precision, Recall = Recall, F1 = F1, Unassigned1 = Unassigned1, Unassigned2 = Unassigned2))
    },
    truth = TruthDF
    ) %>% as.data.frame() %>%  rownames_to_column(var = "metrics") %>%  gather("cell_type", "value", -1) %>%  mutate(methods = predictions) %>%  mutate(data = data) %>%  mutate(value = ifelse(is.na(value), 0, value))
} %>% select(data, methods, cell_type, everything()) %>%  mutate(methods = str_remove(methods, "pred_")) %>% mutate(methods = str_replace(methods, "CellID_G", "CellID(G)")) %>% mutate(methods = str_replace(methods, "CellID_C", "CellID(C)"), value = round(value, digits =3))
PancreasBenchOverall <- PancreasBenchCell  %>%  group_by(methods, data, metrics) %>%  summarise(value = mean(value)) %>%  select(data, everything())

PancreasBenchOverall <- PancreasBenchCell %>%  
    filter(metrics != "Unassigned2") %>%  
    mutate(metrics = ifelse(metrics == "Unassigned1", "Unassigned", metrics))  %>%  
    group_by(methods, data, metrics) %>%  
    summarise(value = mean(value)) %>%  
    select(data, everything())

PancreasBenchCell <- PancreasBenchCell %>%  
    filter(metrics != "Unassigned1") %>%  
    mutate(metrics = ifelse(metrics == "Unassigned2", "Unassigned", metrics))

write_rds(PancreasBenchCell, path = "data/PancreasBenchCell.rds")
write_rds(PancreasBenchOverall, path = "data/PancreasBenchOverall.rds")

# Epithelial --------------------------------------------------------------

map_plas_plas <-
    list(
        c("Basal", "Krt4.13"),
        c("Brush.Tuft", "PNEC"),
        c("Ciliated"),
        c("removed"),
        c("removed"),
        c("removed"),
        c("Ionocytes"),
        c("Secretory","Krt4.13"),
        c("removed")
    ) %>% set_names("Basal","Brush+PNEC","Ciliated","FOXN4+","Interm. basal>secr.", "Interm. secr.>cil.", "Ionocytes", "Secretory", "SLC16A7+" )

map_plas_mon <-
    list(
        c("Basal", "Krt4.13"),
        c("Brush.Tuft"),
        c("Ciliated"),
        c("removed"),
        c("Ionocytes"),
        c("PNEC"),
        c("Secretory", "Krt4.13")
    ) %>% set_names(sort(unique(SeuratEpithelial[[3]]$cell_type1)))
write.xlsx(x = data.frame(`Plasschaert Human` = names(map_plas_plas), `Plasschaert Mouse`  = as.vector(sapply(map_plas_plas, function(x) paste0(x, collapse =", " )))), file = "../FinalTable/SupTable10.xlsx", sheetName = "PlasschaertMouse_Human", append = T)
write.xlsx(x = data.frame(`Montoro` = names(map_plas_mon), `Plasschaert Mouse`  = as.vector(sapply(map_plas_mon, function(x) paste0(x, collapse =", " )))), file = "../FinalTable/SupTable10.xlsx", sheetName = "PlasschaertMouse_Montoro", append = T)

EpithelialBenchCell <- foreach(Seurat = SeuratEpithelial[2:3], mapping = list(map_plas_plas, map_plas_mon), data = c("Plasschaert Human", "Montoro"), .combine = rbind) %:% foreach(predictions = pred_func_name, .combine = rbind) %do% {
    TruthDF <- tibble(Author = as.vector(Seurat$cell_type1), Predictions = as.vector(unlist(Seurat[[predictions]])), map_label = as.vector(sapply(mapping[Seurat$cell_type1], function(x) paste0(x, collapse =", " ))))
    TruthDF <- TruthDF %>% filter(map_label != "removed")
    Unassigned1 <- sum(TruthDF$Predictions == "unassigned")/length(TruthDF$Predictions)
    TruthDF <- TruthDF %>% mutate(positive = mapply(x = .$map_label, y = .$Predictions, FUN =function(x,y){x %like% y}))
    TruthDF <- TruthDF %>% mutate(final_map = ifelse(positive, Author, Predictions))
    PredDF <- sapply(unique(TruthDF$Author), function(x, truth){
        TruthDFpos <- TruthDF %>% filter(Author == x)
        Unassigned2 <- sum(TruthDFpos$Predictions == "unassigned")/length(TruthDFpos$Predictions)
        TruthDFneg <- TruthDF %>% filter(Author != x)
        TP <- sum(TruthDFpos$Author == TruthDFpos$final_map)
        FN <- sum(TruthDFpos$Author != TruthDFpos$final_map)
        FP <- sum(TruthDFneg$final_map %in% mapping[[x]])
        Recall <- TP/ (TP + FN)
        Precision <- TP/ (FP + TP)
        F1 <- 2*(Recall*Precision)/(Recall + Precision) 
        return(c(Precision = Precision, Recall = Recall, F1 = F1, Unassigned1 = Unassigned1, Unassigned2 = Unassigned2))
    },
    truth = TruthDF
    ) %>% as.data.frame() %>%  rownames_to_column(var = "metrics") %>%  gather("cell_type", "value", -1) %>%  mutate(methods = predictions) %>%  mutate(data = data) %>%  mutate(value = ifelse(is.na(value), 0, value))
} %>% select(data, methods, cell_type, everything()) %>% mutate(methods = str_replace(methods, "CellID_G", "CellID(G)")) %>% mutate(methods = str_replace(methods, "CellID_C", "CellID(C)")) %>% mutate(methods = str_remove(methods, "pred_"), value = round(value, digits =3))

EpithelialBenchOverall <- EpithelialBenchCell %>%  
    filter(metrics != "Unassigned2") %>%  
    mutate(metrics = ifelse(metrics == "Unassigned1", "Unassigned", metrics),value = round(value, digits =3))  %>%  
    group_by(methods, data, metrics) %>%  
    summarise(value = mean(value)) %>%  
    select(data, everything())

EpithelialBenchCell <- EpithelialBenchCell %>%  
    filter(metrics != "Unassigned1") %>%  
    mutate(metrics = ifelse(metrics == "Unassigned2", "Unassigned", metrics))

write_rds(EpithelialBenchCell, path = "data/EpithelialBenchCell.rds")
write_rds(EpithelialBenchOverall, path = "data/EpithelialBenchOverall.rds")
