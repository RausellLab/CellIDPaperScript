source("R/utilitary_function.R")
SeuratEpithelial <- readRDS("../SupNote3/data/SeuratEpithelial.rds")

# Ref Montoro -------------------------------------------------------------
montoro_pred <-
  lapply(
    X = pred_func,
    FUN = function(pred) {
      prediction <- pred(SeuratEpithelial[[4]], SeuratEpithelial[[3]])
      return(as.vector(prediction))
    }
  ) %>%
  as.data.frame() %>%
  set_names(paste0(pred_func_name, "__mon")) %>%
  mutate(
    cell_type1 =
      as.vector(SeuratEpithelial[[4]]$cell_type1)
  )

# Ref Plasschaert ---------------------------------------------------------
plasschaert_pred <-
  lapply(
    X = pred_func,
    FUN = function(pred) {
      prediction <- pred(SeuratEpithelial[[4]], SeuratEpithelial[[1]])
      return(as.vector(prediction))
    }
  ) %>%
  set_names(paste0(pred_func_name, "__pla")) %>%
  as.data.frame() %>%
  mutate(
    cell_type1 =
      as.vector(SeuratEpithelial[[4]]$cell_type1)
  )


# benchmark ---------------------------------------------------------------

SeuratEpithelial[[4]] <- AddMetaData(SeuratEpithelial[[4]], plasschaert_pred[, -13] %>% set_rownames(colnames(SeuratEpithelial[[4]])))
SeuratEpithelial[[4]] <- AddMetaData(SeuratEpithelial[[4]], montoro_pred[, -13] %>% set_rownames(colnames(SeuratEpithelial[[4]])))
SeuratEpithelial[[4]]$cell_type2 <- ifelse(SeuratEpithelial[[4]]$cell_type1 %in% c("Endocrine", "Goblet", "Brush.Tuft"), SeuratEpithelial[[4]]$cell_type1, "Rejection")
SeuratIntestinal <- SeuratEpithelial[[4]]
write_rds(x = SeuratIntestinal, path = "../SupNote4/data/SeuratIntestinal.rds")

# mapping -----------------------------------------------------------------
map_airway_intestinal <-
  list(
    c("unassigned"),
    c("PNEC"),
    c("Goblet", "Secretory"),
    c("Brush.Tuft")
  ) %>% set_names(c(
    "Rejection",
    "Endocrine",
    "Goblet",
    "Brush.Tuft"
  ))
MapAirwayIntestinalDF <- tibble(names(map_airway_intestinal), sapply(map_airway_intestinal, function(x) paste0(x,  collapse = ", "))) %>%  set_colnames(c("Intestinal", "Airway"))
xlsx::write.xlsx(x = MapAirwayIntestinalDF, file = "../FinalTable/SupTable10.xlsx", sheetName = "Airway_Intestinal", append = T)


SeuratIntestinal <- read_rds("data/SeuratIntestinal.rds")
mapping <- map_airway_intestinal

IntestinalBenchCell <- foreach(predictions = pred_func_name, .combine = rbind) %:% foreach(ref = c("__pla", "__mon"), .combine = rbind) %do%
  {
    prediction <- paste0(predictions, ref)
    TruthDF <- tibble(Author = as.vector(Seurat$cell_type2), Predictions = as.vector(unlist(Seurat[[prediction]])), map_label = as.vector(sapply(mapping[Seurat$cell_type2], function(x) paste0(x, collapse = ", "))))
    TruthDF <- TruthDF %>% filter(map_label != "removed")
    TruthDF <- TruthDF %>% mutate(positive = mapply(x = .$map_label, y = .$Predictions, FUN = function(x, y) {
      x %like% y
    }))
    TruthDF <- TruthDF %>% mutate(final_map = ifelse(positive, Author, Predictions))
    PredDF <- sapply(unique(TruthDF$Author), function(x, truth) {
      TruthDFpos <- TruthDF %>% filter(Author == x)
      TruthDFneg <- TruthDF %>% filter(Author != x)
      TP <- sum(TruthDFpos$Author == TruthDFpos$final_map)
      FN <- sum(TruthDFpos$Author != TruthDFpos$final_map)
      FP <- sum(TruthDFneg$final_map %in% mapping[[x]])
      Recall <- TP / (TP + FN)
      Precision <- TP / (FP + TP)
      F1 <- 2 * (Recall * Precision) / (Recall + Precision)
      return(c(Precision = Precision, Recall = Recall, F1 = F1))
    },
    truth = TruthDF
    ) %>%
      as.data.frame() %>%
      rownames_to_column(var = "metrics") %>%
      gather("cell_type", "value", -1) %>%
      mutate(methods = prediction) %>%
      mutate(value = ifelse(is.na(value), 0, value))
  } %>%
  separate(col = "methods", into = c("methods", "data"), sep = "__") %>%
  mutate(data = ifelse(data == "pla", "Plasschaert Mouse", "Montoro")) %>%
  dplyr::select(data, methods, cell_type, everything()) %>%
  mutate(methods = str_remove(methods, "pred_")) %>%
  mutate(methods = str_replace(methods, "CellID_G", "CellID(G)")) %>%
  mutate(methods = str_replace(methods, "CellID_C", "CellID(C)"))

FreqDF <- as.data.frame(table(SeuratIntestinal$cell_type2)) %>%  set_colnames(c("cell_type","Freq"))
IntestinalBenchOverall <- IntestinalBenchCell %>%  filter(cell_type != "Rejection") %>% inner_join(FreqDF) %>% 
  group_by(methods, data, metrics) %>%
  summarise(value = mean(value)) %>%
  dplyr::select(data, everything()) %>%
  ungroup() %>%
  mutate(methods = factor(methods, c("CellID(G)", "CellID(C)", "scmap_cluster", "scmap_cell", "Seurat", "MNN", "SCN", "scPred", "CHETAH", "CaSTLe", "scID", "SingleR"))) %>%
  mutate(metrics = factor(metrics, c("Precision", "Recall", "F1")))


write_rds(IntestinalBenchCell, "data/IntestinalBenchCell.rds")
write_rds(IntestinalBenchOverall, "data/IntestinalBenchOverall.rds")

