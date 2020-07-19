source("R/utilitary_function.R")
SeuratPancreas <- readRDS("data/SeuratPancreas.rds")
ref_seurat <- SeuratPancreas[[1]]
bench <- lapply(X = SeuratPancreas[2:3], function(x, ref_seurat) {
        res1 <- lapply(
            X = pred_func,
            FUN = function(pred, ref, seurat) {
                Time <- system.time(predictions <- pred(seurat, ref))
                return(list(Time, predictions))
            },
            ref = ref_seurat,
            seurat = x
        )
        res2 <- as.data.frame(setNames(lapply(res1, function(x){x[[2]]}),pred_func_name))
        Time <- lapply(res1, function(x){x[[1]][3]})
        return(list(res2, Time))
    },
    ref_seurat = ref_seurat
    )

SeuratPancreas[2:3] <- mapply(x = SeuratPancreas[2:3], y = lapply(bench, function(x) as.data.frame(x[[1]])), FUN = function(x,y) AddMetaData(x,y))
SeuratPancreas[2:3] <- mapply(x = SeuratPancreas[2:3], y = lapply(bench, function(x) as.data.frame(x[[2]]) %>% set_colnames(pred_func_name)), FUN = function(x,y) {x@misc$Time <- y; return(x)})

write_rds(SeuratPancreas, "data/SeuratPancreas.rds")

SeuratEpithelial <- readRDS("data/SeuratEpithelial.rds")
ref_seurat <- SeuratEpithelial[[1]]
bench <- lapply(X = SeuratEpithelial[2:3], function(x, ref_seurat) {
    res1 <- lapply(
        X = pred_func,
        FUN = function(pred, ref, seurat) {
            Time <- system.time(predictions <- pred(seurat, ref))
            return(list(Time, predictions))
        },
        ref = ref_seurat,
        seurat = x
    )
    res2 <- as.data.frame(setNames(lapply(res1, function(x){x[[2]]}),pred_func_name))
    Time <- lapply(res1, function(x){x[[1]][3]})
    return(list(res2, Time))
},
ref_seurat = ref_seurat
)

SeuratEpithelial[2:3] <- mapply(x = SeuratEpithelial[2:3], y = lapply(bench, function(x) as.data.frame(x[[1]])), FUN = function(x,y) AddMetaData(x,y))
SeuratEpithelial[2:3] <- mapply(x = SeuratEpithelial[2:3], y = lapply(bench, function(x) as.data.frame(x[[2]]) %>% set_colnames(pred_func_name)), FUN = function(x,y) {x@misc$Time <- y; return(x)})

write_rds(SeuratEpithelial, "data/SeuratEpithelial.rds")


