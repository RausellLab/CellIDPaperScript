source("R/utilitary_function.R")

if(!dir.exists("data")){
    dir.create("data")  
}


# Pancreas ------------------------------------

if(!dir.exists("data/baron")){
    dir.create("data/baron")  
}
if(!dir.exists("data/muraro")){
    dir.create("data/muraro")  
}
if(!dir.exists("data/segerstolpe")){
    dir.create("data/segerstolpe")  
}

if(!dir.exists("data/baronM")){
    dir.create("data/baronM")  
}


# Download Data ---------------------------------------
baron       <- scRNAseq::BaronPancreasData()
baronM      <- scRNAseq::BaronPancreasData(which = "mouse")
muraro      <- scRNAseq::MuraroPancreasData()
segerstolpe <- scRNAseq::SegerstolpePancreasData()
write_rds(baron,"data/baron/baron.rds")
write_rds(baronM,"data/baronM/baronM.rds")
write_rds(muraro,"data/muraro/muraro.rds")
write_rds(segerstolpe,"data/segerstolpe/segerstolpe.rds")


# Epithelial --------------------------------------------------------------

# Create Directory ------------------------------------

create_dir <- function(x){
    if(!dir.exists(x)){
        dir.create(x)  
    }
}

dir_to_create <- c("data",
                   "data/plasschaert.mouse",
                   "data/plasschaert.human",
                   "data/haber",
                   "data/montoro")

lapply(X = dir_to_create, FUN = create_dir)

library(BiocFileCache)
bfc <- BiocFileCache("/home/agarcia/.cache/ExperimentHub")
res <- bfcquery(bfc, "2698", field="rpath", exact=FALSE)
bfcinfo(bfc, rid="BFC16") %>% dplyr::select(rpath)
fl <- file.path(hubCache(eh), "18f183d72124f_2698")
unlink(fl)
# Download data and metadata --------------------------


data_to_download <- c("data/haber/haber.txt.gz",
                      "data/plasschaert.human/plasschaert.human.txt.gz",
                      "data/plasschaert.mouse/plasschaert.mouse.txt.gz",
                      "data/montoro/montoro.txt.gz")

download_data <- function(x,y){
    if(!file.exists(x)){
        download.file(y,
                      destfile = x, method = "wget", quiet = T)  
    }
}
data_url <- c("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332%5Fatlas%5FUMIcounts%2Etxt%2Egz",
              "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102580/suppl/GSE102580%5Ffiltered%5Fnormalized%5Fcounts%5Fhuman%2Etsv%2Egz",
              "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102580/suppl/GSE102580%5Ffiltered%5Fnormalized%5Fcounts%5Fmouse%2Etsv%2Egz",
              "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103354/suppl/GSE103354%5FTrachea%5Fdroplet%5FUMIcounts%2Etxt%2Egz"
)

meta_data_to_download <- c(
    "data/plasschaert.human/plasschaert.human_meta.txt.gz",
    "data/plasschaert.mouse/plasschaert.mouse_meta.txt.gz")


meta_data_url <- c("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102580/suppl/GSE102580%5Fmeta%5Ffiltered%5Fcounts%5Fhuman%2Etsv%2Egz",
                   "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102580/suppl/GSE102580%5Fmeta%5Ffiltered%5Fcounts%5Fmouse%2Etsv%2Egz")


mapply(x=data_to_download,y=data_url,FUN = download_data)
mapply(x=meta_data_to_download,y=meta_data_url,FUN = download_data)
