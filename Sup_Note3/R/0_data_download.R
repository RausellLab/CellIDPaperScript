library(readxl)
library(writexl)
library(data.table)
library(tidyverse)
library(foreach)

#  Create directory -----------------------------------
create_dir <- function(x){
    if(!dir.exists(x)){
        dir.create(x)  
    }
}

dir_to_create <- c("data","figure","data/cite_seq", "data/reap_seq", "data/immune_signature")
lapply(X = dir_to_create, FUN = create_dir)

# Download Data ---------------------------------------

data_url <- c(
              "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100501&format=file",
              "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz",
              "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DADT%5Fumi%2Ecsv%2Egz"
)

download_data <- function(x,y){
    if(!file.exists(x)){
        download.file(y,
                      destfile = x,
                      mode = "wb")  
    }
}

data_to_download <- c("data/reap_seq/reap_seq_raw.tar",
                      "data/cite_seq/cite_seq.csv.gz",
                      "data/cite_seq/cite_seq_protein.csv.gz"
                      )

mapply(x=data_to_download, y=data_url, FUN = download_data)


#Immune Signature
download.file(destfile = "data/immune_signature/xCell.xlsx", url = "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-017-1349-1/MediaObjects/13059_2017_1349_MOESM3_ESM.xlsx", method = "wget")

xCell_signature <- read_xlsx("data/immune_signature/xCell.xlsx")
xCell_signature <- xCell_signature %>% separate(Celltype_Source_ID,c("cell_type","source","number"),"_") %>% gather(value = "genes", key = "gene_num",5:204)
xCell_signature <- xCell_signature %>% dplyr::filter(!is.na(genes)) %>% dplyr::select(-`# of genes`) %>% dplyr::select(-gene_num) %>% dplyr::arrange(cell_type)
xCell_signature <- xCell_signature %>% dplyr::mutate(data=paste0(source,number)) %>% dplyr::group_by(cell_type) %>% dplyr::mutate(ndata=length(unique(data))) %>% dplyr::group_by(cell_type,genes) %>% summarise(occurence=length(genes)/unique(ndata)) %>% dplyr::arrange(occurence)
xCell_signature <- xCell_signature %>% dplyr::group_by(cell_type) %>%  dplyr::arrange(cell_type,-occurence) %>% dplyr::filter(occurence>=0.3)
xCell_signature <- xCell_signature %>% split(xCell_signature$cell_type)
immune_signature <- lapply(xCell_signature,function(x)x$genes)
immune_signature$`CD14+ Monocytes` <- immune_signature$Monocytes
Mono16 <- list(c("FCGR3A","CDKN1C","MTSS1","SIGLEC10","IFITM1","HMOX1","TAGLN","MS4A7","CSF1R","IFITM2","SOD1","CX3CR1","LILRB1","PSCDBP","ITGAL","C6orf187","KLF2","WARS","MAFB","GCH1","CD97","PIK3AP1","MAIL","LYN","BCL2A1","PECAM1"))
names(Mono16) <- "CD16+ Monocytes"
immune_signature <- c(immune_signature, Mono16)
immune_signature$Monocytes <- NULL
names(immune_signature) <- str_remove(names(immune_signature), pattern = "\\+")
immune_signature <- immune_signature[c(
    "HSC",
    "MPP",
    "CMP",
    "GMP",
    "MEP",
    "Erythrocytes",
    "Megakaryocytes",
    "Platelets",
    "B-cells",
    "CD4 T-cells",
    "CD8 T-cells",
    "NK cells",
    "Basophils",
    "Eosinophils",
    "Neutrophils",
    "CD14 Monocytes",
    "CD16 Monocytes",
    "Macrophages",
    "DC",
    "cDC",
    "pDC"
)]
write_rds(immune_signature,"data/immune_signature/immune_signature.rds")


# SupTable1 --------------------------------------------------------------


tabmax <- max(sapply(immune_signature, length))
SupTable1 <- as_tibble(sapply(immune_signature, function(x){
    c(x,as.vector(rep(NA,tabmax -length(x))))
}))

write_xlsx(SupTable1, "../FinalTable/SupTable1.xlsx", col_names = T)
write_tsv(SupTable1, "../FinalTable/SupTable1.tsv", col_names = T)
