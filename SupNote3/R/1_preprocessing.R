set.seed(1)
source("R/utilitary_function.R")
# Load Data -------------------------------------------
list_path_raw_pncreas <- c("data/baron/baron.rds","data/muraro/muraro.rds","data/segerstolpe/segerstolpe.rds", "data/baronM/baronM.rds")
list_sce_pancreas <- lapply(list_path_raw_pncreas,  FUN = read_rds)

# Intersect of Genes ----------------------------------
rownames(list_sce_pancreas[[2]]) <- str_remove(rownames(list_sce_pancreas[[2]]),pattern = "__.*$") #change weird genes annotation
list_sce_pancreas <- lapply(X = list_sce_pancreas,  FUN = function(x){x[names(which(nexprs(counts(x), byrow = T) >= 5)),]})
converter <- convertMouseGeneList(rownames(list_sce_pancreas[[4]]))
list_sce_pancreas[[4]] <- list_sce_pancreas[[4]][rownames(list_sce_pancreas[[4]]) %in% names(converter)]
converted_genes <- converter[rownames(list_sce_pancreas[[4]])]
list_sce_pancreas[[4]] <- list_sce_pancreas[[4]][names(converted_genes)]
rownames(list_sce_pancreas[[4]]) <- converted_genes
list_sce_pancreas <- lapply(X = list_sce_pancreas,  FUN = function(x){x[!duplicated(rownames(x)),]})
list_sce_pancreas <- lapply(list_sce_pancreas, function(x){x[rownames(x) %in% HgProteinCodingGenes,]})
# Remove unbiguous cluster cells ----------------------



list_sce_pancreas[[2]] <- list_sce_pancreas[[2]][,!is.na(list_sce_pancreas[[2]]$label)] 
list_sce_pancreas[[3]] <- list_sce_pancreas[[3]][, (!is.na(list_sce_pancreas[[3]][["cell type"]]))]
list_sce_pancreas[[3]] <- list_sce_pancreas[[3]][, (list_sce_pancreas[[3]][["cell type"]] != "unclassified endocrine cell")]
list_sce_pancreas[[3]][["cell type"]] <- str_remove(list_sce_pancreas[[3]][["cell type"]] ," cell")


list_sce_pancreas[[1]]$cell_type1 <- list_sce_pancreas[[1]]$label
list_sce_pancreas[[2]]$cell_type1 <- list_sce_pancreas[[2]]$label
list_sce_pancreas[[3]]$cell_type1 <- list_sce_pancreas[[3]][["cell type"]]
list_sce_pancreas[[4]]$cell_type1 <- list_sce_pancreas[[4]]$label

list_sce_pancreas[[1]]$cell_type1 <- str_replace(list_sce_pancreas[[1]]$cell_type1,pattern = ".*stellate", "PSC")

list_sce_pancreas[[2]]$cell_type1 <- str_replace(list_sce_pancreas[[2]]$cell_type1, "mesenchymal", "PSC")
list_sce_pancreas[[2]]$cell_type1 <- str_replace(list_sce_pancreas[[2]]$cell_type1, "duct", "ductal")
list_sce_pancreas[[2]]$cell_type1 <- str_replace(list_sce_pancreas[[2]]$cell_type1, "unclear", "unclassified")
list_sce_pancreas[[2]]$cell_type1 <- str_replace(list_sce_pancreas[[2]]$cell_type1, "pp", "gamma")

list_sce_pancreas[[3]]$cell_type1 <- str_replace(list_sce_pancreas[[3]]$cell_type1, "MHC class II", "macrophage")

# Change name for readability -------------------------
colnames(list_sce_pancreas[[1]]) <- paste0("baron",seq(ncol(list_sce_pancreas[[1]])))
colnames(list_sce_pancreas[[2]]) <- paste0("muraro",seq(ncol(list_sce_pancreas[[2]])))
colnames(list_sce_pancreas[[3]]) <- paste0("segerstolpe",seq(ncol(list_sce_pancreas[[3]])))
colnames(list_sce_pancreas[[4]]) <- paste0("baronM",seq(ncol(list_sce_pancreas[[4]])))
list_seurat_pancreas <- pbmcmapply(x = list_sce_pancreas, y = c("baron","muraro","segerstolpe","baronM"),mc.cores = 4, function(x,y) {
    logcounts(x) <- counts(x)
    x <- as.Seurat(x)
    x <- RunSeuratBasics(x)
    x <- RunMCA(x, nmcs =50)
    x@project.name <- y
    return(x)
    })
write_rds(list_seurat_pancreas, path = "data/SeuratPancreas.rds")
rm(list = ls())
# Epithelial --------------------------------------------------------------


haber <- fread("data/haber/haber.txt.gz")
haber <- as.matrix(haber,rownames = "V1")
montoro <- fread("data/montoro/montoro.txt.gz")
montoro <- as.matrix(montoro, rownames = "V1")

plasschaert.mouse <- fread("gzip -dc data/plasschaert.mouse/plasschaert.mouse.txt.gz",skip=2, strip.white = T,blank.lines.skip = F, fill =  TRUE)
plasschaert.mouse_meta <- fread("gzip -dc data/plasschaert.mouse/plasschaert.mouse_meta.txt.gz",skip=8)
plasschaert.mouse_cell_name <- paste0("plasschaert.mouse_",seq(ncol(plasschaert.mouse)-1))
colnames(plasschaert.mouse)<- c("Genes",plasschaert.mouse_cell_name)
plasschaert.mouse_meta[,cells:=plasschaert.mouse_cell_name]
plasschaert.mouse_uninjured <- c("Genes",plasschaert.mouse_meta[clusters_Fig1 != "", cells])
plasschaert.mouse <- plasschaert.mouse[,..plasschaert.mouse_uninjured]
plasschaert.mouse <- as.matrix(data.frame(plasschaert.mouse,row.names = 1))



plasschaert.human <- fread("gzip -dc data/plasschaert.human/plasschaert.human.txt.gz",skip=2, strip.white=T)
plasschaert.human_meta <- fread("gzip -dc data/plasschaert.human/plasschaert.human_meta.txt.gz",skip=5)
plasschaert.human_cell_name <- paste0("plasschaert.human_",seq(ncol(plasschaert.human)-1))
colnames(plasschaert.human)<- c("Genes", plasschaert.human_cell_name)
plasschaert.human_meta[,cells:=plasschaert.human_cell_name]
plasschaert.human <- as.matrix(data.frame(plasschaert.human,row.names = 1))
montoro_cluster        <- stringr::str_remove(colnames(montoro),"^.*_.*_")
haber_cluster        <- stringr::str_remove(colnames(haber),"^.*_.*_")
colnames(haber)      <- paste0("haber","_", seq(length(colnames(haber))))
colnames(montoro)      <- paste0("montoro","_",seq(length(colnames(montoro))))


# ortholog -----------------------------------------------------------------
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
orth.mouse <-  getBM(c("ensembl_gene_id", "mmusculus_homolog_associated_gene_name"),  values=TRUE, mart = human, uniqueRows=TRUE)
ensembl_symbol <- getBM(c("ensembl_gene_id","hgnc_symbol"),  values=TRUE, mart = human, uniqueRows=TRUE)
ortholog_table <- data.table(inner_join(ensembl_symbol,orth.mouse))
ortholog_table <- ortholog_table[(hgnc_symbol %in% rownames(plasschaert.human)) & mmusculus_homolog_associated_gene_name !=""]
ortholog_table <- ortholog_table[mmusculus_homolog_associated_gene_name %in% MgProteinCodingGenes]
unique_hgnc <- names(table(ortholog_table$hgnc_symbol)[table(ortholog_table$hgnc_symbol) ==1])
unique_mmusculus <- names(table(ortholog_table$mmusculus_homolog_associated_gene_name)[table(ortholog_table$mmusculus_homolog_associated_gene_name) ==1])
ortholog_table <- ortholog_table[(hgnc_symbol %in% unique_hgnc) & (mmusculus_homolog_associated_gene_name %in% unique_mmusculus)]
ortholog_vec <- ortholog_table$mmusculus_homolog_associated_gene_name
names(ortholog_vec) <- ortholog_table$hgnc_symbol

# Normalise -------------------------------------------
haber_counts <- haber
montoro_counts <- montoro
plasschaert.mouse_counts <- plasschaert.mouse
plasschaert.human_counts <- plasschaert.human[names(ortholog_vec),]
rownames(plasschaert.human_counts) <- unname(ortholog_vec[rownames(plasschaert.human_counts)])
# Genes intersect -------------------------------------

haber_counts <- as.sparse(haber_counts[rownames(haber_counts) %in% MgProteinCodingGenes,])
montoro_counts <- as.sparse(montoro_counts[rownames(montoro_counts) %in% MgProteinCodingGenes,])
plasschaert.mouse_counts <- as.sparse(plasschaert.mouse_counts[rownames(plasschaert.mouse_counts) %in% MgProteinCodingGenes,])
plasschaert.human_counts <- as.sparse(plasschaert.human_counts)

haber_counts <- haber_counts[nexprs(haber_counts,byrow = T)>=5,]
montoro_counts <- montoro_counts[nexprs(montoro_counts,byrow = T)>=5,]
plasschaert.mouse_counts <- plasschaert.mouse_counts[nexprs(plasschaert.mouse_counts,byrow = T)>=5,]
plasschaert.human_counts <- plasschaert.human_counts[nexprs(plasschaert.human_counts,byrow = T)>=5,]
gc()

#  Get cluster ----------------------------------------

haber_cluster[str_detect(haber_cluster,"Enterocyte")] <- "Enterocyte"
haber_cluster[str_detect(haber_cluster,"TA")] <- "TA"
haber_cluster[str_detect(haber_cluster,"Tuft")] <- "Brush.Tuft"


montoro_cluster[montoro_cluster=="Club"] <- "Secretory"
montoro_cluster[montoro_cluster=="Neuroendocrine"] <- "PNEC"
montoro_cluster[montoro_cluster=="Ionocyte"] <- "Ionocytes"
montoro_cluster[montoro_cluster=="Tuft"] <- "Brush.Tuft"
names(montoro_cluster) <- colnames(montoro)


plasschaert.mouse_cluster       <- plasschaert.mouse_meta[clusters_Fig1!="",clusters_Fig1]
plasschaert.mouse_cluster[plasschaert.mouse_cluster=="Brush"] <- "Brush.Tuft"
plasschaert.mouse_cluster[plasschaert.mouse_cluster=="Cycling Basal (homeostasis)"] <- "Basal"
plasschaert.mouse_cluster[plasschaert.mouse_cluster=="Krt4/13+"] <- "Krt4.13"
names(plasschaert.mouse_cluster)<-  plasschaert.mouse_meta[clusters_Fig1!="",cells]

plasschaert.human_cluster       <- plasschaert.human_meta[,clusters_Fig1]
plasschaert.human_cluster       <- str_replace(plasschaert.human_cluster, "Brush+PNEC", "Brush.PNEC")
names(plasschaert.human_cluster)<- plasschaert.human_meta[,cells]


# Create Seurat -----------------------------------------------------------

haber_seurat <- CreateSeuratObject(haber_counts,min.cells = 5)
montoro_seurat <- CreateSeuratObject(montoro_counts,min.cells = 5)
plasschaert.mouse_seurat <- CreateSeuratObject(plasschaert.mouse_counts,min.cells = 5)
plasschaert.human_seurat <- CreateSeuratObject(plasschaert.human_counts,min.cells = 5)


# Add Celltype  -----------------------------------------------------------

haber_seurat <- AddMetaData(haber_seurat, haber_cluster,col.name = "cell_type1")
montoro_seurat <- AddMetaData(montoro_seurat,montoro_cluster,col.name = "cell_type1")
plasschaert.mouse_seurat <- AddMetaData(plasschaert.mouse_seurat,metadata = plasschaert.mouse_cluster,col.name = "cell_type1")
plasschaert.human_seurat <- AddMetaData(plasschaert.human_seurat,metadata = plasschaert.human_cluster,col.name = "cell_type1")

# Seurat Basics -----------------------------------------------------------
list_seurat_epithelial <- list(plasschaert.mouse_seurat,plasschaert.human_seurat, montoro_seurat, haber_seurat)
list_seurat_epithelial <- lapply(list_seurat_epithelial,  function(x) {
    x <- RunSeuratBasics(x)
    x <- RunMCA(x, nmcs =50)
    return(x)
    })

write_rds(x = list_seurat_epithelial, path = "data/SeuratEpithelial.rds")

# Save --------------------------------------------------------------------


