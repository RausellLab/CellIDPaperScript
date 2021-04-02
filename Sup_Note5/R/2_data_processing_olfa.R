source("R/utilitary_function.R")
library(egg)
dir.create("data/")

# Wu olfa
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120199&format=file", destfile = "data/OlfaWu.tar")
untar("data/OlfaWu.tar", exdir = "data/")

# Rusell olfa
download.file("https://raw.githubusercontent.com/rufletch/p63-HBC-diff/master/ref/oeHBCdiff_clusterLabels.txt", destfile = "data/clusterLabels.txt")
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95601/suppl/GSE95601%5FoeHBCdiff%5FCufflinks%5FeSet%5Fcounts%5Ftable%2Etxt%2Egz", destfile = "data/OlfaRusell.txt.gz")
clusFletcher <- tibble::tribble(
  ~cluster, ~cell_type1,
  1, "Resting.Horizontal.Basal.Cells.(HBCs)",
  2, "Immediate Neuronal Precursor 1 (INP1)",
  3, "Globose Basal Cells (GBCs)",
  4, "Mature Sustentacular Cells",
  5, "Transitional HBC 2",
  7, "Immature Sustentacular Cells",
  8, "Transitional HBC 1",
  9, "Immature Olfactory Sensory Neurons (iOSNs)",
  10, "Immediate Neuronal Precursor 3 (INP3)",
  11, "Microvillous Cells type 1",
  12, "Mature Olfactory Sensory Neurons (mOSNs)",
  14, "Immediate Neuronal Precursor 2 (INP2)",
  15, "Microvillous Cells type 2"
)
clusterCor <- read_tsv("data/clusterLabels.txt", col_names = F) %>% set_colnames(c("cell", "cluster"))
metadata <- column_to_rownames(as.data.frame(inner_join(clusFletcher, clusterCor)), var = "cell")

# Rusell to Seurat --------------------------------------------------------

OlfaRusellCounts <- as.sparse(read.table("data/OlfaRusell.txt.gz"))
OlfaRusellCounts <- OlfaRusellCounts[rownames(OlfaRusellCounts) %in% MgProteinCodingGenes, ]
SeuratFletcher <- CreateSeuratObject(OlfaRusellCounts, min.cells = 5, meta.data = metadata)
SeuratFletcher <- RunSeuratBasics(SeuratFletcher)
SeuratFletcher <- RunMCA(SeuratFletcher, nmcs = 50)
SeuratFletcher$cell_type1[is.na(SeuratFletcher$cell_type1)] <- "unknown"
write_rds(SeuratFletcher, "data/SeuratFletcher.rds")

# Wu to Seurat --------------------------------------------------------
Wucounts <- do.call(cbind, lapply(str_subset(list.files("data/", full.names = T), "GSM"), FUN = read.table))
Wucounts <- as.sparse(Wucounts)
Wucounts <- Wucounts[rownames(Wucounts) %in% MgProteinCodingGenes, ]
SeuratWu <- CreateSeuratObject(Wucounts, min.cells = 5)
SeuratWu <- RunSeuratBasics(SeuratWu)
SeuratWu <- RunMCA(SeuratWu, nmcs = 50)
write_rds(SeuratWu, "data/SeuratWu.rds")

# GS extraction -----------------------------------------------------------
Airway <- read_rds("../SupNote3/data/SeuratEpithelial.rds")[[1]]
Intestinal <- read_rds("../SupNote3/data/SeuratEpithelial.rds")[[4]]
GSAirway <- GetGroupGeneSet(RunMCA(Airway, nmcs = 50), group.by = "cell_type1", dims = 1:50, reduction = "mca", n.features = 200)
GSIntestinal <- GetGroupGeneSet(RunMCA(Intestinal, nmcs = 50), group.by = "cell_type1", dims = 1:50, reduction = "mca", n.features = 200)


# Enrichment with GS ------------------------------------------------------
WuHGT <- RunCellHGT(SeuratWu, pathways = c(Airway = GSAirway["Brush.Tuft"], Intestinal = GSIntestinal["Brush.Tuft"]), dims = 1:50, log.trans = T, n.features = 200)
SeuratWu@assays[["CellID"]] <- CreateAssayObject(WuHGT)
FletcherHGT <- RunCellHGT(SeuratFletcher, pathways = c(Airway = GSAirway["Brush.Tuft"], Intestinal = GSIntestinal["Brush.Tuft"]), dims = 1:50, log.trans = T, n.features = 200)
SeuratFletcher@assays[["CellID"]] <- CreateAssayObject(FletcherHGT)


FletcherCell <- colnames(FletcherHGT)[((FletcherHGT[1, ] > 20) & (FletcherHGT[2, ] > 20))]
WuCell <- colnames(WuHGT)[((WuHGT[1, ] > 20) & (WuHGT[2, ] > 20))]

nSCCWu <- length(WuCell)
nSCCFletcher <- length(FletcherCell)

write_rds(SeuratWu, "data/SeuratWu.rds")
write_rds(SeuratFletcher, "data/SeuratFletcher.rds")
