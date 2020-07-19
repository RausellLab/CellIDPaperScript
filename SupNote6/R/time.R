source("R/utilitary_function.R")
library(bench)
# CellID(G) ---------------------------------------------------------------
ref <- read_rds("data/ref.rds")
query <- read_rds("data/query.rds")
subs <- read_rds("data/subs.rds")



SCN <- mark(
     SCN1 = pred_SCN(subs$subs_1 , ref),
     SCN2 = pred_SCN(subs$subs_2 , ref),
     SCN3 = pred_SCN(subs$subs_3 , ref),
     SCN4 = pred_SCN(subs$subs_4 , ref),
     SCN5 = pred_SCN(subs$subs_5 , ref),
     SCN6 = pred_SCN(subs$subs_6 , ref),
     SCN7 = pred_SCN(subs$subs_7 , ref), 
     check = F, 
     iterations = 1)
write_rds(SCN, path = "data/SCN.rds")


cellidg <- mark(
            cellidg1 = pred_cellidg(subs$subs_1 , ref),
            cellidg2 = pred_cellidg(subs$subs_2 , ref),
            cellidg3 = pred_cellidg(subs$subs_3 , ref),
            cellidg4 = pred_cellidg(subs$subs_4 , ref),
            cellidg5 = pred_cellidg(subs$subs_5 , ref),
            cellidg6 = pred_cellidg(subs$subs_6 , ref),
            cellidg7 = pred_cellidg(subs$subs_7 , ref), 
            check = F, 
            iterations = 1)
write_rds(cellidg, path = "data/cellidg.rds")


cellidc <- mark(cellidc1 = pred_cellidc(subs$subs_1 , ref),
                cellidc2 = pred_cellidc(subs$subs_2 , ref),
                cellidc3 = pred_cellidc(subs$subs_3 , ref),
                cellidc4 = pred_cellidc(subs$subs_4 , ref),
                cellidc5 = pred_cellidc(subs$subs_5 , ref),
                cellidc6 = pred_cellidc(subs$subs_6 , ref),
                cellidc7 = pred_cellidc(subs$subs_7 , ref), 
                check = F, 
                iterations = 1)
write_rds(cellidc, path = "data/cellidc.rds")

scPred <- mark( scPred1 = pred_scPred(subs$subs_1 , ref),
                scPred2 = pred_scPred(subs$subs_2 , ref),
                scPred3 = pred_scPred(subs$subs_3 , ref),
                scPred4 = pred_scPred(subs$subs_4 , ref),
                scPred5 = pred_scPred(subs$subs_5 , ref),
                scPred6 = pred_scPred(subs$subs_6 , ref),
                scPred7 = pred_scPred(subs$subs_7 , ref), 
                check = F, 
                iterations = 1)
write_rds(scPred, path = "data/scPred.rds")


scmap_cluster <- mark(
               scmap_cluster1 = pred_scmap_cluster(subs$subs_1 , ref),
               scmap_cluster2 = pred_scmap_cluster(subs$subs_2 , ref),
               scmap_cluster3 = pred_scmap_cluster(subs$subs_3 , ref),
               scmap_cluster4 = pred_scmap_cluster(subs$subs_4 , ref),
               scmap_cluster5 = pred_scmap_cluster(subs$subs_5 , ref),
               scmap_cluster6 = pred_scmap_cluster(subs$subs_6 , ref),
               scmap_cluster7 = pred_scmap_cluster(subs$subs_7 , ref), 
               check = F, 
               iterations = 1)
write_rds(scmap_cluster, path = "data/scmap_cluster.rds")

scmap_cell <- mark(
  scmap_cell1 = pred_scmap_cell(subs$subs_1 , ref),
  scmap_cell2 = pred_scmap_cell(subs$subs_2 , ref),
  scmap_cell3 = pred_scmap_cell(subs$subs_3 , ref),
  scmap_cell4 = pred_scmap_cell(subs$subs_4 , ref),
  scmap_cell5 = pred_scmap_cell(subs$subs_5 , ref),
  scmap_cell6 = pred_scmap_cell(subs$subs_6 , ref),
  scmap_cell7 = pred_scmap_cell(subs$subs_7 , ref), 
  check = F, 
  iterations = 1)
write_rds(scmap_cell, path = "data/scmap_cell.rds")

Seurat <- mark(
  Seurat1 = pred_Seurat(subs$subs_1 , ref),
  Seurat2 = pred_Seurat(subs$subs_2 , ref),
  Seurat3 = pred_Seurat(subs$subs_3 , ref),
  Seurat4 = pred_Seurat(subs$subs_4 , ref),
  Seurat5 = pred_Seurat(subs$subs_5 , ref),
  Seurat6 = pred_Seurat(subs$subs_6 , ref),
  Seurat7 = pred_Seurat(subs$subs_7 , ref), 
  check = F, 
  iterations = 1)
write_rds(Seurat, path = "data/scmap_cell.rds")

MNN <- mark(
  MNN1 = pred_MNN(subs$subs_1 , ref),
  MNN2 = pred_MNN(subs$subs_2 , ref),
  MNN3 = pred_MNN(subs$subs_3 , ref),
  MNN4 = pred_MNN(subs$subs_4 , ref),
  MNN5 = pred_MNN(subs$subs_5 , ref),
  MNN6 = pred_MNN(subs$subs_6 , ref),
  MNN7 = pred_MNN(subs$subs_7 , ref), 
  check = F, 
  iterations = 1)
write_rds(MNN, path = "data/MNN.rds")

CHETAH <- mark(
  CHETAH = pred_CHETAH(subs$subs_1 , ref),
  CHETAH = pred_CHETAH(subs$subs_2 , ref),
  CHETAH = pred_CHETAH(subs$subs_3 , ref),
  CHETAH = pred_CHETAH(subs$subs_4 , ref),
  CHETAH = pred_CHETAH(subs$subs_5 , ref),
  check = F, 
  iterations = 1)
write_rds(CHETAH, path = "data/CHETAH.rds")

CaSTLe <- mark(
  CaSTLe1 = pred_CaSTLe(subs$subs_1 , ref),
  CaSTLe2 = pred_CaSTLe(subs$subs_2 , ref),
  CaSTLe3 = pred_CaSTLe(subs$subs_3 , ref),
  CaSTLe4 = pred_CaSTLe(subs$subs_4 , ref),
  CaSTLe5 = pred_CaSTLe(subs$subs_5 , ref),
  CaSTLe6 = pred_CaSTLe(subs$subs_6 , ref),
  CaSTLe7 = pred_CaSTLe(subs$subs_7 , ref), 
  check = F, 
  iterations = 1)
write_rds(CaSTLe, path = "data/CaSTLe.rds")

SingleR <- mark(
  SingleR1 = pred_SingleR(subs$subs_1 , ref),
  SingleR2 = pred_SingleR(subs$subs_2 , ref),
  SingleR3 = pred_SingleR(subs$subs_3 , ref),
  SingleR4 = pred_SingleR(subs$subs_4 , ref),
  SingleR5 = pred_SingleR(subs$subs_5 , ref),
  SingleR6 = pred_SingleR(subs$subs_6 , ref),
  SingleR7 = pred_SingleR(subs$subs_7 , ref), 
  check = F, 
  iterations = 1)
write_rds(Seurat, path = "data/SingleR.rds")

scID <- mark(
  scID1 = pred_scID(subs$subs_1 , ref),
  scID2 = pred_scID(subs$subs_2 , ref),
  scID3 = pred_scID(subs$subs_3 , ref),
  scID4 = pred_scID(subs$subs_4 , ref),
  scID5 = pred_scID(subs$subs_5 , ref),
  check = F, 
  iterations = 1)
write_rds(scID, path = "data/scID.rds")

Seurat <- readRDS("data/MNN_query_var.rds")
Seurat <- readRDS("data/CellID_C_query_var.rds")
