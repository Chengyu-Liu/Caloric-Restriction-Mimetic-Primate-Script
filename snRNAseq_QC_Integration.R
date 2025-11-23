#Monkey snRNAseq QC & integration-----by LCY
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(patchwork)
library(ggalluvial)
library(ggplot2)
library(ggridges)
library(SoupX)
library(cowplot)
library(reshape)
library(reshape2)
library(ggpubr)
library(ComplexHeatmap)
library(ggrepel)
set.seed(5)
#----------------------------------------------------------------------------------------------#
plan("multicore", workers = 25)
options(future.globals.maxSize = 300000 * 1024^2)#100000MB~=100G

#ts
ts <- "liver"
#Set parameters
mapping_dir <- 'liuchengyu/mapping/'
tissue_dir <- 'liuchengyu/'

sp_info <- read.csv(paste0(tissue_dir,'sample_info.csv'))

sample_list <- sp_info$sample
name_list <- sample_list

##auto_set ###############-
sample_number <- as.numeric(length(sample_list))
res <- rep(0.8,times=sample_number)
qc_total <- paste0(tissue_dir,'qc_total/')
sample_dir <- paste0(tissue_dir,'analysis/')
#sample_dir <- 'F:/scRNA-seq/Data and Consequences/06_human_liver/02_analysis/'
analysis_dir <-sample_dir
qc_dir <- paste0(analysis_dir, 'qc/')
clusters_dir <- paste0(analysis_dir, 'Clusters/')
#Make dirctionary
dir.create(sample_dir)
dir.create(analysis_dir)
dir.create(qc_dir)
dir.create(qc_total)
dir.create(clusters_dir)
dir.create(DEGs_dir)

#01 creat seurat object------------------------------------------------------------
##01.1 load data & rem bkgd---------------------------------
for (i in seq(1:sample_number)) {
  #i=1
  seurat_Obj <- Read10X_h5(paste0(mapping_dir[i],sample_list[i],"/outs/filtered_feature_bc_matrix.h5"), use.names = TRUE)
  seurat_Obj <- CreateSeuratObject(counts = seurat_Obj, project = sample_list[i], min.cells = 3, min.features = 200)
  seurat_Obj <- NormalizeData(seurat_Obj, normalization.method = "LogNormalize")
  seurat_Obj <- FindVariableFeatures(seurat_Obj, selection.method = "vst")
  seurat_Obj <- ScaleData(seurat_Obj)
  seurat_Obj <- RunPCA(seurat_Obj,  npcs = 30, verbose = F)
  seurat_Obj <- FindNeighbors(seurat_Obj, dims = 1:30)
  seurat_Obj <- FindClusters(seurat_Obj, resolution = 2)
  seurat_Obj <- RunUMAP(seurat_Obj, dims = 1:30)
  #rem bkgd
  #data.data <- load10X(paste0(mapping_dir[i],sample_list[i],'/outs/'))
  data.data <- SoupChannel(LayerData(seurat_Obj,assay="RNA",layer="counts"),LayerData(seurat_Obj,assay="RNA",layer="counts"),calcSoupProfile=F)
  matx <- seurat_Obj@meta.data
  data.data = setClusters(data.data, setNames(matx$seurat_clusters, rownames(matx)))
  umap <- seurat_Obj@reductions$umap@cell.embeddings
  data.data <- setDR(data.data,umap)

  soupProf <- data.frame(row.names = rownames(LayerData(seurat_Obj,assay="RNA",layer="counts")),
                         est=Matrix::rowSums(LayerData(seurat_Obj,assay="RNA",layer="counts"))/sum(LayerData(seurat_Obj,assay="RNA",layer="counts")),
                         counts=Matrix::rowSums(LayerData(seurat_Obj,assay="RNA",layer="counts"))   )
  data.data <- setSoupProfile(data.data,soupProf)
  #
  data.data <- autoEstCont(sc=data.data,forceAccept=T)
  tmp <- adjustCounts(data.data)
  tmp <- CreateSeuratObject(counts = tmp, min.cells = 5, min.features = 200, project = sample_list[i])


  tmp$sample <- sp_info[sp_info$read_name==sample_list[i],]$sample
  tmp$age<- sp_info[sp_info$read_name==sample_list[i],]$age_group
  tmp$gender<- sp_info[sp_info$read_name==sample_list[i],]$gender
  tmp$age_number <- sp_info[sp_info$read_name==sample_list[i],]$age_number


  ##monkey MT_Name
  #rownames(tmp)[rownames(tmp)%in%c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")]

  MT<-c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")
  tmp[["percent.MT"]] <- PercentageFeatureSet(tmp, features = MT)
  #tmp[["percent.MT"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
  ggsave(paste0(qc_dir,name_list[i],"_mt_vln.pdf"), plot = p, width = 10, height = 6)
  ggsave(paste0(qc_total,name_list[i],"_mt_vln.png"), plot = p, width = 10, height = 6)

  plot1 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "percent.MT")
  plot2 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(tmp, feature1 = "percent.MT", feature2 = "nFeature_RNA")
  plot4 <- plot1 + plot2+ plot3
  ggsave(paste0(qc_dir,name_list[i],"_mt_Scatter.pdf"),plot = plot4, width = 12, height = 4)
  ggsave(paste0(qc_total,name_list[i],"_mt_Scatter.png"),plot = plot4, width = 12, height = 4)

  #
  tmp[["percent.RBL"]] <- PercentageFeatureSet(tmp, pattern = "^RP[SL]")
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.RBL"), ncol = 3)
  ggsave(paste0(qc_dir,name_list[i],"_rbl_vln.pdf"), plot = p, width = 10, height = 6)
  ggsave(paste0(qc_total,name_list[i],"_rbl_vln.png"), plot = p, width = 10, height = 6)

  plot1 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "percent.RBL")
  plot2 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- plot1 + plot2
  ggsave(paste0(qc_dir,name_list[i],"_rbl_Scatter.pdf"),plot = plot3, width = 8, height = 4)
  ggsave(paste0(qc_total,name_list[i],"_rbl_Scatter.png"),plot = plot3, width = 8, height = 4)

  #
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_before_qc.rds"))
  assign(sample_list[i], tmp)
  print(paste0(sample_list[i],': ',ncol(tmp),' cells'))
}

##01.2 set mito_percent---------------------------------
for (i in seq(1:sample_number)) {
  tmp <- get(sample_list[i])
  tmp <- subset(tmp, subset = percent.MT < 5 & nFeature_RNA > 200 & nFeature_RNA < 6000)
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
  ggsave(paste0(qc_dir,name_list[i],"_qc_vln.pdf"), plot = p, width = 10, height = 6,limitsize = F)
  ggsave(paste0(qc_total,name_list[i],"_qc_vln.png"),plot = p, width = 10, height = 6,limitsize = F)


  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_after_qc.rds"))
  tmp <- RenameCells(tmp, add.cell.id = name_list[i])
  assign(sample_list[i], tmp)
  print(paste0(sample_list[i],': ',ncol(tmp),' cells'))
}
##01.3 normalization --------------------------------
for (i in seq(1:sample_number)) {
  tmp <- get(sample_list[i])
  print(date())
  print(paste0(sample_list[i], ': SCTransform started'))
  tmp <- SCTransform(tmp, vars.to.regress = "percent.MT", verbose = FALSE,do.scale = T)
  print(date())
  print(paste0(sample_list[i], ': SCTransform finished'))
  #
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_SCT.rds"))
  #
  tmp <- RunPCA(tmp, verbose=F)
  #
  png(paste0(qc_dir,name_list[i],"_pca_heatmap.png"), width=1000,height=2000)
  p=DimHeatmap(tmp, dims=1:30, cells=500, balanced=T)
  print(p)
  dev.off()
  #
  png(paste0(qc_dir,name_list[i],"_ElbowPlot.png"), height = 600, width = 700)
  p<- ElbowPlot(tmp, ndims = 30)
  print(p)
  dev.off()
  #
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_bfPCR.rds"))
  print(paste0("--------------", sample_list[i], " completed (", i, "/" , sample_number ,")--------------"))
  tmp <- RunUMAP(tmp, dims = 1:20, verbose=F)
  tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:20)
  tmp <- FindClusters(tmp, res=res[i])
  tmp[["cluster"]] <- Idents(tmp)
  UMAP <- DimPlot(object = tmp, reduction = "umap", label = TRUE)
  ggsave(paste0(qc_dir,name_list[i],"_umap.pdf"), plot = UMAP, width = 8, height = 6,limitsize = F)
  ggsave(paste0(qc_dir,name_list[i],"_umap.png"), plot = UMAP, width = 8, height = 6,limitsize = F)
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_PCR.rds"))
  assign(sample_list[i], tmp)
  print(paste0("--------------", sample_list[i], " completed (", i, "/" , sample_number ,")--------------"))
}
##01.4 Identify doublets -----------------------------------
pK.df <- data.frame(matrix(nrow=0, ncol=2))
colnames(pK.df) <- c("Sample", "Optimal_pK")

for (i in c(1:sample_number)){
  tmp <- readRDS(paste0(qc_dir,name_list[i],"_PCR.rds"))
  #sweep.res.list <- paramSweep(get(sample_list[i]), PCs = 1:20, sct = T, num.cores=16)
  sweep.res.list <- paramSweep(tmp, PCs = 1:20, sct = T, num.cores=16)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- arrange(bcmvn, desc(BCmetric))$pK[1]
  tmp <- data.frame(Sample=name_list[i], Optimal_pK=pK)
  pK.df <- rbind(pK.df, tmp)
  print(bcmvn)
  print(paste0("--------------", sample_list[i], " completed (", i, "/" , sample_number ,")--------------"))
}
write.csv(pK.df,file =paste0(qc_dir,ts,"_Optimal_pK.csv"), sep = "" )

sample_inforation <- data.frame(matrix(nrow=0, ncol=6))
colnames(sample_inforation) <- c("Sample", "Number","Doublet_prop","AvailableCellNumber","Cell_num_pre","Ratio_use")

for (i in seq(1:sample_number)) {
  #tmp <- get(sample_list[i])
  tmp <- readRDS(paste0(qc_dir,name_list[i],"_PCR.rds"))
  pK.use <- as.numeric(as.character(pK.df$Optimal_pK[i]))
  homotypic.prop <- modelHomotypic(tmp@meta.data$cluster)

  cell_num_pre <- as.numeric(length(tmp@meta.data$orig.ident))
  ratio <- cell_num_pre/1000*0.008

  nExp_poi <- round(ratio*length(tmp@meta.data$orig.ident))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder(tmp, PCs = 1:20, pN = 0.25, pK = pK.use, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
  tmp[["doublet"]] <- tmp[[paste("DF.classifications_0.25", pK.use, nExp_poi.adj, sep="_")]]
  prop <- nExp_poi.adj/length(tmp@meta.data$cluster)
  prop.tmp <- data.frame(Sample=name_list[i], Number=nExp_poi.adj, Doublet_prop=prop,
                         AvailableCellNumber=length(tmp@meta.data$doublet[tmp@meta.data$doublet=='Singlet']),
                         Cell_num_pre=cell_num_pre,Ratio_use=ratio)
  #
  sample_inforation <- rbind(sample_inforation, prop.tmp)
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_doublets.rds"))
  assign(sample_list[i], tmp)
  print(paste0("--------------", sample_list[i], " completed (", i, "/" , sample_number ,")--------------"))
}

write.csv(sample_inforation,file =paste0(qc_dir,ts,"_sample_info.csv"), sep = "" )

for (i in seq(1:sample_number)){
  doub <- DimPlot(get(sample_list[i]), group.by='doublet', cols=c('firebrick', 'grey90'))
  pdf(paste0(qc_dir,name_list[i],"_UMAP_doublet.pdf"), height = 6, width = 8)
  print(doub)
  dev.off()
  png(paste0(qc_dir,name_list[i],"_UMAP_doublet.png"), height = 600, width = 800)
  print(doub)
  dev.off()
}
##01.5 filter out doublet -------------------------
for (i in seq(1:sample_number)){
  tmp <- get(sample_list[i])
  tmp <- subset(tmp, doublet=='Singlet')
  tmp <- SCTransform(tmp, verbose = FALSE, do.scale = T)
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_final.rds"))
  print(paste0(i, ' completed', " (", i, "/" , sample_number ,")"))
}




##01.6 Intergration -------------------------
sample_list <- sp_info$sample
sample_number <- length(sample_list)

for (i in c(1:sample_number)) {
  tmp <- readRDS(file = paste0(qc_dir,sample_list[i],"_final.rds"))
  DefaultAssay(tmp) <- 'RNA'
  assign(sample_list[i],tmp)
}

int.list <- c()
for (i in sample_list) {
  tmp <- get(i)
  tmp2 <- list(tmp)
  int.list <- c(int.list,tmp2)
}
saveRDS(int.list, paste0(analysis_dir,ts,'_int.list.rds'))

merged_obj <- merge(int.list[[1]],
                    y = c(unlist(int.list[2:length(int.list)])),
                    project = "Monkey_liver")
DefaultAssay(merged_obj) <- 'RNA'
merged_obj <- NormalizeData(merged_obj,assay='RNA')
merged_obj <- FindVariableFeatures(merged_obj,assay='RNA')
merged_obj <- ScaleData(merged_obj,assay='RNA')
#merged_obj <- SCTransform(merged_obj)
merged_obj <- RunPCA(merged_obj)
merged_obj <- IntegrateLayers(
  object = merged_obj,
  method = CCAIntegration,
  #normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = F
)
merged_obj <- FindNeighbors(merged_obj, dims = 1:30, reduction = "integrated.cca")
merged_obj <- FindClusters(merged_obj, resolution = 2)
merged_obj <- RunUMAP(merged_obj, dims = 1:30, reduction = "integrated.cca")
saveRDS(merged_obj, paste0(analysis_dir,ts,'_int_v5.rds'))


png(paste0(analysis_dir,'pca_heatmap_v5.png'), height=2000, width=1000)
DimHeatmap(merged_obj, dims = 1:50, cells = 3000, balanced = TRUE)
dev.off()

png(paste0(analysis_dir,'pca_elbow_v5.png'), height=500, width=800)
ElbowPlot(merged_obj, ndims=50)
dev.off()

merged_obj <- JoinLayers(merged_obj,assay='RNA')

saveRDS(merged_obj, file = paste0(analysis_dir,ts,"_beforePC_v5.rds"))