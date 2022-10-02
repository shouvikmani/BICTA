library(Seurat)
library(scCATCH)
Pancreas_LogAndFiltered=read_delim("Pancreas_Processed_X.txt", delim = "\t", col_names = F)
Pancreas_LogAndFiltered=as.data.frame(Pancreas_LogAndFiltered)
Pancreas_LogAndFiltered=t(Pancreas_LogAndFiltered)
#read the file Pancreas_LogAndFiltered and store it in the variable pancreas
dim(pancreas)

#same for the file PancreasMetadata
pancreas_meta=PancreasMetadata
#colnames(pancreas)=pancreas_meta$Cell.ID
#Read the genelist
genes=Pancreas_Genelist
rownames(Pancreas_LogAndFiltered) = make.names(genes$V1, unique=TRUE)
colnames(Pancreas_LogAndFiltered)=Pancreas_Processed_Barcodes$V1
Pancreas_LogAndFiltered[1:10,1:10]
pancreas=as.matrix(Pancreas_LogAndFiltered)
ctrl <- CreateSeuratObject(counts = Pancreas_LogAndFiltered)
ctrl <- FindVariableFeatures(ctrl, 
                             selection.method = "vst",
                             nfeatures = 2000, 
                             verbose = FALSE)
top2000 <- head(x = VariableFeatures(object = ctrl), 
                n =2000)
write.table(top2000, "top2000_variable_genes.txt", col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE)
all.genes <- rownames(ctrl)
#PCA
new.assay <- ctrl[["RNA"]]
ctrl_filtered <- SetAssayData(ctrl, assay = 'RNA', slot = "scale.data", new.data = pancreas)
ctrl_filtered=ctrl

ctrl_filtered <- RunPCA(ctrl_filtered, features = VariableFeatures(object = ctrl_filtered))
VizDimLoadings(ctrl_filtered, dims = 1:2, reduction = "pca")

#determine number of PCs
ElbowPlot(object = ctrl_filtered, 
          ndims = 50)
#clustering
ctrl_filtered <- FindNeighbors(ctrl_filtered, dims = 1:12)
ctrl_filtered <- FindClusters(ctrl_filtered, resolution = c(0.4,0.6,0.8,1,1.2))
ctrl_filtered <- RunUMAP(ctrl_filtered, dims = 1:12,res=0.4)
ctrl_filtered <- RunTSNE(ctrl_filtered, dims = 1:12)
DimPlot(ctrl_filtered, reduction = "umap")




#cell type annotation using scCATCH
clu_markers <- scCATCH::findmarkergenes(object = ctrl_filtered,
                                        species = 'Human',
                                        cluster = 'All',
                                        match_CellMatch = TRUE,
                                        cancer = NULL,
                                        tissue = c('Pancreas'),
                                        cell_min_pct = 0.25,
                                        logfc = 0.25,

                                                                                pvalue = 0.05)
#final list of cell types using scCATCH
library(scCATCH)
clu_ann <- scCATCH(clu_markers$clu_markers,
                   species = 'Human',
                   tissue = c('Pancreas'))
#run singleR
ref_1=HumanPrimaryCellAtlasData()
pred <- SingleR(test=GetAssayData(ctrl_filtered), ref=ref_1, labels=ref_1$label.main)
table(pred$labels,ctrl_filtered$seurat_clusters)
#assigning cluster identities using scCATCH - only 6 clusters could be labelled
new.cluster.ids <-c("Bone Marrow","Bone Marrow","Beta","Bone Marrow","Epithelial","Bone Marrow","Epithelial","Bone Marrow","Unknown", "Unknown","Unknown","Unknown","Unknown","Unknown","Unknown") 
length(new.cluster.ids)
ctrl_copy=ctrl_filtered 
names(new.cluster.ids) <- levels(ctrl_copy)
ctrl_copy <- RenameIdents(ctrl_copy, new.cluster.ids)
DimPlot(ctrl_copy, reduction = "umap", pt.size = 0.5)

#assigning cluster identities using singleR
#assigning cluster identities using scCATCH - only 6 clusters could be labelled
new.cluster.ids <- c(rep("Unknown",3),"Alpha",rep("Unknown",4),"PP","Beta","Delta","Unknown","Beta",rep("Unknown",2))
length(new.cluster.ids)
ctrl_copy_1=ctrl_filtered 
names(new.cluster.ids) <- levels(ctrl_copy_1)
ctrl_copy_1 <- RenameIdents(ctrl_copy_1, new.cluster.ids)
DimPlot(ctrl_copy_1, reduction = "umap", pt.size = 0.5)
