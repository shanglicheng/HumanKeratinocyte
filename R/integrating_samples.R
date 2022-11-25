library(Seurat)
load("~/P_2/rmQCP10H.Rd")

table(sampleInf.qc$Patients)
plot(c(1:1261), match(rownames(sampleInf.qc), colnames(reads.qc)))


reads.DCW11 = reads.qc[,grep("DCW11", sampleInf.qc$Patients)]
reads.P001 = reads.qc[,grep("P001", sampleInf.qc$Patients)]
reads.P004 = reads.qc[,grep("P004", sampleInf.qc$Patients)]
reads.P034 = reads.qc[,grep("P034", sampleInf.qc$Patients)]
reads.P10 = reads.qc[,grep("P10", sampleInf.qc$Patients)]
reads.P35 = reads.qc[,grep("P35", sampleInf.qc$Patients)]
reads.P37 = reads.qc[,grep("P37", sampleInf.qc$Patients)]
reads.P38 = reads.qc[,grep("P38", sampleInf.qc$Patients)]
reads.P8 = reads.qc[,grep("P8", sampleInf.qc$Patients)]

myCreateSeuratObjectFromReadMatrix = function(ReadMatrix, projectName) {
  
  ctrl.data = ReadMatrix
  
  ctrl <- CreateSeuratObject(raw.data = ctrl.data, project = projectName, min.cells = 2)
  ctrl@meta.data$stim <- projectName
  #ctrl <- FilterCells(ctrl, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
  ctrl <- NormalizeData(ctrl)
  ctrl <- ScaleData(ctrl, display.progress = F)
  
  return(ctrl)
}

s.DCW11 = myCreateSeuratObjectFromReadMatrix(ReadMatrix = reads.DCW11, projectName = "DCW11")
s.P001 = myCreateSeuratObjectFromReadMatrix(ReadMatrix = reads.P001, projectName = "P001")
s.P004 = myCreateSeuratObjectFromReadMatrix(ReadMatrix = reads.P004, projectName = "P004")
s.P034 = myCreateSeuratObjectFromReadMatrix(ReadMatrix = reads.P034, projectName = "P034")
s.P10 = myCreateSeuratObjectFromReadMatrix(ReadMatrix = reads.P10, projectName = "P10")
s.P35 = myCreateSeuratObjectFromReadMatrix(ReadMatrix = reads.P35, projectName = "P35")
s.P37 = myCreateSeuratObjectFromReadMatrix(ReadMatrix = reads.P37, projectName = "P37")
s.P38 = myCreateSeuratObjectFromReadMatrix(ReadMatrix = reads.P38, projectName = "P38")
s.P8 = myCreateSeuratObjectFromReadMatrix(ReadMatrix = reads.P8, projectName = "P8")

s.DCW11 = FindVariableGenes(s.DCW11, do.plot = T)
s.P001 = FindVariableGenes(s.P001, do.plot = T)
s.P004 = FindVariableGenes(s.P004, do.plot = T)
s.P034 = FindVariableGenes(s.P034, do.plot = T)
s.P10 = FindVariableGenes(s.P10, do.plot = T)
s.P35 = FindVariableGenes(s.P35, do.plot = T)
s.P37 = FindVariableGenes(s.P37, do.plot = T)
s.P38 = FindVariableGenes(s.P38, do.plot = T)
s.P8 = FindVariableGenes(s.P8, do.plot = T)

g.DCW11 = head(rownames(s.DCW11@hvg.info), 1000)
g.P001 = head(rownames(s.P001@hvg.info), 1000)
g.P004 = head(rownames(s.P004@hvg.info), 1000)
g.P034 = head(rownames(s.P034@hvg.info), 1000)
g.P10 = head(rownames(s.P10@hvg.info), 1000)
g.P35 = head(rownames(s.P35@hvg.info), 1000)
g.P37 = head(rownames(s.P37@hvg.info), 1000)
g.P38 = head(rownames(s.P38@hvg.info), 1000)
g.P8 = head(rownames(s.P8@hvg.info), 1000)

genes.use <- unique(c(g.DCW11, g.P001, g.P004, g.P034, g.P10, g.P35, g.P37, g.P38, g.P8))
genes.use <- intersect(genes.use, rownames(s.DCW11@scale.data))
genes.use <- intersect(genes.use, rownames(s.P001@scale.data))
genes.use <- intersect(genes.use, rownames(s.P004@scale.data))
genes.use <- intersect(genes.use, rownames(s.P034@scale.data))
genes.use <- intersect(genes.use, rownames(s.P10@scale.data))
genes.use <- intersect(genes.use, rownames(s.P35@scale.data))
genes.use <- intersect(genes.use, rownames(s.P37@scale.data))
genes.use <- intersect(genes.use, rownames(s.P38@scale.data))
genes.use <- intersect(genes.use, rownames(s.P8@scale.data))


immune.combined <- RunCCA(s.DCW11, s.P001, genes.use = genes.use, num.cc = 30)
immune.combined <- RunCCA(immune.combined, s.P004, genes.use = genes.use, num.cc = 30)
immune.combined <- RunCCA(immune.combined, s.P034, genes.use = genes.use, num.cc = 30)
immune.combined <- RunCCA(immune.combined, s.P10, genes.use = genes.use, num.cc = 30)
immune.combined <- RunCCA(immune.combined, s.P35, genes.use = genes.use, num.cc = 30)
immune.combined <- RunCCA(immune.combined, s.P37, genes.use = genes.use, num.cc = 30)
immune.combined <- RunCCA(immune.combined, s.P38, genes.use = genes.use, num.cc = 30)
immune.combined <- RunCCA(immune.combined, s.P8, genes.use = genes.use, num.cc = 30)

p1 <- DimPlot(object = immune.combined, reduction.use = "cca", group.by = "stim", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "CC1", group.by = "stim", 
              do.return = TRUE)
pdf(file = paste0(outIndex, "visualize results of CCA plot CC1.pdf"), width = 10, height = 4)
plot_grid(p1, p2)
dev.off()


PrintDim(object = immune.combined, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)

pdf(file = paste0(outIndex, "MetageneBicorPlot.pdf"))
p3 <- MetageneBicorPlot(immune.combined, grouping.var = "stim", dims.eval = 1:30, 
                        display.progress = FALSE)
dev.off()

pdf(file = paste0(outIndex, "DimHeatmap.pdf"))
DimHeatmap(object = immune.combined, reduction.type = "cca", cells.use = 500, 
           dim.use = 1:9, do.balanced = TRUE)
dev.off()

immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim", 
                                 dims.align = 1:20)

pdf(file = paste0(outIndex, "aligned CCA.pdf"), width = 10, height = 4)
p1 <- VlnPlot(object = immune.combined, features.plot = "ACC1", group.by = "stim", 
              do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "ACC2", group.by = "stim", 
              do.return = TRUE)
plot_grid(p1, p2)
dev.off()

# t-SNE and Clustering
immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, 
                           do.fast = T, seed.use = 15555)
immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned", 
                                resolution = 0.6, dims.use = 1:20)
immune.combined@meta.data$Type = sampleInf.qc$Type

cluster = GetClusters(immune.combined)

pdf(file = paste0(outIndex, "TSNEPlot.pdf"), width = 10, height = 8)
p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p3 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "Type")
p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2, p3)
dev.off()

# PCA
immune.combined <- RunPCA(immune.combined, seed.use = 15555)

pdf(file = paste0(outIndex, "PCAPlot.pdf"), width = 10, height = 8)
p1 <- PCAPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p3 <- PCAPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "Type")
p2 <- PCAPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2, p3)
dev.off()
