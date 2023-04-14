library(Seurat)
library(Matrix)
library(dplyr)
library(patchwork)
library(rtracklayer)
library(ggplot2)

# Set file paths
matrix_path <- "1M_v2_20201029_RNA.mtx"
cells_path <- "1M_v2_20201029_RNA_colnames.txt"
features_path <- "1M_v2_20201029_RNA_rownames.txt"
cell_types_path <- "1M_cell_types.tsv"
assignment_conditions <- "1M_assignments_conditions_expid.tsv"
#GTF file
gtf_file <- "Homo_sapiens.GRCh38.109.gtf"


# Read in counts matrix
counts <- readMM(matrix_path)

# Read in cells file and set column names of counts matrix
cells <- read.table(cells_path, header = FALSE)
colnames(counts) <- cells$V1

# Read in features file and set row names of counts matrix
features <- read.table(features_path, header = FALSE)
rownames(counts) <- features[,1]

#Read the annotation file
cell_types <- read.table(cell_types_path, header = TRUE, sep = "\t")
treatment <- read.table(assignment_conditions,header =TRUE ,sep ="\t")
# Create Seurat object using subsetted counts matrix
blood_healthy <- CreateSeuratObject(counts = counts, project = "blood_healthy", min.cells = 3, min.features = 200,Idents =cell_types$cell_type)



# Set barcode and cell type metadata using cell_types data frame
cell_types_ordered <- cell_types[match(colnames(blood_healthy), cell_types$barcode), ]

Idents(blood_healthy) <- cell_types_ordered$cell_type

blood_healthy[["barcode"]] <- cell_types_ordered$barcode
#Add cell type to meta data
blood_healthy[["cell_type"]] <- cell_types_ordered$cell_type
#Add treatment to meta data
blood_healthy[["treatment"]] <- treatment$chem
#Create a new column for percentage of mitochondrial genes
blood_healthy[["percent.mt"]] <- PercentageFeatureSet(blood_healthy, pattern = "^MT-")


#Subset the data with the following QC conditions
blood_healthy <- subset(blood_healthy, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#Normalize the data
blood_healthy <- NormalizeData(blood_healthy, normalization.method = "LogNormalize", scale.factor = 10000)
#Find variable features
blood_healthy <- FindVariableFeatures(blood_healthy, selection.method = "vst", nfeatures = 2000)

#Variable feature plot
top10<-head(VariableFeatures(blood_healthy),10)
#Adding small value to all the features to avoid log(0)
#blood_healthy@assays$RNA@data[top10,]<-blood_healthy@assays$RNA@data[top10,]+0.001

plot1<-VariableFeaturePlot(blood_healthy)
#ggsave("VariableFeaturePlot.png",plot1, width = 10, height = 10)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

#plot1+plot2


blood_healthy <- ScaleData(blood_healthy)


blood_healthy <- RunPCA(blood_healthy, features = VariableFeatures(object = blood_healthy))

#ElbowPlot(blood_healthy)
plan(strategy = "multisession", workers = 20)
blood_healthy <- FindNeighbors(blood_healthy, dims = 1:20)

blood_healthy <- FindClusters(blood_healthy, resolution = 0.5)

blood_healthy <- RunUMAP(blood_healthy, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean",reduction.key = "umap",reduction.name = "umap")
#Save the object file for future use
#saveRDS(blood_healthy, "blood_healthy.rds")
DimPlot(blood_healthy,reduction = "umap",group.by = "cell_type",label = TRUE,pt.size = 0.5)

#Idents() <- pbmc$Subject
##pbmc.markers <- FindAllMarkers(pbmc, assay = 'RNA', only.pos = TRUE, min.pct = 0.25, group_by(Subject) %>% top_n(n = #human-systems-immunology 10, wt = avg_log2FC)
#cluster.averages <- AverageExpression(pbmc, assays = "RNA", return.seurat = TRUE)
#DoHeatmap(cluster.averages, features = top10$gene, size = 2, draw.lines = FALSE)


#Select pDC
pDC <- subset(blood_healthy, subset = cell_type == "pDC")
saveRDS(pDC, "pDC.rds")
###Load the object file to save time
pDC<-readRDS("pDC.rds")
#Find markers
gtf <- rtracklayer::import(gtf_file)
gtf_df <- as.data.frame(gtf)

gene_rows <- gtf_df[gtf_df$seqnames %in% c("X", "Y") & gtf_df$type == "gene", ]

gene_info <- data.frame(gene_name = gene_rows$gene_name, chromosome = gene_rows$seqnames)

x_genes <- gene_info[gene_info$chromosome=="X",]
y_genes=gene_info[gene_info$chromosome=="Y",]

#Use x_genes and y_genes to filter the pDC object
gene_names <- rownames(pDC)
x_genes_present <- intersect(x_genes$gene_name, gene_names)
y_genes_present <- intersect(y_genes$gene_name, gene_names)

subset = gene_names %in% c(x_genes_present, y_genes_present)
matching_genes <- gene_names[subset]
#Selecting only the x amd y genes
pDC <- subset(pDC, features = matching_genes)

#Track Individuals with barcodes
pdc_barcode <-pDC@meta.data$barcode
replacement_barcode <-
Idents(pDC) <- pDC@meta.data$barcode

#Normalize
pDC <- NormalizeData(pDC, normalization.method = "LogNormalize", scale.factor = 10000)
#Find variable features
pDC <- FindVariableFeatures(pDC, selection.method = "vst", nfeatures = 2000)
top10<-head(VariableFeatures(pDC),10)
plot1<-VariableFeaturePlot(pDC)
#Visualize High Cell to cell variation
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Scale the data
pDC <- ScaleData(pDC,features = gene_names)

#Run PCA
pDC <- RunPCA(pDC, features = VariableFeatures(object = pDC))
#Cluster the cells
plan(strategy = "multisession", workers = 20)
pDC <- FindNeighbors(pDC, dims = 1:20)
pDC <- FindClusters(pDC, resolution = 0.5)
pDC <-RunUMAP(pDC, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean",reduction.key = "umap",reduction.name = "umap")
DimPlot(pDC,reduction = "umap",group.by="",pt.size = 0.5)





#Average Expression
#cluster.averages <- AverageExpression(pDC, assays = "RNA", return.seurat = TRUE)

#ymarkers=FindAllMarkers(pDC,ident="all",genes.use=y_genes$gene_name)