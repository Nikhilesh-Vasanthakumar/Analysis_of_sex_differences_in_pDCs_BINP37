library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rtracklayer)
library(methods)
library(Matrix)
library(monocle3)

#Open the seurat object
disco_blood <- readRDS("disco_blood_v01.rds")
#Calculate the percentage of mitochondrial genes
disco_blood@assays$RNA@data <- disco_blood@assays$RNA@counts
disco_blood[["percent.mt"]] <- PercentageFeatureSet(disco_blood, pattern = "^MT-")
#Filter based on Mt content
disco_blood <- subset(disco_blood, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalize the data
disco_blood <- NormalizeData(disco_blood, normalization.method = "LogNormalize", scale.factor = 10000)
#Find variable features
disco_blood <- FindVariableFeatures(disco_blood, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(disco_blood), 10)
plot1 <- VariableFeaturePlot(disco_blood)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes <- rownames(disco_blood)
#Scale the data
disco_blood <- ScaleData(disco_blood, features = all.genes)
#Run PCA
disco_blood <- RunPCA(disco_blood, features = VariableFeatures(object = disco_blood))
plan(strategy = "multicore", workers = 4)
#Fimd clusters
disco_blood <- FindNeighbors(disco_blood, dims = 1:30)
disco_blood <- FindClusters(disco_blood, resolution = 0.5)
#run UMAP
disco_blood <- RunUMAP(disco_blood, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean",reduction.key = "umap",reduction.name = "umap")
DimPlot(disco_blood, reduction = "umap", label = TRUE, label.size = 3, pt.size = 0.5,group.by = "ct")
ggsave("Disco_DimPlot.png", width = 10, height = 10)
#Save the seurat object
saveRDS(disco_blood, "disco_blood_v02.rds")
disco_blood<- readRDS("disco_blood_v02.rds")
#Subset the data based on the cell type pDC
pDC <- subset(disco_blood, subset = ct == "pDC")
gtf_file <- "Homo_sapiens.GRCh38.109.gtf"
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

pDC_meta <- subset(pDC, features = y_genes_present)
pDC_meta <- subset(pDC_meta, subset = nFeature_RNA > 0)
Idents(pDC_meta) <- "sample"
y_genes_expr <- AverageExpression(pDC_meta, features = y_genes_present)
y_genes_mean_expr <- rowMeans(y_genes_expr$RNA)

y_genes_mean_expr <- y_genes_mean_expr > 0.5

y_genes_present_fil <- y_genes_present[y_genes_mean_expr]

qvals <- vector(length = length(y_genes_present_fil))
for (i in seq_along(y_genes_present_fil)) {
  qvals[i] <- quantile(pDC_meta@assays$RNA@data[y_genes_present_fil[i], ], probs = 0.85)
}
Male <-subset(pDC_meta,subset = USP9Y >qvals[1]|EIF1AY >qvals[2]|KDM5D> qvals[3]|RPS4Y1 > qvals[4]|DDX3Y > qvals[5]|UTY > qvals[6])
Male <-subset(pDC_meta,subset = RPS4Y1 >qvals[1])
#Add sex to the meta data
Male <- AddMetaData(Male, metadata = data.frame(sex = rep("male", nrow(Male@meta.data))))
Male@meta.data$sex <- "Male"
pDC[["sex"]]<-NA

#Now mark the subjects as Males in the pDC dataset
#Logic used If the subject values match in Male@meta.data$subject and pDC@meta.data$subject then add male to the corresponding sex column
# Find matching subjects in Male dataset
matching_subjects <- intersect(Male@meta.data$sample, pDC@meta.data$sample)

# Loop through the matching subjects and mark them as males in pDC dataset
for (sample in matching_subjects) {
  pDC@meta.data[pDC@meta.data$sample == sample, "sex"] <- "Male"
}
#Assigning the NA columns as female
pDC@meta.data$sex <- ifelse(is.na(pDC@meta.data$sex), "Female",
                             ifelse(pDC@meta.data$sex == "", "Female", pDC@meta.data$sex))

pDC <- NormalizeData(pDC)
pDC <-ScaleData(pDC)
pDC <- RunPCA(pDC, features = VariableFeatures(object = pDC))

#Perform UMAP
pDC <- FindNeighbors(pDC, dims = 1:20)
pDC <- FindClusters(pDC, resolution = 0.5)
pDC <- RunUMAP(pDC, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean",reduction.key = "umap",reduction.name = "umap")

DimPlot(pDC,group.by = "disease",label = TRUE)
ggsave("pDC_DimPlot_disease.png")
cluster.average <- AverageExpression(pDC,group.by="sample",return.seurat = TRUE)

# Add the new idents to the cluster.average object meta data
cluster.average[["sample"]] <- colnames(cluster.average)

# Create a lookup table of subject and sex and treatment values from the pDC dataset
subject_sex_lookup <- pDC@meta.data[, c("sample", "sex")]
# Remove rows with NAs
subject_sex_lookup <- na.omit(subject_sex_lookup)
# Remove duplicates
subject_sex_lookup <- unique(subject_sex_lookup)
# Order by subject
#subject_sex_lookup <- subject_sex_lookup[order(subject_sex_lookup$subject),]

# Match the subject values between the two datasets
matching_subjects <- intersect(cluster.average@meta.data$sample, subject_sex_lookup$sample)

# Fill in the "sex" column in the cluster.average object with the corresponding values from the pDC dataset
# Find the indices of matching rows in the cluster.average object
matching_rows <- match(matching_subjects, cluster.average@meta.data$sample)
# Extract the corresponding sex values from the subject_sex_lookup table
sex_values <- subject_sex_lookup[match(matching_subjects, subject_sex_lookup[["sample"]]), "sex"]
#treatment_values <-subject_sex_lookup[match(matching_subjects, subject_sex_lookup[["subject"]]), "treatment"]
# Initialize the "sex" column in the cluster.average object with NAs
cluster.average@meta.data$sex <- rep(NA, nrow(cluster.average@meta.data))
#cluster.average@meta.data$treatment <- rep(NA, nrow(cluster.average@meta.data))
# Assign the sex values to matching rows in the "sex" column of cluster.average
cluster.average@meta.data$sex[!is.na(matching_rows)] <- sex_values[!is.na(matching_rows)]
#Assign the treatment values to matching rows in the "treatment" column of cluster.average
#cluster.average@meta.data$treatment[!is.na(matching_rows)] <- treatment_values[!is.na(matching_rows)]


#Normalize the cluster average
cluster.average <- NormalizeData(cluster.average)
#Scale the data
cluster.average <- ScaleData(cluster.average,features =rownames(cluster.average))
#Find variable features
cluster.average <- FindVariableFeatures(object=cluster.average)
#Run PCA
cluster.average <- RunPCA(cluster.average,features = VariableFeatures(object = cluster.average),npcs = 20)
#Run UMAP
cluster.average <- FindNeighbors(cluster.average,dims=1:20)
#Find clusters
cluster.average <- FindClusters(cluster.average, resolution = 0.5)
#Run UMAP
cluster.average <- RunUMAP(cluster.average, dims = 1:20, n.neighbors = 40, min.dist = 0.25, metric = "euclidean",reduction.key = "umap",reduction.name = "umap")



#Check for IFNA1
Idents(cluster.average) <-sex
#DoHeatmap
DoHeatmap(cluster.average,features =c("RPS4Y1","DDX3Y","UTY","USP9Y","EIF1AY","KDM5D"),group.by = "sex")
ggsave("Cluster_heatmap.png")



