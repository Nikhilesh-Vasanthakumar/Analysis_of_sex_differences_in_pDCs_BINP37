
# Understanding sex differences in Immune Response of Plasmacytoid Dendritic Cells (pDCs)

This repository contains the code and analysis for studying the sex differences in the immune response of plasmacytoid dendritic cells (pDCs) under different treatment conditions. The project aims to uncover the molecular characteristics, gene expression patterns, and regulatory mechanisms that drive sex-specific responses in pDCs. The analysis pipeline includes differential gene expression analysis, gene set enrichment analysis (GSEA), and signature score calculation to identify significant expression changes and pathway activations in pDCs across sexes and treatment conditions. The findings contribute to our understanding of immune system biology and have potential implications for developing targeted therapeutic strategies.


## Data Used

The data used for this analysis was collected from the Northern Netherlands population cohort Lifelines and can be accessed from [here](https://eqtlgen.org/sc/datasets/1m-scbloodnl-dataset.html). The dataset includes information from 120 individuals.

For each individual, the following data was collected:

    Unstimulated samples
    3h and 24h samples stimulated with C. albicans
    3h and 24h samples stimulated with M. tuberculosis
    3h and 24h samples stimulated with P. aeruginosa

Two versions of the data are available: v2 and v3. We have used the v2 version of the data for our analysis, which contains 907 genes per cell.

The data is provided in the following files:

    1M_assignments_conditions_expid.tsv: Contains information about the experimental conditions and sample assignments.
    1M_cell_types.tsv: Provides details about the cell types.
    1M_v2_20201029_RNA.mtx: The RNA expression matrix in matrix market format.
    1M_v2_20201029_RNA_colnames.tsv: Contains the column names for the RNA expression matrix.
    1M_v2_20201029_RNA_rownames.tsv: Contains the row names (gene names) for the RNA expression matrix.

Please refer to these files in the Data directory for accessing and working with the provided data.
## Folder Structure

The repository contains the necessary R scripts and files to reproduce the analysis. The folder structure is organized as follows:

```bash
- data/              (contains the dataset files)
- scripts/           (contains R scripts for data processing, analysis, and visualization)
- results/           (stores the generated plots and analysis results)
- README.md          (the main readme file)
- Requirements.R     (script to install required R packages)
```

## Installation

Clone the repository to your local machine using the following command:

```bash
git clone https://github.com/Nikhilesh-Vasanthakumar/BINP37
```
Install R version 4.2.1 on your machine. R can be downloaded from the official R website [here](https://cran.r-project.org/bin/windows/base/old/4.2.1/): 

Open an R environment or an integrated development environment (IDE) that supports R.

Install the required packages by running the provided install_packages.R script. This script will automatically install the necessary packages and their specific versions used in the analysis. Execute the script using the following command in your R environment:

```bash
source("Requirements.R")
```

By following these steps, you will have the necessary packages installed and set up to run the code successfully.

    
## Workflow
The analysis workflow involves several steps to explore and understand the sex-specific differences in pDC expression to different treatments. Here is an overview of the workflow:

1.  **Data Loading and Preprocessing**: The dataset files, including assignments, cell types, RNA expression matrix, column names, and row names, are loaded into the R environment. Quality control and preprocessing steps are performed to filter out low-quality cells and normalize the gene expression data.

2. **Sex Identification**: Sex differences are determined by analyzing the expression patterns of sex marker genes. Heatmaps are used to assign sexes to the samples accurately.

3. **Differential Expression Analysis**: Differential expression analysis is conducted to identify genes that exhibit significant expression differences between males and females at the baseline and different treatment conditions. The Wilcoxon rank-sum test is employed to assess the significance of expression changes.

4. **Gene Set Enrichment Analysis**: Gene Set Enrichment Analysis (GSEA) is performed to identify gene sets that are significantly enriched in different treatment conditions and timepoints. This analysis provides insights into the biological pathways and processes that are differentially regulated between males and females.

5. **Signatue Score Analysis**: For the gene sets enriched in each timepoint and treatment condition, pathway activation differences between sexes are assessed. The average expression per subject is calculated, and Z scores are used to compare the activation differences across different pathways at each timepoint.

## Data Loading and Preprocessing

This section of the script focuses on loading the necessary libraries and performing quality control (QC) measures on the dataset before subsetting the plasmacytoid dendritic cells (pDCs). 

```{r}
# Setting working directory
setwd("SET YOUR DIRECTORY HERE")

# Load the required packages
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(rtracklayer))
suppressMessages(library(methods))
suppressMessages(library(Matrix))
suppressMessages(library(monocle3))

# Load the data

# Set file paths
matrix_path <- "data/1M_v2_20201029_RNA.mtx"
cells_path <- "data/1M_v2_20201029_RNA_colnames.txt"
features_path <- "data/1M_v2_20201029_RNA_rownames.txt"
cell_types_path <- "data/1M_cell_types.tsv"
assignment_conditions <- "data/1M_assignments_conditions_expid.tsv"
gtf_file <- "data/Homo_sapiens.GRCh38.109.gtf"

# Read in counts matrix
counts <- readMM(matrix_path)

# Read in cells file and set column names of counts matrix
cells <- read.table(cells_path, header = FALSE)
colnames(counts) <- cells$V1

# Read in features file and set row names of counts matrix
features <- read.table(features_path, header = FALSE)
rownames(counts) <- features[, 1]

# Read the cell types annotation file
cell_types <- read.table(cell_types_path, header = TRUE, sep = "\t")

# Read the treatment assignment conditions file
treatment <- read.table(assignment_conditions, header = TRUE, sep = "\t")

```
After the data has been loaded we have to create a Seurat object and visualise it to decide the qc parameters for further analysis.

```{r}
blood_healthy <- CreateSeuratObject(counts = counts, project = "blood_healthy", min.cells = 3, min.features = 200,Idents =cell_types$cell_type)
#saveRDS(blood_healthy,file="blood_healthy.rds")
blood_healthy<- readRDS("blood_healthy.rds")

# Match and reorder cell types based on the barcode mapping
cell_types_ordered <- cell_types[match(colnames(blood_healthy), cell_types$barcode), ]

# Set the identity of the cells based on the reordered cell types
Idents(blood_healthy) <- cell_types_ordered$cell_type

# Upload the metadata to the Seurat object
blood_healthy[["barcode"]] <- cell_types_ordered$barcode
blood_healthy[["cell_type"]] <- cell_types_ordered$cell_type

# Add treatment information to the metadata
treatment_ordered <- treatment[match(colnames(blood_healthy), treatment$barcode), ]
blood_healthy[["treatment"]] <- treatment_ordered$timepoint
blood_healthy[["subject"]] <- treatment_ordered$assignment
blood_healthy[["percent.mt"]] <- PercentageFeatureSet(blood_healthy, pattern = "^MT-")

# Visualization of the blood healthy dataset 
metadata <- blood_healthy@meta.data
# Visualize the number of cell counts per cell type via a barplot (Uncomment if needed)
#metadata %>%
#  	ggplot(aes(x=cell_type, fill=cell_type)) +
#  	geom_bar() +
#  	theme_classic() +
#  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
#  	ggtitle("Cell Count per Cell Type")
#CODE FOR FIG 1 IN THE REPORT
#FIG 1A
# Visualize the distribution of genes detected per cell type via a histogram
metadata %>%
  	ggplot(aes(color=cell_type, x=nFeature_RNA, fill=cell_type)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	scale_x_log10() +
  	geom_vline(xintercept = c(200, 2500), linetype = "dashed", color = "red")
ggsave("genes_celltype.png")
#Fig 1B
# Visualize the percentage of mitochondrial genes per subject
metadata %>%
  	ggplot(aes(color=cell_type, x=percent.mt, fill=cell_type)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	scale_x_log10() +
  	geom_vline(xintercept = c(5), linetype = "dashed", color = "red")
ggsave("mt_gene.png")

#FIG 1C
#Idents(blood_healthy) <- "RNA"
#plot_qc <- VlnPlot(blood_healthy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
#ggsave("results/Quality_control.png",plot_qc)

#Subset the data with the following QC conditions
blood_healthy <- subset(blood_healthy, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalize the data
blood_healthy <- NormalizeData(blood_healthy, normalization.method = "LogNormalize", scale.factor = 10000)
#Find variable features
blood_healthy <- FindVariableFeatures(blood_healthy, selection.method = "vst", nfeatures = 2000)

#We will now visualise the variable features in the cells using the and Variable Feature Plot function.

#FIG 1D
#Visualise the top 10 variable features in the cells using VariableFeaturePlot
top10<-head(VariableFeatures(blood_healthy),10)
plot1 <- VariableFeaturePlot(blood_healthy)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
#ggsave("results/VariableFeaturePlot.png",plot2)
```
We will now perform PCA analysis to see the clustering of cells based on cell types.
```
#Scale the data
all.genes <- rownames(blood_healthy)
blood_healthy <- ScaleData(blood_healthy, features = all.genes)
#Perform PCA
blood_healthy <- RunPCA(blood_healthy, features = VariableFeatures(object = blood_healthy))

# Visualize the PCA analysis using PCAPlot
pcaplot <- PCAPlot(blood_healthy, label = TRUE)
pcaplot

# Generate an Elbow plot to identify the optimal number of principal components
ElbowPlot(blood_healthy)

# Save the PCA plot as a PNG file (Uncomment if needed)
#ggsave("PCAPlot_bycell_type.png", pcaplot)

# Visualize the PCA analysis using VizDimLoadings for the first four dimensions
featuresdefpca <- VizDimLoadings(blood_healthy, dims = 1:4, reduction = "pca")
featuresdefpca

# Save the VizDimLoadings plot as a PNG file (Uncomment if needed)
#ggsave("Defining_features.png", featuresdefpca)

# Visualize the combined features for the first 20 dimensions using VizDimLoadings
Combined_features <- VizDimLoadings(blood_healthy, dims = 1:20, reduction = "pca", combine = TRUE)
Combined_features

# Save the combined features plot as a PNG file (Uncomment if needed)
#ggsave("Combined_features_plot.png", Combined_features)
```

Clustering analysis using UMAP
```
# Set the parallel computing strategy to use multiple sessions with 50 workers
#Use multicore instead if running on local system and reduce workers

plan(strategy = "multisession", workers = 50)

# Find the nearest neighbors of cells in the dataset using the first 20 dimensions
blood_healthy <- FindNeighbors(blood_healthy, dims = 1:20)

# Perform clustering on the dataset with a resolution of 0.5
blood_healthy <- FindClusters(blood_healthy, resolution = 0.5, future.seed = TRUE)

# Run UMAP dimensionality reduction on the dataset using the first 20 dimensions, with specific parameters
blood_healthy <- RunUMAP(blood_healthy, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean",
                         reduction.key = "umap", reduction.name = "umap")

#FIG 2 in the report
# Visualize the clustering using DimPlot
DimPlot(blood_healthy, label = TRUE, group.by = "cell_type")
ggsave("results/DimPlot_cell_type.png")

# Visualize the clustering based on treatment groups using DimPlot(Uncomment if necessary)
#DimPlot(blood_healthy, label = TRUE, group.by = "treatment")
#ggsave("DimPlot_treatment.png")

# Visualize the clustering based on both treatment and cell type using DimPlot
#(Uncomment if necessary)
#DimPlot(blood_healthy, label = TRUE, group.by = c("treatment", "cell_type"))
#ggsave("DimPlot_treatment&cell_type.png")
```
Subsetting the pDC cells.
```
#Subset for pDC cells using the subset function
pDC <- subset(blood_healthy, subset = cell_type == "pDC")
```

## Sex Identification

We will now select the X and Y specific genes in the pDC cells to analyse the sex based differential expression.

```{r}
#Select the X and Y specific genes in the pDC cells
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

pDC_meta <- subset(pDC, features = y_genes_present)
pDC_meta <- subset(pDC_meta, subset = nFeature_RNA > 0)
Idents(pDC_meta) <- "subject"


#Find genes in y_genes present whose Mean expression is >0
y_genes_expr <- AverageExpression(pDC_meta, features = y_genes_present,idents = "subject")
y_genes_mean_expr <- rowMeans(y_genes_expr$RNA)
#Get a vector of genes names with mean expression >1
y_genes_mean_expr <- y_genes_mean_expr > 1
#Select only the genes names with True value
y_genes_present_fil <- y_genes_present[y_genes_mean_expr]
#Trying out the filter
qvals <- vector(length = length(y_genes_present_fil))
for (i in seq_along(y_genes_present_fil)) {
  qvals[i] <- quantile(pDC_meta@assays$RNA@data[y_genes_present_fil[i], ], probs = 0.95)
}

Male <-subset(pDC_meta,subset = EIF1AY > qvals[1] | RPS4Y1 >qvals[2])

#Add sex to the meta data
Male <- AddMetaData(Male, metadata = data.frame(sex = rep("male", nrow(Male@meta.data))))
Male@meta.data$sex <- "Male"
pDC[["sex"]]<-NA

#Now mark the subjects as Males in the pDC dataset
#Logic used If the subject values match in Male@meta.data$subject and pDC@meta.data$subject then add male to the corresponding sex column
# Find matching subjects in Male dataset
matching_subjects <- intersect(Male@meta.data$subject, pDC@meta.data$subject)

# Loop through the matching subjects and mark them as males in pDC dataset
for (subject in matching_subjects) {
  pDC@meta.data[pDC@meta.data$subject == subject, "sex"] <- "Male"
}
#Assigning the NA columns as female
pDC@meta.data$sex <- ifelse(is.na(pDC@meta.data$sex), "Female", pDC@meta.data$sex)
#saveRDS(pDC,file="pDC_sex.rds")
#pDC <- readRDS(file = "pDC_sex.rds")
```
## Confirmation of sex assignment and clustering for pDC Dataset

To ensure accurate analysis, we proceeded with confirming the sex assignment and performed clustering on the normalized pDC dataset. This step aimed to visualize the distinct clusters within the pDC population.

```{r}
# Split pDC dataset by sexes
Male <- subset(pDC, subset = sex == "Male")
Female <- subset(pDC, subset = sex == "Female")

# Normalize data for both sexes
Male <- NormalizeData(Male)
Female <- NormalizeData(Female)

# Perform PCA on both sexes
Male <- RunPCA(Male, features = VariableFeatures(object = Male))
Female <- RunPCA(Female, features = VariableFeatures(object = Female))

# Visualize PCA for Male
DimPlot(Male, reduction = "pca", group.by = "treatment")
ggsave("Male_pca_treat.png")

# Visualize PCA for Female
DimPlot(Female, reduction = "pca", group.by = "treatment")
ggsave("Female_pca_treat.png")

# Perform UMAP for Male
plan(strategy = "multisession", workers = 50)
Male <- FindNeighbors(Male, dims = 1:20)
Male <- FindClusters(Male, resolution = 0.5, future.seed = TRUE)
Male <- RunUMAP(Male, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean", reduction.key = "umap", reduction.name = "umap")

# Visualize UMAP for Male
DimPlot(Male, group.by = "treatment", label = TRUE)
ggsave("Male_treat_UMAP.png")

# Perform UMAP for Female
Female <- FindNeighbors(Female, dims = 1:20)
Female <- FindClusters(Female, resolution = 0.5)
Female <- RunUMAP(Female, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean", reduction.key = "umap", reduction.name = "umap")

# Visualize UMAP for Female
DimPlot(Female, group.by = "treatment", label = TRUE)
ggsave("Female_treat_UMAP.png")

# Perform PCA on the combined dataset
pDC <- NormalizeData(pDC)
pDC <- ScaleData(pDC)
pDC <- RunPCA(pDC, features = rownames(pDC))

# Visualize PCA for the combined dataset
DimPlot(pDC, reduction = "pca", group.by = "treatment", label = TRUE)
ggsave("pDC_treat_pca.png")
DimPlot(pDC, reduction = "pca", group.by = c("sex", "treatment"), label = TRUE)
ggsave("pDC_sex&treat_pca.png")

# Perform UMAP on the combined dataset
plan(strategy = "multicore", workers = 4)
pDC <- FindNeighbors(pDC, dims = 1:20)
pDC <- FindClusters(pDC, resolution = 0.5)
pDC <- RunUMAP(pDC, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean", reduction.key = "umap", reduction.name = "umap")

# Visualize UMAP for the combined dataset
#FIG 5
DimPlot(pDC, reduction = "umap", group.by = c("sex", "treatment"), label = TRUE)
ggsave("pDC_sex&treat_UMAP.png")
DimPlot(pDC)
ggsave("pDC_UMAP_ungrouped.png")

# Assign subject identities to pDC dataset
Idents(pDC) <- "subject"

# Visualize the pDC dataset
metadata <- pDC@meta.data

# Create a bar plot of the cells per sex
#FIG 4
metadata %>%
  ggplot(aes(x = sex, fill = sex)) +
  geom_bar() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  ggtitle("pDC Cells per Sex")
ggsave("Cells_per_sex.png")

# Calculate the average expression per subject in the pDC dataset
cluster.average <- AverageExpression(pDC, group.by = "subject", return.seurat = TRUE)

# Extract the numbers from the current idents as new subject identifiers
new_idents <- as.character(gsub("[^0-9]", "", names(Idents(cluster.average))))

# Add the new subject identifiers to the cluster.average object metadata
cluster.average[["subject"]] <- new_idents

# Create a lookup table of subject, sex, and treatment values from the pDC dataset
subject_sex_lookup <- pDC@meta.data[, c("subject", "sex", "treatment")]

# Remove rows with NAs
subject_sex_lookup <- na.omit(subject_sex_lookup)

# Remove duplicates
subject_sex_lookup <- unique(subject_sex_lookup)

# Order the lookup table by subject
subject_sex_lookup <- subject_sex_lookup[order(subject_sex_lookup$subject), ]

# Match the subject values between the cluster.average and subject_sex_lookup datasets
matching_subjects <- intersect(cluster.average@meta.data$subject, subject_sex_lookup$subject)

# Fill in the "sex" column in the cluster.average object with the corresponding values from the pDC dataset
# Find the indices of matching rows in the cluster.average object
matching_rows <- match(matching_subjects, cluster.average@meta.data$subject)

# Extract the corresponding sex and treatment values from the subject_sex_lookup table
sex_values <- subject_sex_lookup[match(matching_subjects, subject_sex_lookup[["subject"]]), "sex"]
treatment_values <- subject_sex_lookup[match(matching_subjects, subject_sex_lookup[["subject"]]), "treatment"]

# Initialize the "sex" column in the cluster.average object with NAs
cluster.average@meta.data$sex <- rep(NA, nrow(cluster.average@meta.data))

# Assign the sex values to matching rows in the "sex" column of cluster.average
cluster.average@meta.data$sex[!is.na(matching_rows)] <- sex_values[!is.na(matching_rows)]

# Normalize the cluster average data
cluster.average <- NormalizeData(cluster.average)

# Scale the data
cluster.average <- ScaleData(cluster.average, features = rownames(cluster.average))

# Find variable features
cluster.average <- FindVariableFeatures(object = cluster.average)

# Perform PCA
cluster.average <- RunPCA(cluster.average, features = VariableFeatures(object = cluster.average))

# Find neighbors
cluster.average <- FindNeighbors(cluster.average, dims = 1:20)

# Find clusters
cluster.average <- FindClusters(cluster.average, resolution = 0.5)

# Perform UMAP
cluster.average <- RunUMAP(cluster.average, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean", reduction.key = "umap", reduction.name = "umap")

# Add subject information to the cluster average object
cluster.average@meta.data$subject <- colnames(cluster.average)

# Generate a heatmap for selected features, grouped by sex
#FIG 3
DoHeatmap(cluster.average, features = c("EIF1AY", "RPS4Y1", "XIST"), group.by = c('sex'), size = 5, draw.lines = FALSE, disp.min = -1, disp.max = 1)
ggsave("pDC_average&sex_heatmap.png")

# Save the cluster.average object as an RDS file (If needed)
#saveRDS(cluster.average, file = "cluster.average_pDC_bloodhealthy.rds")

```
## Differential Expression and GSEA