
# Investigating Sex Differences in the Function of Plasmacytoid Dendritic Cells during Pseudomonas Infection

This repository contains the code and analysis for studying the sex differences in the immune response of plasmacytoid dendritic cells (pDCs) under different treatment conditions. The project aims to uncover the molecular characteristics, gene expression patterns, and regulatory mechanisms that drive sex-specific responses in pDCs. The analysis pipeline includes differential gene expression analysis, gene set enrichment analysis (GSEA), and signature score calculation to identify significant expression changes and pathway activations in pDCs across sexes and treatment conditions. The findings contribute to our understanding of immune system biology and have potential implications for developing targeted therapeutic strategies.

## Data Used

The data used for this analysis was collected from the Northern Netherlands population cohort Lifelines and can be accessed from [here](https://eqtlgen.org/sc/datasets/1m-scbloodnl-dataset.html). The dataset includes information from 120 individuals.

For each individual, the following data was collected:

    Unstimulated samples
    3h and 24h samples stimulated with C. albicans
    3h and 24h samples stimulated with M. tuberculosis
    3h and 24h samples stimulated with P. aeruginosa

Two versions of the data are available: v2 and v3. We have used the v2 version of the data for our analysis, which contains 907 genes per cell.

The data is required for the analysis is provided in the following files:

    1M_assignments_conditions_expid.tsv: Contains information about the experimental conditions and sample assignments.
    1M_cell_types.tsv: Provides details about the cell types.
    1M_v2_20201029_RNA.mtx: The RNA expression matrix in matrix market format.
    1M_v2_20201029_RNA_colnames.tsv: Contains the column names for the RNA expression matrix.
    1M_v2_20201029_RNA_rownames.tsv: Contains the row names (gene names) for the RNA expression matrix.

Please refer to the Data directory for information on accessing and working with the provided data.

The gtf file used to identify the y-chromosome genes for sex identification was taken from [here](https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/)


## Folder Structure

The [Github](https://github.com/Nikhilesh-Vasanthakumar/BINP37) repository contains the necessary R scripts and files to reproduce the analysis. The folder structure is organized as follows:

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
#Create a stacked bar plot of the cells per subject coloured by cell type
# Create a stacked bar plot using ggplot2
library(viridis)
mono <- c("mono 1", "mono 2","mono 4")
tcells <- c("memory CD4T", "naive CD4T", "memory CD8T", "naive CD8T","th1 CD4T","th2 CD4T","reg CD4T" )
nkcells <-c("NK","NKbright","NKdim")

metadata$cell_type1 <- ifelse(metadata$cell_type %in% mono, "Monocytes",
						ifelse(metadata$cell_type %in% tcells, "T cells",
						ifelse(metadata$cell_type %in% nkcells, "NK cells", metadata$cell_type)))
stacked <- as.data.frame(table(metadata$subject, metadata$cell_type1))
colnames(stacked) <- c("subject", "cell_type1", "count")

# Define custom color palette
my_colors <- brewer.pal(10, "Set3")

ggplot(stacked, aes(x = subject, y = count, fill = cell_type1)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Cell frequency per subject")

#Add sex using pDC data

UT_bloodhealthy <- subset(blood_healthy,subset=treatment == "UT")
metadata <- UT_bloodhealthy@meta.data
metadata %>%
  group_by(sex, subject, cell_type) %>%
  summarise(count = n()) %>%
  mutate(pct = count / sum(count)) %>%
		filter(cell_type == "pDC") %>%
  ggplot(aes(x = sex, y = pct, fill = sex)) +
  geom_boxplot() +
  geom_jitter(width = 0.2,aes(colour = sex)) +
  ylab("Proportion of pDCs over all cell types for UT")

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

DimPlot(pDC, reduction = "umap", group.by = c("sex", "treatment"), label = TRUE)
ggsave("pDC_sex&treat_UMAP.png")
DimPlot(pDC)
ggsave("pDC_UMAP_ungrouped.png")

# Assign subject identities to pDC dataset
Idents(pDC) <- "subject"

# Visualize the pDC dataset
metadata <- pDC@meta.data

# Create a bar plot of the cells per sex
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

```{r}

#Now perform DE analysis
library(EnhancedVolcano)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(plotly)
library(stringr)
library(tibble)
library(ggrepel)

# Combine sex and treatment columns
pDC_combined <- pDC
sexsubject <- paste(pDC_combined@meta.data$sex, pDC_combined@meta.data$subject, sep = "_")
sextreatment <- paste(pDC_combined@meta.data$sex, pDC_combined@meta.data$treatment, sep = "_")
pDC_combined[["sex.treatment"]] <- sextreatment

# Split treatment column to create a 'stim' column
pDC_combined[['stim']] <- str_split_fixed(pDC_combined@meta.data$treatment, "h", 2)[,2]
pDC_combined$stim[pDC_combined$stim == ""] <- "UT"

# Visualize data using DimPlot
DimPlot(pDC_combined, group.by = "sex.treatment", label = T, pt.size = 0.5)
ggsave("Dimplot using combined column of sex and treatment.png")
#Add a Disease column to the meta data by splitting the treatment column
Idents(pDC_combined) <- pDC_combined@meta.data$sex.treatment
#Perform differential expression analysis on this combined column for each time point for each disease
unt_diffs <- FindMarkers(pDC_combined, ident.1 = "Male_UT", ident.2 = "Female_UT", min.pct = 0.25, logfc.threshold = 0)

#Volcano plot for Untreated Male vs Female
#FIG 4
EnhancedVolcano(unt_diffs,x='avg_log2FC',y='p_val',lab=rownames(unt_diffs),pCutoff = 0.05,FCcutoff = 0.25,title = "Differentially expressed genes between Males and Females in Untreated Samples")
ggsave("Volcano Plot Male_Female_UT.png")
#Gene Ontology Untreated
#Convert to Ensembl IDs
gene_list <- unt_diffs$avg_log2FC
names(gene_list) <- rownames(unt_diffs)
gene_list <- sort(gene_list,decreasing = TRUE)
#Write gene list to tab-delimited file
ego_UT<-  gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "none")
#Plot GS enriched
#dotplot(ego,x="enrichmentScore",showCategory = 30)
#ggsave("GO_enrichment_UT.png")
#Convert to Entez IDs
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(unt_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
#Perform KEGG enrichment analysis
gene_list <- unt_diffs$avg_log2FC
names(gene_list) <- Ids_entrez
gene_list <- sort(gene_list,decreasing = TRUE)
kegg <- gseKEGG(gene = gene_list, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
#dotplot(kegg,x="enrichmentScore" ,showCategory = 20)
#ggsave("KEGG_enrichment_UT.png")
#Upset plot of enriched KEGG terms
#upsetplot(kegg,y="enrichmentScore",showCategory = 20)

#3H Candida
CA3h_male_diffs <- FindMarkers(pDC_combined, ident.1 = "Male_3hCA", ident.2 = "Male_UT", min.pct = 0.25, logfc.threshold = 0)
colnames(CA3h_male_diffs) <- paste(colnames(CA3h_male_diffs),"Male",sep = "_")
CA3h_female_diffs <- FindMarkers(pDC_combined, ident.1 = "Female_3hCA", ident.2 = "Female_UT", min.pct = 0.25, logfc.threshold = 0)
colnames(CA3h_female_diffs) <- paste(colnames(CA3h_female_diffs),"Female",sep = "_")
#Combine the two dataframes
merged_data <- merge(CA3h_male_diffs, CA3h_female_diffs, by="row.names", all=TRUE)
merged_data$Log2FC <- ifelse(abs(merged_data$avg_log2FC_Male) > 1 & abs(merged_data$avg_log2FC_Female) > 1, ">1",
                            ifelse(abs(merged_data$avg_log2FC_Male) > 1, ">1 in males",
                                   ifelse(abs(merged_data$avg_log2FC_Female) > 1, ">1 in females", "<1")))

# Create the plot
#FIG 5A
plot1_3h <- ggplot(merged_data, aes(x = avg_log2FC_Male, y = avg_log2FC_Female, color = Log2FC)) +
  geom_point() +
  geom_text_repel(data = subset(merged_data, abs(avg_log2FC_Female) > 1),
                  aes(label = Row.names), size = 3, point.padding = unit(0.1, "lines"), color = "black") +
  xlim(-4, 6) +
  ylim(-5, 5) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  ggtitle("DEG after 3h of Candida Treatment") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(values = c(">1" = "red", ">1 in males" = "lightblue", ">1 in females" = "pink", "<1" = "grey"))
ggsave(filename = "CA_3H.png",plot = plot1_3h)


#Male_gene_expr
gene_list_male <-CA3h_male_diffs$avg_log2FC
names(gene_list_male) <- rownames(CA3h_male_diffs)
gene_list_male <- sort(gene_list_male,decreasing = T)
#Female_gene_expr
gene_list_female <-CA3h_female_diffs$avg_log2FC
names(gene_list_female) <- rownames(CA3h_female_diffs)
gene_list_female <- sort(gene_list_female,decreasing = T)
#Perform GO enrichment analysis for male
ego_3hCA_M <- gseGO(gene = gene_list_male , ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
#Plot dotplot
#dotplot(ego,x='enrichmentScore',showCategory=30)
#ggsave("Male_GO_enrichment_CA3h.png")
#upsetplot(ego,y="enrichmentScore",showCategory = 20)

#Perform GO enrichment analysis for female
ego_3hCA_F <- gseGO(gene = gene_list_female , ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
#dotplot(ego,x='enrichmentScore',showCategory=20)
#save the data
#ggsave("Female_GO_enrichment_CA3h.png")
#upsetplot(ego,y="enrichmentScore",showCategory = 20)
#Convert to Entez IDs Male
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(CA3h_male_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_male <-CA3h_male_diffs$avg_log2FC
names(gene_list_male) <- Ids_entrez
gene_list_male <-sort(gene_list_male,decreasing = T)
#Perform KEGG enrichment analysis male
kegg_3hCA_M <- gseKEGG(gene = gene_list_male, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
#dotplot(kegg,x="enrichmentScore",showCategory = 20)
#ggsave("Male_KEGG_enrichment_CA3h.png")
#Convert to Enterez IDs Female
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(CA3h_female_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_female <-CA3h_female_diffs$avg_log2FC
names(gene_list_female) <- Ids_entrez
gene_list_female <-sort(gene_list_female,decreasing = T)
#Perform KEGG enrichment analysis male
kegg_3hCA_F <- gseKEGG(gene = gene_list_female, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
#dotplot(kegg,x="enrichmentScore",showCategory = 20)
#ggsave("Female_KEGG_enrichment_CA3h.png")
#24H Candida
CA24h_male_diffs <- FindMarkers(pDC_combined, ident.1 = "Male_24hCA", ident.2 = "Male_UT", min.pct = 0.25, logfc.threshold = 0)
colnames(CA24h_male_diffs) <- paste(colnames(CA24h_male_diffs),"Male",sep = "_")
CA24h_female_diffs <-FindMarkers(pDC_combined,ident.1 = "Female_24hCA",ident.2 = "Female_UT",min.pct = 0.25,logfc.threshold = 0)
colnames(CA24h_female_diffs) <- paste(colnames(CA24h_female_diffs),"Female",sep = "_")
merged_data <- merge(CA24h_male_diffs, CA24h_female_diffs, by="row.names", all=TRUE)

merged_data$Log2FC <- ifelse(abs(merged_data$avg_log2FC_Male) > 1 & abs(merged_data$avg_log2FC_Female) > 1, ">1",
                            ifelse(abs(merged_data$avg_log2FC_Male) > 1, ">1 in males",
                                   ifelse(abs(merged_data$avg_log2FC_Female) > 1, ">1 in females", "<1")))

#FIG 5B
plot1_24h <-ggplot(merged_data, aes(x = avg_log2FC_Male, y = avg_log2FC_Female, color = Log2FC)) +
  geom_point() +
  geom_text_repel(data = subset(merged_data, abs(avg_log2FC_Female) > 1),
                  aes(label = Row.names), size = 3, point.padding = unit(0.1, "lines"), color = "black") +
  xlim(-4, 6) +
  ylim(-5, 5) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  ggtitle("DEG after 24h of Candida Treatment") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(values = c(">1" = "red", ">1 in males" = "lightblue", ">1 in females" = "pink", "<1" = "grey"))
plot1_24h
ggsave(filename = "CA_24H.png",plot = plot1_24h)

#Male_gene_expr
gene_list_male <-CA24h_male_diffs$avg_log2FC
names(gene_list_male) <- rownames(CA24h_male_diffs)
gene_list_male <- sort(gene_list_male,decreasing = T)
#Female_gene_expr
gene_list_female <-CA24h_female_diffs$avg_log2FC
names(gene_list_female) <- rownames(CA24h_female_diffs)
gene_list_female <- sort(gene_list_female,decreasing = T)
#Perform GO enrichment analysis for male
ego_24hCA_M <- gseGO(gene = gene_list_male , ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
#Plot dotplot
#dotplot(ego,x='enrichmentScore',showCategory=30)
#ggsave("Male_GO_enrichment_CA24h.png")
#Perform GO enrichment analysis for female
ego_24hCA_F <- gseGO(gene = gene_list_female , ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
#dotplot(ego,x='enrichmentScore',showCategory=20)
#ggsave("Female_GO_enrichment_CA24h.png")
#Convert to Entez IDs Male
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(CA24h_male_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_male <-CA24h_male_diffs$avg_log2FC
names(gene_list_male) <- Ids_entrez
gene_list_male <-sort(gene_list_male,decreasing = T)
#Perform KEGG enrichment analysis male
kegg_24hCA_M <- gseKEGG(gene = gene_list_male, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
#dotplot(kegg,x="enrichmentScore",showCategory = 30)
#ggsave("Male_KEGG_enrichment_CA24h.png")
#Convert to Enterez IDs Female
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(CA24h_female_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_female <-CA24h_female_diffs$avg_log2FC
names(gene_list_female) <- Ids_entrez
gene_list_female <-sort(gene_list_female,decreasing = T)
#Perform KEGG enrichment analysis male
kegg_24hCA_F <- gseKEGG(gene = gene_list_female, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
#dotplot(kegg,x="enrichmentScore",showCategory = 30)
#ggsave("Female_KEGG_enrichment_CA24h.png")

#3h Pseudomonas
PA3h_male_diffs <- FindMarkers(pDC_combined, ident.1 = "Male_3hPA", ident.2 = "Male_UT", min.pct = 0.25, logfc.threshold = 0)
colnames(PA3h_male_diffs) <- paste(colnames(PA3h_male_diffs),"Male",sep = "_")
PA3h_female_diffs <- FindMarkers(pDC_combined, ident.1 = "Female_3hPA", ident.2 = "Female_UT", min.pct = 0.25, logfc.threshold = 0)
colnames(PA3h_female_diffs) <- paste(colnames(PA3h_female_diffs),"Female",sep = "_")
#Combine the two dataframes
merged_data <- merge(PA3h_male_diffs, PA3h_female_diffs, by="row.names", all=TRUE)
merged_data$Log2FC <- ifelse(abs(merged_data$avg_log2FC_Male) > 1 & abs(merged_data$avg_log2FC_Female) > 1, ">1",
                            ifelse(abs(merged_data$avg_log2FC_Male) > 1, ">1 in males",
                                   ifelse(abs(merged_data$avg_log2FC_Female) > 1, ">1 in females", "<1")))
#FIG 6A
# Create the plot
plot2_3h <- ggplot(merged_data, aes(x = avg_log2FC_Male, y = avg_log2FC_Female, color = Log2FC)) +
  geom_point() +
  geom_text_repel(data = subset(merged_data, abs(avg_log2FC_Female) > 1),
                  aes(label = Row.names), size = 3, point.padding = unit(0.1, "lines"), color = "black") +
  xlim(-4, 6) +
  ylim(-5, 5) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  ggtitle("DEG after 3h of PA Treatment") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(values = c(">1" = "red", ">1 in males" = "lightblue", ">1 in females" = "pink", "<1" = "grey"))
ggsave("3hPA_DEG.png", plot2_3h, width = 10, height = 10, units = "cm")

#Male_gene_expr
gene_list_male <-PA3h_male_diffs$avg_log2FC
names(gene_list_male) <- rownames(PA3h_male_diffs)
gene_list_male <- sort(gene_list_male,decreasing = T)
#Female_gene_expr
gene_list_female <-PA3h_female_diffs$avg_log2FC
names(gene_list_female) <- rownames(PA3h_female_diffs)
gene_list_female <- sort(gene_list_female,decreasing = T)
#Perform GO enrichment analysis for male
ego_3hPA_M <- gseGO(gene = gene_list_male , ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
#Plot dotplot
#Fig 9A
#dotplot(ego,x='enrichmentScore',showCategory=30)
#ggsave("Male_GO_enrichment_PA3h.png")
#Perform GO enrichment analysis for female
ego_3hPA_F <- gseGO(gene = gene_list_female , ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
#Fig 9B
#dotplot(ego,x='enrichmentScore',showCategory=30)
#ggsave("Female_GO_enrichment_PA3h.png")
#Convert to Entez IDs Male
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(PA3h_male_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_male <-PA3h_male_diffs$avg_log2FC
names(gene_list_male) <- Ids_entrez
gene_list_male <-sort(gene_list_male,decreasing = T)
#Perform KEGG enrichment analysis male
kegg_3hPA_M <- gseKEGG(gene = gene_list_male, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
dotplot(kegg,x="enrichmentScore",showCategory = 30)
#ggsave("Male_KEGG_enrichment_PA3h.png")
#Convert to Enterez IDs Female
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(PA3h_female_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_female <-PA3h_female_diffs$avg_log2FC
names(gene_list_female) <- Ids_entrez
gene_list_female <-sort(gene_list_female,decreasing = T)
#Perform KEGG enrichment analysis male
kegg_3hPA_F <- gseKEGG(gene = gene_list_female, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
#dotplot(kegg,x="enrichmentScore",showCategory = 30)
#ggsave("Female_KEGG_enrichment_PA3h.png")

#24h Pseudomonas
PA24h_male_diffs <- FindMarkers(pDC_combined, ident.1 = "Male_24hPA", ident.2 = "Male_UT", min.pct = 0.25, logfc.threshold = 0)
colnames(PA24h_male_diffs) <- paste(colnames(PA24h_male_diffs),"Male",sep = "_")
PA24h_female_diffs <-FindMarkers(pDC_combined,ident.1 = "Female_24hPA",ident.2 = "Female_UT",min.pct = 0.25,logfc.threshold = 0)
colnames(PA24h_female_diffs) <- paste(colnames(PA24h_female_diffs),"Female",sep = "_")
merged_data <- merge(PA24h_male_diffs, PA24h_female_diffs, by="row.names", all=TRUE)

merged_data$Log2FC <- ifelse(abs(merged_data$avg_log2FC_Male) > 1 & abs(merged_data$avg_log2FC_Female) > 1, ">1",
                            ifelse(abs(merged_data$avg_log2FC_Male) > 1, ">1 in males",
                                   ifelse(abs(merged_data$avg_log2FC_Female) > 1, ">1 in females", "<1")))

#FIG 6B
plot2_24h <-ggplot(merged_data, aes(x = avg_log2FC_Male, y = avg_log2FC_Female, color = Log2FC)) +
  geom_point() +
  geom_text_repel(data = subset(merged_data, abs(avg_log2FC_Female) > 1),
                  aes(label = Row.names), size = 3, point.padding = unit(0.1, "lines"), color = "black") +
  xlim(-4, 6) +
  ylim(-5, 5) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  ggtitle("DEG after 24h of Pseudomonas Treatment") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(values = c(">1" = "red", ">1 in males" = "lightblue", ">1 in females" = "pink", "<1" = "grey"))
plot2_24h
ggsave(filename = "PA_24H.png",plot = plot2_24h)


#Male_gene_expr
gene_list_male <-PA24h_male_diffs$avg_log2FC
names(gene_list_male) <- rownames(PA24h_male_diffs)
gene_list_male <- sort(gene_list_male,decreasing = T)
#Female_gene_expr
gene_list_female <-PA24h_female_diffs$avg_log2FC
names(gene_list_female) <- rownames(PA24h_female_diffs)
gene_list_female <- sort(gene_list_female,decreasing = T)
#Perform GO enrichment analysis for male
ego_24hPA_M <- gseGO(gene = gene_list_male , ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
#Plot dotplot
#dotplot(ego,x='enrichmentScore',showCategory=30)
#ggsave("Male_GO_enrichment_PA24h.png")
#Perform GO enrichment analysis for female
ego_24hPA_F <- gseGO(gene = gene_list_female , ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
#dotplot(ego,x='enrichmentScore',showCategory=30)
#ggsave('Female_GO_enrichment_PA24h.png')
#Convert to Entez IDs Male
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(PA24h_male_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_male <-PA24h_male_diffs$avg_log2FC
names(gene_list_male) <- Ids_entrez
gene_list_male <-sort(gene_list_male,decreasing = T)
#Perform KEGG enrichment analysis male
kegg_24hPA_M <- gseKEGG(gene = gene_list_male, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
#dotplot(kegg,x="enrichmentScore",showCategory = 30)
#ggsave("Male_KEGG_enrichment_PA24h.png")
#Convert to Enterez IDs Female
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(PA24h_female_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_female <-PA24h_female_diffs$avg_log2FC
names(gene_list_female) <- Ids_entrez
gene_list_female <-sort(gene_list_female,decreasing = T)
#Perform KEGG enrichment analysis male
kegg_24hPA_F <- gseKEGG(gene = gene_list_female, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
#dotplot(kegg,x="enrichmentScore",showCategory = 30)
#ggsave("Female_KEGG_enrichment_PA24h.png")
#3H M.Tuberculosis
MTB3h_male_diffs <- FindMarkers(pDC_combined, ident.1 = "Male_3hMTB", ident.2 = "Male_UT", min.pct = 0.25, logfc.threshold = 0)
colnames(MTB3h_male_diffs) <- paste(colnames(MTB3h_male_diffs),"Male",sep = "_")
MTB3h_female_diffs <- FindMarkers(pDC_combined, ident.1 = "Female_3hMTB", ident.2 = "Female_UT", min.pct = 0.25, logfc.threshold = 0)
colnames(MTB3h_female_diffs) <- paste(colnames(MTB3h_female_diffs),"Female",sep = "_")
#Combine the two dataframes
merged_data <- merge(MTB3h_male_diffs, MTB3h_female_diffs, by="row.names", all=TRUE)
merged_data$Log2FC <- ifelse(abs(merged_data$avg_log2FC_Male) > 1 & abs(merged_data$avg_log2FC_Female) > 1, ">1",
                            ifelse(abs(merged_data$avg_log2FC_Male) > 1, ">1 in males",
                                   ifelse(abs(merged_data$avg_log2FC_Female) > 1, ">1 in females", "<1")))
#FIG 7A
plot3_3h <- ggplot(merged_data, aes(x = avg_log2FC_Male, y = avg_log2FC_Female, color = Log2FC)) +
  geom_point() +
  geom_text_repel(data = subset(merged_data, abs(avg_log2FC_Female) > 1),
                  aes(label = Row.names), size = 3, point.padding = unit(0.1, "lines"), color = "black") +
  xlim(-4, 6) +
  ylim(-5, 5) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  ggtitle("DEG after 3h of Mycobacterium Tuberculosis") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(values = c(">1" = "red", ">1 in males" = "lightblue", ">1 in females" = "pink", "<1" = "grey"))
plot3_3h
ggsave(filename = "MTB_3H.png",plot = plot3_3h)

#Male_gene_expr
gene_list_male <-MTB3h_male_diffs$avg_log2FC
names(gene_list_male) <- rownames(MTB3h_male_diffs)
gene_list_male <- sort(gene_list_male,decreasing = T)
#Female_gene_expr
gene_list_female <-MTB3h_female_diffs$avg_log2FC
names(gene_list_female) <- rownames(MTB3h_female_diffs)
gene_list_female <- sort(gene_list_female,decreasing = T)
#Perform GO enrichment analysis for male
ego_3hMTB_M <- gseGO(gene = gene_list_male , ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
#Plot dotplot
#dotplot(ego,x='enrichmentScore',showCategory=30)
#ggsave("Male_GO_enrichment_MTB3h.png")
#Perform GO enrichment analysis for female
ego_3hMTB_F <- gseGO(gene = gene_list_female , ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
#dotplot(ego,x='enrichmentScore',showCategory=30)
#ggsave("Female_GO_enrichment_MTB3h.png")
#Convert to Entez IDs Male
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(MTB3h_male_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_male <-MTB3h_male_diffs$avg_log2FC
names(gene_list_male) <- Ids_entrez
gene_list_male <-sort(gene_list_male,decreasing = T)
#Perform KEGG enrichment analysis male
kegg_3hMTB_M <- gseKEGG(gene = gene_list_male, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
#dotplot(kegg,x="enrichmentScore",showCategory = 30)
#ggsave("Male_KEGG_enrichment_MTB3h.png")
#Convert to Enterez IDs Female
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(MTB3h_female_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_female <-MTB3h_female_diffs$avg_log2FC
names(gene_list_female) <- Ids_entrez
gene_list_female <-sort(gene_list_female,decreasing = T)
#Perform KEGG enrichment analysis male
kegg_3hMTB_F <- gseKEGG(gene = gene_list_female, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
#dotplot(kegg,x="enrichmentScore",showCategory = 20)
#ggsave("Female_KEGG_enrichment_MTB3h.png")
#24H M.Tuberculosis
MTB24h_male_diffs <- FindMarkers(pDC_combined, ident.1 = "Male_24hMTB", ident.2 = "Male_UT", min.pct = 0.25, logfc.threshold = 0)
colnames(MTB24h_male_diffs) <- paste(colnames(MTB24h_male_diffs),"Male",sep = "_")
MTB24h_female_diffs <-FindMarkers(pDC_combined,ident.1 = "Female_24hMTB",ident.2 = "Female_UT",min.pct = 0.25,logfc.threshold = 0)
colnames(MTB24h_female_diffs) <- paste(colnames(MTB24h_female_diffs),"Female",sep = "_")
merged_data <- merge(MTB24h_male_diffs, MTB24h_female_diffs, by="row.names", all=TRUE)
merged_data$Log2FC <- ifelse(abs(merged_data$avg_log2FC_Male) > 1 & abs(merged_data$avg_log2FC_Female) > 1, ">1",
                            ifelse(abs(merged_data$avg_log2FC_Male) > 1, ">1 in males",
                                   ifelse(abs(merged_data$avg_log2FC_Female) > 1, ">1 in females", "<1")))
#FIG 7B
plot3_24h <-ggplot(merged_data, aes(x = avg_log2FC_Male, y = avg_log2FC_Female, color = Log2FC)) +
  geom_point() +
  geom_text_repel(data = subset(merged_data, abs(avg_log2FC_Female) > 1),
                  aes(label = Row.names), size = 3, point.padding = unit(0.1, "lines"), color = "black") +
  xlim(-4, 6) +
  ylim(-5, 5) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  ggtitle("DEG after 24h of Mycobacterium Tuberculosis") +
  theme(text = element_text(size = 8)) +
  scale_color_manual(values = c(">1" = "red", ">1 in males" = "lightblue", ">1 in females" = "pink", "<1" = "grey"))
plot3_24h
ggsave(filename = "MTB_24H.png",plot = plot3_24h)


#Male_gene_expr
gene_list_male <-MTB24h_male_diffs$avg_log2FC
names(gene_list_male) <- rownames(MTB24h_male_diffs)
gene_list_male <- sort(gene_list_male,decreasing = T)
#Female_gene_expr
gene_list_female <-MTB24h_female_diffs$avg_log2FC
names(gene_list_female) <- rownames(MTB24h_female_diffs)
gene_list_female <- sort(gene_list_female,decreasing = T)
#Perform GO enrichment analysis for male
ego_24hMTB_M <- gseGO(gene = gene_list_male , ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
#Plot dotplot
#dotplot(ego,x='enrichmentScore',showCategory=30)
#ggsave("Male_GO_enrichment_MTB24h.png")
#Perform GO enrichment analysis for female
ego_24hMTB_F <- gseGO(gene = gene_list_female , ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
#dotplot(ego,x='enrichmentScore',showCategory=30)
#ggsave("Female_GO_enrichment_MTB24h.png")
#Convert to Entez IDs Male
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(MTB24h_male_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_male <-MTB24h_male_diffs$avg_log2FC
names(gene_list_male) <- Ids_entrez
gene_list_male <-sort(gene_list_male,decreasing = T)
#Perform KEGG enrichment analysis male
kegg_24hMTB_M <- gseKEGG(gene = gene_list_male, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
#dotplot(kegg,x="enrichmentScore",showCategory = 30)
#ggsave("Male_KEGG_enrichment_MTB24h.png")
#Convert to Enterez IDs Female
Ids_entrez <- mapIds(org.Hs.eg.db,keys=rownames(MTB24h_female_diffs),column="ENTREZID", keytype="SYMBOL", multiVals="first")
gene_list_female <-MTB24h_female_diffs$avg_log2FC
names(gene_list_female) <- Ids_entrez
gene_list_female <-sort(gene_list_female,decreasing = T)
#Perform KEGG enrichment analysis male
kegg_24hMTB_F <- gseKEGG(gene = gene_list_female, organism = "hsa", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
#dotplot(kegg,x="enrichmentScore",showCategory = 30)
#ggsave("Female_KEGG_enrichment_MTB24h.png")
```
Now we are comparing between the sexes for each treatment at 3h and 24h
```
#Further subsetting to facilitate comparisions
#Subset the cluster.average to select only the UT individuals
UT <- subset(pDC,subset = treatment == 'UT')
cluster.average_UT <-AverageExpression(UT,group.by = 'subject',return.seurat = T)
#Normalize the data
cluster.average_UT <- NormalizeData(cluster.average_UT)
#Scale the data
cluster.average_UT <- ScaleData(cluster.average_UT)
#Find Variable Genes
cluster.average_UT <- FindVariableFeatures(cluster.average_UT)
#Perform a PCA
cluster.average_UT <- RunPCA(cluster.average_UT,dims=1:20)
#Find Neighbors
cluster.average_UT <- FindNeighbors(cluster.average_UT)
#Find Clusters
cluster.average_UT <- FindClusters(cluster.average_UT,resolution = 0.5)

#Run UMAP
cluster.average_UT <- RunUMAP(cluster.average_UT, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean",reduction.key = "umap",reduction.name = "umap")

subject_sex_lookup <- pDC@meta.data[, c("subject", "sex")]
# Remove rows with NAs
subject_sex_lookup <- na.omit(subject_sex_lookup)
# Remove duplicates
subject_sex_lookup <- unique(subject_sex_lookup)
# Order by subject
subject_sex_lookup <- subject_sex_lookup[order(subject_sex_lookup$subject),]
#Add subject column to cluster.average_UT
cluster.average_UT@meta.data$subject <-names(Idents(cluster.average_UT))
# Match the subject values between the two datasets
matching_subjects <- intersect(cluster.average_UT@meta.data$subject, subject_sex_lookup$subject)

# Fill in the "sex" column in the cluster.average object with the corresponding values from the pDC dataset
# Find the indices of matching rows in the cluster.average object
matching_rows <- match(matching_subjects, cluster.average_UT@meta.data$subject)
# Extract the corresponding sex values from the subject_sex_lookup table
sex_values <- subject_sex_lookup[match(matching_subjects, subject_sex_lookup[["subject"]]), "sex"]
# Initialize the "sex" column in the cluster.average object with NAs
cluster.average_UT@meta.data$sex <- rep(NA, nrow(cluster.average_UT@meta.data))

# Assign the sex values to matching rows in the "sex" column of cluster.average
cluster.average_UT@meta.data$sex[!is.na(matching_rows)] <- sex_values[!is.na(matching_rows)]


#CA
# Subset pDC_CA
pDC_CA <-subset(pDC_combined,subset= treatment == '3hCA')
cluster.average_3hCA <-AverageExpression(pDC_CA,group.by = 'subject' ,return.seurat= T)
#Normalize the data
cluster.average_3hCA <- NormalizeData(cluster.average_3hCA)
#Scale the data
cluster.average_3hCA <- ScaleData(cluster.average_3hCA)
#Find variable genes
cluster.average_3hCA <- FindVariableFeatures(cluster.average_3hCA)
#Run PCA
cluster.average_3hCA <- RunPCA(cluster.average_3hCA,dims=1:20)
#Find Neighbors
cluster.average_3hCA <- FindNeighbors(cluster.average_3hCA)
#Find Clusters
cluster.average_3hCA <- FindClusters(cluster.average_3hCA, resolution = 0.5)
#Run UMAP
cluster.average_3hCA <- RunUMAP(cluster.average_3hCA, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean",reduction.key = "umap",reduction.name = "umap")
#Add the sex and subject information.
subject_sex_lookup <- pDC@meta.data[, c("subject", "sex")]
# Remove rows with NAs
subject_sex_lookup <- na.omit(subject_sex_lookup)
# Remove duplicates
subject_sex_lookup <- unique(subject_sex_lookup)
# Order by subject
subject_sex_lookup <- subject_sex_lookup[order(subject_sex_lookup$subject),]
#Add subject column to cluster.average_3hCA
cluster.average_3hCA@meta.data$subject <-names(Idents(cluster.average_3hCA))
# Match the subject values between the two datasets
matching_subjects <- intersect(cluster.average_3hCA@meta.data$subject, subject_sex_lookup$subject)

# Fill in the "sex" column in the cluster.average object with the corresponding values from the pDC dataset
# Find the indices of matching rows in the cluster.average object
matching_rows <- match(matching_subjects, cluster.average_3hCA@meta.data$subject)
# Extract the corresponding sex values from the subject_sex_lookup table
sex_values <- subject_sex_lookup[match(matching_subjects, subject_sex_lookup[["subject"]]), "sex"]
# Initialize the "sex" column in the cluster.average object with NAs
cluster.average_3hCA@meta.data$sex <- rep(NA, nrow(cluster.average_3hCA@meta.data))

# Assign the sex values to matching rows in the "sex" column of cluster.average
cluster.average_3hCA@meta.data$sex[!is.na(matching_rows)] <- sex_values[!is.na(matching_rows)]

#CA24h
#Subset from pDC_combined
pDC_combined_CA24h <- subset(pDC_combined, subset = treatment=='24hCA')
#Write the same code as the previous block
cluster.average_24hCA <- AverageExpression(pDC_combined_CA24h,group.by = 'subject',return.seurat = T)
#Normalize the data
cluster.average_24hCA <- NormalizeData(cluster.average_24hCA)
#Find variable genes
cluster.average_24hCA <- FindVariableFeatures(cluster.average_24hCA,selection.method = 'vst',nfeatures = 2000)
#Find PCA
cluster.average_24hCA <- RunPCA(cluster.average_24hCA)
#Find Neighbors
cluster.average_24hCA <- FindNeighbors(cluster.average_24hCA, dims = 1:20)
#Find clusters
cluster.average_24hCA <- FindNeighbors(cluster.average_24hCA, dims = 1:20)
#Run UMAP
cluster.average_24hCA <- cluster.average_3hCA <- RunUMAP(cluster.average_3hCA, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean",reduction.key = "umap",reduction.name = "umap")

#Add the sex information and subject information.
#Add subject column to cluster.average_24hCA
cluster.average_24hCA@meta.data$subject <-names(Idents(cluster.average_24hCA))
# Match the subject values between the two datasets
matching_subjects <- intersect(cluster.average_24hCA@meta.data$subject, subject_sex_lookup$subject)

# Fill in the "sex" column in the cluster.average object with the corresponding values from the pDC dataset
# Find the indices of matching rows in the cluster.average object
matching_rows <- match(matching_subjects, cluster.average_24hCA@meta.data$subject)
# Extract the corresponding sex values from the subject_sex_lookup table
sex_values <- subject_sex_lookup[match(matching_subjects, subject_sex_lookup[["subject"]]), "sex"]
# Initialize the "sex" column in the cluster.average object with NAs
cluster.average_24hCA@meta.data$sex <- rep(NA, nrow(cluster.average_24hCA@meta.data))

# Assign the sex values to matching rows in the "sex" column of cluster.average
cluster.average_24hCA@meta.data$sex[!is.na(matching_rows)] <- sex_values[!is.na(matching_rows)]

#Pathways enriched at 3h for Pseudomonas when compared with UT
#Subset the 3hPA data from pDC
pDC_combined_PA3h <- subset(pDC_combined, subset = treatment=='3hPA')


cluster.average_3hPA <- AverageExpression(pDC_combined_PA3h, group.by = 'subject',return.seurat = TRUE)
#Normalize the data
cluster.average_3hPA <- NormalizeData(cluster.average_3hPA)
#Find variable genes
cluster.average_3hPA <- FindVariableFeatures(cluster.average_3hPA, selection.method = "vst", nfeatures = 2000)
#Find PCA
cluster.average_3hPA <- RunPCA(cluster.average_3hPA,features =rownames(cluster.average_3hPA))
#Find Neighbors
cluster.average_3hPA <- FindNeighbors(cluster.average_3hPA, dims = 1:20)
#Find Clust
cluster.average_3hPA <- FindClusters(cluster.average_3hPA, resolution = 0.5)
#Find UMAP
cluster.average_3hPA <- RunUMAP(cluster.average_3hPA, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean",reduction.key = "umap",reduction.name = "umap")
#Add subject column to cluster.average_3hhPA
cluster.average_3hPA@meta.data$subject <-names(Idents(cluster.average_3hPA))
# Match the subject values between the two datasets
matching_subjects <- intersect(cluster.average_3hPA@meta.data$subject, subject_sex_lookup$subject)

# Fill in the "sex" column in the cluster.average object with the corresponding values from the pDC dataset
# Find the indices of matching rows in the cluster.average object
matching_rows <- match(matching_subjects, cluster.average_3hPA@meta.data$subject)
# Extract the corresponding sex values from the subject_sex_lookup table
sex_values <- subject_sex_lookup[match(matching_subjects, subject_sex_lookup[["subject"]]), "sex"]
# Initialize the "sex" column in the cluster.average object with NAs
cluster.average_3hPA@meta.data$sex <- rep(NA, nrow(cluster.average_3hPA@meta.data))

# Assign the sex values to matching rows in the "sex" column of cluster.average
cluster.average_3hPA@meta.data$sex[!is.na(matching_rows)] <- sex_values[!is.na(matching_rows)]

#PA 24h
pDC_combined_PA24h <- subset(pDC_combined, subset = treatment=='24hPA')
cluster.average_24hPA <- AverageExpression(pDC_combined_PA3h, group.by = 'subject',return.seurat = TRUE)
#Normalize the data
cluster.average_24hPA <- NormalizeData(cluster.average_24hPA)
#Find variable genes
cluster.average_24hPA <- FindVariableFeatures(cluster.average_24hPA, selection.method = "vst", nfeatures = 2000)
#Find PCA
cluster.average_24hPA <- RunPCA(cluster.average_24hPA,features =rownames(cluster.average_24hPA))
#Find Neighbors
cluster.average_24hPA <- FindNeighbors(cluster.average_24hPA, dims = 1:20)
#Find Clust
cluster.average_24hPA <- FindClusters(cluster.average_24hPA, resolution = 0.5)
#Find UMAP
cluster.average_24hPA <- RunUMAP(cluster.average_24hPA, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean",reduction.key = "umap",reduction.name = "umap")
#Add subject column to cluster.average_3hhPA
cluster.average_24hPA@meta.data$subject <-names(Idents(cluster.average_24hPA))
# Match the subject values between the two datasets
matching_subjects <- intersect(cluster.average_24hPA@meta.data$subject, subject_sex_lookup$subject)

# Fill in the "sex" column in the cluster.average object with the corresponding values from the pDC dataset
# Find the indices of matching rows in the cluster.average object
matching_rows <- match(matching_subjects, cluster.average_24hPA@meta.data$subject)
# Extract the corresponding sex values from the subject_sex_lookup table
sex_values <- subject_sex_lookup[match(matching_subjects, subject_sex_lookup[["subject"]]), "sex"]
# Initialize the "sex" column in the cluster.average object with NAs
cluster.average_24hPA@meta.data$sex <- rep(NA, nrow(cluster.average_24hPA@meta.data))

# Assign the sex values to matching rows in the "sex" column of cluster.average
cluster.average_24hPA@meta.data$sex[!is.na(matching_rows)] <- sex_values[!is.na(matching_rows)]

#Pathways enriched in MTB compared with UT
#MTB 3h
pDC_combined_MTB3h <- subset(pDC_combined, subset = treatment=='3hMTB')
cluster.average_3hMTB <- AverageExpression(pDC_combined_MTB3h, group.by = 'subject',return.seurat = TRUE)
#Normalize the data
cluster.average_3hMTB <- NormalizeData(cluster.average_3hMTB)
#Find variable genes
cluster.average_3hMTB <- FindVariableFeatures(cluster.average_3hMTB, selection.method = "vst", nfeatures = 2000)
#Find PCA
cluster.average_3hMTB <- RunPCA(cluster.average_3hMTB,features =rownames(cluster.average_3hMTB))
#Find Neighbors
cluster.average_3hMTB <- FindNeighbors(cluster.average_3hMTB, dims = 1:20)
#Find Clust
cluster.average_3hMTB <- FindClusters(cluster.average_3hMTB, resolution = 0.5)
#Find UMAP
cluster.average_3hMTB <- RunUMAP(cluster.average_3hMTB, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean",reduction.key = "umap",reduction.name = "umap")

#Add subject column to cluster.average_3hhPA
cluster.average_3hMTB@meta.data$subject <-names(Idents(cluster.average_3hMTB))
# Match the subject values between the two datasets
matching_subjects <- intersect(cluster.average_3hMTB@meta.data$subject, subject_sex_lookup$subject)

# Fill in the "sex" column in the cluster.average object with the corresponding values from the pDC dataset
# Find the indices of matching rows in the cluster.average object
matching_rows <- match(matching_subjects, cluster.average_3hMTB@meta.data$subject)
# Extract the corresponding sex values from the subject_sex_lookup table
sex_values <- subject_sex_lookup[match(matching_subjects, subject_sex_lookup[["subject"]]), "sex"]
# Initialize the "sex" column in the cluster.average object with NAs
cluster.average_3hMTB@meta.data$sex <- rep(NA, nrow(cluster.average_3hMTB@meta.data))

# Assign the sex values to matching rows in the "sex" column of cluster.average
cluster.average_3hMTB@meta.data$sex[!is.na(matching_rows)] <- sex_values[!is.na(matching_rows)]

#24h MTB

pDC_combined_MTB24h <- subset(pDC_combined, subset = treatment=='24hMTB')

cluster.average_24hMTB <- AverageExpression(pDC_combined_MTB3h, group.by = 'subject',return.seurat = TRUE)
#Normalize the data
cluster.average_24hMTB <- NormalizeData(cluster.average_24hMTB)
#Find variable genes
cluster.average_24hMTB <- FindVariableFeatures(cluster.average_24hMTB, selection.method = "vst", nfeatures = 2000)
#Find PCA
cluster.average_24hMTB <- RunPCA(cluster.average_24hMTB,features =rownames(cluster.average_24hMTB))
#Find Neighbors
cluster.average_24hMTB <- FindNeighbors(cluster.average_24hMTB, dims = 1:20)
#Find Clust
cluster.average_24hMTB <- FindClusters(cluster.average_24hMTB, resolution = 0.5)
#Find UMAP
cluster.average_24hMTB <- RunUMAP(cluster.average_24hMTB, dims = 1:20, n.neighbors = 50, min.dist = 0.25, metric = "euclidean",reduction.key = "umap",reduction.name = "umap")

#Add subject column to cluster.average_3hhPA
cluster.average_24hMTB@meta.data$subject <-names(Idents(cluster.average_24hMTB))
# Match the subject values between the two datasets
matching_subjects <- intersect(cluster.average_24hMTB@meta.data$subject, subject_sex_lookup$subject)

matching_rows <- match(matching_subjects, cluster.average_24hMTB@meta.data$subject)
# Extract the corresponding sex values from the subject_sex_lookup table
sex_values <- subject_sex_lookup[match(matching_subjects, subject_sex_lookup[["subject"]]), "sex"]
# Initialize the "sex" column in the cluster.average object with NAs
cluster.average_24hMTB@meta.data$sex <- rep(NA, nrow(cluster.average_24hMTB@meta.data))

# Assign the sex values to matching rows in the "sex" column of cluster.average
cluster.average_24hMTB@meta.data$sex[!is.na(matching_rows)] <- sex_values[!is.na(matching_rows)]


#Perform gsea on the subset datasets to be used for comparisions
#Comparing the gene expression differences at each time point
Idents(pDC_CA) <- pDC_CA@meta.data$sex
#Finding DE genes between the sexes
CA_3h <- FindMarkers(pDC_CA,ident.1="Male",ident.2="Female",min.pct = 0.25,logfc.threshold = 0)
EnhancedVolcano(CA_3h,lab = rownames(CA_3h),x = 'avg_log2FC',y = 'pct.1',title = '3h CA Males vs Females',pCutoff = 0.5,FCcutoff = 0.25)
gene_list <- CA_3h$avg_log2FC
names(gene_list) <- rownames(CA_3h)
gene_list <- sort(gene_list,decreasing = TRUE)
#Performing GSEA for 3h CA
gsea_3hCA  <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "none")
dotplot(gsea_3hCA,x='enrichmentScore',showCategory=30)
p_threshold <- 0.01
enrichment_score_threshold <- 0.5
#Selecting the common enriched gene sets between the different analyses (3-hour time point):
Common_enriched <- Reduce(intersect, list(ego_3hCA_M$Description,ego_3hCA_F$Description,gsea_3hCA$Description))
#Filter the common enriched gene sets based on significance thresholds:
Common_enriched_sig <- Common_enriched[gsea_3hCA$p.adjust[gsea_3hCA$Description %in% Common_enriched] < p_threshold &
                                      abs(gsea_3hCA$enrichmentScore[gsea_3hCA$Description %in% Common_enriched]) > enrichment_score_threshold]
#Now compare with 24h
Idents(pDC_combined_CA24h) <-pDC_combined_CA24h@meta.data$sex
CA_24h <- FindMarkers(pDC_combined_CA24h,ident.1 = 'Male',ident.2 = 'Female',min.pct = 0.25,logfc.threshold = 0)
EnhancedVolcano(CA_24h,lab = rownames(CA_24h),x = 'avg_log2FC',y = 'pct.1',title = '3h CA Males vs Females',pCutoff = 0.5,FCcutoff = 0.25)
gene_list <- CA_24h$avg_log2FC
names(gene_list) <- rownames(CA_24h)
gene_list <- sort(gene_list,decreasing = TRUE)
#Performing GSEA for 24h CA
gsea_24hCA  <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "none")
#Plotting the results of the GSEA
dotplot(gsea_24hCA,x='enrichmentScore',showCategory=30)
#Selecting the common enriched gene sets between the different analyses
Common_enriched <- intersect(gsea_3hCA$Description, gsea_24hCA$Description)
Common_enriched_sig <- Common_enriched[
  gsea_3hCA$p.adjust[gsea_3hCA$Description %in% Common_enriched] < p_threshold &
  abs(gsea_3hCA$enrichmentScore[gsea_3hCA$Description %in% Common_enriched]) > enrichment_score_threshold
]
#Selecting common differentially expressed genesets  between all the timepoints
Common_enriched_all_CA <- Reduce(intersect, list(ego_UT$Description,ego_3hCA_M$Description,ego_3hCA_F$Description, gsea_3hCA$Description, gsea_24hCA$Description))

#Common_enriched_sig_all <- Common_enriched_all[
#  gsea_3hCA$p.adjust[gsea_3hCA$Description %in% Common_enriched_all] < p_threshold &
#  abs(gsea_3hCA$enrichmentScore[gsea_3hCA$Description %in% Common_enriched_all]) > enrichment_score_threshold &
#  gsea_24hCA$p.adjust[gsea_24hCA$Description %in% Common_enriched_all] < p_threshold &
#  abs(gsea_24hCA$enrichmentScore[gsea_24hCA$Description %in% Common_enriched_all]) > enrichment_score_threshold
#]
```
```{r}
#Comparing for PA
Idents(pDC_combined_PA3h) <- pDC_combined_PA3h@meta.data$sex
#Finding DE genes between the sexes
PA_3h <- FindMarkers(pDC_combined_PA3h,ident.1="Male",ident.2="Female",min.pct = 0.25,logfc.threshold = 0)
EnhancedVolcano(PA_3h,lab = rownames(PA_3h),x = 'avg_log2FC',y = 'pct.1',title = '3h CA Males vs Females',pCutoff = 0.5,FCcutoff = 0.25)
gene_list <- PA_3h$avg_log2FC
names(gene_list) <- rownames(PA_3h)
gene_list <- sort(gene_list,decreasing = TRUE)
#Performing GSEA
gsea_3hPA  <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "none")
#Fig 9C
dotplot(gsea_3hPA,x='enrichmentScore',showCategory=30)
p_threshold <- 0.01
enrichment_score_threshold <- 0.5
#Selecting the common enriched gene sets between the different analyses (3-hour time point):
Common_enriched <- Reduce(intersect, list(ego_3hPA_M$Description,ego_3hPA_F$Description,gsea_3hPA$Description))
Common_enriched_sig <- Common_enriched[gsea_3hPA$p.adjust[gsea_3hPA$Description %in% Common_enriched] < p_threshold &
                                      abs(gsea_3hPA$enrichmentScore[gsea_3hPA$Description %in% Common_enriched]) > enrichment_score_threshold]

Idents(pDC_combined_PA24h) <-pDC_combined_PA24h@meta.data$sex
#Finding DE genes between the sexes
PA_24h <- FindMarkers(pDC_combined_PA24h,ident.1 = 'Male',ident.2 = 'Female',min.pct = 0.25,logfc.threshold = 0)
gene_list <- PA_24h$avg_log2FC
names(gene_list) <- rownames(PA_24h)
gene_list <- sort(gene_list,decreasing = TRUE)
gsea_24hPA  <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "none")
#Fig 9D
dotplot(gsea_24hPA,x='enrichmentScore',showCategory=30)

Common_enriched <- intersect(gsea_24hPA$Description, gsea_24hPA$Description)
#Finding the common enriched gene sets between the different analyses (24-hour time point):
Common_enriched_sig <- Common_enriched[
  gsea_24hPA$p.adjust[gsea_24hPA$Description %in% Common_enriched] < p_threshold &
  abs(gsea_24hPA$enrichmentScore[gsea_24hPA$Description %in% Common_enriched]) > enrichment_score_threshold
]
#Selecting common differentially expressed genesets  between all the timepoints
Common_enriched_all_PA <- Reduce(intersect, list(ego_UT$Description,ego_3hPA_M$Description,ego_3hPA_F$Description, gsea_3hPA$Description, gsea_24hPA$Description))
```
```{r}
#Comparing for MTB
Idents(pDC_combined_MTB3h) <- pDC_combined_MTB3h@meta.data$sex
#Finding DE genes between the sexes
MTB_3h <-FindMarkers(pDC_combined_MTB3h,ident.1 = "Male",ident.2 = "Female",min.pct = 0.25,logfc.threshold = 0)
EnhancedVolcano(MTB_3h,lab = rownames(MTB_3h),x = 'avg_log2FC',y = 'pct.1',title = '3h CA Males vs Females',pCutoff = 0.5,FCcutoff = 0.25)
gene_list <- MTB_3h$avg_log2FC
names(gene_list) <- rownames(MTB_3h)
gene_list <- sort(gene_list,decreasing = TRUE)
gsea_3hMTB  <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "none")
dotplot(gsea_3hMTB,x='enrichmentScore',showCategory=30)
p_threshold <- 0.01
enrichment_score_threshold <- 0.5
#Selecting the common enriched gene sets between the different analyses
Common_enriched <- Reduce(intersect, list(ego_3hMTB_M$Description,ego_3hMTB_F$Description,gsea_3hMTB$Description))
Common_enriched_sig <- Common_enriched[gsea_3hMTB$p.adjust[gsea_3hMTB$Description %in% Common_enriched] < p_threshold &
                                      abs(gsea_3hMTB$enrichmentScore[gsea_3hMTB$Description %in% Common_enriched]) > enrichment_score_threshold]
Idents(pDC_combined_MTB24h) <- pDC_combined_MTB24h@meta.data$sex
#Finding DE genes between the sexes
gsea_24MTB <- FindMarkers(pDC_combined_MTB24h,ident.1 = "Male",ident.2 = "Female",min.pct = 0.25,logfc.threshold = 0)
EnhancedVolcano(gsea_24MTB,lab = rownames(gsea_24MTB),x = 'avg_log2FC',y = 'pct.1',title = '3h CA Males vs Females',pCutoff = 0.5,FCcutoff = 0.25)
gene_list <- gsea_24MTB$avg_log2FC
names(gene_list) <- rownames(gsea_24MTB)
gene_list <- sort(gene_list,decreasing = TRUE)
#Finding enriched gene sets
gsea_24hMTB  <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "none")
#Plotting the enriched gene sets
dotplot(gsea_24hMTB,x='enrichmentScore',showCategory=30)
Common_enriched_3h_24h <- intersect(gsea_3hMTB$Description, gsea_24hMTB$Description)
#Selecting common differentially expressed genesets  between all the timepoints
Common_enriched_all_MTB <- Reduce(intersect, list(ego_UT$Description,ego_3hMTB_M$Description,ego_3hMTB_F$Description, gsea_3hMTB$Description, gsea_24hMTB$Description))
```

## Visualisation of Metabolic and Immune Pathways at baseline
Metabolic pathways play a crucial role in various biological processes and are responsible for essential cellular functions. By analyzing the differential expression of genes within specific metabolic pathways, we can gain insights into potential sex-specific differences in metabolic regulation between the sexes.

```{r}
#Analysing the metabolic pathways between sexe at baseline
my_palette <- paletteer::paletteer_c("grDevices::Blue-Red", n = 20, direction = -1)
#Fig 8A
#KEGG Carbon metabolism
i <- DoHeatmap(cluster.average_UT,features = c("GRSF1","OCIAD1","RAB21","SOCS1","PTGDS","CSF2RB","MARCH1"),group.by = 'sex')+ scale_fill_gradientn(colors = rev(brewer.pal(10, "RdBu")), oob = scales::oob_squish_any, limits = c(-2, 2))

#Fig 8B
#Nitric Oxide metabolic pathway
DoHeatmap(cluster.average_UT,features=c("TSPO","KLF4","HSP90AA1","CD36","IL1B","AIF1","SOD2"),group.by = 'sex',draw.lines = T,)+scale_fill_gradientn(colors = rev(brewer.pal(1000, "RdBu")),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))

#Fig 8C
#Biosynthesis of Amino acids
DoHeatmap(cluster.average_UT,features = c("GRSF1","OCIAD1","RAB21","CSF2RB","MARCH1"),group.by = 'sex',disp.min = -1,disp.max = 1)+ scale_fill_gradientn(colors = rev(brewer.pal(10, "RdBu")), oob = scales::oob_squish_any, limits = c(-2, 2))

#Fig 8D
#IL-6 Production
DoHeatmap(cluster.average_UT,features =c("ARRB2","TNFAIP3","PYCARD","TYROBP","CD36","IL1B","AIF1"),group.by = 'sex',disp.min = -1,disp.max = 1)+scale_fill_gradientn(colors = rev(brewer.pal(1000, "RdBu")),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
```

## Signature Score Calculation and Visualisation for Pseudomonas
The signature score analysis involves calculating scores that represent the activity or enrichment of predefined gene sets or pathways within gene expression datasets. These scores provide insights into the biological processes or molecular mechanisms that are associated with the response to Pseudomonas infection.

The genes used at each timepoint are the differentially expressed found.

### Response to Bacterium
```{r}
#SIGNATURE SCORE ANALYSIS
library(paletteer)
library(RColorBrewer)
library(gridExtra)
library(stringr)
library(tidyr)
library(gridExtra)
source ("scripts/Gene_Formatting.R")
#Analysing the Gene sets for Common Pseudomonas enriched pathways and plotting the signature Scores
#PA
#Subset according to stim
pDC_PA_stim <- subset(pDC_combined, stim == "PA" | stim == "UT")
#Response to Bacterium
genes_UT <-format_gene_list("TNFAIP3/GPX1/GNLY/NR4A1/TSPO/PYCARD/C15orf48/ST13/CD36/IL1B/CCL5/S100A8/SOD2/FOS/S100A9/PPBP/LYZ")

genes_3h<-format_gene_list("CXCL11/CXCL10/TNF/ISG15/CCL2/CMPK2/S100A9/C15orf48/IRF8/SOD2/TNFRSF1B/RGS1/ARID5A/WDR83/LYN/PLAC8/XBP1/FOS/PTGIR/LY96/HLA-E/RNASE6/SPN/NR4A1/TMF1/GSDMD/USP18/ADH5/LYZ/RNF213/FCER1G/CFD/JAK2/MAP2K7/CEBPB/HNRNPA0/IL6R/CARD16/HLA-B/CHMP5/RAB14/HERC6/HLA-A/CD96/CD274/CASP4/MAPK14/IFNAR1/GCH1/PYCARD/ERAP1/HADHB/PPP1R11")

genes_24h <- format_gene_list("PLAC8/IL1B/S100A8/GNLY/RNF213/CEBPB/S100A9/FCER1G")

genes <- unique(c(genes_24h, genes_3h, genes_UT))

heatmap_resp_bac__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))
saveRDS(heatmap_resp_bac__UT, file = "1.rds")
heatmap_resp_bac__3hPA <- DoHeatmap(cluster.average_3hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hPA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))
saveRDS(heatmap_resp_bac__3hPA, file = "2.rds")
heatmap_resp_bac__24hPA <- DoHeatmap(cluster.average_24hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hPA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))
saveRDS(heatmap_resp_bac__24hPA, file = "3.rds")
grid.arrange(heatmap_resp_bac__UT, heatmap_resp_bac__3hPA, heatmap_resp_bac__24hPA, ncol = 3)
#Find Signature Scores
genes_resp_bac <- genes

#Subset the data based on the genes
#Subset the data based on the genes
resp_bacma <- subset(pDC_combined, features = genes_resp_bac)
average_resp_bac <- AverageExpression(pDC_combined,group.by = c('subject','treatment'),slot = "scale.data",features = genes_resp_bac,return.seurat = F)

average_resp_bac_PA <- AverageExpression(pDC_PA_stim,group.by = c('subject','treatment'),slot = "scale.data",features = genes_resp_bac,return.seurat = F)
meta <- pDC_combined@meta.data
#Transpose the data
average_resp_bac <- t(average_resp_bac$RNA)

#Convert the data to dataframe
average_resp_bac <- as.data.frame(average_resp_bac)

#Find the row means
average_resp_bac <- cbind(rownames(average_resp_bac),rowMeans(average_resp_bac))

#Convert to Dataframe
average_resp_bac <- as.data.frame(average_resp_bac)

#Rename  columns into subject and timepoints
average_resp_bac <- separate(average_resp_bac,col = "V1",into =c("subject","timepoints"),sep = "_")

colnames(average_resp_bac) <- c("subject","timepoints","avg_exprerssion")

#Convert to dataframe
average_resp_bac <- as.data.frame(average_resp_bac)

average_resp_bac$avg_exprerssion <- as.numeric(average_resp_bac$avg_exprerssion)

#Find the z-score
average_resp_bac$z_score <- zscore(average_resp_bac$avg_exprerssion)



average_resp_bac$timepoints <-factor(average_resp_bac$timepoints)


average_resp_bac$sex <- meta$sex[match(average_resp_bac$subject,meta$subject)]
ggplot(average_resp_bac, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Response to Bacterium ") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))



#Plotting the signature scores for Pseudomonas only
average_resp_bac_PA <- t(average_resp_bac_PA$RNA)
average_resp_bac_PA <- as.data.frame(average_resp_bac_PA)
average_resp_bac_PA <- cbind(rownames(average_resp_bac_PA),rowMeans(average_resp_bac_PA))
average_resp_bac_PA <- as.data.frame(average_resp_bac_PA)
average_resp_bac_PA <- separate(average_resp_bac_PA,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_resp_bac_PA) <- c("subject","timepoints","avg_exprerssion")
average_resp_bac_PA <- as.data.frame(average_resp_bac_PA)
average_resp_bac_PA$avg_exprerssion <- as.numeric(average_resp_bac_PA$avg_exprerssion)
average_resp_bac_PA$z_score <- zscore(average_resp_bac_PA$avg_exprerssion)
average_resp_bac_PA$timepoints <-factor(average_resp_bac_PA$timepoints,levels = c("UT","3hPA","24hPA"))
average_resp_bac_PA$sex <- meta$sex[match(average_resp_bac_PA$subject,meta$subject)]
#FIg 10B
ggplot(average_resp_bac_PA, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Response to Bacterium in Pseudomonas Treatment") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
  ```
### Inflammatory response
```{r}
#inflammatory response
genes_UT <- format_gene_list("NEAT1/TNFAIP3/GPX1/IFNGR2/CD44/KLF4/PYCARD/TYROBP/LGALS1/TIMP1/CD36/IL1B/CCL5/AIF1/S100A8/FOS/S100A9/PPBP/LYZ")

genes_3h <-format_gene_list("CXCL11/CXCL10/TNF/CTSC/CCL2/S100A9/JUN/TNFRSF1B/AIF1/KLF4/THBS1/SUCNR1/WDR83/PTPRC/P2RX1/DDT/LYN/FOS/IFI35/PTGIR/CD40/PLSCR1/RAC1/PTPN2/LY96/HLA-E/IGFBP4/NKG7/IL10RB/GSDMD/USP18/LILRA5/TRADD/LYZ/MGLL/TOLLIP/ETS1/TLR7/GPSM3/FYN/ANXA1/PSMA6/SCYL1/IFI16/CD47/JAK2/CEBPB/CXCR3/HNRNPA0/IL6R/LIPA/CXCR4/CD96/IL17RA/APOL2/CASP4/MAPK14")

genes_24h <- format_gene_list("CCL4/IL1B/S100A8/IL2RA/CXCR4/CCL5/CEBPB/TNFRSF1B/S100A9/CXCL2/MIF/NEAT1/CST7/TNFRSF4/PTPN2/RB1/CD44/TNFAIP3/PTPRC/NFKBIZ/CXCL1/RAC1/CCL2/FURIN/HLA-E/IL10RB/CTSC/MAPKAPK2/REL/NDUFC2/NINJ1")

genes <- unique(c(genes_24h, genes_3h, genes_UT))

heatmap_infla__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))
saveRDS(heatmap_infla__UT, file = "1.rds")
heatmap_infla__3hPA <- DoHeatmap(cluster.average_3hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hPA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))
saveRDS(heatmap_infla__3hPA, file = "2.rds")
heatmap_infla__24hPA <- DoHeatmap(cluster.average_24hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hPA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))
saveRDS(heatmap_infla__24hPA, file = "3.rds")
grid.arrange(heatmap_infla__UT, heatmap_infla__3hPA, heatmap_infla__24hPA, ncol = 3)

#Find Signature Scores
genes_infla <- genes


#Subset the data based on the genes
inflama <- subset(pDC_combined, features = genes_infla)
average_infla <- AverageExpression(pDC_combined,group.by = c('subject','treatment'),slot = "scale.data",features = genes_infla,return.seurat = F)
#Transpose the data
average_infla <- t(average_infla$RNA)
#Convert the data to dataframe
average_infla <- as.data.frame(average_infla)
#Find the row means
average_infla <- cbind(rownames(average_infla),rowMeans(average_infla))
#Convert to Dataframe
average_infla <- as.data.frame(average_infla)
#Rename  columns into subject and timepoints
average_infla <- separate(average_infla,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_infla) <- c("subject","timepoints","avg_exprerssion")
#Convert to dataframe
average_infla <- as.data.frame(average_infla)
average_infla$avg_exprerssion <- as.numeric(average_infla$avg_exprerssion)
#Find the z-score
average_infla$z_score <- zscore(average_infla$avg_exprerssion)

meta <- pDC_combined@meta.data
average_infla$timepoints <-factor(average_infla$timepoints)
average_infla$sex <- meta$sex[match(average_infla$subject,meta$subject)]

#Remove NA data rows
average_infla <- average_infla[complete.cases(average_infla),]

ggplot(average_infla, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Inflamatory Response") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#For only PA stimulated samples
#Average expression of PA stimulated samples
#Fig 10A
average_infla_PA <- AverageExpression(pDC_PA_stim,group.by = c('subject','treatment'),slot = "scale.data",features = genes_infla,return.seurat = F)
average_infla_PA <- t(average_infla_PA$RNA)
average_infla_PA <- as.data.frame(average_infla_PA)
average_infla_PA <- cbind(rownames(average_infla_PA),rowMeans(average_infla_PA))
average_infla_PA <- as.data.frame(average_infla_PA)
average_infla_PA <- separate(average_infla_PA,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_infla_PA) <- c("subject","timepoints","avg_exprerssion")
average_infla_PA <- as.data.frame(average_infla_PA)
average_infla_PA$avg_exprerssion <- as.numeric(average_infla_PA$avg_exprerssion)
average_infla_PA$z_score <- zscore(average_infla_PA$avg_exprerssion)
average_infla_PA$timepoints <-factor(average_infla_PA$timepoints,levels = c("UT","3hPA","24hPA"))
average_infla_PA$sex <- meta$sex[match(average_infla_PA$subject,meta$subject)]

ggplot(average_infla_PA, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Inflamatory Response for Pseudomonas treatment") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
  ```
### Blood vessel morphogenesis
  ```{r}
  genes_UT<- format_gene_list("ID1/CYBB/MTDH/PRKCB/MYH9/CXCR4/GTF2I/B4GALT1/SPI1/ACTG1/TNFAIP3/GPX1/ANXA2/PTEN/SOX4/NR4A1/KLF4/BAX/ZFP36L1/RGCC/FLNA/PKM/IL1B")

genes_3h <- format_gene_list("CXCL10/TNF/CCL2/EREG/RGCC/JUN/KLF4/THBS1/CITED2/GADD45A/XBP1/CD40/JUNB/NR4A1/PRKCB/PANK2/STAT1/STK4/MECP2/RIN2/WASF2/MMP19/RNF213/GNA13/NAA15/AAMP/ETS1/ANXA1/FMNL3/HIF1AN/TGFBR1/CXCR3/CXCR4/ITGA5/UBP1/MAPK14/ADIPOR2/SEMA4A/HIF1A/ERAP1/PDCL3/KLF2/FOXJ2/CYP1B1/ADAM8/SYK/SRF/PGK1/PRDM1/ACTG1/EGLN1/RAP1A/JAK1/HMOX1/SEC24B/SIRT1/NOTCH1/BTG1/PIK3C2A/WARS2")

genes_24h <- format_gene_list("IL1B/RNF213/CXCR4/SAT1/TGFBR2/CCM2/SOX4/GNA13/JAK1/TNFAIP3/CCL2/ANXA2/BTG1/TNFAIP2/NRP2/PGK1/RAP1A/NINJ1/STK4/ADIPOR2/RBPJ/PPP1R16B/SP100/ITGB1/NFE2L2/HIF1A/SPHK1/SEMA4A/HTATIP2/AIMP1/EMC10/NOTCH2/GTF2I/CD40/ZMIZ1/ZFP36L1")

genes <- unique(c(genes_24h, genes_3h, genes_UT))
heatmap_blood_morph__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
     scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))
saveRDS(heatmap_blood_morph__UT, file = "1.rds")
heatmap_blood_morph__3hPA <- DoHeatmap(cluster.average_3hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
      scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
      guides(fill = FALSE, color = FALSE) + ggtitle("3hPA") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))
saveRDS(heatmap_blood_morph__3hPA, file = "2.rds")
heatmap_blood_morph__24hPA <- DoHeatmap(cluster.average_24hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
      scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
      guides(color = FALSE) + ggtitle("24hPA") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))
saveRDS(heatmap_blood_morph__24hPA, file = "3.rds")
    grid.arrange(heatmap_blood_morph__UT, heatmap_blood_morph__3hPA, heatmap_blood_morph__24hPA, ncol = 3)
#Subset the data based on the genes

genes_blood_morph <- genes

#Subset the data based on the genes
blood_morphma <- subset(pDC_combined, features = genes_blood_morph)
average_blood_morph <- AverageExpression(pDC_combined,group.by = c('subject','treatment'),slot = "scale.data",features = genes_blood_morph,return.seurat = F)

#Transpose the data
average_blood_morph <- t(average_blood_morph$RNA)
#Convert the data to dataframe
average_blood_morph <- as.data.frame(average_blood_morph)
#Find the row means
average_blood_morph <- cbind(rownames(average_blood_morph),rowMeans(average_blood_morph))
#Convert to Dataframe
average_blood_morph <- as.data.frame(average_blood_morph)
#Rename  columns into subject and timepoints
average_blood_morph <- separate(average_blood_morph,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_blood_morph) <- c("subject","timepoints","avg_exprerssion")
#Convert to dataframe
average_blood_morph <- as.data.frame(average_blood_morph)
average_blood_morph$avg_exprerssion <- as.numeric(average_blood_morph$avg_exprerssion)
#Find the z-score
average_blood_morph$z_score <- zscore(average_blood_morph$avg_exprerssion)

meta <- pDC_combined@meta.data
average_blood_morph$timepoints <-factor(average_blood_morph$timepoints)
average_blood_morph$sex <- meta$sex[match(average_blood_morph$subject,meta$subject)]

#Remove NA data rows
average_blood_morph <- average_blood_morph[complete.cases(average_blood_morph),]

ggplot(average_blood_morph, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Blood vessel Morphogenesis") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


#For PA treated samples
#Average expression of PA stimulated samples
average_blood_morph_PA <- AverageExpression(pDC_PA_stim,group.by = c('subject','treatment'),slot = "scale.data",features = genes_blood_morph,return.seurat = F)
average_blood_morph_PA <- t(average_blood_morph_PA$RNA)
average_blood_morph_PA <- as.data.frame(average_blood_morph_PA)
average_blood_morph_PA <- cbind(rownames(average_blood_morph_PA),rowMeans(average_blood_morph_PA))
average_blood_morph_PA <- as.data.frame(average_blood_morph_PA)
average_blood_morph_PA <- separate(average_blood_morph_PA,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_blood_morph_PA) <- c("subject","timepoints","avg_exprerssion")
average_blood_morph_PA <- as.data.frame(average_blood_morph_PA)
average_blood_morph_PA$avg_exprerssion <- as.numeric(average_blood_morph_PA$avg_exprerssion)
average_blood_morph_PA$z_score <- zscore(average_blood_morph_PA$avg_exprerssion)
average_blood_morph_PA$timepoints <-factor(average_blood_morph_PA$timepoints,levels = c("UT","3hPA","24hPA"))
average_blood_morph_PA$sex <- meta$sex[match(average_blood_morph_PA$subject,meta$subject)]
ggplot(average_blood_morph_PA, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Blood vessel Morphogenesis for Pseudomonas Treatment") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
```
### Immune Response
```{r}
genes_UT <- format_gene_list("ISG15/PTGDS/IFIT2/IFI44L/MX1/CTSC/SEMA7A/CSF2RB/NFKBID/ISG20/GZMB/NINJ1/TNFRSF4/CCL3/HLA-DQA1/PLSCR1/ADAR/TAP1/SOCS1/IRF1/CD47/HLA-F/IL6ST/LTB/RBCK1/YTHDF2/GPR183/HLA-A/HLA-DQB1/RFTN1/LPXN/BIRC3/ICOSLG/XRCC5/TNFSF13B/REL/SLAMF7/PSMA1/CDC42/HLA-B/STAT6/PPP1R14B/TRIM22/MAD2L2/NFKBIZ/CLEC4C/DGKZ/MX2/STAT1/NR4A3/ERCC1/MPEG1/ACTR3/BIRC2/PSMB10/LITAF/TNIP1/LILRB4/SLC7A5/TRAF4/SP100/TAPBP/NFKB2/IFI16/PTPN2")

genes_3h <- format_gene_list("CXCL11/CXCL10/IFIT2/TNF/CTSC/TNFSF10/NMB/RSAD2/ISG15/IRF1/LST1/MX1/EXOSC9/PARP14/CCL2/IFI6/IFIH1/IFIT3/S100A9/IRF8/EREG/IFI44L/GZMB/HLA-DQA2/PHF14/GBP1/IRF7/ISG20/KDM5D/RGCC/TNFRSF1B/AIF1/OASL/USP15/RGS1/STX7/LTB/THBS1/ARID5A/HSPA1B/CTSL/STOML2/GAPT/SUCNR1/IFITM3/STAT2/PTPRC/PARP9/LYN/FKBP1A/XBP1/LSM14A/IFI35/CD40/EIF2AK2/ST3GAL1/PLSCR1/PTPN2/KLRB1/TAP1/PARP1/LY96/AZI2/CD164/HLA-E/RNASE6/APBB1IP/LCP2/JUNB/DENND1B/ZBTB1/SPN/PHPT1/NKG7/IL10RB/PRKCB/TNFRSF21/GSDMD/USP18/STAT1/SOCS1/LILRA5/LYZ/HLA-C/DOCK10/TOLLIP/PTPN1/FCER1G/ETS1/TLR7/METTL3/FYN/ANXA1/PTPN11/KAT5/CFD/IFI16/CD47/JAK2/MAP2K7/BANF1/CEBPB/TRIM25/SAMSN1/CXCR3/IL6R/DAPK1/CD3D/HLA-B/EXOSC6/UBQLN1/CXCR4/PAK1/MORC3/SMPD1/HLA-A/STX8/SAMHD1/CD96/IL17RA/DDX1/IFIT1/CD274/TRAF2/YTHDF2/CASP4/SLAMF7/MAPK14/POLR3D/CMTM3/ADAR/IFNAR1/MYO1G/ATP7A/XRCC5/SEMA4A/GCH1/CREBBP/CDC37/CD48/ZC3HAV1/PYCARD/ERAP1/RHBDF2/ATG5/ITGAL")

genes_24h <- format_gene_list("HLA-DRB5/CCL4/IL1B/S100A8/GNLY/IL2RA/CXCR4/CCL5/CEBPB/TNFRSF1B/S100A9/LST1/TNFSF10/IFITM2/EIF2AK4/CXCL2/GZMB/HLA-DOA/MIF/CORO1A/IRF4/CD58/ISG20/IL32/CD59/RIF1/JAK1/FCER1G/CST7/TNFRSF4/UBE2J1/KYNU/PTPN2/TRIM38/LGALS3/TNFAIP3/PTPRC/NFKBIZ/WAS/TRAF4/CXCL1/HLA-B/MEF2C/CCL2/ACTR2/UBE2N/NR1H2/PQBP1/DHX36/FES/OPTN/FURIN/PHPT1/FTH1/PKN1/HLA-E/IL4I1/IL10RB/RNF19B/CTSC/XRCC5/HSPD1/MAPKAPK2/REL/NINJ1/GAPDH/IL2RG/TYK2/ICAM1/RFTN1/NKG7/LITAF/UBQLN1/RBPJ/CYLD/HLA-DQA2/IRF1/PRDX1/NFKB2/FCGRT/WNK1/NPLOC4/SP100")

genes <- unique(c(genes_24h, genes_3h, genes_UT))
 heatmap_immune_res__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
     scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

    heatmap_immune_res__3hPA <- DoHeatmap(cluster.average_3hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
      scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
      guides(fill = FALSE, color = FALSE) + ggtitle("3hPA") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

    heatmap_immune_res__24hPA <- DoHeatmap(cluster.average_24hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
      scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
      guides(color = FALSE) + ggtitle("24hPA") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

    grid.arrange(heatmap_immune_res__UT, heatmap_immune_res__3hPA, heatmap_immune_res__24hPA, ncol = 3)
#Subset the data based on the genes

genes_immune_res <- genes

#Subset the data based on the genes
immune_resma <- subset(pDC_combined, features = genes_immune_res)
average_immune_res <- AverageExpression(pDC_combined,group.by = c('subject','treatment'),slot = "scale.data",features = genes_immune_res,return.seurat = F)

#Transpose the data
average_immune_res <- t(average_immune_res$RNA)
#Convert the data to dataframe
average_immune_res <- as.data.frame(average_immune_res)
#Find the row means
average_immune_res <- cbind(rownames(average_immune_res),rowMeans(average_immune_res))
#Convert to Dataframe
average_immune_res <- as.data.frame(average_immune_res)
#Rename  columns into subject and timepoints
average_immune_res <- separate(average_immune_res,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_immune_res) <- c("subject","timepoints","avg_exprerssion")
#Convert to dataframe
average_immune_res <- as.data.frame(average_immune_res)
average_immune_res$avg_exprerssion <- as.numeric(average_immune_res$avg_exprerssion)
#Find the z-score
average_immune_res$z_score <- zscore(average_immune_res$avg_exprerssion)

meta <- pDC_combined@meta.data
average_immune_res$timepoints <-factor(average_immune_res$timepoints)
average_immune_res$sex <- meta$sex[match(average_immune_res$subject,meta$subject)]

#Remove NA data rows
average_immune_res <- average_immune_res[complete.cases(average_immune_res),]

ggplot(average_immune_res, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Immune Response") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#For PA treated samples
#Average expression of PA stimulated samples
average_immune_res_PA <- AverageExpression(pDC_PA_stim,group.by = c('subject','treatment'),slot = "scale.data",features = genes_immune_res,return.seurat = F)
average_immune_res_PA <- t(average_immune_res_PA$RNA)
average_immune_res_PA <- as.data.frame(average_immune_res_PA)
average_immune_res_PA <- cbind(rownames(average_immune_res_PA),rowMeans(average_immune_res_PA))
average_immune_res_PA <- as.data.frame(average_immune_res_PA)
average_immune_res_PA <- separate(average_immune_res_PA,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_immune_res_PA) <- c("subject","timepoints","avg_exprerssion")
average_immune_res_PA <- as.data.frame(average_immune_res_PA)
average_immune_res_PA$avg_exprerssion <- as.numeric(average_immune_res_PA$avg_exprerssion)
average_immune_res_PA$z_score <- zscore(average_immune_res_PA$avg_exprerssion)
average_immune_res_PA$timepoints <-factor(average_immune_res_PA$timepoints,levels = c("UT","3hPA","24hPA"))
average_immune_res_PA$sex <- meta$sex[match(average_immune_res_PA$subject,meta$subject)]
ggplot(average_immune_res_PA, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Immune Response for Pseudomonas Treatment") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
```
### Response to Cytokine
```{r}
#response to cytokine
genes_UT <- format_gene_list("ISG15/IFIT2/TCL1A/AGPAT2/MX1/CSF2RB/MT2A/TNFRSF4/CCL3/PLSCR1/ADAR/SOCS1/IRF1/SKIL/CD47/EIF5A/IL6ST/UGCG/MIR142/YTHDF2/CREB1/BIRC3/HSPA5/TRIM44/IL10RA/XRCC5/TNFSF13B/GAS6/REL/CDC42/STAT6/ZFAND6/MX2/STAT1/CSF2RA/ACTR3/BIRC2/HCLS1/LILRB4/NFKB1/SMPD3/SP100/IFI16/PTPN2/MBP/HAX1/CYLD/HSP90AB1/NFE2L2/NCL/IL3RA/TPR/SOCS3/RPLP0/SRSF3/KIF5B")

genes_3h <- format_gene_list("CXCL11/CXCL10/IFIT2/TNF/ISG15/IRF1/MX1/PARP14/CCL2/IFIH1/IFIT3/IRF8/EREG/GBP1/IRF7/TNFRSF1B/AIF1/OASL/TNFRSF18/KLF4/THBS1/HSPA1B/IFITM3/STAT2/CTBP2/PTPRC/PARP9/LYN/XBP1/FOS/LSM14A/TCL1A/CD40/EIF2AK2/SRM/STAT4/UGCG/PLSCR1/PTPN2/PDGFB/YY1/AZI2/DDOST/IL10RB/TNFRSF21/USP18/STAT1/TXNDC17/SOCS1/LILRA5/TRADD/TPR/COMMD7/ADAM10/TOLLIP/AGPAT2/PTPN1/FCER1G/METTL3/ANXA1/PTPN11/IFI16/CD47/JAK2/MAP2K7/CEBPB/TRIM25/CXCR3/IL6R/DAPK1/CARD16/CXCR4/SH2B3/SMPD1/STX8/SAMHD1/IL17RA/IFIT1/PDIA3/CD274/TRAF2/YTHDF2/CASP4/MAPK14/ADAR/IFNAR1/XRCC5/ADIPOR2/GCH1/ZFP36L2/CDC37/HIF1A/PYCARD")

genes_24h <- format_gene_list("CCL4/IL1B/MT2A/IL2RA/CXCR4/CCL5/CEBPB/TNFRSF1B/IFITM2/TNFRSF18/CXCL2/CORO1A/PLP2/CD58/RIF1/JAK1/FCER1G/TNFRSF4/KYNU/PTPN2/BAD/LIMS1/CD44/TNFAIP3/PTPRC/WAS/PPP2CB/CXCL1/RAD23B/CCL2/ACTR2/NR1H2/CREB1/NRP2/HAX1/HSPA5/IL10RB/XRCC5/P4HB/MAPKAPK2/REL/GAPDH/ADIPOR2/EEF1E1/IL2RG/TYK2/ICAM1/PDIA3")

genes <- unique(c(genes_24h, genes_3h, genes_UT))
 heatmap_res_cyto__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
     scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

    heatmap_res_cyto__3hPA <- DoHeatmap(cluster.average_3hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
      scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
      guides(fill = FALSE, color = FALSE) + ggtitle("3hPA") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

    heatmap_res_cyto__24hPA <- DoHeatmap(cluster.average_24hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
      scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
      guides(color = FALSE) + ggtitle("24hPA") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

    grid.arrange(heatmap_res_cyto__UT, heatmap_res_cyto__3hPA, heatmap_res_cyto__24hPA, ncol = 3)
#Subset the data based on the genes

genes_res_cyto <- genes

#Subset the data based on the genes
res_cytoma <- subset(pDC_combined, features = genes_res_cyto)
average_res_cyto <- AverageExpression(pDC_combined,group.by = c('subject','treatment'),slot = "scale.data",features = genes_res_cyto,return.seurat = F)

#Transpose the data
average_res_cyto <- t(average_res_cyto$RNA)
#Convert the data to dataframe
average_res_cyto <- as.data.frame(average_res_cyto)
#Find the row means
average_res_cyto <- cbind(rownames(average_res_cyto),rowMeans(average_res_cyto))
#Convert to Dataframe
average_res_cyto <- as.data.frame(average_res_cyto)
#Rename  columns into subject and timepoints
average_res_cyto <- separate(average_res_cyto,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_res_cyto) <- c("subject","timepoints","avg_exprerssion")
#Convert to dataframe
average_res_cyto <- as.data.frame(average_res_cyto)
average_res_cyto$avg_exprerssion <- as.numeric(average_res_cyto$avg_exprerssion)
#Find the z-score
average_res_cyto$z_score <- zscore(average_res_cyto$avg_exprerssion)

meta <- pDC_combined@meta.data
average_res_cyto$timepoints <-factor(average_res_cyto$timepoints)
average_res_cyto$sex <- meta$sex[match(average_res_cyto$subject,meta$subject)]

#Remove NA data rows
average_res_cyto <- average_res_cyto[complete.cases(average_res_cyto),]

ggplot(average_res_cyto, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Response to Cytokine") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#For PA treated samples
#Average expression of PA stimulated samples
average_res_cyto_PA <- AverageExpression(pDC_PA_stim,group.by = c('subject','treatment'),slot = "scale.data",features = genes_res_cyto,return.seurat = F)
average_res_cyto_PA <- t(average_res_cyto_PA$RNA)
average_res_cyto_PA <- as.data.frame(average_res_cyto_PA)
average_res_cyto_PA <- cbind(rownames(average_res_cyto_PA),rowMeans(average_res_cyto_PA))
average_res_cyto_PA <- as.data.frame(average_res_cyto_PA)
average_res_cyto_PA <- separate(average_res_cyto_PA,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_res_cyto_PA) <- c("subject","timepoints","avg_exprerssion")
average_res_cyto_PA <- as.data.frame(average_res_cyto_PA)
average_res_cyto_PA$avg_exprerssion <- as.numeric(average_res_cyto_PA$avg_exprerssion)
average_res_cyto_PA$z_score <- zscore(average_res_cyto_PA$avg_exprerssion)
average_res_cyto_PA$timepoints <-factor(average_res_cyto_PA$timepoints,levels = c("UT","3hPA","24hPA"))
average_res_cyto_PA$sex <- meta$sex[match(average_res_cyto_PA$subject,meta$subject)]
ggplot(average_res_cyto_PA, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Response to cytokine for Pseudomonas Treatment") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
```
### Cytokine-mediated signaling pathway
```{r}
#cytokine-mediated signaling pathway
genes_UT <- format_gene_list("ISG15/AGPAT2/MX1/CSF2RB/TNFRSF4/CCL3/ADAR/SOCS1/IRF1/EIF5A/IL6ST/UGCG/YTHDF2/BIRC3/TRIM44/IL10RA/TNFSF13B/GAS6")

genes_3h <- format_gene_list("CXCL11/CXCL10/TNF/ISG15/IRF1/MX1/PARP14/CCL2/EREG/IRF7/TNFRSF1B/OASL/TNFRSF18/HSPA1B/IFITM3/STAT2/PTPRC/PARP9/LYN/LSM14A/STAT4/UGCG/PTPN2/PDGFB/AZI2/IL10RB/USP18/STAT1/TXNDC17/SOCS1/LILRA5/TRADD/COMMD7/TOLLIP/AGPAT2/PTPN1/FCER1G/METTL3/PTPN11/JAK2/CXCR3/IL6R/CARD16/CXCR4/SH2B3/SAMHD1/IL17RA/TRAF2/YTHDF2/CASP4/ADAR/IFNAR1/ADIPOR2/CDC37/HIF1A/PYCARD")

genes_24h <- format_gene_list("CCL4/IL1B/IL2RA/CXCR4/CCL5/TNFRSF1B/IFITM2/TNFRSF18/CXCL2/PLP2/JAK1/FCER1G/TNFRSF4/PTPN2/BAD/LIMS1/CD44/TNFAIP3/PTPRC/PPP2CB/CXCL1/CCL2/NR1H2/IL10RB/P4HB/ADIPOR2/IL2RG/TYK2/PADI2/CYLD/TRAF1/IRF1/WNK1/TRADD/SP100")

genes <- unique(c(genes_24h, genes_3h, genes_UT))
 heatmap_cyto_sig__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
     scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

    heatmap_cyto_sig__3hPA <- DoHeatmap(cluster.average_3hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
      scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
      guides(fill = FALSE, color = FALSE) + ggtitle("3hPA") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

    heatmap_cyto_sig__24hPA <- DoHeatmap(cluster.average_24hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
      scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
      guides(color = FALSE) + ggtitle("24hPA") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

    grid.arrange(heatmap_cyto_sig__UT, heatmap_cyto_sig__3hPA, heatmap_cyto_sig__24hPA, ncol = 3)
#Subset the data based on the genes

genes_cyto_sig <- genes

#Subset the data based on the genes
cyto_sigma <- subset(pDC_combined, features = genes_cyto_sig)
average_cyto_sig <- AverageExpression(pDC_combined,group.by = c('subject','treatment'),slot = "scale.data",features = genes_cyto_sig,return.seurat = F)

#Transpose the data
average_cyto_sig <- t(average_cyto_sig$RNA)
#Convert the data to dataframe
average_cyto_sig <- as.data.frame(average_cyto_sig)
#Find the row means
average_cyto_sig <- cbind(rownames(average_cyto_sig),rowMeans(average_cyto_sig))
#Convert to Dataframe
average_cyto_sig <- as.data.frame(average_cyto_sig)
#Rename  columns into subject and timepoints
average_cyto_sig <- separate(average_cyto_sig,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_cyto_sig) <- c("subject","timepoints","avg_exprerssion")
#Convert to dataframe
average_cyto_sig <- as.data.frame(average_cyto_sig)
average_cyto_sig$avg_exprerssion <- as.numeric(average_cyto_sig$avg_exprerssion)
#Find the z-score
average_cyto_sig$z_score <- zscore(average_cyto_sig$avg_exprerssion)

meta <- pDC_combined@meta.data
average_cyto_sig$timepoints <-factor(average_cyto_sig$timepoints)
average_cyto_sig$sex <- meta$sex[match(average_cyto_sig$subject,meta$subject)]

#Remove NA data rows
average_cyto_sig <- average_cyto_sig[complete.cases(average_cyto_sig),]

ggplot(average_cyto_sig, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Cytokine-Mediated Signalling Pathway") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#Average expression of PA stimulated samples
average_cyto_sig_PA <- AverageExpression(pDC_PA_stim,group.by = c('subject','treatment'),slot = "scale.data",features = genes_cyto_sig,return.seurat = F)
average_cyto_sig_PA <- t(average_cyto_sig_PA$RNA)
average_cyto_sig_PA <- as.data.frame(average_cyto_sig_PA)
average_cyto_sig_PA <- cbind(rownames(average_cyto_sig_PA),rowMeans(average_cyto_sig_PA))
average_cyto_sig_PA <- as.data.frame(average_cyto_sig_PA)
average_cyto_sig_PA <- separate(average_cyto_sig_PA,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_cyto_sig_PA) <- c("subject","timepoints","avg_exprerssion")
average_cyto_sig_PA <- as.data.frame(average_cyto_sig_PA)
average_cyto_sig_PA$avg_exprerssion <- as.numeric(average_cyto_sig_PA$avg_exprerssion)
average_cyto_sig_PA$z_score <- zscore(average_cyto_sig_PA$avg_exprerssion)
average_cyto_sig_PA$timepoints <-factor(average_cyto_sig_PA$timepoints,levels = c("UT","3hPA","24hPA"))
average_cyto_sig_PA$sex <- meta$sex[match(average_cyto_sig_PA$subject,meta$subject)]
ggplot(average_cyto_sig_PA, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Cytokine Mediated signalling pathway for Pseudomonas Treatment") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
```
### Respnse to Biotic stimulus
```
genes_UT <-format_gene_list("NR1H2/TUBB/TNFAIP3/IFITM2/GPX1/GNLY/IFNAR2/IFNGR2/NR4A1/TSPO/COTL1/KLF4/BAX/PYCARD/C15orf48/TYROBP/TRIM38/ST13/GAPDH/HSP90AA1/CD36/IL1B/CCL5/DDX3X/AIF1/S100A8/SOD2/FOS/S100A9/PPBP/LYZ")

genes_3h <- format_gene_list("CXCL11/CXCL10/IFIT2/TNF/NMB/RSAD2/ISG15/IRF1/MX1/PARP14/CCL2/IFI6/IFIH1/IFIT3/CMPK2/S100A9/C15orf48/IRF8/EREG/IFI44L/GZMB/SOD2/GBP1/IRF7/ISG20/TNFRSF1B/AIF1/OASL/USP15/RGS1/KLF4/ARID5A/IRF9/ZMPSTE24/NT5C3A/IFITM3/GTF2F1/WDR83/STAT2/PTPRC/PARP9/LYN/PLAC8/XBP1/FOS/LSM14A/IFI35/PTGIR/CD40/EIF2AK2/PLSCR1/PTPN2/SLFN11/KLRB1/PARP1/LY96/AZI2/HLA-E/RNASE6/ZBTB1/SPN/NR4A1/TMF1/NKG7/IL10RB/GSDMD/USP18/STAT1/MECP2/COTL1/SOCS1/LILRA5/ADH5/LYZ/HLA-C/RNF213/EXOC1/TOLLIP/PTPN1/FCER1G/TLR7/METTL3/FYN/ANXA1/PTPN11/KAT5/CFD/IFI16/CD47/JAK2/MAP2K7/BANF1/CEBPB/TRIM25/HNRNPA0/IL6R/DAPK1/CARD16/HLA-B/CHMP5/RAB14/HERC6/CXCR4/PAK1/MORC3/SMPD1/HLA-A/STX8/SAMHD1/CD96/IL17RA/DDX1/IFIT1/CD274/YTHDF2/CASP4/SLAMF7/MAPK14/STMN1/POLR3D/ARMC5/ADAR/IFNAR1/XRCC5/GCH1/CREBBP/CDC37/HIF1A/ZC3HAV1/PYCARD/ERAP1/ATG5/HADHB/PPP1R11/ILF3/ZYX")

genes_24h <- format_gene_list("C15orf48/SOD2/PLAC8/CCL4/IL1B/S100A8/GNLY/RNF213/CXCR4/CCL5/CEBPB/TNFRSF1B/S100A9/IFITM2/EIF2AK4/PMAIP1/CXCL2/GZMB/MIF/CORO1A/GTF2F1/IRF4/CD58/OTUD5/ISG20/JAK1/FCER1G/DDX17/KYNU/PTPN2/TRIM38/LGALS3/TNFAIP3/PTPRC/WAS/TRAF4/CXCL1/HLA-B/MEF2C/CCL2/ACTR2/NR1H2/PQBP1/DHX36/AUP1/FES/OPTN/ODC1/TSPO/HSPA5/HLA-E/IL4I1/CCT5/IL10RB/RNF19B/BNIP3L/FKBP5/XRCC5/HSPD1/MAPKAPK2/REL/PPM1B/GAPDH/JUND")

genes <- unique(c(genes_UT, genes_3h, genes_24h))

heatmap_biotic__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
     scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

    heatmap_biotic__3hPA <- DoHeatmap(cluster.average_3hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
      scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
      guides(fill = FALSE, color = FALSE) + ggtitle("3hPA") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

    heatmap_biotic__24hPA <- DoHeatmap(cluster.average_24hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
      scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
      guides(color = FALSE) + ggtitle("24hPA") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

    grid.arrange(heatmap_biotic__UT, heatmap_biotic__3hPA, heatmap_biotic__24hPA, ncol = 3)
#Subset the data based on the genes

genes_biotic <- genes

#Subset the data based on the genes
bioticma <- subset(pDC_combined, features = genes_biotic)
average_biotic <- AverageExpression(pDC_combined,group.by = c('subject','treatment'),slot = "scale.data",features = genes_biotic,return.seurat = F)

#Transpose the data
average_biotic <- t(average_biotic$RNA)
#Convert the data to dataframe
average_biotic <- as.data.frame(average_biotic)
#Find the row means
average_biotic <- cbind(rownames(average_biotic),rowMeans(average_biotic))
#Convert to Dataframe
average_biotic <- as.data.frame(average_biotic)
#Rename  columns into subject and timepoints
average_biotic <- separate(average_biotic,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_biotic) <- c("subject","timepoints","avg_exprerssion")
#Convert to dataframe
average_biotic <- as.data.frame(average_biotic)
average_biotic$avg_exprerssion <- as.numeric(average_biotic$avg_exprerssion)
#Find the z-score
average_biotic$z_score <- zscore(average_biotic$avg_exprerssion)

meta <- pDC_combined@meta.data
average_biotic$timepoints <-factor(average_biotic$timepoints)
average_biotic$sex <- meta$sex[match(average_biotic$subject,meta$subject)]

#Remove NA data rows
average_biotic <- average_biotic[complete.cases(average_biotic),]

ggplot(average_biotic, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Response to Biotic Stimuls") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#Average expression of PA stimulated samples
average_biotic_PA <- AverageExpression(pDC_PA_stim,group.by = c('subject','treatment'),slot = "scale.data",features = genes_biotic,return.seurat = F)
average_biotic_PA <- t(average_biotic_PA$RNA)
average_biotic_PA <- as.data.frame(average_biotic_PA)
average_biotic_PA <- cbind(rownames(average_biotic_PA),rowMeans(average_biotic_PA))
average_biotic_PA <- as.data.frame(average_biotic_PA)
average_biotic_PA <- separate(average_biotic_PA,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_biotic_PA) <- c("subject","timepoints","avg_exprerssion")
average_biotic_PA <- as.data.frame(average_biotic_PA)
average_biotic_PA$avg_exprerssion <- as.numeric(average_biotic_PA$avg_exprerssion)
average_biotic_PA$z_score <- zscore(average_biotic_PA$avg_exprerssion)
average_biotic_PA$timepoints <-factor(average_biotic_PA$timepoints,levels = c("UT","3hPA","24hPA"))
average_biotic_PA$sex <- meta$sex[match(average_biotic_PA$subject,meta$subject)]
ggplot(average_biotic_PA, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Response to biotic stimulus for Pseudomonas Treatment") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
```
### Cell-Cell Adhesion
```{r}
genes_UT <- format_gene_list("SPI1/ICAM3/CYFIP2/MSN/ACTG1/MYADM/DNAJB6/ANXA2/ACTB/PLAUR/PTEN/MAP3K8/CD44/SOX4/FXYD5/ITGAE/KLF4/PYCARD/ZFP36L1/RGCC/FLNA/LGALS1/IL1B/CCL5/AIF1/S100A8/S100A9")

genes_3h <- format_gene_list("CCL2/S100A9/IRF1/ANXA1/S100A11/S100A8/SOCS1/CLIC1/RUNX3/HLA-E/LYN/CD44/ADAM19/ANXA2/HLA-A/PTPRC/IL1RN/SPN/LAPTM5/FERMT3/SMARCC1/PTPN2/XBP1/FXYD5/EMB/MSN/MYL12A/ACTB/CD99/TNFSF13B/CD164")

genes_24h <- format_gene_list("HLA-DRB5/IL1B/S100A8/IL2RA/CCL5/CEBPB/S100A9/VMP1/HLA-DOA/CORO1A/TGFBR2/CD58/SOX4/JAK1/ADAM19/PTPN2/BAD/LIMS1/LGALS3/CD44/PTPRC/NFKBIZ/CCL2/ANXA2/HLA-E/IL4I1/BRD7/HSPD1/MSN/NINJ1/FXYD5/CD2AP/IL2RG/TYK2/ICAM1/CTNNA1")

genes <- unique(c(genes_UT,genes_3h,genes_24h))

heatmap_cell_cell__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
     scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

    heatmap_cell_cell__3hPA <- DoHeatmap(cluster.average_3hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
      scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
      guides(fill = FALSE, color = FALSE) + ggtitle("3hPA") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

    heatmap_cell_cell__24hPA <- DoHeatmap(cluster.average_24hPA, features = genes, group.by = 'sex', draw.lines = TRUE) +
      scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
      guides(color = FALSE) + ggtitle("24hPA") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

    grid.arrange(heatmap_cell_cell__UT, heatmap_cell_cell__3hPA, heatmap_cell_cell__24hPA, ncol = 3)
#Subset the data based on the genes

genes_cell_cell <- genes

#Subset the data based on the genes
cell_cellma <- subset(pDC_combined, features = genes_cell_cell)
average_cell_cell <- AverageExpression(pDC_combined,group.by = c('subject','treatment'),slot = "scale.data",features = genes_cell_cell,return.seurat = F)

#Transpose the data
average_cell_cell <- t(average_cell_cell$RNA)
#Convert the data to dataframe
average_cell_cell <- as.data.frame(average_cell_cell)
#Find the row means
average_cell_cell <- cbind(rownames(average_cell_cell),rowMeans(average_cell_cell))
#Convert to Dataframe
average_cell_cell <- as.data.frame(average_cell_cell)
#Rename  columns into subject and timepoints
average_cell_cell <- separate(average_cell_cell,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_cell_cell) <- c("subject","timepoints","avg_exprerssion")
#Convert to dataframe
average_cell_cell <- as.data.frame(average_cell_cell)
average_cell_cell$avg_exprerssion <- as.numeric(average_cell_cell$avg_exprerssion)
#Find the z-score
average_cell_cell$z_score <- zscore(average_cell_cell$avg_exprerssion)

meta <- pDC_combined@meta.data
average_cell_cell$timepoints <-factor(average_cell_cell$timepoints)
average_cell_cell$sex <- meta$sex[match(average_cell_cell$subject,meta$subject)]

#Remove NA data rows
average_cell_cell <- average_cell_cell[complete.cases(average_cell_cell),]

ggplot(average_cell_cell, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Cell-Cell Adhesion") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#Average expression of PA stimulated samples
average_cell_cell_PA <- AverageExpression(pDC_PA_stim,group.by = c('subject','treatment'),slot = "scale.data",features = genes_cell_cell,return.seurat = F)
average_cell_cell_PA <- t(average_cell_cell_PA$RNA)
average_cell_cell_PA <- as.data.frame(average_cell_cell_PA)
average_cell_cell_PA <- cbind(rownames(average_cell_cell_PA),rowMeans(average_cell_cell_PA))
average_cell_cell_PA <- as.data.frame(average_cell_cell_PA)
average_cell_cell_PA <- separate(average_cell_cell_PA,col = "V1",into =c("subject","timepoints"),sep = "_")
colnames(average_cell_cell_PA) <- c("subject","timepoints","avg_exprerssion")
average_cell_cell_PA <- as.data.frame(average_cell_cell_PA)
average_cell_cell_PA$avg_exprerssion <- as.numeric(average_cell_cell_PA$avg_exprerssion)
average_cell_cell_PA$z_score <- zscore(average_cell_cell_PA$avg_exprerssion)
average_cell_cell_PA$timepoints <-factor(average_cell_cell_PA$timepoints,levels = c("UT","3hPA","24hPA"))
average_cell_cell_PA$sex <- meta$sex[match(average_cell_cell_PA$subject,meta$subject)]
ggplot(average_cell_cell_PA, aes(x = timepoints, y = z_score, color = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # add a dashed red line at y=0
  xlab("Timepoint") +
  ylab("Zscore") +
  ggtitle("Cell-Cell Adhesion for Pseudomonas Treatment") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
```

## Supplementary Analysis
Apart from the analysis on Pseudomonas we visulaised the difference in sexes for their immune response against the treatments MTB,CA and even other pathways in PA.
```
library(paletteer)
library(RColorBrewer)
library(gridExtra)
library(stringr)
library(tidyr)
library(gridExtra)
source ("scripts/Gene_Formatting.R")
#Analysing the metabolic pathways between sexe at baseline
my_palette <- paletteer::paletteer_c("grDevices::Blue-Red", n = 20, direction = -1)

#GO-nitric oxide metabolic process(Gene list from the core_enrichment genes)
j <-DoHeatmap(cluster.average_UT,features=c("TSPO","KLF4","HSP90AA1","CD36","IL1B","AIF1","SOD2"),group.by = 'sex',draw.lines = T,)+scale_fill_gradientn(colors = rev(brewer.pal(1000, "RdBu")),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
j<-saveRDS(j, file = "nitric.rds")

#Using all the genes from the pathway
DoHeatmap(cluster.average_UT,features = c("TLR2", "NOS1", "CYP1B1", "CYB5R3", "DDAH2", "NQO1", "GCH1", "CYB5B", "NOS3", "ARG2", "MTARC1", "TLR4", "TICAM1", "SLC7A6", "KLRK1", "AKT1", "POR", "CPS1", "MTARC2", "RORA", "GCHFR", "ACP5", "SPR", "NOS2", "CAV1"),group.by = 'sex',disp.min = -1,disp.max = 1)

#KEGG Carbon metabolism
i <- DoHeatmap(cluster.average_UT,features = c("GRSF1","OCIAD1","RAB21","SOCS1","PTGDS","CSF2RB","MARCH1"),group.by = 'sex')+ scale_fill_gradientn(colors = rev(brewer.pal(10, "RdBu")), oob = scales::oob_squish_any, limits = c(-2, 2))
saveRDS(i, file = "Carbon.rds")
#Glycolysis
DoHeatmap(cluster.average_UT,features = c("GAPDH","ADH5","ADPGK","PGK1","LDHA","ALDH3B1","MINPP1","PGAM1","PFKM","PFKL","ALDH9A1","ALDH2","PKM","GPI","ALDH7A1","PGAM4","FBP1","ALDOA","ENO2","LDHB","ALDH1B1","HK3"),group.by = 'sex',disp.min = -1,disp.max = 1)+ scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-1, 1))

#TCA cycle
DoHeatmap(cluster.average_UT,features = c("PDHB","MDH2","OGDH","SDHC","FH","DLAT","ACO2","SDHA","IDH3G","ACLY","CS","SUCLG1","DLD","SDHD","IDH3A"),group.by = 'sex',disp.min = -2,disp.max = 2)+scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"), oob = scales::oob_squish_any,
                       limits = c(-2, 2))

# Oxidative phosphorylation
DoHeatmap(cluster.average_UT,features = c("NDUFA11","COX5B","ATP5ME","COX8A","UQCRQ","NDUFB1","NDUFA1","ATP5MC2","NDUFS6","COX7C","UQCR11","COX6B1","NDUFB3","ATP5MG","UQCRH","COX7B2"),group.by = 'sex')+scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"), oob = scales::oob_squish_any,
                       limits = c(-2, 2))

#Upregulation of Microphilial/Macrophage activation_GO
DoHeatmap(cluster.average_UT,features = c("CCL3", "STAP1","CTSC","CCL3","MIR142"),group.by = 'sex') + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))

#IL-6 Production
a<-DoHeatmap(cluster.average_UT,features =c("ARRB2","TNFAIP3","PYCARD","TYROBP","CD36","IL1B","AIF1"),group.by = 'sex',disp.min = -1,disp.max = 1)+scale_fill_gradientn(colors = rev(brewer.pal(1000, "RdBu")),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
saveRDS(a, file = "IL6.rds")
#IL-17 Production
DoHeatmap(cluster.average_UT,features=c("RFTN1","SLC7A5","IFNG" ,"LILRB2", "AIF1",  "AKIRIN2" ,"CD36", "EREG", "ZBTB20", "HYAL2" ,"APP" ,"TSLP" ,"F2R","PTAFR","SETD4","NOS2","TLR3","HMGB1"),group.by='sex')+scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Positive regulation of monocyte chemotaxis_GO
DoHeatmap(cluster.average_UT,features =c("LGMN","CCL5","AIF1"),group.by = 'sex',disp.min = -1,disp.max = 1)+scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))

#Biosynthesis of amino acids-KEGG
k <-DoHeatmap(cluster.average_UT,features = c("GRSF1","OCIAD1","RAB21","CSF2RB","MARCH1"),group.by = 'sex',disp.min = -1,disp.max = 1)+ scale_fill_gradientn(colors = rev(brewer.pal(10, "RdBu")), oob = scales::oob_squish_any, limits = c(-2, 2))
saveRDS(k, file = "Amino.rds")


#Cholesterol Metabolism-KEGG
DoHeatmap(cluster.average_UT,features = c("CCT3","HNRNPA2B1","LIMD2","RPL28"),group.by = 'sex',disp.min = -1.5,disp.max = 1.5)+  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2))

#Pentose Phosphate Pathway-KEGG

DoHeatmap(cluster.average_UT, features = c("OCIAD1", "RAB21", "PTGDS"),
          group.by = 'sex', disp.min = -1, disp.max = 1) +
   scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2))
```

position the caret at any line or the code chunk, then click "+".

The code chunk appears:
```{r}
#Pathways enriched at 3h CA when compared with UT

# Male
#Regulation of Viral genome replication-GO
DoHeatmap(cluster.average_3hCA,features = c("ISG15","ISG20","MX1","IFI16","PLSCR1","OAS2","EIF2AK2","IFIH1","ADAR","OAS1"),group.by = 'sex',disp.min = -1,disp.max = 1)+scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Using all the genes for the gene ontology
DoHeatmap(cluster.average_3hCA, features = c("DDB1", "PKN2", "PPIE", "AICDA", "TRIM6", "OAS2", "IFI16", "APOBEC3B", "PPIA", "ZNFX1", "HMGA2", "CNOT7", "CD28", "BANF1", "VAPB", "STAU1", "BCL2", "ISG20", "CXCL8", "TASOR", "MX1", "N4BP1", "IFIT5", "OAS1", "FAM111A", "TMEM39A", "EIF2AK2", "IFNB1", "HACD3", "LARP1", "TRIM28", "SRPK1", "GBP7", "RAD23A", "APOBEC3C", "APOBEC3F", "FKBP6", "PPARA", "CCL5", "RESF1", "SHFL", "TOP2A", "TRIM38", "MPHOSPH8", "SETDB1", "TNIP1", "IFIH1", "ZC3HAV1", "APOBEC3D", "ADARB1", "ADAR", "PABPC1", "APOBEC3G", "IFITM1", "NOTCH1", "IFIT1", "APOBEC3H", "SRPK2", "DDX3X", "RSAD2", "BST2", "RNASEL", "PDE12", "BTBD17", "MORC2", "IFNL3", "TNF", "OAS3", "PPIH", "MAVS", "APOBEC3A", "PROX1", "OASL", "TARBP2", "LTF", "ISG15", "SLPI", "TOP2B", "ILF3", "INPP5K", "PPID", "PLSCR1", "IFITM2", "IFITM3"), group.by = 'sex', disp.min = -2, disp.max = 2)

#Defense response to virus-GO
DoHeatmap(cluster.average_3hCA,features = c("CXCL10","ISG15","PMAIP1","BIRC3","ISG20","IFIT2","MX1","IFIT3","IFI44L","IFI16","PLSCR1","STAT1","OAS2","EIF2AK2","NT5C3A","IFIH1","ADAR","HERC5","IRF1","IFI6","TANK","OAS1","MX2"),group.by = 'sex',disp.min = -1,disp.max = 1)+scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))

#All the genes
gene_list <- c("IFNL4", "DDX60L", "TLR2", "CARD9", "CASP1", "IL12B", "DDX21", "DHX15", "FADD", "TRAF3", "DHX58", "FCN3", "TBK1", "TRIM6", "OAS2", "TRIM52", "DHX36", "IFI16", "LILRB1", "AGBL4", "TBKBP1", "APOBEC3B", "DMBT1", "DEFA3", "DEFA1", "GBP2", "GBP1", "NT5C3A", "ATG7", "RIGI", "RTP4", "ISG20", "HERC5", "UBE2W", "CD207", "RELA", "BECN1", "IL6", "MLKL", "AZU1", "PYCARD", "DDIT4", "RIPK3", "STAT1", "ZBP1", "CD40", "IRF3", "IFIT5", "TLR7", "TRIM27", "UNC13D", "SLFN13", "DHX16", "OAS1", "DDX60", "IL10RB", "ABCF3", "EIF2AK2", "TRIM56", "IFNA2", "IFNA10", "IFNA21", "IFNA5", "IFNA17", "IFNB1", "TRIM28", "GBP7", "UBE2N", "PRF1", "PML", "IKBKE", "APOBEC3C", "AZI2", "APOBEC3F", "GARIN5A", "UNC93B1", "TRIM8", "STAT2", "BCL2L1", "STING1", "URS0000338542_9606", "CLPB", "sting_human-1", "alt-cbc_human", "MYD88", "MAP3K14", "TRIM22", "IFNE", "USP44", "NLRC5", "IFNL1", "IFIH1", "SETD2", "EIF2AK4", "IRF9", "ZC3HAV1", "NLRP9", "USP29", "APOBEC3D", "TANK", "IRF5", "ADARB1", "ADAR", "RNF185", "APOBEC3G", "GBP3", "TRIM41", "IFNGR2", "MBL2", "RSAD2", "BST2", "GPAM", "TLR8", "IFNK", "IRF7", "IFNL2", "IFNL3", "TNF", "OAS3", "NLRP6", "IFNG", "CARD8", "F2RL1", "MAVS", "USP20", "PARP9", "APOBEC3A", "OASL", "IFNAR2", "TRIM26", "IFNA16", "IFNW1", "EXOSC5", "POLR3A", "GBP5", "IFI27", "IRF2")

DoHeatmap(cluster.average_3hCA, features = gene_list, group.by = 'sex', disp.min = -1, disp.max = 1)

#Negative Regulation of Immune System Process-GO
DoHeatmap(cluster.average_3hCA,features=c("IL2RA","IL7R","ISG15","AHR","IFI16","PARP14","IL4I1","ADTRP","LPXN","PAG1","CCL3","HLA-A","ADAR","IRF1","HLA-F","RUNX3","CD47","USP18","OAS1","HLA-E","CASP3","CCL2","HLA-DRB1","SERPINB9","SAMSN1","BTN2A2","NMI","TNFAIP3","ELF1"),group.by = 'sex',disp.min = -1,disp.max = 1)+scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))

#Cell Death-GO
DoHeatmap(cluster.average_3hCA, features = c("CXCL10","OPN3","CTSC","TXN","BCL2A1","IL2RA","ENDOG","IL7R","CFLAR","BTG1","GCLM","SLC7A11","NFKB1","CDKN1A","PMAIP1","BIRC3","IFIT2","MX1","IFIT3","SQSTM1","AHR","IFI16","PLSCR1","TNFSF10","PTPN1","BID","ATF5","STAT1","TNFAIP8","EIF2AK2","AKR1C3","KMO","MARCKS","CCR7","ARL6IP5","GZMB","CCL3","NINJ1","SRGN","XAF1","STK4","ADAR","EMP3","IRF1","IFI6","HLA-F","TAF9","CASP3","TRAF1","PKM","SP100","DNAJA1","CCL2","CD38","ACTN4","FAS","NUP62","ATF4","SERPINB9","RTN4","RBCK1","HSP90B1"), group.by = 'sex', disp.min = -1, disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Cell Migration-GO
DoHeatmap(cluster.average_3hCA, features = c("SMPD3", "CCL5", "AIF1", "B4GALT1", "CORO1A", "MAP2K3", "RAC2", "SELPLG", "PHACTR1", "SERPINF1", "GPX1", "LDLRAD4", "DUSP1", "PPBP", "PYCARD", "ACTG1", "KLF4", "ARHGDIB", "S100A8", "NR4A2", "NR4A1", "S100A9", "MYADM", "ID1", "RGCC", "CXCR4", "FCER1G"), group.by = "sex", disp.min = -1, disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Skeletal Muscle cell differentitaion-GO
DoHeatmap(cluster.average_3hCA,features = c("NR4A1","FOS"),group.by="sex",disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Using All the genes
DoHeatmap(cluster.average_3hCA,features =c("HIVEP3","MED20","SAP30","LEMD2","GTF3C5","NR4A1","SMYD1","MAFF","MEF2D","MEGF10","HEYL","PHOX2B","VAX1","SELENON","BCL9L","BTG2","ANKRD33","FOXN2","KRAS","CDON","FOS","DMRTA2","MYF6","KLF5","SCX","NOTCH1","MYF5","WNT3A","CITED2","SHH","ATF3","HMG20B","KLHL40","BCL9","ZNF689","EGR1","EGR2","RB1","EMD","PAX5","BARX2","SOX11","SOX8","MSTN","ASB2","MYOD1","COPS2"),group.by ='sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#CYTOPLASMIC TRANSLATION-GO
DoHeatmap(cluster.average_3hCA,features = c("RPS27A","FAU","EIF3G","RPS20","RPL38","RPS21","RPS27","RPL22L1","RPL36","RPS19","RPL36A","YBX1","RPS11","EIF2S3","RPL18A","RPS10","RPL17","RPS28","RPS6","RPL23","RPL23A","RPS16","RPS12","RPL37A","RPL3","RPL35","RPLP2","RPL9","RPL35A","RPL8","RPL27","RPL7","RPLP0","RPL19","RPL37","RPS23","RPS8","RPS24","RPS4X","UBA52","RPS14","RPS25","EIF3E","RPS13","RPL32","RPS5","RPL30","RPL29","RPL11","RPL13","RPS3","RPL28","RPS3A","RPS15A","EIF4B","RPL18","RPS18","RPL14","RPL22","RPL27A","RPL31","RPS2","EIF3D","RPS15","RPL10","RPS7","RPL6","RPL12","RPS9","RPL26","RPL15","EIF3F","RPL4","RPL10A","MIF4GD","EIF3K","RPL13A","EIF3L","EIF3H","RPL5","RPSA","RPL7A"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))

#KEGG-RIG 1 LIKE RECEPTOR SIGNALLING
DoHeatmap(cluster.average_3hCA,features=c("ADAR" ,"CXCL10" ,"IFIH1" ,"ISG15" ,"NFKB1" ,"TANK"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Ferroptosis-Kegg
DoHeatmap(cluster.average_3hCA,features = c("FTH1" ,"FTL" ,"GCLM" ,"SAT1" ,"SLC7A11"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#AutoImmune Thyroid Disease
DoHeatmap(cluster.average_3hCA,features=c("FAS", "GZMB", "HLA-A", "HLA-DQA1", "HLA-DQA2", "HLA-DRB1", "HLA-E", "HLA-F"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Apoptosis-kegg
DoHeatmap(cluster.average_3hCA,features = c("ATF4", "BCL2A1", "BID", "BIRC3", "CASP3", "CFLAR", "CTSC", "ENDOG", "FAS", "GZMB", "IL3RA", "NFKB1", "PMAIP1", "TNFSF10", "TRAF1"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Ribosome-KEGG
DoHeatmap(cluster.average_3hCA,features=c("MRPL23", "RPL10", "RPL10A", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL26", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL3", "RPL30", "RPL31", "RPL32", "RPL35", "RPL35A", "RPL36", "RPL37", "RPL37A", "RPL38", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPLP0", "RPLP2", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS18", "RPS19", "RPS2", "RPS20", "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS28", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPSA", "UBA52"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Fc Epsilon RI Signaling Pathway-kegg
DoHeatmap(cluster.average_3hCA,features = c("FCER1G", "MAP2K3" ,"RAC2"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
```

Type any R code in the chunk, for example:
```{r}
#Pathways enriched at 24h CA when compared with UT
#MALE-GO
#Farensyl diphosphate metabolic process - GO
DoHeatmap(cluster.average_24hCA,features=c("FDPS","HMGCS1","FDFT1"),group.by = 'sex',disp.min = -2,disp.max = 2) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Nitric oxide biosynthetic process - GO
DoHeatmap(cluster.average_24hCA,features = c("EEF1D","EIF4A1","DNAJC3","EEF1A1","EIF3K","MYCL","PABPC1","RPLP1","RPL41","CCND3","MIF4GD","RPL24","TCF4","RPS25","RRBP1","PIM2","IL1B","IFI16","RPL18","RPS29","RPL12","JMJD1C","CCNL1","RPL10","RPLP0","FAU","RPL11","CXXC5","EEF2","NFKBIA","RPL38","RPS12","RPS27A","EIF3F","POLB","CSRNP1","RPL7A","ZFAT","RPS5","ZNF791","RBM3","RPSA","IRF4","CDKN2D","SERP1","USF2","PRKCB","RPL15","RPS21","RPL19","SP110","RPL13A","RPL31","GAS6","RPL27A","RPS9","RPL35A","RPL39","RPL36A","FOSL2","RPS3","JUND","S100A8","RPL18A","RPL21","RPL37","ARID5A","RPS11","NR4A3","RPL37A","EIF1","RPL4","RPS23","RPLP2","RPL7","EIF3L","RPS2","IER2","RPS4X","RPL23A","RPL30","S100A9","SIK1","RPS27","HES4","RPL29","RPL9","RPL35","RPL22","RPL26","RPS16","RPS28","MYBL2","RPL34","RPL27","RPS10","RPS24","VIM","RPL23","RPS7","RPS15A","RPS13","CD4","RPS6","RPS14","SPIB","RPL17","ODC1","KLF10","RPL32","SMPD3","NR4A1","RPL10A","RPL5","RPL13","EZR","GUK1","MAP2K3","APP","NR3C1","RPL3","EEF1B2","SKIL","TSC22D3","JUNB","RGCC","ARID3A","RPS8","RPS18","COMMD6","RPL6","XBP1","NR4A2","TGFB1","NCF1","RPS3A","PLAC8","ZFP36L2","IRF7","FOSB","IRF8","KDM6B","ZNF331","BTG2","PPP1R15A","RGS2","ZFP36"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Gene Expression- GO
DoHeatmap(cluster.average_24hCA,features = c("TCF4","SRSF6","RPS25","MIR142","RRBP1","PIM2","IL1B","SAMHD1","IFI16","RPL18","SRSF3","RPS29","SPON2","RPL12","JMJD1C","CCNL1","RPL10","ISCA1","RPLP0","FAU","RPL11","CXXC5","EEF2","SRGN","NFKBIA","RPL38","RPS12","RPS27A","EIF3F","POLB","CSRNP1","RPL7A","ZFAT","RPS5","ZNF791","RBM3","RPSA","SPCS1","IRF4","TNFRSF21","SERP1","USF2","LAPTM5","FBL","PRKCB","RPL15","RPS21","RPL19","MBP","SP110","RPL13A","RPL31","AHNAK","GAS6","RPL27A","ZC3HAV1","RPS9","RPL35A","GAPT","RPL39","RPL36A","FOSL2","RPS3","JUND","S100A8","RPL18A","RPL21","RPL37","ARID5A","RPS11","NR4A3","RPL37A","ISG15","SLC7A5","EIF1","RPL4","CCL3","RPS23","SRSF7","RPLP2","RPL7","EIF3L","RPS2","IER2","RPS4X","RPL23A","RPL30","S100A9","SIK1","RPS27","HES4","RPL29","RPL9","RPL35","RPL22","RPL26","RPS16","RPS28","DNAJB9","MYBL2","RPL34","RPL27","RPS10","RPS24","VIM","RPL23","TYROBP","RPS7","RPS15A","RPS13","CD4","RPS6","RPS14","SPIB","RPL17","KLF10","PTPRS","RPL32","GPX1","SERPINF1","NR4A1","RPL10A","RPL5","RPL13","EZR","MAP2K3","MIR24-2","MYADM","APP","NR3C1","WNT10A","RPL3","EEF1B2","AGPAT2","LMNA","SKIL","SNHG7","TSC22D3","ERP29","JUNB","RGCC","ARID3A","RBM38","MZB1","RPS8","RPS18","THBD","COMMD6","RPL6","XBP1","NR4A2","TGFB1","NCF1","RPS3A","PLAC8","NFKBID","ZFP36L2","CYBA","CRIP1","IRF7","FOSB","IRF8","KDM6B","PLD4","ZNF331","FCER1G","LILRA4","BTG2","PPP1R15A","RGS2","ZFP36","PTGDS"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Female Pregnanacy- GO
DoHeatmap(cluster.average_24hCA,features = c("UCP2","TIMP1","IL1B","ARHGDIB","CTSB","RPL29","JUNB","THBD","FOS","FOSB","RGS2"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#rRNA processing- GO
DoHeatmap(cluster.average_24hCA,features =c ("RPS19", "RPS15", "RPL14", "RPS25", "FBL", "RPL11", "RPS21", "RPL7A", "RPL35A", "RPS27", "RPL7", "RPS16", "RPS7", "RPS28", "RPL35", "RPL27", "RPS6", "RPS14", "RPL26", "RPS24", "RPL5", "RPL10A", "RPS8"),group.by='sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#KEGG
#Glycolysis / Gluconeogenesis-KEGG
DoHeatmap(cluster.average_24hCA,features = c("ALDH2", "ALDOA", "ENO1" ,"GAPDH" ,"LDHA" ,"PFKL" ,"PGAM1","PGK1","PKM","TPI1"),group.by='sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Trepenoid backbone biosynthesis-KEGG
DoHeatmap(cluster.average_24hCA,features=c("ACAT2","FDPS","HMGCS1", "IDI1" ,"MVD"),group.by = 'sex',disp.min = -1,disp.max = 1) +scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Covid 19 -Kegg
DoHeatmap(cluster.average_24hCA,features=c("FAU", "IL1B", "ISG15", "MX1", "MX2", "NFKBIA", "PRKCB", "RPL10", "RPL10A", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL3", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36A", "RPL37", "RPL37A", "RPL38", "RPL39", "RPL4", "RPL41", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL9", "RPLP0", "RPLP1", "RPLP2", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS18", "RPS19", "RPS2", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS27", "RPS27A", "RPS28", "RPS29", "RPS3", "RPS3A", "RPS4X","RPS4Y1","RPS5","RPS6", "RPS7" ,"RPS8", "RPS9", "RPSA"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Valine, leucine and isoleucine degradation-KEGG
DoHeatmap(cluster.average_24hCA,features=c("ACAT2", "ALDH2", "HMGCS1" ,"IL4I1"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
```

```{r}
#Pathways enriched at 3h for Pseudomonas when compared with UT
#regulation of lymphocyte chemotaxis
DoHeatmap(cluster.average_3hPA,features = c("CXCL10","CCL2","CCL4","CCL3"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))

#NK cell chemotaxix -GO
DoHeatmap(cluster.average_3hPA,features =c ("CCL2","CCL4","CCL3"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))

#Calcium Ion-GO
DoHeatmap(cluster.average_3hPA,features = c("CXCL10","CCL2","VMP1","CCL4","CCL3"),group.by = 'sex')

#Cytoplasmic Translation-GO
DoHeatmap(cluster.average_3hPA,features=c("SH3BGRL", "YBX1", "RPS26", "RPL35", "RPL41", "RPL9", "RPS25", "RPS5", "RPL18A", "RPL22", "RPL30", "RPS4X", "RPL11", "RPS6", "RPL21", "RPL35A", "RPL37A", "RPS3", "RPL14", "RPL34", "RPS27", "RPL13", "RPL31", "RPL32", "RPS14", "RPL27", "RPL23", "RPS15", "RPS18", "RPS9", "EIF2S3", "RPS8", "RPS13", "RPL12", "RPS15A", "RPS24", "RPL18", "RPS28", "RPS23", "EIF3D", "RPL27A", "RPL10", "RPS2", "RPL3", "RPL26", "RPS7", "RPL7", "EIF3E", "RPLP0", "RPL15", "RPL13A", "RPL23A", "RPS3A", "RPL7A", "EIF4B", "RPL10A", "RPL4", "RPL6", "RPL5", "EIF3F", "EIF3H", "EIF3L", "RPSA", "MIF4GD"),group.by='sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))

#Cold Induced thermogenesis-GO
DoHeatmap(cluster.average_3hPA,features =c("IRF4","CXCR4","KDM6B","UCP2","PLAC8"),group.by = 'sex',disp.min = -1,disp.max = 1)

#Dicarboxylic acid metabolic Process-GO
DoHeatmap(cluster.average_3hPA,features = c("SLC7A11","GCLM","MTHFD2","KMO"	),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Neuron differentiation-GO
DoHeatmap(cluster.average_3hPA,features = c("ITM2C", "RAC2", "CLN8", "BCL11A", "LST1", "PHACTR1", "BTG2", "IER2", "S100A9", "SERPINF1", "KLF4", "NR4A2", "SOX4", "ID1", "CXCR4", "RGS2"),group.by = 'sex',disp.min = -1,disp.max = 1) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))

#Toll-like receptor signaling pathway-KEGG
DoHeatmap(cluster.average_3hPA,features = c("CCL3", "CCL4", "CXCL10", "IL1B", "NFKB1", "STAT1"),group.by='sex') + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Viral protein interaction with cytokine and cytokine receptor-KEGG
DoHeatmap(cluster.average_3hPA,features = c("CCL2", "CCL3", "CCL4", "CCR7", "CXCL10", "CXCL3", "IL2RA" ,"TNFSF10"),group.by = 'sex') + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#TNF signaling pathway-KEGG
DoHeatmap(cluster.average_3hPA,features = c("CCL2", "CCL3", "CCL4", "CCR7", "CXCL10","CXCL3", "IL2RA", "TNFSF10"),group.by='sex',disp.min = -2,disp.max = 2) + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))

```
```{r}
#Pathways enriched at 24h for Pseudomonas aeruginosa compared with UT


#GO
#MHC class II protein complex-GO
DoHeatmap(cluster.average_24hPA,features = c("HLA-DQA1", "HLA-DPA1", "HLA-DRB1", "HLA-DRA", "HLA-DPB1", "HLA-DQB1", "HLA-DQA2", "HLA-DOA", "HLA-DRB5"),group.by = 'sex') + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#Sterol biosynthetic process-GO
DoHeatmap(cluster.average_24hPA,features =c("SQLE", "FDPS", "HMGCS1", "FDFT1", "MSMO1", "INSIG1", "MVD", "IDI1", "LSS", "HMGCR", "EBP", "DHCR24", "ACLY"),group.by = 'sex') + scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))
#peptide biosynthetic process-GO
DoHeatmap(cluster.average_24hPA,features=c("RPS29", "FAU", "MGST2", "EIF4B", "RPL38", "RPS15", "RPS20", "DNAJC3", "CIRBP", "RPS4Y1", "TCOF1", "RPS25", "EEF1D", "RPL12", "RPS21", "RPS5", "RPL18", "RPLP2", "RPL36A", "RPL37A", "MKNK2", "RPL7A", "EEF1A1", "RPS11", "MIF4GD", "PABPC1", "RPL37", "RPS27A", "RPLP0", "RBM3", "RPL11", "RPL23", "RPS16", "SERP1", "RPSA", "RPL27A", "RRBP1", "RPS9", "TMED2", "RPL39", "RPL13A", "RPL35A", "RPS23", "RPS10", "EIF1", "RPL10", "RPL19", "VIM", "RPL31", "RPL30", "RPS3", "RPL22", "RPL15", "RPL17", "RPL18A", "RPS28", "RPL29", "RPL27", "RPS15A", "RPL23A", "RPS2", "RPL21", "RPL35", "RPL4", "RPS27", "EEF2", "RPS7", "RPL26", "RPS4X", "RPL9", "RPS14", "EIF3F", "RPL32", "RPS6", "RPS24", "RPL7", "RPL5", "EEF1B2", "RPS13", "RPL10A", "RPL34", "RPL13", "EIF3L", "RPS8", "APP", "RPS18", "RPL3", "RPL6", "RPS3A", "ZFP36L2", "BTG2", "PPP1R15A", "RGS2", "ZFP36"),group.by = 'sex')+ scale_fill_gradientn(colors = brewer.pal(1000, "RdBu"),
                       oob = scales::oob_squish_any,
                       limits = c(-2, 2))

```


```{r}
#Comparing GO between UT, 3hCA, 24hCA
#First pathway
#Nitric Oxide Metabolic Process
heatmap_NO_UT <- DoHeatmap(cluster.average_UT, features = c("TSPO", "KLF4", "HSP90AA1", "CD36", "IL1B", "AIF1", "SOD2"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_NO_UT
heatmap_NO_3hCA <- DoHeatmap(cluster.average_3hCA, features = c("TSPO", "KLF4", "HSP90AA1", "CD36", "IL1B", "AIF1", "SOD2"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_NO_24hCA <- DoHeatmap(cluster.average_24hCA, features = c("TSPO", "KLF4", "HSP90AA1", "CD36", "IL1B", "AIF1", "SOD2"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_NO_24hCA
grid.arrange(heatmap_NO_UT, heatmap_NO_3hCA, heatmap_NO_24hCA, ncol = 3)
```
```{r}
#Carbon Metabolic Process
heatmap_C_UT <- DoHeatmap(cluster.average_UT, features = c("GRSF1","OCIAD1","RAB21","SOCS1","PTGDS","CSF2RB","MARCH1"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_C_3hCA <- DoHeatmap(cluster.average_3hCA, features = c("GRSF1","OCIAD1","RAB21","SOCS1","PTGDS","CSF2RB","MARCH1"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_C_24hCA <- DoHeatmap(cluster.average_24hCA, features = c("GRSF1","OCIAD1","RAB21","SOCS1","PTGDS","CSF2RB","MARCH1"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_C_UT, heatmap_C_3hCA, heatmap_C_24hCA, ncol = 3)

```
```{r}
#Glycolysis
heatmap_Gly_UT <- DoHeatmap(cluster.average_UT, features = c("GAPDH","ADH5","ADPGK","PGK1","LDHA","ALDH3B1","MINPP1","PGAM1","PFKM","PFKL","ALDH9A1","ALDH2","PKM","GPI","ALDH7A1","PGAM4","FBP1","ALDOA","ENO2","LDHB","ALDH1B1","HK3"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Gly_3hCA <- DoHeatmap(cluster.average_3hCA, features = c("GAPDH","ADH5","ADPGK","PGK1","LDHA","ALDH3B1","MINPP1","PGAM1","PFKM","PFKL","ALDH9A1","ALDH2","PKM","GPI","ALDH7A1","PGAM4","FBP1","ALDOA","ENO2","LDHB","ALDH1B1","HK3"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Gly_24hCA <- DoHeatmap(cluster.average_24hCA, features = c("GAPDH","ADH5","ADPGK","PGK1","LDHA","ALDH3B1","MINPP1","PGAM1","PFKM","PFKL","ALDH9A1","ALDH2","PKM","GPI","ALDH7A1","PGAM4","FBP1","ALDOA","ENO2","LDHB","ALDH1B1","HK3"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Gly_UT, heatmap_Gly_3hCA, heatmap_Gly_24hCA, ncol = 3)

```
```{r}
#TCA Cycle

heatmap_TCA_UT <- DoHeatmap(cluster.average_UT, features = c("PDHB","MDH2","OGDH","SDHC","FH","DLAT","ACO2","SDHA","IDH3G","ACLY","CS","SUCLG1","DLD","SDHD","IDH3A"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_TCA_3hCA <- DoHeatmap(cluster.average_3hCA, features = c("PDHB","MDH2","OGDH","SDHC","FH","DLAT","ACO2","SDHA","IDH3G","ACLY","CS","SUCLG1","DLD","SDHD","IDH3A"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_TCA_24hCA <- DoHeatmap(cluster.average_24hCA, features = c("PDHB","MDH2","OGDH","SDHC","FH","DLAT","ACO2","SDHA","IDH3G","ACLY","CS","SUCLG1","DLD","SDHD","IDH3A"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_TCA_UT, heatmap_TCA_3hCA, heatmap_TCA_24hCA, ncol = 3)
```
```{r}
#Oxidative Phosphorylation

heatmap_OP_UT <- DoHeatmap(cluster.average_UT, features = c("NDUFA11","COX5B","ATP5ME","COX8A","UQCRQ","NDUFB1","NDUFA1","ATP5MC2","NDUFS6","COX7C","UQCR11","COX6B1","NDUFB3","ATP5MG","UQCRH","COX7B2"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_OP_3hCA <- DoHeatmap(cluster.average_3hCA, features = c("NDUFA11","COX5B","ATP5ME","COX8A","UQCRQ","NDUFB1","NDUFA1","ATP5MC2","NDUFS6","COX7C","UQCR11","COX6B1","NDUFB3","ATP5MG","UQCRH","COX7B2"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_OP_24hCA <- DoHeatmap(cluster.average_24hCA, features = c("NDUFA11","COX5B","ATP5ME","COX8A","UQCRQ","NDUFB1","NDUFA1","ATP5MC2","NDUFS6","COX7C","UQCR11","COX6B1","NDUFB3","ATP5MG","UQCRH","COX7B2"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_OP_UT, heatmap_OP_3hCA, heatmap_OP_24hCA, ncol = 3)
```
```{r}
#Farnesyl Diphosphate Metabolic Process
heatmap_FDP_UT <- DoHeatmap(cluster.average_UT, features = c("FDPS","HMGCS1","FDFT1"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_FDP_3hCA <- DoHeatmap(cluster.average_3hCA, features = c("FDPS","HMGCS1","FDFT1"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_FDP_24hCA <- DoHeatmap(cluster.average_24hCA, features = c("FDPS","HMGCS1","FDFT1"), group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_FDP_UT, heatmap_FDP_3hCA, heatmap_FDP_24hCA, ncol = 3)
```
```{r}
#Cellular Nitrogen Compound Biosynthetic Process
genes<-c("RPS4Y1", "TCOF1", "DDX5", "BCL11A", "RPS25", "EEF1D", "RPL12", "S100A9", "RPS21", "ELF2", "RPS5", "RPL18", "RPLP2", "RPL36A", "RPL37A", "RBCK1", "IRF1", "MKNK2", "RPL7A", "EEF1A1", "ADA", "JMJD1C", "PRKCB", "RPS11", "MIF4GD", "PABPC1", "ZFAT", "RPL37", "BAZ1A", "RPS27A", "RPLP0", "CLN8", "RBM3", "NADK", "S100A8", "RPL11", "PIM2", "CSRNP1", "RPL23", "IFI16", "RPS16", "TRIM22", "USF2", "SERP1", "RPSA", "RPL27A", "RRBP1", "RPS9", "TMED2", "ZNF791", "RPL39", "RPL13A", "RPL35A", "RPS23", "POLB", "RPS10", "GAS6", "EIF1", "NR4A3", "RPL10", "RPL19", "VIM", "CCNL1", "RPL31", "RPL30", "RPS3", "RPL22", "RPL15", "SP110", "JUND", "CDKN2D", "RPL17", "RPL18A", "RPS28", "DCK", "RPL29", "ARID5A", "SPIB", "RPL27", "IRF4", "RPS15A", "FOSL2", "RPL23A", "RPS2", "RPL21", "RPL35", "SIK1", "CD4", "RPL4", "RPS27", "EEF2", "RPS7", "RPL26", "NFKBIA", "RPS4X", "RPL9", "CXXC5", "ODC1", "IER2", "RPS14", "EIF3F", "RPL32", "RPS6", "RPS24", "RPL7", "RPL5", "HES4", "EEF1B2", "RPS13", "RPL10A", "RPL34", "RPL13", "SMPD3", "NR3C1", "EZR", "KLF10", "MYBL2", "SKIL", "NR4A1", "COMMD6", "EIF3L", "TSC22D3", "RPS8", "MAP2K3", "APP", "RGCC", "RPS18", "GUK1", "RPL3", "JUNB", "RPL6", "TGFB1", "ARID3A", "NR4A2", "XBP1", "RPS3A", "ZFP36L2", "NCF1", "FOSB", "ZNF331", "BTG2", "KDM6B", "PPP1R15A", "IRF7","IRF8","RGS2","PLAC8","ZFP36")
heatmap_NCBUT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_NCB3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_NCB24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_NCBUT, heatmap_NCB3hCA, heatmap_NCB24hCA, ncol = 3)
```
```{r}
#Cellular amide metabolic process
genes<-c("RPS4Y1","TCOF1","RPS25","EEF1D","RPL12","TAPBP","RPS21","RPS5","RPL18","RPLP2","RPL36A","RPL37A","MKNK2","RPL7A","EEF1A1","ADA","RPS11","MIF4GD","PABPC1","RPL37","RPS27A","RPLP0","CLN8","RBM3","RPL11","RPL23","RPS16","SPCS1","SERP1","RPSA","RPL27A","RRBP1","RPS9","TMED2","RPL39","RPL13A","RPL35A","RPS23","RPS10","EIF1","RPL10","RPL19","VIM","RPL31","RPL30","RPS3","RPL22","RPL15","RPL17","RPL18A","RPS28","RPL29","RPL27","RPS15A","RPL23A","RPS2","RPL21","RPL35","RPL4","RPS27","EEF2","RPS7","RPL26","RPS4X","RPL9","RPS14","EIF3F","RPL32","RPS6","RPS24","RPL7","RPL5","EEF1B2","RPS13","RPL10A","RPL34","RPL13","SMPD3","EIF3L","RPS8","APP","RPS18","RPL3","RPL6","RPS3A","ZFP36L2","BTG2","PPP1R15A","RGS2","ZFP36")

heatmap_AM_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_AM_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_AM_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_AM_UT, heatmap_AM_3hCA, heatmap_AM_24hCA, ncol = 3)
```
```{r}
#Peptide Metabolic Process
genes <- c("RPS29","FAU","MGST2","EIF4B","RPL38","RPS15","RPS20","DNAJC3","CIRBP","RPS4Y1","TCOF1","RPS25","EEF1D","RPL12","TAPBP","RPS21","RPS5","RPL18","RPLP2","RPL36A","RPL37A","MKNK2","RPL7A","EEF1A1","RPS11","MIF4GD","PABPC1","RPL37","RPS27A","RPLP0","RBM3","RPL11","RPL23","RPS16","SPCS1","SERP1","RPSA","RPL27A","RRBP1","RPS9","TMED2","RPL39","RPL13A","RPL35A","RPS23","RPS10","EIF1","RPL10","RPL19","VIM","RPL31","RPL30","RPS3","RPL22","RPL15","RPL17","RPL18A","RPS28","RPL29","RPL27","RPS15A","RPL23A","RPS2","RPL21","RPL35","RPL4","RPS27","EEF2","RPS7","RPL26","RPS4X","RPL9","RPS14","EIF3F","RPL32","RPS6","RPS24","RPL7","RPL5","EEF1B2","RPS13","RPL10A","RPL34","RPL13","EIF3L","RPS8","APP","RPS18","RPL3","RPL6","RPS3A","ZFP36L2","BTG2","PPP1R15A","RGS2","ZFP36")
heatmap_Pep_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Pep_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Pep_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Pep_UT, heatmap_Pep_3hCA, heatmap_Pep_24hCA, ncol = 3)
```
```{r}
#rRNA metabolic process
genes <- c("RPL14", "GTF3A", "RPS15", "TCOF1", "RPS25", "RPS21", "RPL7A", "RPL11", "FBL", "RPS16", "RPL35A", "RPS28", "RPL27", "RPL35", "RPS27", "RPS7", "RPL26", "ISG20", "RPS14", "RPS6", "RPS24", "RPL7", "RPL5", "RPL10A", "RPS8")
heatmap_rRNA_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_rRNA_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_rRNA_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_rRNA_UT, heatmap_rRNA_3hCA, heatmap_rRNA_24hCA, ncol = 3)

```
```{r}
#Apoptotic process
genes <- c("ADAR", "SLC25A6", "POLB", "SIVA1", "GAS6", "TPT1", "RPL10", "TNFRSF21", "CCNL1", "RPS3", "CDKN2D", "CTSB", "PMAIP1", "PIM3", "SIK1", "RPS7", "RPL26", "NFKBIA", "CCL3", "TRAF4", "LMNA", "RNF130", "RPS6", "UCP2", "B4GALT1", "NR3C1", "MYBL2", "SKIL", "NR4A1", "GADD45B", "CXCR4", "TSC22D3", "HSP90B1", "APP", "RGCC", "DNASE1L3", "TGFB1", "ERP29", "NR4A2", "MCL1", "XBP1", "IFIT2", "HERPUD1", "RPS3A", "MX1", "ITM2C", "NCF1", "MZB1", "CRIP1", "TCL1A", "BTG2", "PTCRA", "PPP1R15A", "IRF7", "PLAC8", "ZFP36")

heatmap_apop_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_apop_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_apop_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_apop_UT, heatmap_apop_3hCA, heatmap_apop_24hCA, ncol = 3)
```
```{r}
#Ferroptosis
genes <- c("FTH1", "FTL", "GCLM" ,"SAT1", "SLC7A11")
heatmap_ferrop_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_ferrop_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_ferrop_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_ferrop_UT, heatmap_ferrop_3hCA, heatmap_ferrop_24hCA, ncol = 3)
```
```{r}
#Carboxylic acid metabolic process
genes <- c("PKM", "GSTP1", "ENO1", "TPI1", "GAPDH", "GPX4", "ACAT2", "IL4I1", "SCD", "MSMO1", "INSIG1", "LDHA", "MTHFD2", "DBI", "ALDOA", "ADTRP", "KYNU", "PGD", "FADS2", "IDH2", "PGAM1", "VDAC1", "FADS1", "PGK1", "ELOVL5", "PFKL", "ACLY", "STARD4", "ACSL4", "ARL2")

heatmap_Caboxy_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Caboxy_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Caboxy_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Caboxy_UT, heatmap_Caboxy_3hCA, heatmap_Caboxy_24hCA, ncol = 3)
```
```{r}
#Cholestrol biosynthetic process
genes <- c("FDPS","HMGCS1","FDFT1","MSMO1","INSIG1","MVD","IDI1","LSS","HMGCR","EBP","DHCR24","ACLY")
heatmap_Chole_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))
heatmap_Chole_UT
heatmap_Chole_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Chole_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))
heatmap_Chole_24hCA
grid.arrange(heatmap_Chole_UT, heatmap_Chole_3hCA, heatmap_Chole_24hCA, ncol = 3)

```
```{r}
#Cell Migration-GO
genes <- c("SMPD3", "CCL5", "AIF1", "B4GALT1", "CORO1A", "MAP2K3", "RAC2", "SELPLG", "PHACTR1", "SERPINF1", "GPX1", "LDLRAD4", "DUSP1", "PPBP", "PYCARD", "ACTG1", "KLF4", "ARHGDIB", "S100A8", "NR4A2", "NR4A1", "S100A9", "MYADM", "ID1", "RGCC", "CXCR4", "FCER1G")

heatmap_Migra_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Migra_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Migra_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Migra_UT, heatmap_Migra_3hCA, heatmap_Migra_24hCA, ncol = 3)
```
```{r}
#Valine, leucine and isoleucine degradation
genes <- c("ACAT2", "ALDH2", "HMGCS1" ,"IL4I1")
heatmap_Stero_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Stero_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Stero_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Stero_UT, heatmap_Stero_3hCA, heatmap_Stero_24hCA, ncol = 3)
```
```{r}
#Steroid biosynthesis Pathway
genes <- c("SQLE","FDPS","HMGCS1","FDFT1","MSMO1","INSIG1","NFKB1","MVD","IDI1","LSS","HMGCR","EBP","DHCR24","ACLY","STARD4")
heatmap_Stero_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Stero_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Stero_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Stero_UT, heatmap_Stero_3hCA, heatmap_Stero_24hCA, ncol = 3)
```
```{r}
#Calcium signalling Pathway
genes <- c("CALM2", "CXCR4", "GNA15" ,"ORAI2" ,"PRKCB" ,"SLC25A5" ,"SLC25A6", "VEGFB")
heatmap_Cal_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Cal_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Cal_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(200, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Cal_UT, heatmap_Cal_3hCA, heatmap_Cal_24hCA, ncol = 3)
```
```{r}
#Negative Regulation of Immune System Process-GO
genes <- c("IL2RA","IL7R","ISG15","AHR","IFI16","PARP14","IL4I1","ADTRP","LPXN","PAG1","CCL3","HLA-A","ADAR","IRF1","HLA-F","RUNX3","CD47","USP18","OAS1","HLA-E","CASP3","CCL2","HLA-DRB1","SERPINB9","SAMSN1","BTN2A2","NMI","TNFAIP3","ELF1")
heatmap_Negati_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Negati_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Negati_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Negati_UT, heatmap_Negati_3hCA, heatmap_Negati_24hCA, ncol = 3)

```
```{r}
#Defense Response to Virus-GO
genes <- c("CXCL10","ISG15","PMAIP1","BIRC3","ISG20","IFIT2","MX1","IFIT3","IFI44L","IFI16","PLSCR1","STAT1","OAS2","EIF2AK2","NT5C3A","IFIH1","ADAR","HERC5","IRF1","IFI6","TANK","OAS1","MX2")

heatmap_Def_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Def_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Def_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Def_UT, heatmap_Def_3hCA, heatmap_Def_24hCA, ncol = 3)

```
```{r}
#Regulation of Lympocyte Chemotaaxis-GO
genes <- c("CXCL10","CCL2","CCL4","CCL3")

heatmap_Lym_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Lym_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Lym_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Lym_UT, heatmap_Lym_3hCA, heatmap_Lym_24hCA, ncol = 3)
```
```{r}
#IL-6 production-GO
genes <- c("ARRB2","TNFAIP3","PYCARD","TYROBP","CD36","IL1B","AIF1")
heatmap_IL6_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_IL6_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_IL6_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_IL6_UT, heatmap_IL6_3hCA, heatmap_IL6_24hCA, ncol = 3)
```
```{r}
#IL-17 production-GO
genes <-c("RFTN1","SLC7A5","IFNG" ,"LILRB2", "AIF1",  "AKIRIN2" ,"CD36", "EREG", "ZBTB20", "HYAL2" ,"APP" ,"TSLP" ,"F2R")
heatmap_IL17_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_IL17_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_IL17_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_IL17_UT, heatmap_IL17_3hCA, heatmap_IL17_24hCA, ncol = 3)
```
```{r}
#Cell death-GO
genes <- c("CXCL10","OPN3","CTSC","TXN","BCL2A1","IL2RA","ENDOG","IL7R","CFLAR","BTG1","GCLM","SLC7A11","NFKB1","CDKN1A","PMAIP1","BIRC3","IFIT2","MX1","IFIT3","SQSTM1","AHR","IFI16","PLSCR1","TNFSF10","PTPN1","BID","ATF5","STAT1","TNFAIP8","EIF2AK2","AKR1C3","KMO","MARCKS","CCR7","ARL6IP5","GZMB","CCL3","NINJ1","SRGN","XAF1","STK4","ADAR","EMP3","IRF1","IFI6","HLA-F","TAF9","CASP3","TRAF1","PKM","SP100","DNAJA1","CCL2","CD38","ACTN4","FAS","NUP62","ATF4","SERPINB9","RTN4","RBCK1","HSP90B1")
heatmap_death__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_death__3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_death__24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_death__UT, heatmap_death__3hCA, heatmap_death__24hCA, ncol = 3)
```
```{r}
#Auto immune thyroid desease-GO
genes <- c("FAS", "GZMB", "HLA-A", "HLA-DQA1", "HLA-DQA2", "HLA-DRB1", "HLA-E", "HLA-F")
heatmap_auto_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_auto_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_auto_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_auto_UT, heatmap_auto_3hCA, heatmap_auto_24hCA, ncol = 3)
```
```{r}
#MHC class II peptide complex assembly-GO
genes <- c("HLA-DQA1","HLA-DPA1","HLA-DRB1","HLA-DRA","HLA-DPB1","HLA-DQB1","HLA-DQA2","HLA-DOA","HLA-DRB5")
heatmap_MHC_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_MHC_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_MHC_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_MHC_UT, heatmap_MHC_3hCA, heatmap_MHC_24hCA, ncol = 3)
```
```{r}
#T cell mediated immunity-GO
genes <- c("IL7R","MYO1G","IL4I1","CD70","RFTN1","HLA-DRB1","HLA-DRA","CTSH","HLA-A","MALT1")
heatmap_T_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_T_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_T_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_T_UT, heatmap_T_3hCA, heatmap_T_24hCA, ncol = 3)
```
```{r}
#RIG -I like receptor signaling pathway-GO
genes <- c("ADAR", "CXCL10", "IFIH1", "ISG15", "NFKB1", "TANK")
heatmap_Rig_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Rig_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Rig_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Rig_UT, heatmap_Rig_3hCA, heatmap_Rig_24hCA, ncol = 3)
```
```{r}
#NOD like receptor signaling pathway-GO
genes <- c("BIRC3", "CCL2" ,"CXCL3" ,"IFI16" ,"IL1B" ,"NFKB1", "OAS2", "STAT1", "TXN")
heatmap_NOD_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_NOD_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_NOD_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_NOD_UT, heatmap_NOD_3hCA, heatmap_NOD_24hCA, ncol = 3)
```
```{r}
#Necropotosis-GO
genes <- c("BID", "BIRC3", "CFLAR", "EIF2AK2", "FTH1", "FTL", "IL1B", "SQSTM1", "STAT1", "TNFSF10")
heatmap_Necro_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Necro_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Necro_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Necro_UT, heatmap_Necro_3hCA, heatmap_Necro_24hCA, ncol = 3)
```
```{r}
#Cytokine-Cytokine receptor Interaction
genes <-c("CCL2","CCL3","CCL4","CCR7","CXCL10","CXCL3","IL1B","IL1RN","IL2RA","IL7R","TNFRSF4","TNFSF10")
heatmap_cytocyto_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_cytocyto_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_cytocyto_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_cytocyto_UT, heatmap_cytocyto_3hCA, heatmap_cytocyto_24hCA, ncol = 3)
```
```{r}
#Intestinal immune network for IgA production
genes <-c("BCL6","CD79A","CD79B","CD83","CXCL12","CXCL13","CXCR4","IL10","IL21","IL21R","IL6","LTB","LTBR","TNFRSF13B")
heatmap_intestine_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_intestine_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_intestine_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_intestine_UT, heatmap_intestine_3hCA, heatmap_intestine_24hCA, ncol = 3)
```
```{r}
#Ribosome
genes <-c("MRPL23", "RPL10", "RPL10A", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL26", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL3", "RPL30", "RPL31", "RPL32", "RPL35", "RPL35A", "RPL36", "RPL37", "RPL37A", "RPL38", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPLP0", "RPLP2", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS18", "RPS19", "RPS2", "RPS20", "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS28", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPSA", "UBA52")

heatmap_Ribosome_UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Ribosome_3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Ribosome_24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Ribosome_UT, heatmap_Ribosome_3hCA, heatmap_Ribosome_24hCA, ncol = 3)
```
```{r}
#Visualizing the common pathways in the three time points
#response to cytokine-GO
source("C:/Users/vnikh/PycharmProjects/BINP37/Gene_Formatting.R")
genes_24h <- c("BCLAF1", "SIRPA", "IL10RB", "IL10RA", "CREB1", "PRPF8", "TMSB4X", "SPI1", "LGALS9", "NDUFA13", "SH2B3", "ZFP36", "TAF9", "KIF5B", "TNFRSF18", "ENDOG", "MBP", "CD44", "CXCL1", "TPR", "ETV3", "CXCL3", "PYCARD", "CD40", "CD58", "UBE2G2", "TFRC", "PTPN2", "BIRC3", "PARP9", "RAD23B", "ACTR3", "TRIM44", "RNMT", "HCLS1", "IL13RA1", "TP53", "CYLD", "LAMP3", "CCL22", "PTPRC", "IFI16", "TANK", "ZFP36L1", "CSF2RA", "IFNGR1", "HCK", "IL7R", "CDC42", "TNFRSF4", "NRP2", "IFNAR1", "SETD2", "LAPTM5", "HIF1A", "EBI3", "GSN", "CSF2RB", "SAMHD1", "CD70")

genes_3h <- format_gene_list("CXCL10/CCL3/ISG15/IFIT2/IFIT3/IFIH1/CXCL2/MX1/OAS1/IL1B/SOCS1/RNMT/USP18/TNFSF13B/XAF1/PARP14/ZFP36/IRF1/EIF2AK2/MX2/NMI/LGALS9/XRCC5/XBP1/PNPT1/SP100/B3GNT2/MIR142/IL3RA/TNFRSF4/PTPRC/GCLM/FCER1G/STAT3/BST2/CD38/ENDOG/NDUFA13/TRIM44/IFI16/CCL4/MTF2/LILRA4/MYD88/IRF8/ADAR/APPL1/NUB1/PML/DHX9")

genes_UT <- format_gene_list("CXCL10/IL2RA/ENDOG/IL7R/ISG15/GCLM/NFKB1/TNFRSF4/BIRC3/IFIT2/MX1/IFIT3/CCL4/IFI16/PLSCR1/PTPN1/STAT1/PARP14/NFAT5/OAS2/EIF2AK2/KMO/IFIH1/CCR7/CCL3/XAF1/B3GNT2/ADAR/RBMX/IRF1/TAF9/TANK/CD47/USP18/OAS1/MX2/CASP3/TRAF1/SP100/CCL2/CD38/ACTN4/FAS")

#Select the unique genes from the three time points
genes <- unique(c(genes_24h, genes_3h, genes_UT))

heatmap_Cytokine__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Cytokine__3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Cytokine__24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Cytokine__UT, heatmap_Cytokine__3hCA, heatmap_Cytokine__24hCA, ncol = 3)

```
```{r}
#Negative Regulation of Immune system process-GO

genes_UT <- format_gene_list("IL2RA/IL7R/ISG15/AHR/IFI16/PARP14/IL4I1/ADTRP/LPXN/PAG1/CCL3/HLA-A/ADAR/IRF1/HLA-F/RUNX3/CD47/USP18/OAS1/HLA-E/CASP3/CCL2/HLA-DRB1/SERPINB9/SAMSN1/BTN2A2/NMI/TNFAIP3/ELF1")

genes_3h <- format_gene_list("CCL3/ISG15/OAS1/SOCS1/USP18/PARP14/IRF1/HLA-E/NMI/LGALS9/CD55/HLA-A/HLA-B/PTPRC/BST2/BTN2A2/IFI16/IL4I1/ADAR/APPL1")

genes_24h <- format_gene_list("CD46/PTEN/BTN2A2/TNFAIP3/PARP1/RUNX3/CUEDC2/TYROBP/PRKAR1A/ARRB2/USP15/TRIB1/SPI1/LGALS9/CNN2/DUSP3/SPN/ID2/HLA-DRB1/AHR/PTPN2/PTGER4/TMBIM6/CD86/PTPRC/IFI16/GPX1/HCK/IL7R/IL4I1/CD55/RASSF5/LAPTM5/HLA-DOA/CST7/SAMHD1")

genes <- unique(c(genes_24h, genes_3h, genes_UT))

heatmap_Negative_Immune__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Negative_Immune__3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Negative_Immune__24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Negative_Immune__UT, heatmap_Negative_Immune__3hCA, heatmap_Negative_Immune__24hCA, ncol = 3)
```
```{r}
#[3] "regulation of cold-induced thermogenesis"
genes_UT <- format_gene_list("IRF4/KDM6B/UCP2/CXCR4/PLAC8")

genes_3h <- format_gene_list("ADIPOR1/ATF4/RHEB/NR1H2/CXCR4/UCP2/FABP5/KDM6B/IRF4")

genes_24h <- format_gene_list("PLAC8/IRF4/FABP5/CXCR4/GNAS/DYNC1H1")

genes <- unique(c(genes_24h, genes_3h, genes_UT))

heatmap_Cold__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Cold__3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_Cold__24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_Cold__UT, heatmap_Cold__3hCA, heatmap_Cold__24hCA, ncol = 3)
```
```{r}
#regulation of immune system process
genes_UT<- format_gene_list("CXCL10/CTSC/IL2RA/IL7R/ISG15/SLC7A11/CDKN1A/TNFRSF4/BIRC3/AHR/CCL4/IFI16/PLSCR1/PTPN1/BID/STAT1/PARP14/IL4I1/ADTRP/MYO1G/CD83/LPXN/IFIH1/PAG1/CCR7/CCL3/NINJ1/HLA-A/ADAR/IRF1/HLA-F/RUNX3/CD47/USP18/OAS1/HLA-E/CASP3/HLA-DQA1/CCL2/CD38/HLA-DRB1/SERPINB9/RTN4/RBCK1/SAMSN1/HLA-DQA2/ZNFX1/SLC15A4/BTN2A2/NMI/TNFAIP3/RFTN1")

genes_3h <- format_gene_list("CXCL10/CCL3/ISG15/IFIH1/CTSC/OAS1/IL1B/SOCS1/USP18/SLC15A4/TNFSF13B/RBCK1/PARP14/ZFP36/IRF1/HLA-E/NMI/LGALS9/XRCC5/XBP1/IFI35/CD55/MIR142/LGMN/HLA-A/HLA-B/ARID5A/PQBP1/TNFRSF4/PTPRC/NFKBIZ/FCER1G/STAT3/CDC73/BST2/BTN2A2/CD38/JUNB/IFI16/CCL4/CD99/CALR/IL4I1/LILRA4/MYD88/ADAR/APPL1/CDKN1A/B2M/KHDRBS1/ZNFX1/DHX9/TAPBP/SMARCC1/IRF7/ELF1/HLA-F/MSN")

genes_24h <- format_gene_list("PRKDC/NOTCH2/PQBP1/RPS6KA3/HLA-DPA1/CD46/LGALS1/PTEN/BTN2A2/HLA-DPB1/SUPT6H/SH3KBP1/TNFAIP3/PARP1/HLA-DQB1/HSPH1/PRMT1/RUNX3/CUEDC2/HAX1/TYROBP/PRKAR1A/PLSCR1/NR4A3/CALR/ARRB2/SIRPA/USP15/TRIB1/RB1/GPSM3/JUNB/CDKN1A/CREB1/RASSF2/SPI1/CLPTM1/LGALS9/CNN2/SWAP70/BID/DHPS/ZFP36/DUSP3/TNFRSF18/NPLOC4/PPP2R3C/PHB2/EP300/PYCARD/HLA-DMB/EIF6/SPN/CD40/ID2/HLA-DRB1/TFRC/AHR/PTPN2/BIRC3/PARP9/CNOT4/PTGER4/ACTB/PJA2/TMBIM6/CD86/HCLS1/CYLD/PTPRC/SASH3/IFI16/ZFP36L1/GPX1/CTSH/HCK/ICOSLG/UBE2N/HSP90AA1/IL7R/LY96/HLA-DQA2/IL4I1/CD83/CD55/BAX/TNFRSF4/RASSF5/LAPTM5/DDX3X/HIF1A/EBI3/HLA-DOA/CST7/SAMHD1/CD70")


genes <- unique(c(genes_24h, genes_3h, genes_UT))

heatmap_immune_reg__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_immune_reg__3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_immune_reg__24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_immune_reg__UT, heatmap_immune_reg__3hCA, heatmap_immune_reg__24hCA, ncol = 3)
```
```{r}
#CA-
#1-cell surface receptor signaling pathway"

genes_UT <- format_gene_list("ISG15/AGPAT2/MX1/SEMA7A/CSF2RB/NFKBID/TNFRSF4/CCL3/PLSCR1/ADAR/WNT10A/SOCS1/IRF1/SKIL/CD47/EIF5A/PLEK/ATP6AP2/P2RY6/IL6ST/UGCG/RBCK1/STK4/SULF2/YTHDF2/HLA-A/ARHGDIA/HLA-DQB1/RFTN1/LPXN/BIRC3/HSPA5/TRIM44/IL10RA/ICOSLG/NR4A2/ARID4B/BAIAP2/PMAIP1/TNFSF13B/GAS6/RBM4/CDC42/RNF213/STAT6/RAP1A/IRF9/P2RY14/APH1A/LEPROTL1/MAD2L2/NFKBIZ/DGKZ/PMEPA1/STAT1/NR4A3/RB1CC1/CSF2RA/BIRC2/VPS35/DDX5/HCLS1/SUB1/BID/LILRB4/NFKB1/SMPD3/SP100/ZFAND5/SNX6/PTPN2/DDB1")

genes_3h <- format_gene_list("CXCL10/CCL3/ISG15/TNFSF10/IFI6/CXCL2/MX1/ATF3/OAS1/IL1B/SOCS1/RNF213/USP18/TNFSF13B/RBCK1/PARP14/IRF1/PILRB/NMI/SP100/CD55/LGMN/HLA-A/IL3RA/SUB1/IFT57/BAG1/TNFRSF4/USP47/PTPRC/GCLM/NFKBIZ/FCER1G/STAT3/STMN1/CDC73/BTN2A2/CD38/OCIAD1/SUDS3/USP34/NDUFA13/TRIM44/CCL4")

genes_24h <- format_gene_list("TMSB4X/RASSF2/SPI1/NDUFA13/CSNK1D/BID/TRIM33/SH2B3/PFDN5/SIAH2/AKT1S1/OSBPL8/DUSP3/CPEB4/TNFRSF18/CSNK1A1/CD44/CXCL1/PHB2/TCTN3/CXCL3/BAIAP2/EP300/PYCARD/USP34/SPN/CD40/HLA-DRB1/ECE1/PTPN2/RAB14/BIRC3/PARP9/MBD2/MTMR4/PSME3/RAB7A/PDCD6/TRIM44/CD86/HCLS1/IL13RA1/SPRED2/TP53/HSPB1/CYLD/CCL22/PTPRC/TANK/PEA15/CSF2RA/GPX1/IFNGR1/HCK/ICOSLG/UBE2N/IL7R/HIVEP1/LY96/CDC42/CD55/BAX/VPS35/TNFRSF4/NRP2/BCL2A1/IFNAR1/CCNC/CALCRL/LAPTM5/MTSS1/DDX3X/HIF1A/ARHGDIA/EBI3/CSF2RB/PTGIR/SAMHD1/CD70")

genes <- unique(c(genes_24h, genes_3h, genes_UT))

heatmap_cell_sur__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_cell_sur__3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_cell_sur__24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_cell_sur__UT, heatmap_cell_sur__3hCA, heatmap_cell_sur__24hCA, ncol = 3)

```
```{r}
#2 Immune Response
genes_UT <- format_gene_list("ISG15/PTGDS/IFIT2/IFI44L/MX1/CTSC/SEMA7A/CSF2RB/NFKBID/ISG20/GZMB/NINJ1/TNFRSF4/CCL3/HLA-DQA1/PLSCR1/ADAR/TAP1/SOCS1/IRF1/CD47/HLA-F/IL6ST/LTB/RBCK1/YTHDF2/GPR183/HLA-A/HLA-DQB1/RFTN1/LPXN/BIRC3/ICOSLG/XRCC5/TNFSF13B/REL/SLAMF7/PSMA1/CDC42/HLA-B/STAT6/PPP1R14B/TRIM22/MAD2L2/NFKBIZ/CLEC4C/DGKZ/MX2/STAT1/NR4A3/ERCC1/MPEG1/ACTR3/BIRC2/PSMB10/LITAF/TNIP1/LILRB4/SLC7A5/TRAF4/SP100/TAPBP/NFKB2/IFI16/PTPN2")

genes_3h <- format_gene_list("CXCL10/CCL3/ISG15/IFIT2/IFIT3/TNFSF10/LTB/IFIH1/HERC5/CTSC/IFI6/CXCL2/MX1/RGS1/OAS1/IL1B/SOCS1/IFI44L/USP18/SLC15A4/NKG7/TNFSF13B/RBCK1/EXOSC9/PARP14/IRF1/EIF2AK2/MX2/HLA-E/IFI44/NMI/LGALS9/XRCC5/XBP1/IFI35/SP100/ZC3HAV1/CD55/GNLY/TAP1/HLA-A/HLA-B/ARID5A/PQBP1/TNFRSF4/PTPRC/NFKBIZ/FCER1G/STAT3/BST2/MPEG1/BTN2A2/CD38/JUNB/IFI16/CCL4/HLA-C/CD164/IL4I1/LILRA4/MYD88/ISG20/SLAMF7/IRF8/ADAR/APPL1/B2M/NUB1/KHDRBS1/ZNFX1/PML/CNPY3/DHX9/GPR65/TAPBP")

genes_24h <- format_gene_list("DUSP3/KIF5B/ARL8B/GNLY/MBP/NPLOC4/CXCL1/PPP2R3C/PHB2/CXCL3/EP300/DHX15/PYCARD/HLA-DMB/SPN/CD40/CD58/HLA-DRB1/MCOLN2/TFRC/PRDX1/AHR/PTPN2/BIRC3/PARP9/TUBB4B/PTGER4/ACTR3/PJA2/CD86/TP53/CYLD/LAMP3/DBNL/CCL22/PTPRC/SASH3/IFI16/TANK/GPX1/IFNGR1/CTSH/HCK/ICOSLG/UBE2N/HSP90AA1/IL7R/LY96/ALCAM/HLA-DQA2/CDC42/IL4I1/CD83/CD55/BAX/TNFRSF4/IFNAR1/SETD2/LAPTM5/DDX3X/EBI3/GSN/HLA-DOA/CSF2RB/CST7/SAMHD1/CD70/LTB")

genes <- unique(c(genes_24h, genes_3h, genes_UT))

heatmap_immu__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_immu__3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_immu__24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_immu__UT, heatmap_immu__3hCA, heatmap_immu__24hCA, ncol = 3)

```
```{r}
#response to cytokine

genes_UT <- format_gene_list("ISG15/IFIT2/TCL1A/AGPAT2/MX1/CSF2RB/MT2A/TNFRSF4/CCL3/PLSCR1/ADAR/SOCS1/IRF1/SKIL/CD47/EIF5A/IL6ST/UGCG/MIR142/YTHDF2/CREB1/BIRC3/HSPA5/TRIM44/IL10RA/XRCC5/TNFSF13B/GAS6/REL/CDC42/STAT6/ZFAND6/MX2/STAT1/CSF2RA/ACTR3/BIRC2/HCLS1/LILRB4/NFKB1/SMPD3/SP100/IFI16/PTPN2/MBP/HAX1/CYLD/HSP90AB1/NFE2L2/NCL/IL3RA/TPR/SOCS3/RPLP0/SRSF3/KIF5B")

genes_3h <- format_gene_list("CXCL10/CCL3/ISG15/IFIT2/IFIT3/IFIH1/CXCL2/MX1/OAS1/IL1B/SOCS1/RNMT/USP18/TNFSF13B/XAF1/PARP14/ZFP36/IRF1/EIF2AK2/MX2/NMI/LGALS9/XRCC5/XBP1/PNPT1/SP100/B3GNT2/MIR142/IL3RA/TNFRSF4/PTPRC/GCLM/FCER1G/STAT3/BST2/CD38/ENDOG/NDUFA13/TRIM44/IFI16/CCL4/MTF2/LILRA4/MYD88/IRF8/ADAR/APPL1/NUB1/PML/DHX9")

genes_24h <- format_gene_list("BCLAF1/SIRPA/IL10RB/IL10RA/CREB1/PRPF8/TMSB4X/SPI1/LGALS9/NDUFA13/SH2B3/ZFP36/TAF9/KIF5B/TNFRSF18/ENDOG/MBP/CD44/CXCL1/TPR/ETV3/CXCL3/PYCARD/CD40/CD58/UBE2G2/TFRC/PTPN2/BIRC3/PARP9/RAD23B/ACTR3/TRIM44/RNMT/HCLS1/IL13RA1/TP53/CYLD/LAMP3/CCL22/PTPRC/IFI16/TANK/ZFP36L1/CSF2RA/IFNGR1/HCK/IL7R/CDC42/TNFRSF4/NRP2/IFNAR1/SETD2/LAPTM5/HIF1A/EBI3/GSN/CSF2RB/SAMHD1/CD70")

genes <- unique(c(genes_24h, genes_3h, genes_UT))

heatmap_cyto__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_cyto__3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_cyto__24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_cyto__UT, heatmap_cyto__3hCA, heatmap_cyto__24hCA, ncol = 3)
```
```{r}
#Cytokine Mediated Signaling Pathway
genes_UT <- format_gene_list("ISG15/AGPAT2/MX1/CSF2RB/TNFRSF4/CCL3/ADAR/SOCS1/IRF1/EIF5A/IL6ST/UGCG/YTHDF2/BIRC3/TRIM44/IL10RA/TNFSF13B/GAS6")
genes_3h <- format_gene_list("CXCL10/CCL3/ISG15/CXCL2/MX1/OAS1/IL1B/SOCS1/USP18/TNFSF13B/PARP14/IRF1/NMI/SP100/IL3RA/TNFRSF4/PTPRC/FCER1G/STAT3/TRIM44/CCL4/LILRA4/MYD88/ADAR/APPL1")
genes_24h <- format_gene_list("TNFRSF18/CD44/CXCL1/CXCL3/PYCARD/PTPN2/BIRC3/PARP9/TRIM44/IL13RA1/TP53/CYLD/CCL22/PTPRC/TANK/CSF2RA/IFNGR1/HCK/IL7R/TNFRSF4/IFNAR1/LAPTM5/HIF1A/EBI3/CSF2RB/SAMHD1/CD70")

genes <- unique(c(genes_24h, genes_3h, genes_UT))

heatmap_sig__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_sig__3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_sig__24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_sig__UT, heatmap_sig__3hCA, heatmap_sig__24hCA, ncol = 3)
```
```{r}
#cellular response to cytokine stimulus
genes_UT <-format_gene_list("ISG15/IFIT2/TCL1A/AGPAT2/MX1/CSF2RB/MT2A/TNFRSF4/CCL3/ADAR/SOCS1/IRF1/CD47/EIF5A/IL6ST/UGCG/YTHDF2/CREB1/BIRC3/HSPA5/TRIM44/IL10RA/XRCC5/TNFSF13B/GAS6/CDC42/STAT6/ZFAND6/STAT1/CSF2RA/ACTR3/BIRC2/HCLS1/LILRB4/NFKB1/SMPD3/SP100/IFI16/PTPN2/HAX1/CYLD/HSP90AB1/NFE2L2/NCL/IL3RA/TPR/SOCS3/RPLP0/SRSF3/KIF5B")
genes_3h <- format_gene_list("CXCL10/CCL3/ISG15/IFIT2/IFIT3/CXCL2/MX1/OAS1/IL1B/SOCS1/RNMT/USP18/TNFSF13B/PARP14/ZFP36/IRF1/NMI/LGALS9/XRCC5/XBP1/PNPT1/SP100/B3GNT2/IL3RA/TNFRSF4/PTPRC/GCLM/FCER1G/STAT3/NDUFA13/TRIM44/IFI16/CCL4/MTF2/LILRA4/MYD88/IRF8/ADAR/APPL1")
genes_24h <- format_gene_list("BCLAF1/SIRPA/IL10RB/IL10RA/CREB1/PRPF8/TMSB4X/SPI1/LGALS9/NDUFA13/SH2B3/ZFP36/KIF5B/TNFRSF18/CD44/CXCL1/TPR/ETV3/CXCL3/PYCARD/CD40/CD58/UBE2G2/TFRC/PTPN2/BIRC3/PARP9/RAD23B/ACTR3/TRIM44/RNMT/HCLS1/IL13RA1/TP53/CYLD/CCL22/PTPRC/IFI16/TANK/ZFP36L1/CSF2RA/IFNGR1/HCK/IL7R/CDC42/TNFRSF4/NRP2/IFNAR1/LAPTM5/HIF1A/EBI3/GSN/CSF2RB/SAMHD1/CD70")

genes <- unique(c(genes_24h, genes_3h, genes_UT))

heatmap_res_cyto__UT <- DoHeatmap(cluster.average_UT, features = genes, group.by = 'sex', draw.lines = TRUE) +
 scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) + guides(fill = FALSE, color = FALSE) + ggtitle("UT") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y =element_text(size =10, face = "bold", vjust = 0),plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_res_cyto__3hCA <- DoHeatmap(cluster.average_3hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(fill = FALSE, color = FALSE) + ggtitle("3hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

heatmap_res_cyto__24hCA <- DoHeatmap(cluster.average_24hCA, features = genes, group.by = 'sex', draw.lines = TRUE) +
  scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), oob = scales::oob_squish_any, limits = c(-2, 2)) +
  guides(color = FALSE) + ggtitle("24hCA") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), plot.margin = unit(c(1, 0, 1, 1), "cm"))

grid.arrange(heatmap_res_cyto__UT, heatmap_res_cyto__3hCA, heatmap_res_cyto__24hCA, ncol = 3)
```

## Conclusion

In conclusion, this study explored the role of plasmacytoid dendritic cells (pDCs) in the context of Pseudomonas infection and investigated potential sex differences in pDC function and gene expression. The findings revealed distinct gene expression patterns and pathway activations between males and females at different timepoints during Pseudomonas treatment. Additionally, the study highlighted the integral role of pDCs in antitumor immune responses, autoimmune responses, and their potential as therapeutic targets. Further investigations into sex differences and the involvement of pDCs in infectious diseases and autoimmune disorders, such as lupus, can provide valuable insights for understanding immune responses and developing targeted therapies.
