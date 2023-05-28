# Define the list of required packages with version details
required_packages <- c(
  "Seurat = 4.3.0",
  "dplyr = 1.1.2",
  "ggplot2 = 3.4.2",
  "patchwork = 1.1.2",
  "rtracklayer = 1.58.0",
  "methods = 4.2.1",
  "Matrix = 1.5-4",
  "monocle3 = 1.3.1",
  "SeuratWrappers = 0.3.1",
  "EnhancedVolcano = 1.16.0",
  "AnnotationDbi = 1.60.2",
  "clusterProfiler = 4.6.2",
  "org.Hs.eg.db = 3.16.0",
  "stringr = 1.5.0",
  "plotly = 4.10.1",
  "tibble = 3.2.1",
  "ggrepel = 0.9.3",
  "paletteer = 1.5.0",
  "RColorBrewer = 1.1-3",
  "gridExtra = 2.3",
  "tidyr = 1.3.0",
  "e1071 = 1.7-13"
)

# Install required packages
install_packages <- function(packages) {
  for (package in packages) {
    install.packages(package, dependencies = TRUE)
  }
}

# Call the install_packages function
install_packages(required_packages)
