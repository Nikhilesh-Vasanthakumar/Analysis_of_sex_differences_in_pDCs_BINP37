#Script to replace the / with "" and seperate with ,

#Eg NIK/HIL/ESH to "NIK","HIL","ESH"

#Input: Gene list with / seperated genes
#Output: Gene Vector

#Read the gene list
format_gene_list <- function(input_string) {
  gene_list <- strsplit(input_string, "/")
  gene_list <- unlist(gene_list)
  return(gene_list)
}

# Calculate Z-score for the selected genes
zscore <- function(x) {
  (x - mean(x)) / sd(x)
}


