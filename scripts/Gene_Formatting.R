#!/usr/bin/env Rscript

# Script containing function to replace the '/' with '' and separate with commas, and a function to calculate the signature score based on expression data

# Function to format gene list
# This function takes an input string and splits it by '/' to separate multiple genes.
# The resulting gene list is returned as a character vector.
format_gene_list <- function(input_string) {
  gene_list <- strsplit(input_string, "/")  # Split the input string by '/'
  gene_list <- unlist(gene_list)  # Convert the split result to a vector
  return(gene_list)  # Return the formatted gene list
}

# Function to calculate Z-score
# This function calculates the Z-score for a given numeric vector.
# The Z-score is obtained by subtracting the mean from each value and dividing by the standard deviation.
# The resulting Z-scored vector is returned.
zscore <- function(x) {
  (x - mean(x)) / sd(x)  # Calculate Z-score for the input vector
}

# Explanation and Comments:
# The 'format_gene_list' function is designed to handle gene lists provided as a string with genes separated by '/'.
# It splits the input string by '/' and returns the resulting genes as a character vector.
# For example, if the input string is "gene1/gene2/gene3", the function will return a character vector c("gene1", "gene2", "gene3").

# The 'zscore' function calculates the Z-score for a given numeric vector.
# It subtracts the mean from each value and divides by the standard deviation to obtain the Z-score.
# This function is useful for standardizing gene expression data or any other numeric data for further analysis.
# The resulting Z-scored vector represents the number of standard deviations a data point is away from the mean.

