#!/usr/bin/Rscript

#Implementation of a Support Vector Machine (SVM) model for binary classification using the e1071 package in R. Here is a brief explanation of the code:
#Data Preparation: The gene expression data for three genes, "EIF1AY", "RPS4Y1" and "XIST" are loaded into the variable X, and the corresponding sex labels ("Male" or "Female") are converted into numeric values (y).
#Data Splitting: The dataset is split into training and test sets using a 80:20 ratio. The training set (train_X and train_y) is used to train the SVM model, while the test set (test_X and test_y) is used for evaluating the model's performance.
#SVM Model Training: The SVM model is trained using the svm() function with a linear kernel, a cost parameter of 1, and setting probability = TRUE to obtain probability estimates for classification.
#Prediction and Evaluation: The trained SVM model is used to make predictions on the test set (test_X) by calling the predict() function. The predicted probabilities and class labels are obtained.
#Performance Evaluation: The predicted class labels are compared with the true labels (test_y) to calculate the accuracy of the SVM model.

# Load the e1071 library for SVM implementation
library(e1071)

# Define the gene names of interest
gene_names <- c("EIF1AY", "RPS4Y1","XIST")

# Extract the gene expression data matrix for the selected genes
X <- t(cluster.average@assays$RNA@data[gene_names, ])

# Convert the sex labels to numeric values (0 for Female, 1 for Male)
y <- as.numeric(cluster.average@meta.data$sex == "Male")

# Split the data into training and test sets
set.seed(123)
train_idx <- sample(nrow(X), 0.8 * nrow(X))  # Select random indices for training set
test_idx <- setdiff(seq_len(nrow(X)), train_idx)  # Remaining indices for test set
train_X <- as.matrix(X[train_idx, ])  # Training set features
train_y <- y[train_idx]  # Training set labels
test_X <- as.matrix(X[test_idx, ])  # Test set features
test_y <- y[test_idx]  # Test set labels

# Train the SVM model with a linear kernel, cost parameter of 1, and probability estimation
model <- svm(train_X, train_y, kernel = "linear", cost = 1, probability = TRUE)

# Make predictions on the test set with probability information
pred_y <- predict(model, test_X, probability = TRUE)

# Classify predictions as male (1) or female (0) based on a threshold of 0.5
pred_y_class <- ifelse(pred_y > 0.5, 1, 0)

# Compare predicted labels with true labels to calculate accuracy
accuracy <- mean(pred_y_class == test_y)

# Save the trained SVM model (optional)
# saveRDS(model, file = "svm_model.rds")
