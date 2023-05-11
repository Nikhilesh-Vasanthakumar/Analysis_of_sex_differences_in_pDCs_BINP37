#Load the dataset
gene_names <- c("EIF1AY", "RPS4Y1")
X <- t(cluster.average@assays$RNA@data[gene_names, ])
y <- as.numeric(cluster.average@meta.data$sex == "Male")
#Split the data into training and test sets
set.seed(123)
train_idx <- sample(nrow(X), 0.8 * nrow(X))
test_idx <- setdiff(seq_len(nrow(X)), train_idx)
train_X <- as.matrix(X[train_idx, ])
train_y <- y[train_idx]
test_X <- as.matrix(X[test_idx, ])
test_y <- y[test_idx]

# Train the SVM model
model <- svm(train_X, train_y, kernel = "linear", cost = 1, probability = TRUE)

# Make predictions on the test set with probability information
pred_y <- predict(model, test_X, probability = TRUE)

#If value >0.5 then classify as male or else female for pred_y
# Classify predictions as male or female
pred_y_class <- ifelse(pred_y > 0.5, 1, 0)
# Compare predicted labels with true labels
accuracy <- mean(pred_y_class == test_y)

#Save the model
saveRDS(model, file = "svm_model.rds")