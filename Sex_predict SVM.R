# Prepare the data
X <- t(pDC@assays$RNA@data) # gene expression data
y <- as.numeric(pDC@meta.data$sex == "Male") # binary labels (0 for female, 1 for male)
#Filter out the NA values and corresponding cells
na_idx <- which(is.na(y))
X <- X[-na_idx, ]
y <- y[-na_idx]

# Split the data into training and test sets
set.seed(123)
train_idx <- sample(nrow(X), 0.8 * nrow(X))
test_idx <- setdiff(seq_len(nrow(X)), train_idx)
train_X <- X[train_idx, ]
train_y <- y[train_idx]
test_X <- X[test_idx, ]
test_y <- y[test_idx]

# Train the SVM model
library(e1071)
model <- svm(train_X, train_y, kernel = "linear", cost = 1, probability = TRUE)

# Make predictions on the test set with probability information
pred_y <- predict(model, test_X, probability = TRUE)

# Extract the probability estimates for the positive class (male)
prob_pos <- attr(pred_y, "probabilities")[, 2]

library(caret)
confusionMatrix(factor(ifelse(prob_pos > 0.5, "Male", "Female")), factor(ifelse(test_y == 1, "Male", "Female")))
