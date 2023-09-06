library(GA) # Assuming the GA package is required for the EQ_optim function

# Load the data
MetaTara <- read.csv('Supp_data2_Ocean.csv')
responses <- read.csv('X.csv', header=FALSE)
y <- as.matrix(responses[1])
M <- as.matrix(MetaTara[3:ncol(MetaTara)])

# Splitting function
split_data <- function(M, y, ratio = 0.5) {
  sample_size <- floor(ratio * nrow(M))
  train_indices <- sample(seq_len(nrow(M)), size = sample_size)
  
  train_M <- M[train_indices, ]
  train_y <- y[train_indices, ]
  
  test_M <- M[-train_indices, ]
  test_y <- y[-train_indices, ]
  
  list(train_M = train_M, train_y = train_y, test_M = test_M, test_y = test_y)

}

# Compute R^2 based on the model and test data
evaluate_performance <- function(model, test_M, test_y) {
  # Compute predicted values for the test set
  predictedOutput <- test_M %*% model$x
  
  # Compute residuals
  residuals <- test_y - predictedOutput
  
  # Compute sum of squared residuals
  SSR <- sum(residuals^2)
  
  # Compute total sum of squares
  TSS <- sum((test_y - mean(test_y))^2)
  
  # Calculate R^2 as 1 minus the ratio of SSR to TSS
  R2 <- 1 - (SSR / TSS)
  
  return(R2)
}

# Function to compute AIC
compute_AIC <- function(trainingData, trainingOutput, coefficients) {
  aicValues <- numeric(length(coefficients))
  sortedTaxaIndices <- order(coefficients, decreasing = TRUE)
  
  for (n in 1:length(coefficients)) {
    idx <- sortedTaxaIndices[1:n]
    groupAssemblage <- rep(0, length(coefficients))
    groupAssemblage[idx] <- 1
    
    r2 <- evaluate_performance(groupAssemblage, trainingData, trainingOutput)
    ssr <- (1 - r2) * sum((trainingOutput - mean(trainingOutput))^2)
    
    k <- ncol(trainingData) + 1
    L <- exp(-ssr/2)
    aic <- 2*k - 2*log(L)
    
    aicValues[n] <- aic
  }
  
  return(aicValues)
}

# Function to compute cumulative R^2
calculateCumulativeR2 <- function(allCoefficients, allOutSampleR2) {
  cumulativeR2 <- rowSums(allCoefficients * allOutSampleR2)
  return(cumulativeR2)
}

# Loop to repeat the process 100 times
results_list <- list()
all_AICs <- matrix(0, 100, ncol(M))

for(i in 1:100) {
  data_split <- split_data(M, y)
  
  # Step (b) 1: Compute EQO on training data
  result <- EQ_optim(pattern = pattern, M = data_split$train_M, y = data_split$train_y)
  
  # Step (b) 2: Calculate the AIC values under different group sizes
  AIC_values <- compute_AIC(data_split$train_M, data_split$train_y, result$members)
  all_AICs[i, ] <- AIC_values
  
  # Step (c): Calculate the cross-validation R^2 on test subset
  performance <- evaluate_performance(result, data_split$test_M, data_split$test_y)
  
  # Store results for this iteration
  results_list[[i]] <- list(result = result, performance = performance, AIC = AIC_values)
}

# Step (d) 1: Calculate the cumulative R^2 for each taxon
cumulativeR2 <- calculateCumulativeR2(results_list$allBinaryCoefficients, results_list$allOutSampleR2)

# Step (d) 2: Get the best coefficients by selecting the top k = optimalGroupSize taxa based on cumulative R^2
optimalGroupSize <- which.min(rowMeans(all_AICs))
top_taxa <- order(cumulativeR2, decreasing = TRUE)[1:optimalGroupSize]
cross_validated_group <- M[, top_taxa]

# Print results
print(results_list)
print(cross_validated_group)
