args <- commandArgs(trailingOnly = TRUE)
group_size <- as.integer(args[1])

# Initialize data_split as a list
data_split <- list()

# Read training data and output from CSV files
M <- as.matrix(read.csv('EQO/data/trainingData.csv', header=FALSE))
y <- as.matrix(read.csv('EQO/data/trainingOutput.csv', header=FALSE))

# Source the required R script
source("EQO/EQOFunctions/EQO_GA.R")

# Execute the EQ_optim function
result <- EQ_optim('c', M, y, Nmax = group_size, amin=0, amax=1, popSize=100, maxIter=200, parallel=TRUE, monitor=FALSE)

# Save the coefficients to a CSV file for MATLAB to read
write.csv(result$x, 'data/coefficients.csv', row.names=FALSE)
