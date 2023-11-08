args <- commandArgs(trailingOnly = TRUE)
group_size <- as.integer(args[1])
training_data_filename <- args[2]
training_output_filename <- args[3]
coefficients_filename <- args[4]

# Read training data and output from CSV files
M <- as.matrix(read.csv(training_data_filename, header=FALSE))
y <- as.matrix(read.csv(training_output_filename, header=FALSE))

# Source the required R script
source("EQO/EQOFunctions/EQO_GA.R")

# Execute the EQ_optim function
result <- EQ_optim('c', M, y, Nmax = group_size, amin=0, amax=1, popSize=100, maxIter=200, parallel=FALSE, monitor=FALSE)

# Save the coefficients to a CSV file for MATLAB to read
write.csv(result$x, coefficients_filename, row.names=FALSE)