#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <numberOfTaxaInAGroup_list> <numSamples_list> <regressionMethod>"
    exit 1
fi

# Extract the numberOfTaxaInAGroup_list, numSamples_list, and regressionMethod from the command-line arguments
numberOfTaxaInAGroup_list=($1)
numSamples_list=($2)
regressionMethod=$3

# Loop over numberOfTaxaInAGroup_list and numSamples_list and run MATLAB script in parallel using screen
for numberOfTaxaInAGroup in ${numberOfTaxaInAGroup_list[@]}; do
    for numSamples in ${numSamples_list[@]}; do
        screen -dmS "screen_${numberOfTaxaInAGroup}_${numSamples}" matlab -nodisplay -r "RunParallel($numberOfTaxaInAGroup, $numSamples, '$regressionMethod'); exit;"
        
        # # Check if regressionMethod is 'OLS' and wait for 30 seconds
        # if [ "$regressionMethod" == "OLS" ]; then
        #     sleep 5
        # fi
    done
done
