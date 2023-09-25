#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <numberOfTaxaInAGroup_list> <numSamples_list>"
    exit 1
fi

# Extract the numberOfTaxaInAGroup_list and numSamples_list from the command-line arguments
numberOfTaxaInAGroup_list=$1
numSamples_list=$2

# Loop over numberOfTaxaInAGroup_list and numSamples_list and run MATLAB script in parallel using screen
for numberOfTaxaInAGroup in $numberOfTaxaInAGroup_list; do
    for numSamples in $numSamples_list; do
        screen -dmS "screen_$numberOfTaxaInAGroup_$numSamples" matlab -nodisplay -r "RunParallel($numberOfTaxaInAGroup, $numSamples); exit;"
    done
done
