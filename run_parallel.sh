#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -lt 6 ]; then
    echo "Usage: $0 <numberOfTaxaInAGroup_list> <numSamples_list> <regressionMethod> <phylogenyDependency> <noiseLevel> <betaEps> [<realAbd>] [<usePhylogeny>]"
    exit 1
fi

# Extract the parameters from the command-line arguments
numberOfTaxaInAGroup_list=($1)
numSamples_list=($2)
regressionMethod=$3
phylogenyDependency=$4
noiseLevel=$5
betaEps=$6

# Set a default value for realAbd and usePhylogeny if it is not provided
realAbd=${7:-"EMPTY_ARRAY"}
index_list=(${8:-""})
usePhylogeny=${9:-"EMPTY_ARRAY"}

# Loop over numberOfTaxaInAGroup_list and numSamples_list and run MATLAB script in parallel using screen
for numberOfTaxaInAGroup in ${numberOfTaxaInAGroup_list[@]}; do
    for numSamples in ${numSamples_list[@]}; do
        # If index_list is not provided, then run the MATLAB script without index
        if [ -z "$index_list" ]; then
            screen -dmS "screen_${numberOfTaxaInAGroup}_${numSamples}" matlab -nodisplay -r "RunParallel($numberOfTaxaInAGroup, $numSamples, '$regressionMethod', $phylogenyDependency, $noiseLevel, $betaEps, '$realAbd', [], '$usePhylogeny'); exit;"
        else
            # Loop over index_list and run MATLAB script in parallel using screen
            for index in ${index_list[@]}; do
                screen -dmS "screen_${numberOfTaxaInAGroup}_${numSamples}_${index}" matlab -nodisplay -r "RunParallel($numberOfTaxaInAGroup, $numSamples, '$regressionMethod', $phylogenyDependency, $noiseLevel, $betaEps, '$realAbd', $index, '$usePhylogeny'); exit;"
            done
        fi
        # Wait for 30 seconds
        # sleep 30
    done
done