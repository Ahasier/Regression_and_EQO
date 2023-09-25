function [taxaAbundance, functionalOutput, extraArgOut] = generateData(treeData, phylogenyDependency, numberOfTaxaInAGroup, noiseLevel, numberOfSamples, settings)
% GENERATEDATA loads real data or generates synthetic data based on the provided parameters for regression analysis.
%
% INPUTS:
%   treeData: data related to the phylogenetic tree
%   phylogenyDependency: parameter to model dependency in the data based on phylogeny
%   numberOfTaxaInAGroup: the number of taxa to be grouped together in synthetic data
%   noiseLevel: the amount of noise to be added to the synthetic data
%   numberOfSamples: the number of samples to generate for the synthetic data
%   settings: a struct containing various settings for data generation/loading
%
% OUTPUTS:
%   taxaAbundance: matrix containing taxa abundance data (samples x taxa)
%   functionalOutput: vector containing the functional outputs (samples x 1)
%   syntheticCoefficients: vector containing the true coefficients used to generate synthetic data (taxa x 1), empty for real data

% Check if using real data based on the provided settings
% If using real data, load them. Otherwise, generate synthetic data.
if usingRealData(settings)
    % Load real data from a file or database
    [taxaAbundance, functionalOutput] = loadRealData(numberOfSamples);
    % For real data, synthetic coefficients are not applicable, so return empty
    extraArgOut = {};
else
    % Generate synthetic data based on the input parameters and settings
    % This include generating synthetic microbial abundance data,
    % synthetic functional output data, and the true coefficients used to generate this synthetic data
    [taxaAbundance, functionalOutput, syntheticCoefficients] = generateSyntheticData(treeData, phylogenyDependency, numberOfTaxaInAGroup, noiseLevel, numberOfSamples, settings);
    extraArgOut = {syntheticCoefficients};
end

% If phylogenetic information is incorporated, use the grouping method to
% generate the grouped abundance to handle it.
if useExtraFeatures(settings)
    extraPhyloVars.numTaxa = size(taxaAbundance, 2);
    [extraPhyloVars.numBranches, extraPhyloVars.addedLeaves, extraPhyloVars.Idx] = groupPhylogeny(get(treeData));
    [extraArgOut{end + 1}] = extraPhyloVars;
end
end
