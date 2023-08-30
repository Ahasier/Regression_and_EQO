function [taxaAbundance, functionalOutput, syntheticCoefficients, selectedLeaves, treeStructure] = generateSyntheticData(treeData, phylogenyDependency, numberOfTaxaInAGroup, noiseLevel, numberOfSamples, settings)
% GENERATESYNTHETICDATA simulates microbial abundance and functional output data based on 
% the provided parameters and settings.
%
% Inputs:
%   treeData: data related to the phylogenetic tree
%   phylogenyDependency: parameter to model dependency in the data based on phylogeny
%   numberOfTaxaInAGroup: the number of taxa to be grouped together in synthetic data
%   noiseLevel: the amount of noise to be added to the synthetic data
%   numberOfSamples: the number of samples to generate for the synthetic data
%   settings: a struct containing various settings for data generation
%
% Outputs:
%   taxaAbundance: matrix containing taxa abundance data (samples x taxa)
%   functionalOutput: vector containing the functional outputs (samples x 1)
%   syntheticCoefficients: vector containing the true coefficients used to generate synthetic data (taxa x 1)
%   selectedLeaves: logical vector indicating which traits (leaves) were selected to generate synthetic data
%   treeStructure: additional information about the tree structure

%% Extracting tree structure for potential further use
treeStructure = get(treeData);

%% 1. Generate covariance matrix from phylogenetic distances
covarianceMatrix = createCovarianceMatrix(treeData, phylogenyDependency);

%% 2. Evolve traits on the phylogenetic tree using the covariance matrix
traitValues = evolveTraitsOnTree(covarianceMatrix);

%% 3. Select traits based on trait values and generate corresponding coefficients
[syntheticCoefficients, selectedLeaves] = selectTraitsAndGenerateCoefficients(traitValues, numberOfTaxaInAGroup, settings);

%% 4. Generate taxa abundance data
taxaAbundance = generateAbundanceData(treeStructure, numberOfSamples, settings);

%% 5. Construct the observable output using the generated data
functionalOutput = constructObservableOutput(taxaAbundance, syntheticCoefficients, noiseLevel);
end

%% Helper Functions
% Computes a covariance matrix based on phylogenetic distances and a given dependency parameter.
function covarianceMatrix = createCovarianceMatrix(treeData, phylogenyDependency)
distanceMatrix = pdist(treeData);
squareDistanceMatrix = squareform(distanceMatrix);

if phylogenyDependency == 0
    covarianceMatrix = eye(size(squareDistanceMatrix));
else
    covarianceMatrix = exp(-squareDistanceMatrix./phylogenyDependency);
end
end

% Simulates the evolution of traits along a phylogenetic tree based on a given covariance matrix.
function traitValues = evolveTraitsOnTree(covarianceMatrix)
matrixLength = length(covarianceMatrix);
traitValues = mvnrnd(ones(matrixLength,1), covarianceMatrix, 1);
end

% Selects a subset of traits based on their values and generates corresponding synthetic coefficients.
function [syntheticCoefficients, selectedLeaves] = selectTraitsAndGenerateCoefficients(traitValues, numberOfTaxaInAGroup, settings)
cap = 0.5;

[~, traitSortIdx] = sort(traitValues, 'descend');
selectedTraitsIndices = traitSortIdx(1:numberOfTaxaInAGroup);

syntheticCoefficients = zeros(length(traitValues), 1);

if isfield(settings, 'Beta0')
    syntheticCoefficients(selectedTraitsIndices) = normrnd(settings.Beta0, settings.BetaEps, 1, numberOfTaxaInAGroup);
else
    % Use default value 1 in case when settings.Beta0 is not defined.
    syntheticCoefficients(selectedTraitsIndices) = normrnd(1, settings.BetaEps, 1, numberOfTaxaInAGroup);
end

while any(syntheticCoefficients(selectedTraitsIndices) < cap)
    idxRege = syntheticCoefficients(selectedTraitsIndices) < cap;
    syntheticCoefficients(selectedTraitsIndices(idxRege)) = normrnd(settings.Beta0, settings.BetaEps, 1, sum(idxRege));
end

selectedLeaves = (syntheticCoefficients ~= 0);
end

% Generates synthetic microbial abundance data using either a model to mimic real data or a normal distribution.
function taxaAbundance = generateAbundanceData(treeData, numberOfSamples, settings)
xm = 87.71;
a = 4.254;

if generateMimicRealAbundance(settings)
    Nseed = rand(treeData.NumLeaves, numberOfSamples)';
    taxaAbundance = (- xm * (-1 + (1 - Nseed).^(1/a))) ./ ((1 - Nseed).^(1/a));
    taxaAbundance = floor(taxaAbundance);
else
    taxaAbundance = normrnd(10, 10, treeData.NumLeaves, numberOfSamples)';
end
end

% Calculates synthetic functional output data based on synthetic abundance data and coefficients, with added noise.
function functionalOutput = constructObservableOutput(taxaAbundance, syntheticCoefficients, noiseLevel)
functionalOutput = taxaAbundance * syntheticCoefficients + normrnd(0, noiseLevel, size(taxaAbundance, 1), 1);
end

% Checks the settings to determine whether to mimic real abundance data or generate synthetic abundance data using a normal distribution.
function isRealAbundance = generateMimicRealAbundance(settings)
isRealAbundance = isfield(settings, 'RealAbd') && strcmp(settings.RealAbd, 'On');
end