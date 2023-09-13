function [taxaAbundance, functionalOutput, syntheticCoefficients, selectedLeaves, treeStructure] = generateSyntheticData(treeData, phylogenyDependency, numberOfTaxaInAGroup, noiseLevel, numberOfSamples, settings)
% GENERATESYNTHETICDATA simulates microbial abundance and functional output data based on 
% the provided parameters and settings.
%
% INPUTS:
%   treeData: data related to the phylogenetic tree
%   phylogenyDependency: parameter to model dependency in the data based on phylogeny
%   numberOfTaxaInAGroup: the number of taxa to be grouped together in synthetic data
%   noiseLevel: the amount of noise to be added to the synthetic data
%   numberOfSamples: the number of samples to generate for the synthetic data
%   settings: a struct containing various settings for data generation
%
% OUTPUTS:
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
% If generating mimic-real abundance data, generate them from the specific long-tail (Pareto) distribution
if isGeneratingMimicRealAbundance(settings)
    realData = loadRealTreeData().N1;
    taxaAbundance = generateMimicRealAbundance(realData, numberOfSamples, treeData.NumLeaves);
else
    % Otherwise, generate abundance data from a random Guassian.
    taxaAbundance = normrnd(10, 10, numberOfSamples, treeData.NumLeaves);
end
end

% Calculates synthetic functional output data based on synthetic abundance data and coefficients, with added noise.
function functionalOutput = constructObservableOutput(taxaAbundance, syntheticCoefficients, noiseLevel)
functionalOutput = taxaAbundance * syntheticCoefficients + normrnd(0, noiseLevel, size(taxaAbundance, 1), 1);
end

% Checks the settings to determine whether to mimic real abundance data or generate synthetic abundance data using a normal distribution.
function isRealAbundance = isGeneratingMimicRealAbundance(settings)
isRealAbundance = isfield(settings, 'RealAbd') && strcmp(settings.RealAbd, 'On');
end

% % Generate mimic-real abundance data from long-tail distribution
% function taxaAbundance = generateMimicRealAbundance(NumLeaves, numberOfSamples)
% % Set parameter values for the distribution of mimic-real taxa abundance
% [xm, a] = setMimicRealAbundanceDistributionParameters();
% 
% % 1. Set the row sums of abundance matrix following a long-tail distribution
% Nseed1 = rand(numberOfSamples)';
% rowSums = (- xm.dist1 * (-1 + (1 - Nseed1).^(1/a.dist1))) ./ ((1 - Nseed1).^(1/a.dist1));
% 
% % 2. Set the values of matrix elements from a long-tail distribution
% % Randomly generate taxa abundance data from the determined distribution
% Nseed2 = rand(NumLeaves, numberOfSamples)';
% taxaAbundance = (- xm.dist2 * (-1 + (1 - Nseed2).^(1/a.dist2))) ./ ((1 - Nseed2).^(1/a.dist2));
% 
% % 3. Normalize the matrix to the determined row sums
% taxaAbundance = taxaAbundance.*rowSums./sum(taxaAbundance, 2);
% 
% % Taxa abundance takes integer values
% taxaAbundance = floor(taxaAbundance);
% end

% % Generate mimic-real abundance data from long-tail distribution
% function taxaAbundance = generateMimicRealAbundance(realData, numSamples, numTaxa)
% % Extract row and column sums from real data
% realRowSums = sum(realData, 2);
% realColSums = sum(realData, 1);
% 
% % Simulate row and column sums using empirical distributions
% syntheticRowSums = datasample(realRowSums, numSamples);
% syntheticColSums = datasample(realColSums, numTaxa);
% 
% % Initialize synthetic matrix with Pareto-distributed values
% [xm, alpha] = setMimicRealAbundanceDistributionParameters();
% Nseed = rand(numSamples, numTaxa);
% taxaAbundance = xm.*(1 - Nseed).^(-1/alpha) - xm;
% 
% % Bi-normalize synthetic matrix to match row and column sums
% for iter = 1:1000
%     % Row normalization
%     taxaAbundance = bsxfun(@rdivide, taxaAbundance, sum(taxaAbundance, 2)) .* syntheticRowSums;
%     
%     % Column normalization
%     taxaAbundance = bsxfun(@rdivide, taxaAbundance, sum(taxaAbundance, 1)) .* syntheticColSums;
%     
%     % Convergence check can be added here if required
% end
% 
% % Introduce sparsity based on real data sparsity
% sparsityRate = sum(sum(realData == 0)) / numel(realData);
% mask = rand(numSamples, numTaxa) < sparsityRate;
% taxaAbundance(mask) = 0;
% end

function taxaAbundance = generateMimicRealAbundance(realData, numSamples, numTaxa)
% Generate a synthetic microbial abundance matrix.
%
% INPUTs:
%   numSamples: Number of samples (rows).
%   numTaxa: Number of taxa (columns).

% Fraction of taxa that are generalists
generalistFraction = 0.2137;
% Fraction of taxa that are dominant in each sample
dominantFraction = 0.1016;

% Initialize output synethetic abundance matrix
taxaAbundance = zeros(numSamples, numTaxa);

% Determine the number of generalists and specialists
numGeneralists = round(numTaxa * generalistFraction);

% Randomly select taxa to be generalists
generalistIndices = randsample(numTaxa, numGeneralists, false);

% Assign abundances for generalists
for i = 1:length(generalistIndices)
    taxaAbundance(:, generalistIndices(i)) = exprnd(0.75, [numSamples, 1]);
end

% Assign abundances for specialists
specialistIndices = setdiff(1:numTaxa, generalistIndices);
for i = 1:length(specialistIndices)
    sampleIndices = randsample(numSamples, randi([1, floor(numSamples/2)]), false);
    taxaAbundance(sampleIndices, specialistIndices(i)) = exprnd(0.75, [length(sampleIndices), 1]);
end

% Introduce diversity within taxa (dominant vs. rare species)
for j = 1:numTaxa
    dominantIndices = randsample(numSamples, round(numSamples * dominantFraction), false);
    taxaAbundance(dominantIndices, j) = taxaAbundance(dominantIndices, j) .* (rand([length(dominantIndices), 1]) * 5 + 5);
end

% Extract row and column sums from real data
realRowSums = sum(realData, 2);

% Simulate row and column sums using empirical distributions
syntheticRowSums = datasample(realRowSums, numSamples);

% Normalize such that each row sums to the empirical distribution
taxaAbundance = taxaAbundance ./ sum(taxaAbundance, 2) .* syntheticRowSums;

% Taxa abundance takes integer values
taxaAbundance = floor(taxaAbundance);
end

% Set parameter values for the distribution of mimic-real abundance matrix elements.
% function [xm, alpha] = setMimicRealAbundanceDistributionParameters()
% % Parameters for the first long-tail distribution (row sums)
% xm = 1.565;
% alpha = 0.6313;
% end