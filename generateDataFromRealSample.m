function [taxaAbundanceStandardized, functionalOutputStandardized, extraArgOut] = generateDataFromRealSample(index, noiseLevel, settings, numSamples)
global paths

filename = 'subsample1taxa50';

load([paths.data, filename, '.mat'], 'abundanceData', 'prunedbranchedTaxa', 'subT');

taxaAbundance = abundanceData(randsample(size(abundanceData, 1), numSamples), :);

selectedLeaves = prunedbranchedTaxa(:, index);

if isfield(settings, 'BetaEps') && settings.BetaEps > 0
    syntheticCoefficients = zeros(length(selectedLeaves), 1);
    syntheticCoefficients(logical(selectedLeaves)) = normrnd(1, settings.BetaEps, 1, sum(selectedLeaves));
else
    syntheticCoefficients = selectedLeaves;
end

functionalOutput = taxaAbundance * syntheticCoefficients + normrnd(0, noiseLevel, size(taxaAbundance, 1), 1);

% save([paths.data, filename, '_samplesize', num2str(numSamples),'.mat'], 'taxaAbundance', 'functionalOutput', 'prunedbranchedTaxa');

% Standarize the taxa abundance matrix and functional outputs
taxaAbundanceStandardized = standarization(taxaAbundance);
% functionalOutputStandardized = standarization(functionalOutput);
% taxaAbundanceStandardized = centerMatrix(taxaAbundance);
functionalOutputStandardized = centerMatrix(functionalOutput);

extraArgOut = {syntheticCoefficients};

% If phylogenetic information is incorporated, use the grouping method to
% generate the grouped abundance to handle it.
if useExtraFeatures(settings)
    extraPhyloVars.numTaxa = size(taxaAbundance, 2);
    [extraPhyloVars.numBranches, extraPhyloVars.addedLeaves, extraPhyloVars.Idx] = groupPhylogeny(subT);
    [extraArgOut{end + 1}] = extraPhyloVars;
end
end

% Function to standardize a matrix
function [standardizedM] = standarization(M)
% Calculate the mean of each column
colMeans = mean(M);

% Calculate the standard devidation of each column
colStds = std(M);

% Subtract the column means from each element in the corresponding column
standardizedM = (M - colMeans)./colStds;

standardizedM(isnan(standardizedM)) = 0;
end