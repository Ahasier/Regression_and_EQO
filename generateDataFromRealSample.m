function [taxaAbundance, functionalOutput, extraArgOut] = generateDataFromRealSample(index, noiseLevel, settings, numSamples, biOrNumbi)
global paths

load([paths.data, 'data/realdata_subsample5.mat'], 'abundanceData', 'prunedbranchedTaxa');

taxaAbundance = abundanceData(randSample(size(abundanceData, 1), numSamples), :);

selectedLeaves = prunedbranchedTaxa(:, index);

if strcmp(biOrNumbi, 'Binary')
    syntheticCoefficients = selectedLeaves;
elseif strcmp(biOrNumbi, 'NonBinary')
    syntheticCoefficients(selectedLeaves) = normrnd(1, settings.BetaEps, 1, sum(selectedLeaves));
else
    error("Flag must be either binary or non-binary.")
end

functionalOutputBeforeCenter = taxaAbundance * syntheticCoefficients + normrnd(0, noiseLevel, size(taxaAbundance, 1), 1);
functionalOutput = centerMatrix(functionalOutputBeforeCenter);

extraArgOut = {syntheticCoefficients};
end