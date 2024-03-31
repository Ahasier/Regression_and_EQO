clc;clear
% taxaAbundance = exprnd(1, [1000 500]);
% taxaAbundance = normrnd(1, 1, [1000 500]);
numTaxa = 500;
numSamples = 1000;

% noiseLevel = 0.1 * mean(stdTaxaAbundance);
% SNR = 1;
SNR = 1;

stdTaxaAbundance = exp(-numTaxa/100:0.02:numTaxa/100 - 0.01);
taxaAbundance = zeros(numSamples, numTaxa);
for n = 1:numTaxa
    taxaAbundance(:, n) = normrnd(0, stdTaxaAbundance(n), [numSamples 1]);
%     taxaAbundance(:, n) = exprnd(stdTaxaAbundance(n), [numSamples 1]);
end

noiseLevel = median(stdTaxaAbundance) * sqrt(numTaxa) / SNR;

selectedLeaves = zeros(numTaxa,1);selectedLeaves(randperm(numTaxa, numTaxa / 5)) = 1;

if isfield(settings, 'BetaEps') && settings.BetaEps > 0
    syntheticCoefficients(selectedLeaves) = normrnd(1, settings.BetaEps, 1, sum(selectedLeaves));
else
    syntheticCoefficients = selectedLeaves;
end

functionalOutput = taxaAbundance * syntheticCoefficients + normrnd(0, noiseLevel, size(taxaAbundance, 1), 1);

taxaAbundanceStandardized = normalize(taxaAbundance);
functionalOutputStandardized = normalize(functionalOutput);
% taxaAbundanceStandardized = centerMatrix(taxaAbundance);
% functionalOutputStandardized = centerMatrix(functionalOutput);

coefficients = taxaAbundanceStandardized \ functionalOutputStandardized;

meanTaxonAbundance = mean(taxaAbundance);
stdTaxonAbundance = std(taxaAbundance);

trueTaxaIndex = selectedLeaves > 0;
irrelevantTaxaIndex = selectedLeaves == 0;

trueTaxaMeanAbundance = meanTaxonAbundance(trueTaxaIndex);
irrelevantTaxaMeanAbundance = meanTaxonAbundance(irrelevantTaxaIndex);
trueTaxaStdAbundance = stdTaxonAbundance(trueTaxaIndex);
irrelevantTaxaStdAbundance = stdTaxonAbundance(irrelevantTaxaIndex);

trueTaxaCoefficients = coefficients(trueTaxaIndex);
irrelevantTaxaCoefficients = coefficients(irrelevantTaxaIndex);

% figure();
% h1 = scatter(log(irrelevantTaxaMeanAbundance), irrelevantTaxaCoefficients, 'filled', 'blue');
% hold on;
% h2 = scatter(log(trueTaxaMeanAbundance), trueTaxaCoefficients, 'filled', 'red');
% box on;
% ax = plotstyle(gca, 1);
% xlabel('mean abundance of a taxon');
% ylabel('coefficients');
% title('Guassian random abundance, standardized')

figure();
h1 = scatter(log(irrelevantTaxaStdAbundance), log(irrelevantTaxaCoefficients), 'filled', 'blue');
hold on;
h2 = scatter(log(trueTaxaStdAbundance), log(trueTaxaCoefficients), 'filled', 'red');
box on;
grid on;
ax = plotstyle(gca, 1);
xlabel('log(std(mean abundance)) of a taxon');
ylabel('log(coefficients)');
title(['Guassian random abundance, standardized, \kappa=', num2str(SNR)])