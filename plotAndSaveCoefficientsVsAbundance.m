function plotAndSaveCoefficientsVsAbundance(coefficients, trainingData, index, flag)
global paths

% Compute parameters
numSamples = size(trainingData, 1) * 2;

load([paths.data, 'subsample1taxa50_samplesize', num2str(numSamples),'.mat'], 'taxaAbundance', 'prunedbranchedTaxa');

stdTaxonAbundance = std(taxaAbundance);

trueTaxaIndex = prunedbranchedTaxa(:,index) > 0;
irrelevantTaxaIndex = prunedbranchedTaxa(:,index) == 0;

trueTaxaStdAbundance = stdTaxonAbundance(trueTaxaIndex);
irrelevantTaxaStdAbundance = stdTaxonAbundance(irrelevantTaxaIndex);

trueTaxaCoefficients = coefficients(trueTaxaIndex);
irrelevantTaxaCoefficients = coefficients(irrelevantTaxaIndex);

if flag == 1
    stringFlag = '_standardized';
elseif flag == 0
    stringFlag = '_nonstandardized';
else
    error('flag must be either 0 or 1');
end

id = 25;
noiseLevel = 0;

figure();
h1 = scatter(log(irrelevantTaxaStdAbundance), irrelevantTaxaCoefficients, 'filled', 'blue');
hold on;
h2 = scatter(log(trueTaxaStdAbundance), trueTaxaCoefficients, 'filled', 'red');
box on;
ax = plotstyle(gca, 1);
grid on
xlabel('log(std(abundance) of a taxon)');
ylabel('coefficients');
title([num2str(numSamples), ' sample, ', stringFlag(2:end), ', id ', num2str(id),', noise level ', num2str(noiseLevel)])

figure();
h1 = scatter(log(irrelevantTaxaStdAbundance), log(irrelevantTaxaCoefficients), 'filled', 'blue');
hold on;
h2 = scatter(log(trueTaxaStdAbundance), log(trueTaxaCoefficients), 'filled', 'red');
box on;
ax = plotstyle(gca, 1);
grid on
xlabel('log(std(abundance) of a taxon)');
ylabel('log(coefficients)');
title([num2str(numSamples), ' sample, ', stringFlag(2:end), ', id ', num2str(id),', noise level ', num2str(noiseLevel)])

trueCoefficients = ones(50,1);
trueCoefficients(irrelevantTaxaIndex) = 0;
deltaX = coefficients - trueCoefficients .* std(taxaAbundance)';
dX2 = deltaX.^2;
varA = var(taxaAbundance);
SNR = varA'./dX2;

figure();scatter(find(trueTaxaIndex),log(SNR(trueTaxaIndex)), 'filled', 'red')
hold on;scatter(find(irrelevantTaxaIndex),log(SNR(irrelevantTaxaIndex)), 'filled', 'blue')
grid on;box on;h = plotstyle(gca, 1);
xlabel('Taxa index');
ylabel('log(Actual SNR)');
title([num2str(numSamples), ' sample, ', stringFlag(2:end), ', id ', num2str(id),', noise level ', num2str(noiseLevel)])

save(['results/coefficientsVsAbundance_samplesize', num2str(numSamples), '_index', num2str(index), stringFlag, '.mat'], 'trueTaxaMeanAbundance', 'trueTaxaCoefficients', 'irrelevantTaxaMeanAbundance', 'irrelevantTaxaCoefficients')
end