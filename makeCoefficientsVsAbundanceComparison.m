function makeCoefficientsVsAbundanceComparison(numSamples1, numSamples2, index)
file1 = ['results/coefficientsVsAbundance_samplesize', num2str(numSamples1), '_index', num2str(index), '_nonstandardized', '.mat'];
file2 = ['results/coefficientsVsAbundance_samplesize', num2str(numSamples1), '_index', num2str(index), '_standardized', '.mat'];
file3 = ['results/coefficientsVsAbundance_samplesize', num2str(numSamples2), '_index', num2str(index), '_nonstandardized', '.mat'];
file4 = ['results/coefficientsVsAbundance_samplesize', num2str(numSamples2), '_index', num2str(index), '_standardized', '.mat'];

results1 = load(file1);
results2 = load(file2);
results3 = load(file3);
results4 = load(file4);

figure();
subplot(2,2,1)
h1.a = scatter(log(results1.irrelevantTaxaMeanAbundance), results1.irrelevantTaxaCoefficients, 'filled', 'blue');
hold on;
h1.b = scatter(log(results1.trueTaxaMeanAbundance), results1.trueTaxaCoefficients, 'filled', 'red');
box on;
ax1 = plotstyle(gca, 1);
xlabel('log(mean abundance of a taxon)');
ylabel('coefficients');
title('non standardized, not enough samples')
subplot(2,2,2)
h2.a = scatter(log(results2.irrelevantTaxaMeanAbundance), results2.irrelevantTaxaCoefficients, 'filled', 'blue');
hold on;
h2.b = scatter(log(results2.trueTaxaMeanAbundance), results2.trueTaxaCoefficients, 'filled', 'red');
box on;
ax2 = plotstyle(gca, 1);
xlabel('log(mean abundance of a taxon)');
ylabel('coefficients');
title('standardized, not enough samples')
subplot(2,2,3)
h3.a = scatter(log(results3.irrelevantTaxaMeanAbundance), results3.irrelevantTaxaCoefficients, 'filled', 'blue');
hold on;
h3.b = scatter(log(results3.trueTaxaMeanAbundance), results3.trueTaxaCoefficients, 'filled', 'red');
box on;
xlabel('log(mean abundance of a taxon)');
ylabel('coefficients');
title('non standardized, enough samples')
ax3 = plotstyle(gca, 1);
subplot(2,2,4)
h4.a = scatter(log(results4.irrelevantTaxaMeanAbundance), results4.irrelevantTaxaCoefficients, 'filled', 'blue');
hold on;
h4.b = scatter(log(results4.trueTaxaMeanAbundance), results4.trueTaxaCoefficients, 'filled', 'red');
box on;
ax4 = plotstyle(gca, 1);
xlabel('log(mean abundance of a taxon)');
ylabel('coefficients');
title('standardized, enough samples')
end