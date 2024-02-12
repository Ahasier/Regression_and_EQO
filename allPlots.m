numberOfTaxaInAGroup_list = 5:5:50;
numSamples_list = 10:10:200;
numSamples_list2 = 10:10:200;

% Initialize all neccessary parametters from configurations files
[numPermutations, phylogenyDependency, noiseLevel, meshGrid, settings] = initializations('OLS');

% % Figure 2: Comparison between LASSO and EQO on binary ground truth
% Figure2(numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% 
% % Figure 3: Accuracy Change on Non-binary Ground Truth
% Figure3(numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
% 
% % Figure 4: Comparison between LASSO and EQO on non-binary ground truth and
% % mimic-real abundance data
% Figure4(numberOfTaxaInAGroup_list, numSamples_list2, meshGrid);
% 
% % Figure 6: Comparison between phylo-regularized LASSO and EQO on 
% % phylogenetcally-dependent non-binary ground truth and mimic-real 
% % abundance data
% Figure6(numberOfTaxaInAGroup_list, numSamples_list2, meshGrid);
% 
% % Figure S2: Comparison among LASSO and EQO with different noise levels.
% FigureS2(numberOfTaxaInAGroup_list, numSamples_list, meshGrid);

groupsize_list = [10 30 50];
samplesize_list = [150 100 50];
figure('Renderer', 'painters', 'Position', [0 0 1600 1600]);
for n = 1:3
    groupsize = groupsize_list(n);
    for m = 1:3
        samplesize = samplesize_list(m);
        
        subplot(3,3,3*(m-1)+n)
        [hScatterLASSO, hScatterEQO] = plotVaryEpsilon(groupsize,samplesize);
        xlim([0 0.5])
    end
end

%% Helper functions
function Map = loadMapFromFile(regressionMethod, fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid)
Map = accessToAccMap('results/Accuracy/', regressionMethod, fullIdentifier, 'Load');
Map = Map(numberOfTaxaInAGroup_list/meshGrid.TaxaGroup, numSamples_list/meshGrid.Samples);
end

function resultMap = accessToAccMap(resultsPath, regressionMethod, fullIdentifier, saveOrLoad)
filename = [resultsPath, 'Acc', regressionMethod, fullIdentifier, '.csv'];
if strcmp(saveOrLoad, 'Load')
    resultMap = csvread(filename);
end
end

function Figure2(numberOfTaxaInAGroup_list, numSamples_list, meshGrid)
fullIdentifierF2 = '_Beta01_BetaEps0_TrNaN';
fullIdentifierF2_2 = '_Beta01_BetaEps0_TrNaN_maxLambda10';
EQOMapF2 = loadMapFromFile('EQO', fullIdentifierF2, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
LASSOMapF2 = loadMapFromFile('LASSO', fullIdentifierF2_2, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
figure('Renderer', 'painters', 'Position', [0 0 2000 1000]);
subplot(2,4,1)
h0 = plotBeta('binary');
subplot(2,4,2)
h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMapF2, 'LASSO', 1);
subplot(2,4,3)
h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMapF2, 'EQO', 1);
subplot(2,4,4)
h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMapF2, 'LASSO', EQOMapF2, 'EQO', 1);
end

function Figure3(numberOfTaxaInAGroup_list, numSamples_list, meshGrid)
fullIdentifierF3 = '_Beta01_BetaEps0.5_TrNaN';
fullIdentifierF3_2 = '_Beta01_BetaEps0.5_TrNaN_maxLambda10';
EQOMapF3 = loadMapFromFile('EQO', fullIdentifierF3, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
LASSOMapF3 = loadMapFromFile('LASSO', fullIdentifierF3_2, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
subplot(2,4,5)
h0 = plotBeta('non-binary');
subplot(2,4,6)
h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMapF3, 'LASSO', 1);
subplot(2,4,7)
h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMapF3, 'EQO', 1);
subplot(2,4,8)
[hScatterOLS, hScatterEQO] = plotVaryEpsilon(30,110);
xlim([0 0.5])
end

function Figure4(numberOfTaxaInAGroup_list, numSamples_list2, meshGrid)
fullIdentifierF4 = '_Beta01_BetaEps0.5_TrNaN_maxLambda10_RealAbdOn_usePhylogenyOn_convexd8'; 
fullIdentifierF4_2 = '_Beta01_BetaEps0.5_TrNaN_maxLambda10_RealAbdOn_usePhylogenyOff_convexd8'; 

EQOMapF4 = loadMapFromFile('LASSO', fullIdentifierF4, numberOfTaxaInAGroup_list, numSamples_list2, meshGrid);
LASSOMapF4 = loadMapFromFile('LASSO', fullIdentifierF4_2, numberOfTaxaInAGroup_list, numSamples_list2, meshGrid);
load('results/RealAbundance.mat');
figure('Renderer', 'painters', 'Position', [0 0 1600 320]);
% subplot(1,4,1)
% h0 = imagesc(taxaAbundance);
% xlabel('Functional group size','fontsize',16);
% ylabel('Number of samples','fontsize',16);
% title('Realistic taxa abundance (centered)')
ax3 = plotstyle(gca, 1);
colormap(ax3, bluewhitered)
subplot(1,3,1)
h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list2, EQOMapF4, 'LASS0 on Phylo-dependent data', 1);
subplot(1,3,2)
h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list2, LASSOMapF4, 'LASSO on Phylo-independent data', 1);
subplot(1,3,3)
h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list2, EQOMapF4, 'LASS0 on Phylo-dependent data', LASSOMapF4, 'Phylo-independent data', 1);
end

function Figure6(numberOfTaxaInAGroup_list, numSamples_list2, meshGrid)
fullIdentifierF6_1 = '_Beta01_BetaEps0.5_TrNaN_RealAbdOn_usePhylogenyOff';
fullIdentifierF6_2 = '_Beta01_BetaEps0.5_TrNaN_maxLambda10_RealAbdOn_usePhylogenyOn_PhyloIndependent';
fullIdentifierF6_3 = '_Beta01_BetaEps0.5_TrNaN_maxLambda10_RealAbdOn_usePhylogenyOn_PhyloD1B';
EQOMapF6 = loadMapFromFile('EQO', fullIdentifierF6_1, numberOfTaxaInAGroup_list, numSamples_list2, meshGrid);
LASSOMapF6_2 = loadMapFromFile('LASSO', fullIdentifierF6_2, numberOfTaxaInAGroup_list, numSamples_list2, meshGrid);
LASSOMapF6_3 = loadMapFromFile('LASSO', fullIdentifierF6_3, numberOfTaxaInAGroup_list, numSamples_list2, meshGrid);
figure('Renderer', 'painters', 'Position', [0 0 1500 640]);
subplot(2,3,1)
h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list2, LASSOMapF6_2, 'LASSO', 2);
subplot(2,3,2)
h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list2, LASSOMapF6_3, 'phylo-regularized LASSO', 2);
subplot(2,3,3)
h3 = plotMap(numberOfTaxaInAGroup_list, numSamples_list2, EQOMapF6, 'EQO', 2);
subplot(2,3,4)
h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list2, LASSOMapF6_3, 'phylo-regularized LASSO', LASSOMapF6_2, 'LASSO', 2);
subplot(2,3,5)
h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list2, LASSOMapF6_2, 'LASSO', EQOMapF6, 'EQO', 2);
subplot(2,3,6)
h3 = plotGapMap(numberOfTaxaInAGroup_list, numSamples_list2, LASSOMapF6_3, 'phylo-regularized LASSO', EQOMapF6, 'EQO', 2);
end

function FigureS2(numberOfTaxaInAGroup_list, numSamples_list, meshGrid)
fullIdentifierFS2_1 = '_Beta01_BetaEps0_TrNaN_maxLambda10_noise01';
fullIdentifierFS2_2 = '_Beta01_BetaEps0_TrNaN_maxLambda10';
fullIdentifierFS2_3 = '_Beta01_BetaEps0_TrNaN_maxLambda10_noise10';
fullIdentifierFS2_4 = '_Beta01_BetaEps0_TrNaN_noise01';
fullIdentifierFS2_5 = '_Beta01_BetaEps0_TrNaN';
fullIdentifierFS2_6 = '_Beta01_BetaEps0_TrNaN_noise10';
LASSOMapNoise01 = loadMapFromFile('LASSO', fullIdentifierFS2_1, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
LASSOMapNoise1 = loadMapFromFile('LASSO', fullIdentifierFS2_2, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
LASSOMapNoise10 = loadMapFromFile('LASSO', fullIdentifierFS2_3, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
EQOMapNoise01 = loadMapFromFile('EQO', fullIdentifierFS2_4, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
EQOMapNoise1 = loadMapFromFile('EQO', fullIdentifierFS2_5, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
EQOMapNoise10 = loadMapFromFile('EQO', fullIdentifierFS2_6, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);
figure('Renderer', 'painters', 'Position', [0 0 1500 640]);
subplot(2,3,1)
h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMapNoise01, 'LASSO (noise level 0.1)', 1);
subplot(2,3,2)
h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMapNoise1, 'LASSO (noise level 1)', 1);
subplot(2,3,3)
h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, LASSOMapNoise10, 'LASSO (noise level 10)', 1);
subplot(2,3,4)
h1 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMapNoise01, 'EQO (noise level 0.1)', 1);
subplot(2,3,5)
h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMapNoise1, 'EQO (noise level 1)', 1);
subplot(2,3,6)
h2 = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMapNoise10, 'EQO (noise level 10)', 1);
end

function h = plotBeta(flag)
dx = 0.001;
x = 0:dx:2;
% lenx = length(x);
if strcmp(flag, 'binary')
    % Define the values and counts
    values = [0, 1];
    counts = [70, 30];
    
    % Create a bar plot for the histogram
    h = bar(values, counts, 'BarWidth', 0.2, 'FaceColor', 'cyan', 'FaceAlpha', 0.3, 'EdgeColor', '#0072BD', 'LineWidth', 2);
    set(gca, 'XTick', [0 1], 'YLim', [0 80]);
    ylabel('Count');
    xlim([-0.5 1.5])
    ylim([0 80]);
%     beta = zeros(lenx,1);
%     beta(1) = 5e7;
%     beta(1/dx + 1) = 5e7;
%     
%     hold on;
%     h = plot(x, beta, 'linewidth', 2);
elseif strcmp(flag, 'non-binary')
    mu = 1;
    sigma = 0.5;
    beta = (0.5/0.8410).*(1/(sigma * sqrt(2 * pi)))*exp(- 0.5*((x - mu) / sigma).^2);
    beta(1) = 5e7;
    beta(2:0.5/dx)  = 0;
    
    hold on;
    h = plot(x, beta, 'linewidth', 2);
    
    % Add lines for mean and variance
    yLimits = ylim();
    line([mu mu], yLimits, 'Color', 'red', 'LineStyle', '--', 'linewidth', 2);
    line([mu-sigma mu-sigma], yLimits, 'Color', 'blue', 'LineStyle', '--', 'linewidth', 2);
    line([mu+sigma mu+sigma], yLimits, 'Color', 'blue', 'LineStyle', '--');
    
    % Fill under curve
    area(x, beta, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    % Legend
    legend('Distribution', 'Mean (\mu)', 'Standard Deviation (\sigma)', 'Location', 'NorthEast');
    
    % Set y-label and y-limit
    ylabel('Probability density');
    ylim([0 1]);
end
% xlim([-0.5, 2])
xlabel('Ground true \beta');
title('Probability distribution of ground truth');
ax = plotstyle(gca, 1);

hold off;
end

function [hScatterLASSO, hScatterEQO] = plotVaryEpsilon(groupsize, samplesize)
numberOfTaxaInAGroup = 20;
numSamples = 120;
betaEps_list = 0:0.01:0.5;

% Get estimation accuracies from ground truth beta with different epsilon
accuraciesLASSO = load(['results/Accuracy/accuraciesVaryingEps_LASSO_Beta01_BetaEps0_TrNaN_maxLambda10_K',num2str(groupsize),'_nSpl',num2str(samplesize),'.mat']).accuracies;
accuraciesEQO = load(['results/Accuracy/accuraciesVaryingEps_EQO_Beta01_BetaEps0_TrNaN_K',num2str(groupsize),'_nSpl',num2str(samplesize),'.mat']).accuracies;

% Define colors for the plots
colorLASSO = [0.2, 0.6, 0.8]; % Shade of blue for LASSO
colorEQO = [0.8, 0.2, 0.2]; % Shade of red for EQO

% --- LASSO Scatter plot and fit ---
hScatterLASSO = scatter(betaEps_list, accuraciesLASSO, 50, colorLASSO, 'filled', 'DisplayName', 'LASSO');
hold on;

% Linear fit for LASSO
pLASSO = polyfit(betaEps_list, accuraciesLASSO, 1);
x_fit = linspace(min(betaEps_list), max(betaEps_list), 100);
y_fit_LASSO = polyval(pLASSO, x_fit);
plot(x_fit, y_fit_LASSO, '-', 'Color', colorLASSO, 'LineWidth', 2);

% Error range for LASSO
std_dev_LASSO = std(accuraciesLASSO - polyval(pLASSO, betaEps_list));
y_upper_LASSO = y_fit_LASSO + std_dev_LASSO;
y_lower_LASSO = y_fit_LASSO - std_dev_LASSO;
fillColorLASSO = colorLASSO + [0.4, 0.4, 0.4];
fillColorLASSO(fillColorLASSO > 1) = 1;
fill([x_fit, fliplr(x_fit)], [y_upper_LASSO, fliplr(y_lower_LASSO)], fillColorLASSO, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% --- EQO Scatter plot and fit ---
hScatterEQO = scatter(betaEps_list, accuraciesEQO, 50, colorEQO, 'filled', 'DisplayName', 'EQO');

% Linear fit for EQO
pEQO = polyfit(betaEps_list, accuraciesEQO, 1);
y_fit_EQO = polyval(pEQO, x_fit);
plot(x_fit, y_fit_EQO, '-', 'Color', colorEQO, 'LineWidth', 2);

% Error range for EQO
std_dev_EQO = std(accuraciesEQO - polyval(pEQO, betaEps_list));
y_upper_EQO = y_fit_EQO + std_dev_EQO;
y_lower_EQO = y_fit_EQO - std_dev_EQO;
fillColorEQO = colorEQO + [0.4, 0.4, 0.4];
fillColorEQO(fillColorEQO > 1) = 1;
fill([x_fit, fliplr(x_fit)], [y_upper_EQO, fliplr(y_lower_EQO)], fillColorEQO, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel('\sigma_\beta');
ylabel('Accuracy');
title(['True group size = ',num2str(groupsize),', ',num2str(samplesize),' samples'], 'fontsize', 24)
box on;
ax = plotstyle(gca, 1);
legend([hScatterLASSO, hScatterEQO], 'Location', 'best'); % Display legends

% hold off;
end