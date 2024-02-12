% Generating figure 1: an illustrative diagrams showing our method.

% First, generating a toy tree with some leaves highlighted indicating 
% selection of the ground true functional group
clear
%% Generating a toy tree
toyTree = generateToyTree();
% treeData = get(tree100);

%% Genreating mock data
[numPermutations, phylogenyDependency, noiseLevel, meshGrid, settings] = initializations('LASSO');
[taxaAbundance, functionalOutput, syntheticCoefficients, selectedLeaves, treeStructure] = generateSyntheticData(toyTree, 100, 4, 1, 12, settings);

%% Plots
view(toyTree, selectedLeaves);

% figure();
% imagesc(covarianceMatrix);
% ax = plotstyle(gca, 1);
% xlabel('Taxa','fontsize',18)
% ylabel('Taxa','fontsize',18)
% title('Covariance matrix','fontsize',20)
% colorbar

figure();
imagesc(taxaAbundance);
ax = plotstyle(gca, 1);
xlabel('Samples','fontsize',18)
ylabel('Taxa','fontsize',18)
title('Mock abundance (centered)','fontsize',20)
colormap(ax, bluewhitered)
colorbar

figure()
imagesc(functionalOutput);
ax = plotstyle(gca, 1);
ylabel('Sample','fontsize',18)
title('Functional ouputs (centered)','fontsize',20)
colormap(ax, bluewhitered)
colorbar

% figure();
% imagesc(N2);
% ax = plotstyle(gca, 1);
% xlabel('Samples','fontsize',18)
% ylabel('Taxa','fontsize',18)
% title('Mock abundance after grouping','fontsize',20)
% colorbar

figure();
imagesc(results.accuracy);
ax = plotstyle(gca, 1);
title('Prediction accuracy','fontsize',20)
colormap(ax, bluewhitered)
colorbar