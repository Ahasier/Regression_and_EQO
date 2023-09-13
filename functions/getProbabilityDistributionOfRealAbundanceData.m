function [countRowSums, countColSums, probRowSums, probColSums, binRowSums, binColSums] = getProbabilityDistributionOfRealAbundanceData(realData, binWidth)
% Extract row and column sums from real data
realRowSums = sum(realData, 2);
realColSums = sum(realData, 1);

[countRowSums, binRowSums] = histcounts(realRowSums, 'BinWidth', binWidth);
probRowSums = countRowSums/sum(countRowSums);
binRowSums(1) = [];

[countColSums, binColSums] = histcounts(realColSums, 'BinWidth', binWidth);
probColSums = countColSums/sum(countColSums);
binColSums(1) = [];
end