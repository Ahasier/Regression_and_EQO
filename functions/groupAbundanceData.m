function [groupedAbundance] = groupAbundanceData(abundanceData, numBranches, addedLeaves, Idx)
% GROUPABUNDANCEDATA groups phylogenetic branches using the indexing method to handle
% extra phylogenetic information. Function groupAbundanceData will then incorporates
% it to form the grouped abundance table.
%
% Inputs:
%   abundanceData: The original abundance data.
%   addedLeaves: The leaves that were added during the grouping process.
%   Idx: Indexes used for grouping.
%
% Outputs:
%   groupedAbundance: The grouped abundance data after incorporating extra phylogenetic information.

% Initialize an array to store the added abundance data.
addedAbundance = zeros(size(abundanceData, 1), numBranches - 1);

% Loop to add abundance data for the added leaves.
for id = 1:length(Idx)
    addedAbundance(:, Idx(id)) = addedAbundance(:, Idx(id)) + abundanceData(:, addedLeaves(id));
end

% Combine the original and added abundance data.
groupedAbundance = [abundanceData, addedAbundance];
end