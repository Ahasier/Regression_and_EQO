function [groupedAbundance, addedLeaves, Idx] = groupAbundanceData(abundanceData, treeData)
% GROUPABUNDANCEDATA groups abundance data using the indexing method to handle
% extra phylogenetic information and incorporates it to form the grouped
% abundance table.
%
% Inputs:
%   abundanceData: The original abundance data.
%   treeData: The phylogenetic tree data.
%
% Outputs:
%   groupedAbundance: The grouped abundance data after incorporating extra phylogenetic information.
%   addedLeaves: The leaves that were added during the grouping process.
%   Idx: Indexes used for grouping.

% Initialize an array to store the added abundance data.
addedAbundance = zeros(size(abundanceData, 1), treeData.NumBranches - 1);

% Extract pointers for branches from the tree data.
addedLeaves = treeData.Pointers(1:treeData.NumBranches - 1, :);
addedLeaves = addedLeaves(:);

% Create an index array.
Idx = [1:treeData.NumBranches - 1; 1:treeData.NumBranches - 1]';
Idx = Idx(:);

% Loop to handle branches that point to other branches instead of leaves.
while max(addedLeaves) > treeData.NumLeaves
    % Find branches that point to other branches.
    SelIdx = addedLeaves > treeData.NumLeaves;
    
    % Extract pointers for those branches.
    tp = treeData.Pointers(addedLeaves(SelIdx) - treeData.NumLeaves, :);
    addedLeaves = [addedLeaves; tp(:)];
    addedLeaves(SelIdx) = [];
    
    % Update the index array.
    Idx = [Idx; Idx(SelIdx); Idx(SelIdx)];
    Idx(SelIdx) = [];
end

% Loop to add abundance data for the added leaves.
for id = 1:length(Idx)
    addedAbundance(:, Idx(id)) = addedAbundance(:, Idx(id)) + abundanceData(:, Idx(id));
end

% Combine the original and added abundance data.
groupedAbundance = [abundanceData, addedAbundance];
end