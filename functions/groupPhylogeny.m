function [numBranches, addedLeaves, Idx] = groupPhylogeny(treeData)
% GROUPPHYLOGENY groups phylogenetic branches using the indexing method to handle
% extra phylogenetic information. Function groupAbundanceData will then incorporates
% it to form the grouped abundance table.
%
% Inputs:
%   abundanceData: The original abundance data.
%   treeData: The phylogenetic tree data.
%
% Outputs:
%   addedLeaves: The leaves that were added during the grouping process.
%   Idx: Indexes used for grouping.

% Retrieve number of branches in the phylogenetic tree
numBranches = treeData.NumBranches;

% Extract pointers for branches from the tree data.
addedLeaves = treeData.Pointers(1:numBranches - 1, :);
addedLeaves = addedLeaves(:);

% Create an index array.
Idx = [1: numBranches - 1; 1:numBranches - 1]';
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
end