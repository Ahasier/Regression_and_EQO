function [groupedAbundance, addedLeaves, Idx] = groupAbundanceData(abundanceData,treeData)
addedAbundance = zeros(treeData.NumBranches - 1, size(abundanceData,2));
%% Indexing method
addedLeaves = treeData.Pointers(1:treeData.NumBranches - 1,:);
addedLeaves = addedLeaves(:);
Idx = [1:treeData.NumBranches - 1;1:treeData.NumBranches - 1]';
Idx = Idx(:);
while max(addedLeaves) > treeData.NumLeaves
    SelIdx = addedLeaves > treeData.NumLeaves;
    tp = treeData.Pointers(addedLeaves(SelIdx) - treeData.NumLeaves,:);
    addedLeaves = [addedLeaves; tp(:)];
    addedLeaves(SelIdx) = [];
    Idx = [Idx; Idx(SelIdx);Idx(SelIdx)];
    Idx(SelIdx) = [];
end
for id = 1:length(Idx)
    addedAbundance(Idx(id),:) = addedAbundance(Idx(id),:) + abundanceData(addedLeaves(id),:);
end
groupedAbundance = [abundanceData;addedAbundance];
end