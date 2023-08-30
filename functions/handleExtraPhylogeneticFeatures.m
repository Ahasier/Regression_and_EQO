% If phylogenetic information is incorporated, update estimated coefficients using extra features. Currently unneeded.
function updatedBeta = handleExtraPhylogeneticFeatures(AllBeta, TaraNames, Idx, Ladd)
if useExtraFeatures
    upAllBeta = AllBeta(1:length(TaraNames));
    for id = 1:length(Idx)
        upAllBeta(Ladd(id)) = upAllBeta(Ladd(id)) + AllBeta(length(TaraNames) + Idx(id));
    end
    updatedBeta = upAllBeta;
else
    updatedBeta = AllBeta;
end
end