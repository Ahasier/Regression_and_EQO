function [index1, index2] = indicesOfAccuracyMatrixElement(numberOfTaxaInAGroup, numSamples, meshGrid)
index1 = numberOfTaxaInAGroup / meshGrid.TaxaGroup;
index2 = numSamples / meshGrid.Samples;
end