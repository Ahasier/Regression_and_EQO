function centeredMatrix = centerMatrix(M)
% centerMatrix centers a matrix such that the mean of each column is zero.
%
% Inputs:
%   M - A matrix to be centered.
%
% Outputs:
%   centeredMatrix - The centered matrix.

% Calculate the mean of each column
colMeans = mean(M);

% Subtract the column means from each element in the corresponding column
centeredMatrix = M - colMeans;
end