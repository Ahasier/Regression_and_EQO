function accuracy = calculateAccuracy(crossValidatedCoefficients, syntheticCoefficients, settings)
% caculateAccuracy Computes the accuracy of the model based on cross-validated coefficients
% and synthetic coefficients (true values).
% INPUT:
%   crossValidatedCoefficients: The coefficients estimated by the cross-validation procedure
%   syntheticCoefficients: The true coefficients (or synthetic ones for simulation)
%   options: A struct that may contain a custom threshold value
% OUTPUT:
%   accuracy: A scalar that indicates the accuracy of the model

% Check if a custom threshold is specified in options
if isfield(settings, 'Threshold') && ~strcmp(settings.Threshold, 'cv')
    accuracy = MCC(crossValidatedCoefficients >= settings.Threshold, syntheticCoefficients >= settings.Threshold);
else
    % If no custom threshold, use a default value of 1e-4
    accuracy = MCC(crossValidatedCoefficients >= 1e-4, syntheticCoefficients >= 1e-4);
end
end

%% MCC
function phi = MCC(B, beta)
% MCC Computes the Matthews correlation coefficient between two binary vectors
% INPUT:
%    B: A binary vector (e.g., model's prediction)
%    beta: A binary vector (e.g., true labels)
% OUTPUT:
%    phi: Matthews correlation coefficient (a scalar)

% Convert continuous values to binary (0 or 1)
bB = B ~= 0;

% Compute the elements of confusion matrix
n11 = sum(bB .* beta);
n10 = sum(bB .* (1 - beta));
n01 = sum((1 - bB) .* beta);
n00 = sum((1 - bB) .* (1 - beta));

% Compute Matthews correlation coefficient
denominator = sqrt(sum(bB) * sum(1 - bB) * sum(beta) * sum(1 - beta));

if denominator == 0
    phi = 0;
else
    phi = (n11 * n00 - n10 * n01) / denominator;
end
end
