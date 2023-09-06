function visualizeResults(allErrorsAtAllThresholds, allR2OutSample, beta0)
% VISUALIZERESULTS Plots squared errors and R^2 correlation coefficients
% against threshold values for different lambda values.
%
% Inputs:
% - allErrorsAtAllThresholds: Matrix of squared errors for different thresholds.
% - allR2OutSample: Matrix of R^2 values for out-of-sample data for different thresholds.
% - beta0: Initial coefficient or threshold value.

numLambdas = size(allErrorsAtAllThresholds, 1);
thresholdValues = 0:0.1:2*beta0; % As per your initialization

for l = 1:numLambdas
    figure;
    
    % Plot squared errors vs. threshold values
    subplot(2,1,1);
    plot(thresholdValues, allErrorsAtAllThresholds(l, :));
    title(['Squared Errors vs. Thresholds for Lambda ' num2str(l)]);
    xlabel('Threshold Value');
    ylabel('Squared Error');
    
    % Plot R^2 vs. threshold values
    subplot(2,1,2);
    plot(thresholdValues, allR2OutSample(l, :));
    title(['R^2 vs. Thresholds for Lambda ' num2str(l)]);
    xlabel('Threshold Value');
    ylabel('R^2 Value');
    
    % Save the plots (optional)
    saveas(gcf, ['results_lambda_' num2str(l) '.png']);
end
end
