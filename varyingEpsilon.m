numberOfTaxaInAGroup = 20;
numSamples = 120;
betaEps_list = 0:0.01:0.5;

% Get estimation accuracies from ground truth beta with different epsilon
% accuraciesOLS = runWithVaryingEpsilon(numberOfTaxaInAGroup, numSamples, 'OLS', betaEps_list);

% Plot accuracy curve
figure();

% Define colors for the plots
colorOLS = [0.2, 0.6, 0.8]; % Shade of blue for OLS
colorEQO = [0.8, 0.2, 0.2]; % Shade of red for EQO

% --- OLS Scatter plot and fit ---
scatter(betaEps_list, accuraciesOLS, 50, colorOLS, 'filled', 'DisplayName', 'OLS');
hold on;

% Linear fit for OLS
pOLS = polyfit(betaEps_list, accuraciesOLS, 1);
x_fit = linspace(min(betaEps_list), max(betaEps_list), 100);
y_fit_OLS = polyval(pOLS, x_fit);
plot(x_fit, y_fit_OLS, '-', 'Color', colorOLS, 'LineWidth', 2);

% Error range for OLS
std_dev_OLS = std(accuraciesOLS - polyval(pOLS, betaEps_list));
y_upper_OLS = y_fit_OLS + std_dev_OLS;
y_lower_OLS = y_fit_OLS - std_dev_OLS;
fillColorOLS = colorOLS + [0.4, 0.4, 0.4];
fillColorOLS(fillColorOLS > 1) = 1;
fill([x_fit, fliplr(x_fit)], [y_upper_OLS, fliplr(y_lower_OLS)], fillColorOLS, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% --- EQO Scatter plot and fit ---
scatter(betaEps_list, accuracies, 50, colorEQO, 'filled', 'DisplayName', 'EQO');

% Linear fit for EQO
pEQO = polyfit(betaEps_list, accuracies, 1);
y_fit_EQO = polyval(pEQO, x_fit);
plot(x_fit, y_fit_EQO, '-', 'Color', colorEQO, 'LineWidth', 2);

% Error range for EQO
std_dev_EQO = std(accuracies - polyval(pEQO, betaEps_list));
y_upper_EQO = y_fit_EQO + std_dev_EQO;
y_lower_EQO = y_fit_EQO - std_dev_EQO;
fillColorEQO = colorEQO + [0.4, 0.4, 0.4];
fillColorEQO(fillColorEQO > 1) = 1;
fill([x_fit, fliplr(x_fit)], [y_upper_EQO, fliplr(y_lower_EQO)], fillColorEQO, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel('\epsilon', 'fontsize', 16);
ylabel('Accuracy', 'fontsize', 16);
box on;
ax = plotstyle(gca);
legend('Location', 'best'); % Display legends

hold off;