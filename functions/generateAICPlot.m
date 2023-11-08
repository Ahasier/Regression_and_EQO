function generateAICPlot(aicValues, optimalGroupSize)
figure();
h = plot(aicValues, 'linewidth', 2);
xlim([0 optimalGroupSize*2])
xlabel('Group size', 'fontsize', 18);
ylabel('AIC value', 'fontsize', 18);
ax = plotstyle(gca);
end