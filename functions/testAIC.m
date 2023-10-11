figure();
plot(aicValues, 'lineWidth', 2)
title('AIC Curve', 'fontsize', 18);
xlabel('Group Size', 'fontsize', 18);
ylabel('AIC Values', 'fontsize', 18);
ax = plotstyle(gca);