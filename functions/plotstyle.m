%% Plotsyle
function ax = plotstyle(gca, flag)
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 2;
ax.YAxis.Color = 'k';

if flag == 1
    ax.FontSize = 16;
    ax.XLabel.FontSize = 20;
    ax.YLabel.FontSize = 20;
    ax.Title.FontSize = 20;
elseif flag == 2
    ax.FontSize = 14;
    ax.XLabel.FontSize = 18;
    ax.YLabel.FontSize = 18;
%     ax.Title.FontSize = 16;
end

end