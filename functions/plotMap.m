function h = plotMap(Ks, nSpls, LASSOMap, flag, plotstyleflag)
h = imagesc(Ks, nSpls,LASSOMap');

% Set the missing data pixels to grey
h.AlphaData = ~isnan(LASSOMap');
ax = gca; ax.Color = 0.8*[1 1 1];

xlabel('Functional group size');
ylabel('Number of samples');
title(['Accuracy of ', flag])
caxis([0 1])
set(gca,'YDir','normal')
ax1 = plotstyle(gca, plotstyleflag);
% colormap(ax, bluewhitered)
colorbar
end