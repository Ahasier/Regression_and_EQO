function h = plotMap(Ks, nSpls, LASSOMap, flag, plotstyleflag)
h = imagesc(Ks, nSpls,LASSOMap');
xlabel('Functional group size');
ylabel('Number of samples');
title(['Accuracy of ', flag])
caxis([0 1])
set(gca,'YDir','normal')
ax1 = plotstyle(gca, plotstyleflag);
% colormap(ax, bluewhitered)
colorbar
end