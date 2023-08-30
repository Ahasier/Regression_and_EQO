function h = plotMap(Ks, nSpls, LASSOMap, flag)
h = imagesc(Ks, nSpls,LASSOMap');
xlabel('Number of Taxa Selected for \beta = 1');
ylabel('Number of Samples');
title(['Accuracy of ', flag])
caxis([0 1])
set(gca,'YDir','normal')
ax1 = plotstyle(gca);
% colormap(ax, bluewhitered)
colorbar
end