function h = plotGapMap(Ks, nSpls, Map1, flag1, Map2, flag2)
h = imagesc(Ks, nSpls,Map1' - Map2');
xlabel('Number of Taxa Selected for \beta = 1');
ylabel('Number of Samples');
title(['Accuracy gap between ',flag1,' and ',flag2])
caxis([-1 1])
set(gca,'YDir','normal')
ax3 = plotstyle(gca);
colormap(ax3, bluewhitered)
colorbar
end