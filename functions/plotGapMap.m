function h = plotGapMap(Ks, nSpls, Map1, flag1, Map2, flag2, plotstyleflag)
h = imagesc(Ks, nSpls,Map1' - Map2');
xlabel('Functional group size');
ylabel('Number of samples');
if strcmp(flag1(1:3), 'phy')
    title({'Accuracy gap between ',[flag1,' and ',flag2], ['(blue means ', flag2, ' is better)']})
else
    title({['Accuracy gap between ',flag1,' and ',flag2], ['(blue means ', flag2, ' is better)']})
end
caxis([-1 1])
set(gca,'YDir','normal')
ax3 = plotstyle(gca, plotstyleflag);
colormap(ax3, bluewhitered)
colorbar
end