figure,imagesc(score(:,:,13))
axis image
axis off
colormap gray

figure,imagesc(conf(:,:,13))
axis image
axis off
colormap jet
colorbar
c=colorbar
ylabel(c,'confidence')

figure,imagesc(score(:,:,13).*conf(:,:,13))
axis image
axis off
colormap gray

figure,imagesc(score(:,:,14))
axis image
axis off
colormap gray

figure,imagesc(conf(:,:,14))
axis image
axis off
colormap jet
caxis([0 1])
colorbar
c=colorbar
ylabel(c,'confidence')

figure,imagesc(score(:,:,14).*conf(:,:,14))
axis image
axis off
colormap gray