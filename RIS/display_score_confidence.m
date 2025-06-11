i = 26;
spectra_names = {'lamp black', 'bone black', 'burnt umber', ...
    'red ochre', 'Yellow Ochre', 'gypsum', 'chalk', 'Madder lake', ...
    'Realgar', 'Malachite', 'Orpiment', 'Indigo', 'Azurite', ...
    'Red lead', 'vermilion', 'Green earth', 'verdigris', 'lead white', ...
    'ultramarine', 'Naples Yellow', 'Smalt', 'indian yellow', ...
    'Copper resinate', 'lead tin yellow', 'Van Dyke brown', ...
    'Carmine lake', 'cobalt blue', 'raw umber', 'Cd Red', ...
    'titanium white', 'Cd yellow', 'zinc white'};

tmp = score(:,:,i);
tmp2 = conf(:,:,i);
tmp3 = tmp.*tmp2;
figure,subplot(121),imagesc(tmp)
title(spectra_names{i})
colormap gray
axis image
axis off
caxis([0 1])
subplot(122),imagesc(tmp3)
colormap gray
axis image
axis off
caxis([0 1])
freezeColors
%{
subplot(122),imagesc(tmp2)
colormap jet
axis image
axis off
caxis([0 1])
%}
unfreezeColors