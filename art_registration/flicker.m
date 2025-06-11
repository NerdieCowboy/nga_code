trial = 'gourmet';
[cube,lambda] = read_cube('Gourmet1_full_rot180crop');

m0 = 1599;   %lines (y)
n0 = 1145;   %samples(x)
p0 = 215;    %bands

fig=figure;
set(fig,'DoubleBuffer','on');
fn = strcat(trial,'_movie');    %AVI filename
aviobj = avifile(fn,'compression','None');
aviobj.fps = 10;    %frame rate (frames/sec)

%forward
for i = 1:215
    imshow(cube(:,:,i))
    F = getframe(fig);
    aviobj = addframe(aviobj,F);
    pause(0.1)
end

%backward
for i = 214:-1:2
    imshow(cube(:,:,i))
    F = getframe(fig);
    aviobj = addframe(aviobj,F);
    pause(0.1)
end

close(fig)
aviobj = close(aviobj);
