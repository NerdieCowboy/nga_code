clear all
close all
clc

fn = 'test2';
fnw = 'white2';
delta = 16;
h_size = 33;
[~,lambda,m0,n0,p0] = read_cube(fn,'bsq',1,1,1,1,1,1);

h0 = ones(1,1,h_size);
h0 = h0/(sum(h0(:)));
cube1 = uint16(zeros(m0,n0,p0));
for m = 1:delta:m0
    m1 = m;
    m2 = m+delta-1;
    if (m2 > m0)
        m2 = m0;
    end
    for n = 1:delta:n0
        n1 = n;
        n2 = n+delta-1;
        if (n2 > n0)
            n2 = n0;
        end
        h = repmat(h0,[(m2-m1+1) (n2-n1+1) 1]);
        [cube,~,~,~,~] = read_cube(fn,'bsq',m1,m2,n1,n2,[],[]);
        cube = single(cube');
        [w,~,~,~,~] = read_cube(fnw,'bsq',m1,m2,n1,n2,[],[]);
        w = single(w');
        cube1 = uint16((cube./w)*(2^14-1));
        clear cube w
    end
end

write_cube(cube1,fn,['rot_ff_' fn '.bsq'],'uint16')

fnh_out = ['rot_ff_' fn '.hdr'];
fid1 = fopen(fnh_out,'w');
fprintf(fid1,'ENVI\n');
fprintf(fid1,'description = {}\n');
fprintf(fid1,'samples = %u\n',n0);
fprintf(fid1,'lines = %u\n',m0);
fprintf(fid1,'bands = %u\n',p0);
fprintf(fid1,'header offset = 0\n');
fprintf(fid1,'file type = ENVI Standard\n');
fprintf(fid1,'data type = 12\n');
fprintf(fid1,'interleave = bsq\n');
fprintf(fid1,'byte order = 0\n');
fclose(fid1);

