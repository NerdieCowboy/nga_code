%% Extract IR subimages
function [cube] = cube_rotate_ff4(trial,w,n0)

%{
trial = 'tp3_CaSO4_83msec_1849';
w = 'white_83msec_128avg_1837';
n0 = 1025;  n0 -> # of scanlines
[cube] = cube_rotate_ff4(trial,w,n0);
%}

fn_full = [trial '.img'];
fn_fullw = [w '.img'];
fid = fopen(fn_full,'r');

header = fread(fid,13107712);  %header
clear header
m0 = 1280;
p0 = 1024;

fidw = fopen(fn_fullw,'r');
header = fread(fidw,13107712);  %header
clear header

tmp2 = fread(fidw,[m0,p0],'uint16',0,'ieee-le');
tmp2 = tmp2';
tmp2 = single(tmp2);
fclose(fidw);
        
n1 = 0;
cube = uint16(zeros(p0,m0,n0));
while ~feof(fid)
    tmp = fread(fid,[m0,p0],'uint16',0,'ieee-le');
    if (~isempty(tmp))
        n1 = n1 + 1;
        tmp = tmp';
        tmp = single(tmp);

        cube(:,:,n1) = uint16((tmp./tmp2)*(2^14));
        clear tmp
    end
end
clear tmp2
fclose(fid);
cube = shiftdim(cube,1);
[m0,n0,p0] = size(cube);
multibandwrite(uint16(cube),[trial '1'],'bsq','machfmt','ieee-le')

fnh_out = [trial '1.hdr'];
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

end
