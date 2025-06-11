%% Extract IR subimages
function [cube] = cal_xswir(trial,n1,n2)

%{
trial = 'xe_83msec_256avg_2';
% spatial crop coordinates
n1 = 300;
n2 = 999;
cal_xswir(trial,n1,n2);
%}

fn_full = [trial '.img'];
fid = fopen(fn_full,'r');

header = fread(fid,13107712);  %header
clear header
m0 = 1280;
p0 = 1024;

cube = fread(fid,[m0,p0],'uint16',0,'ieee-le');
cube = cube';
cube = single(cube);
cube = cube(:,n1:n2);
mu = mean(cube,2);
[m0,n0] = size(mu);
fclose(fid);
multibandwrite(single(mu),[trial '1'],'bsq','machfmt','ieee-le')

fnh_out = [trial '1.hdr'];
fid1 = fopen(fnh_out,'w');
fprintf(fid1,'ENVI\n');
fprintf(fid1,'description = {}\n');
fprintf(fid1,'samples = %u\n',n0);
fprintf(fid1,'lines = %u\n',m0);
fprintf(fid1,'bands = %u\n',1);
fprintf(fid1,'header offset = 0\n');
fprintf(fid1,'file type = ENVI Standard\n');
fprintf(fid1,'data type = 4\n');
fprintf(fid1,'interleave = bsq\n');
fprintf(fid1,'byte order = 0\n');
fclose(fid1);

end
