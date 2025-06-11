%% Extract IR subimages
function [cube] = cube_rotate_ff(ir_fn,w_fn)

tmp = regexp(ir_fn,'\\');
if (isempty(tmp))
    tmp = regexp(ir_fn,'\/');
end
tmp1 = tmp(end);
tmp2 = length(ir_fn);
fn_path = ir_fn(1:tmp1);
fn = ir_fn((tmp1+1):tmp2);
clear tmp tmp1 tmp2
tmp = regexp(fn,'\.');
trial = fn(1:(tmp-1));
clear tmp

tmp = regexp(w_fn,'\\');
if (isempty(tmp))
    tmp = regexp(w_fn,'\/');
end
tmp1 = tmp(end);
tmp2 = length(w_fn);
fn_pathw = w_fn(1:tmp1);
fnw = w_fn((tmp1+1):tmp2);
clear tmp tmp1 tmp2
tmp = regexp(fnw,'\.');
w = fnw(1:(tmp-1));
clear tmp

fn_full = [fn_path trial '.img'];
fn_fullw = [fn_pathw w '.img'];
fid = fopen(fn_full,'r');
%{
header = fread(fid,2);  %header
m0 = fread(fid,[1,1],'uint16',0,'ieee-le');
p0 = fread(fid,[1,1],'uint16',0,'ieee-le');
header = fread(fid,13107706);  %header
%}
header = fread(fid,13107712);  %header
clear header
m0 = 1280;
n0 = 1023;
p0 = 1024;

fidw = fopen(fn_fullw,'r');
header = fread(fidw,13107712);  %header
clear header

n1 = 0;
cube = uint16(zeros(p0,m0,n0));
while ~feof(fid)
    tmp = fread(fid,[m0,p0],'uint16',0,'ieee-le');
    tmp2 = fread(fidw,[m0,p0],'uint16',0,'ieee-le');
    if (~isempty(tmp))
        n1 = n1 + 1;
        tmp = tmp';
        tmp = single(tmp);
        tmp2 = tmp2';
        tmp2 = single(tmp2);
        
        %cube(:,:,n1) = uint16((tmp)*(2^14));
        cube(:,:,n1) = uint16((tmp./tmp2)*(2^14));
        clear tmp tmp2
    end
end
fclose(fid);
%fclose(fidw);
cube = shiftdim(cube,1);
[m0,n0,p0] = size(cube);
multibandwrite(uint16(cube),[fn_path trial '1'],'bsq','machfmt','ieee-le')

fnh_out = [fn_path trial '1.hdr'];
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

%clear cube

end
