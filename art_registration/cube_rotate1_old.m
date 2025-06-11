%% Extract IR subimages
function [] = cube_rotate1(ir_fn)

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

fn_full = [fn_path trial '.img'];
fid = fopen(fn_full,'r');
header = fread(fid,2);  %header
m0 = fread(fid,[1,1],'uint16',0,'ieee-le');
n0 = fread(fid,[1,1],'uint16',0,'ieee-le');
header = fread(fid,506);  %header
%missing 512 Bytes
clear header

cube = uint16(zeros(m0,n0,1024));
p1 = 0;
while ~feof(fid)
    tmp = fread(fid,[m0,n0],'uint16',0,'ieee-le');
    if (~isempty(tmp))
        p1 = p1 + 1;
        
        %tmp = tmp';
        %tmp = tmp(end:-1:1,:);
        
        [~,~,p0] = size(cube);
        if (p1 > p0)
            cube = cat(3,cube,zeros(m0,n0,128));
        end
        cube(:,:,p1) = tmp;
        clear tmp
    end
end
fclose(fid);
[~,~,p0] = size(cube);
cube(:,:,(p1+1):p0) = [];
cube = shiftdim(cube,2);
[m0,n0,p0] = size(cube);

inc = 1;
cube = single(cube);
cube2 = single(zeros(m0,n0,floor(p0/3)));
for p = 3:3:p0
    cube2(:,:,inc) = single(mean(cube(:,:,(p-2):p),3));
    inc = inc + 1;
end
clear cube
multibandwrite(uint16(cube2),[fn_path trial],'bsq','machfmt','ieee-le')
[m0,n0,p0] = size(cube2);
clear cube2

fnh_out = [fn_path trial '.hdr'];
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
