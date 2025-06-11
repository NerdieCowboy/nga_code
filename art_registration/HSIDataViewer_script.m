


% hard-coded
pc = 1; % pc(1), mac(0)

% user selection
rot = 1;    % rotate cube (1), don't rotate (0)
ff = 1;     % FF cube (1), don't FF

img_path = '';
img_fn = ''; % data filename
img_ext = '';   % data file extension
if (pc == 1)
    fn_full = [img_path '\' img_fn '.' img_ext];
elseif (pc ==0)
    fn_full = [img_path '/' img_fn '.' img_ext];
end
w = ''; % white image filename
w_ext = ''; % white file extension
if (pc == 1)
    w_full = [w_path '\' w_fn '.' w_ext];
elseif (pc ==0)
    w_full = [w_path '/' w_fn '.' w_ext];
end
n0 = 1025;    % number of scanlines

header_sz = 13107712;
x0 = 1280;  % number of columns (before rotating)
y0 = 1024;  % number of rows (before rotating)

fid = fopen(fn_full,'r');
header = fread(fid,header_sz);  %header
clear header

fidw = fopen(fn_fullw,'r');
header = fread(fidw,header_sz);  %header
clear header

if (rot == 1)
    tmp2 = fread(fidw,[x0,y0],'uint16',0,'ieee-le');
    tmp2 = tmp2';
    tmp2 = single(tmp2);
    fclose(fidw);
    
    n1 = 0;
    cube = uint16(zeros(y0,x0,n0));
    while ~feof(fid)
        tmp = fread(fid,[x0,y0],'uint16',0,'ieee-le');
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
elseif (rot == 0)
    tmp2 = fread(fidw,[y0,x0],'uint16',0,'ieee-le');
    tmp2 = single(tmp2);
    fclose(fidw);
    
    n1 = 0;
    cube = uint16(zeros(y0,x0,n0));
    while ~feof(fid)
        tmp = fread(fid,[y0,x0],'uint16',0,'ieee-le');
        if (~isempty(tmp))
            n1 = n1 + 1;
            tmp = single(tmp);

            cube(:,:,n1) = uint16((tmp./tmp2)*(2^14));
            clear tmp
        end
    end
    clear tmp2
    fclose(fid);
end

cube = shiftdim(cube,1);
[m0,n0,p0] = size(cube);

if (pc == 1)
    trial_out = [img_path '\'];
    fnh_out = [img_path '\'];
elseif (pc ==0)
    trial_out = [img_path '/'];
    fnh_out = [img_path '/'];
end
if (rot == 1)
    trial_out = [trial_out 'rot_'];
    fnh_out = [fnh_out 'rot_'];
end
if (ff == 1)
    trial_out = [trial_out 'ff_'];
    fnh_out = [fnh_out 'ff_'];
end
if ((ff==1) || (rot == 1))
    trial_out = [trial_out img_fn];
    fnh_out = [fnh_out img_fn '.hdr'];
else
    trial_out = [trial_out img_fn '1'];
    fnh_out = [fnh_out img_fn '1.hdr'];
end
multibandwrite(uint16(cube),trial_out,'bsq','machfmt','ieee-le')

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