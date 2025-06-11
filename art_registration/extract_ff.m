%% Extract IR subimages
function [] = extract_ff(fn,w,scanwidth)

h_sz = 388;
m_vnir = 640;
n_vnir = scanwidth;
p_vnir = 320;
m_swir = 640;
n_swir = scanwidth;
p_swir = 400;

fid = fopen([fn '.vxs'],'r');
header = fread(fid,h_sz);  %header
clear header

vnir = uint16(zeros(p_vnir,m_vnir,n_vnir));
n1 = 1;
while (n1 <= n_vnir)
    tmp = fread(fid,[p_vnir,m_vnir],'ushort',0,'ieee-le');
    vnir(:,:,n1) = uint16(tmp);
    clear tmp
    n1 = n1 + 1;
end

swir = uint16(zeros(p_swir,m_swir,n_swir));
n1 = 1;
while (n1 <= n_swir)
    tmp = fread(fid,[p_swir,m_swir],'ushort',0,'ieee-le');
    swir(:,:,n1) = uint16(tmp);
    clear tmp
    n1 = n1 + 1;
end

vnir_dk = fread(fid,[p_vnir,m_vnir],'ushort',0,'ieee-le');
vnir = vnir - uint16(repmat(vnir_dk,[1,1,n_vnir]));
clear vnir_dk

swir_dk = fread(fid,[p_swir,m_swir],'ushort',0,'ieee-le');
fclose(fid);
swir = swir - uint16(repmat(swir_dk,[1,1,n_swir]));
clear swir_dk

fid = fopen([w '.vxs'],'r');
header = fread(fid,h_sz);  %header
clear header

vnir_w = uint16(zeros(p_vnir,m_vnir,n_vnir));
n1 = 1;
while (n1 <= n_vnir)
    tmp = fread(fid,[p_vnir,m_vnir],'ushort',0,'ieee-le');
    vnir_w(:,:,n1) = uint16(tmp);
    clear tmp
    n1 = n1 + 1;
end

swir_w = uint16(zeros(p_swir,m_swir,n_swir));
n1 = 1;
while (n1 <= n_swir)
    tmp = fread(fid,[p_swir,m_swir],'ushort',0,'ieee-le');
    swir_w(:,:,n1) = uint16(tmp);
    clear tmp
    n1 = n1 + 1;
end

vnir_w_dk = fread(fid,[p_vnir,m_vnir],'ushort',0,'ieee-le');
vnir_w = vnir_w - uint16(repmat(vnir_w_dk,[1,1,n_vnir]));
clear vnir_w_dk
vnir = single(vnir);
vnir = vnir./single(vnir_w);
clear vnir_w
vnir = uint16(vnir*((2^14)-1));
vnir = shiftdim(vnir,1);
multibandwrite(vnir,['ff_vnir_' fn],'bsq');
clear vnir

swir_w_dk = fread(fid,[p_swir,m_swir],'ushort',0,'ieee-le');
fclose(fid);
swir_w = swir_w - uint16(repmat(swir_w_dk,[1,1,n_swir]));
clear swir_w_dk
swir = single(swir);
swir = swir./single(swir_w);
clear swir_w
swir = uint16(swir*((2^14)-1));
swir = shiftdim(swir,1);
multibandwrite(swir,['ff_swir_' fn],'bsq');
clear swir
%{
vnir = shiftdim(vnir,1);
multibandwrite(vnir,['ff_vnir_' fn],'bsq');
swir = shiftdim(swir,1);
multibandwrite(swir,['ff_swir_' fn],'bsq');
%}

fnh_out = ['ff_vnir_' fn '.hdr'];
fid1 = fopen(fnh_out,'w');
fprintf(fid1,'ENVI\n');
fprintf(fid1,'description = {\n');
fprintf(fid1,' File Imported into ENVI.}\n');
fprintf(fid1,'samples = %u\n',n_vnir);
fprintf(fid1,'lines = %u\n',m_vnir);
fprintf(fid1,'bands = %u\n',p_vnir);
fprintf(fid1,'header offset = 0\n');
fprintf(fid1,'file type = ENVI Standard\n');
fprintf(fid1,'data type = 12\n');
fprintf(fid1,'interleave = bsq\n');
fprintf(fid1,'byte order = 0\n');
fclose(fid1);

fnh_out = ['ff_swir_' fn '.hdr'];
fid1 = fopen(fnh_out,'w');
fprintf(fid1,'ENVI\n');
fprintf(fid1,'description = {\n');
fprintf(fid1,' File Imported into ENVI.}\n');
fprintf(fid1,'samples = %u\n',n_swir);
fprintf(fid1,'lines = %u\n',m_swir);
fprintf(fid1,'bands = %u\n',p_swir);
fprintf(fid1,'header offset = 0\n');
fprintf(fid1,'file type = ENVI Standard\n');
fprintf(fid1,'data type = 12\n');
fprintf(fid1,'interleave = bsq\n');
fprintf(fid1,'byte order = 0\n');
fclose(fid1);

end
