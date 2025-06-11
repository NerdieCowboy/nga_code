function [timg] = apply_transform(fn,blknum0,m1,n1,Ntran,init_reg,Xstar,feather_width,xrf_mode,createavi)
%% Transform repeat function
% apply_transform.m
% Given that a transform has been applied to an image (image 1) in order to 
%   register it with another image (image 0). This function applies that same 
%   transform to a second image (image 2). If image 1 and image 2 were 
%   originally registered, then function will register image2 with image 0.
% Author: Damon Conover
% Email: dconover@gwmail.gwu.edu
% Latest Revision: 26 March 2014
%
% Usage: apply_transform(fn,init_reg,Xstar)
%   inputs: fn -> full path to the file being transformed
%           m1 -> full image size (vertical)
%           n1 -> full image size (horizontal)
%           Ntran -> order of the polynomial transformation (0 -> bilinear)
%           init_reg -> full path to init_reg.csv file
%           Xstar -> full path to Xstar.csv file
%           createavi -> create AVI movie of transform steps (0-> off, 1->on)
%   output:    transformed image
%

warning('off', 'all')
clear mex

turnonplot = 0;

% extract path and trial name
tmp = regexp(fn,'\\');
pc = 1;
if (isempty(tmp))
    pc = 0;
    tmp = regexp(fn,'\/');
end
tmp1 = tmp(end);
tmp2 = length(fn);
fn_path = fn(1:tmp1);
fn1 = fn((tmp1+1):tmp2);
clear tmp tmp1 tmp2
tmp = regexp(fn1,'\.');
trial = fn1(1:(tmp-1));
clear tmp
if (pc == 1)
    fn_path2 = [fn_path trial '\'];
elseif (pc == 0)
    fn_path2 = [fn_path '/' trial '/'];
end

if ((turnonplot == 1) && (createavi == 1))
    fig=figure;
    set(fig,'DoubleBuffer','on');
    fnm = strcat(trial,'_movie');    %AVI filename
    aviobj = avifile(fnm,'compression','Cinepak','quality',100);
    aviobj.fps = 10;    %frame rate (frames/sec)
end

% subimage info
if (xrf_mode == 1)
    test = 0;
    i = 1;
    while (test == 0)
        if (pc == 1)
            fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',i)) '.bin'];
        elseif (pc == 0)
            fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03d',i)) '.bin'];
        end
        if (exist(fn_full2,'file') == 2)
            fid = fopen(fn_full2,'r');
            m0 = double(fread(fid,1,'*single'));
            n0 = double(fread(fid,1,'*single'));
            fclose(fid);
            test = 1;
        end
        i = i + 1;
    end 
else
    test = 0;
    i = 1;
    while (test == 0)
        if (pc == 1)
            fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
        elseif (pc == 0)
            fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03d',i)) '.tif'];
        end
        if (exist(fn_full2,'file') == 2)
            info = imfinfo(fn_full2);
            m0 = info.Height;
            n0 = info.Width;
            bd1 = info.BitsPerSample;
            bd1 = bd1(1);
            test = 1;
        end
        i = i + 1;
    end
end

% segment full image into overlapping blocks
if (xrf_mode == 1)
    listOftiffs = dir(fullfile(fn_path2,'*.bin'));
else
    listOftiffs = dir(fullfile(fn_path2,'*.tif'));
end
p2 = numel(listOftiffs);

% apply the same transform to image 2 as was applied to image 1
Xstar00 = csvread(Xstar);
init_reg00 = csvread(init_reg);

fn_out = [trial '1.tif'];
img2_reg = zeros(m1+2*(m0-1),n1+2*(n0-1));
%overlap_msk = zeros(m1+2*(m0-1),n1+2*(n0-1));
blknum = uint16(zeros(m1+2*(m0-1),n1+2*(n0-1)));
blknum(m0:(m0+m1-1),n0:(n0+n1-1)) = blknum0;
if (exist(fn_out,'file') == 2)
    tmp = double(imread(fn_out));
    img2_reg(m0:(m0+m1-1),n0:(n0+n1-1)) = tmp;
    clear tmp
    %overlap_msk = logical(img2_reg>0);
end
[m10 n10] = size(img2_reg);

cnt = 0;
i = 1;
while (cnt < p2)
    if (xrf_mode == 1)
        fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',i)) '.bin'];
        if (exist(fn_full2,'file') == 2)
            cnt = cnt + 1;
            fid = fopen(fn_full2,'r');
            m0 = double(fread(fid,1,'*single'));
            n0 = double(fread(fid,1,'*single'));
            img21 = double(fread(fid,[m0,n0],'*single'));
            fclose(fid);
        end
    else
        fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
        if (exist(fn_full2,'file') == 2)
            cnt = cnt + 1;
            img21 = double(imread(fn_full2,'tif'));
        end
    end
    if (exist(fn_full2,'file') == 2)
        [img20,m21,n21] = apply_mapping(img21,init_reg00(i,:),Xstar00(i,:),Ntran,m10,n10);
        clear img21
        
        img2_reg = mosaic1(img2_reg,blknum,i,img20,m21,n21);
        clear img20
        delete(fn_full2);

        if (turnonplot == 1)
            tmp2 = double(img2_reg);
            tmp2 = tmp2/max(tmp2(:));
            tmp2 = 255*tmp2;
            imshow(uint8(tmp2))
            clear tmp2
            pause(0.1);

            if (createavi == 1)
                F = getframe(fig);
                aviobj = addframe(aviobj,F);
            end
        end
    end
    clear img21
    i = i + 1;
end
clear init_reg00 Xstar00

if (createavi == 1)
    close(fig)
    aviobj = close(aviobj);
end

% write file to disk
if (pc == 1)
    fn_out = [fn_path '\' trial '1.tif'];
elseif (pc == 0)
    fn_out = [fn_path '/' trial '1.tif'];
end 
if (xrf_mode == 1)
    timg = single(img2_reg(m0:(m0+m1-1),n0:(n0+n1-1)));
else
    if (bd1 == 16)
        imwrite(uint16(img2_reg(m0:(m0+m1-1),n0:(n0+n1-1))),fn_out,'tif','Compression','None');
        timg = [];
    elseif (bd1 == 8)
        imwrite(uint8(img2_reg(m0:(m0+m1-1),n0:(n0+n1-1))),fn_out,'tif','Compression','None');
        timg = [];
    end
end
clear img2_reg
end

function [z,m21,n21] = apply_mapping(I2,init_reg,X,Ntran,mfull,nfull)
%% apply_mapping
% Transform second image so that it is registered with the first
%   Input: I2 -> image
%          X -> transform coefficients
%   Outputs: z -> registered image

[m0 n0] = size(I2);
[x0,y0] = meshgrid(0:(n0-1),0:(m0-1));
x0 = x0 + init_reg(2);
y0 = y0 + init_reg(1);
if (Ntran == 0)
    v2 = X(1)*x0 + X(2)*y0 + X(3)*x0.*y0 + X(4);
    u2 = X(5)*x0 + X(6)*y0 + X(7)*x0.*y0 + X(8);
else
    v2 = zeros(m0,n0);
    u2 = zeros(m0,n0);
    pnt = 1;
    for i = 0:Ntran
        for j = 0:(Ntran-i)
            v2 = v2 + X(pnt)*((x0.^i).*(y0.^j));
            u2 = u2 + X(pnt+sum((1:(Ntran+1))))*((x0.^i).*(y0.^j));
            pnt = pnt + 1;
        end
    end
end
clear x0 y0

m21 = ceil(min(u2(:)));
if (m21 < 1)
    m21 = 1;
end
m22 = floor(max(u2(:)));
if (m22 > mfull)
    m22 = mfull;
end
n21 = ceil(min(v2(:)));
if (n21 < 1)
    n21 = 1;
end
n22 = floor(max(v2(:)));
if (n22 > nfull)
    n22 = nfull;
end
F = TriScatteredInterp(double(v2(:)),double(u2(:)),double(I2(:)),'natural');
clear u2 v2 I2
[v1,u1] = meshgrid(n21:n22,m21:m22);
z = F(v1,u1);
z(isnan(z)) = 0;
% removed 20140306
%{
msk21 = imfill((z>0),'holes');
se = strel('square',3);
msk21 = imerode(msk21,se);
z = z.*msk21;
clear msk21
%}
clear F u1 v1
end

function [img2_reg] = mosaic1(img2_reg,blknum,blk,img20,m21,n21)
%% mosaicing

[m0,n0] = size(img20);
m22 = m21 + m0 - 1;
n22 = n21 + n0 - 1;
blknum0 = blknum(m21:m22,n21:n22);
msk = (blknum0==blk);
img2_reg(m21:m22,n21:n22) = img2_reg(m21:m22,n21:n22).*single(~msk) + img20.*single(msk);
clear msk img20
end