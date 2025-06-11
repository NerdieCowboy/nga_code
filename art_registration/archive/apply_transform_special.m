function [] = apply_transform_special(fn,m1,n1,Ntran,init_reg,Xstar,createavi)
%% Transform repeat function
% apply_transform.m
% Given that a transform has been applied to an image (image 1) in order to 
%   register it with another image (image 0). This function applies that same 
%   transform to a second image (image 2). If image 1 and image 2 were 
%   originally registered, then function will register image2 with image 0.
% Author: Damon Conover
% Email: dconover@gwmail.gwu.edu
% Latest Revision: 6 February 2013
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

feather_width = 64;
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
    %fn_path2 = [fn_path '\' trial '\'];
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
test = 0;
i = 1;
while (test == 0)
    if (pc == 1)
        fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',i)) '.tif']
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

% segment full image into overlapping blocks
listOftiffs = dir(fullfile(fn_path2,'*.tif'));
p2 = numel(listOftiffs);

% apply the same transform to image 2 as was applied to image 1
Xstar00 = csvread(Xstar);
init_reg00 = csvread(init_reg);

fn_out = [trial '1.tif'];
img2_reg = zeros(m1+2*(m0-1),n1+2*(n0-1));
overlap_msk = zeros(m1+2*(m0-1),n1+2*(n0-1));
if (exist(fn_out,'file') == 2)
    tmp = double(imread(fn_out));
    img2_reg(m0:(m0+m1-1),n0:(n0+n1-1)) = tmp;
    clear tmp
    overlap_msk = logical(img2_reg>0);
end
[m10 n10] = size(img2_reg);

i = 130;
while (i < p2)
    fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',i)) '.tif']
    if (exist(fn_full2,'file') == 2)
        img21 = double(imread(fn_full2,'tif'));
        [img20,m21,n21] = apply_mapping(img21,init_reg00(i,:),Xstar00(i,:),Ntran,m10,n10);
        clear img21
        
        [img2_reg,overlap_msk] = mosaic(img2_reg,overlap_msk,img20,m21,n21,init_reg00(i,:),feather_width);
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
    i = i + 1
end
i = 129;
while (i >= 1)
    fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
    if (exist(fn_full2,'file') == 2)
        img21 = double(imread(fn_full2,'tif'));
        [img20,m21,n21] = apply_mapping(img21,init_reg00(i,:),Xstar00(i,:),Ntran,m10,n10);
        clear img21
        
        [img2_reg,overlap_msk] = mosaic(img2_reg,overlap_msk,img20,m21,n21,init_reg00(i,:),feather_width);
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
    i = i - 1
end
clear overlap_msk init_reg00 Xstar00

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
if (bd1 == 16)
    imwrite(uint16(img2_reg(m0:(m0+m1-1),n0:(n0+n1-1))),fn_out,'tif','Compression','None');
elseif (bd1 == 8)
    imwrite(uint8(img2_reg(m0:(m0+m1-1),n0:(n0+n1-1))),fn_out,'tif','Compression','None');
end
clear img2_reg
end

function [img2_reg,overlap_msk] = mosaic(img2_reg,overlap_msk,img20,m21,n21,init_reg,feather_width)
%% mosaicing

msk21 = logical(imfill((img20>0),'holes'));
%offm = round(m1*0.1);
%offn = round(n1*0.1);
%m21 = 1 + init_reg(1) - offm - 1;
%m22 = m1 + 2*offm + init_reg(1) - offm - 1;
%n21 = 1 + init_reg(2) - offn - 1;
%n22 = n1 + 2*offn + init_reg(2) - offn - 1;
[m0,n0] = size(img20);
m22 = m21 + m0 - 1;
n22 = n21 + n0 - 1;
ov_msk_sub = overlap_msk(m21:m22,n21:n22);
intersect_msk = and(ov_msk_sub,msk21);
img2_reg(m21:m22,n21:n22) = img2_reg(m21:m22,n21:n22) + img20.*(msk21-intersect_msk);

se = strel('square',3);
C1 = imdilate(intersect_msk,se)-intersect_msk;
A1 = and(C1,msk21);
[A1y A1x] = find(A1==1);
clear A1
B1 = and(C1,ov_msk_sub);
[B1y B1x] = find(B1==1);
clear C1 B1

se = strel('square',feather_width);
tmp = intersect_msk - imerode(intersect_msk,se);
[Cy Cx] = find(tmp==1);
clear tmp intersect_msk

if (~isempty([A1y A1x]) && ~isempty([B1y B1x]) && ~isempty([Cy Cx]))
    [~,dA] = dsearchn([A1y A1x],[Cy Cx]);
    [~,dB] = dsearchn([B1y B1x],[Cy Cx]);

    if (~isempty(dA))
        for i = 1:length(Cy)
            if (dB(i) >= dA(i))
                Cy_sub = Cy(i) + m21 - 1;
                Cx_sub = Cx(i) + n21 - 1;
                %img2_reg(Cy_sub,Cx_sub) = (img20(Cy(i),Cx(i))*dB(i) + img2_reg(Cy_sub,Cx_sub)*dA(i))/(dA(i) + dB(i));
                img2_reg(Cy_sub,Cx_sub) = (img20(Cy(i),Cx(i))*abs(feather_width-dA(i)) + img2_reg(Cy_sub,Cx_sub)*dA(i))/(abs(feather_width-dA(i))+dA(i));
            end
        end
    end
    clear kA kB dA dB
end
clear Cx Cy A1x A1y B1x B1y

overlap_msk(m21:m22,n21:n22) = or(ov_msk_sub,msk21);
clear msk21 img20
end


function [z,m21,n21] = apply_mapping(I2,init_reg,X,Ntran,mfull,nfull)
%% apply_mapping
% Transform second image so that it is registered with the first
%   Input: I2 -> image
%          X -> transform coefficients
%   Outputs: z -> registered image

[m0 n0] = size(I2);
%{
offm = round(m0*0.1);
offn = round(n0*0.1);
z = single(zeros((m0+2*offm),(n0+2*offn)));
mend = init_reg(1) + m0 - 1;
nend = init_reg(2) + n0 - 1;

stopm = 0;
m11 = init_reg(1);
m12 = m11 + gridsz - 1;
if (m12 >= mend)
    m12 = mend;
    stopm = stopm + 1;
end
while (stopm < 2)
    stopn = 0;
    n11 = init_reg(2);
    n12 = n11 + gridsz - 1;
    if (n12 >= nend)
        n12 = nend;
        stopn = stopn + 1;
    end
    while (stopn < 2)
        [x0,y0] = meshgrid(n11:n12,m11:m12);
        if (Ntran == 0)
            u2 = X(1)*x0 + X(2)*y0 + X(3)*x0.*y0 + X(4);
            v2 = X(5)*x0 + X(6)*y0 + X(7)*x0.*y0 + X(8);
        else
            u2 = zeros(size(y0));
            v2 = zeros(size(x0));
            pnt = 1;
            for i = 0:Ntran
                for j = 0:(Ntran-i)
                    u2 = u2 + X(pnt)*((x0.^i).*(y0.^j));
                    v2 = v2 + X(pnt+sum((1:(Ntran+1))))*((x0.^i).*(y0.^j));
                    pnt = pnt + 1;
                end
            end
        end
        clear x0 y0
        msk = (u2>=1).*(u2<=nfull).*(v2>=1).*(v2<=mfull);
        u2(msk==0) = [];
        v2(msk==0) = [];
        tmp = I2((m11-init_reg(1)+1):(m12-init_reg(1)+1),(n11-init_reg(2)+1):(n12-init_reg(2)+1));
        tmp(msk==0) = [];
        clear msk
        F = TriScatteredInterp(double(u2(:)),double(v2(:)),double(tmp(:)),'natural');
        clear tmp
        n21 = ceil(min(u2(:)));
        n22 = floor(max(u2(:)));
        m21 = ceil(min(v2(:)));
        m22 = floor(max(v2(:)));
        clear u2 v2
        [u1,v1] = meshgrid(n21:n22,m21:m22);
        z0 = F(u1,v1);
        z0(isnan(z0)) = 0;
        clear F u1 v1
        
        msk = (z0~=0);      
        m21 = m21 - init_reg(1) + offm + 1;
        m22 = m22 - init_reg(1) + offm + 1;
        n21 = n21 - init_reg(2) + offn + 1;
        n22 = n22 - init_reg(2) + offn + 1;
        z(m21:m22,n21:n22) = z(m21:m22,n21:n22).*(1-msk) + z0.*msk;
        clear z0 msk

        n11 = n12 - round(gridsz/8) + 1;
        n12 = n11 + gridsz - 1;
        if (n12 >= nend)
            n12 = nend;
            stopn = stopn + 1;
        end
    end
    m11 = m12 - round(gridsz/8) + 1;
    m12 = m11 + gridsz - 1;
    if (m12 >= mend)
        m12 = mend;
        stopm = stopm + 1;
    end
end
%}
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
msk21 = imfill((z>0),'holes');
se = strel('square',3);
msk21 = imerode(msk21,se);
z = z.*msk21;
clear F u1 v1 msk21
end
