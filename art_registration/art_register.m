function [] = art_register(rgb_fn,sc,ir_fn,j01,j02,repcnt_max,max_shift,minpts,Npoly,Ntran,memory_limited,morepts,useresults,filteron,addsharp,debugmode,usegui,feather_width,mxd)
%% Registration/Mosaic algorithm
% art_register.m
% Author: Damon Conover
% Email: dconover@gwu.edu
% Latest Revision: 15 July 2015
%

% settings
gpu = 0;    %no gpu=0, gpu=1
turnonoutput = 1;   %print maximum disparity and number of features points to screen -> 1, otherwise -> 0
turnonplot = 1;     %display mosiac after each block is added -> 1, otherwise -> 0
if (usegui ~= 0)
    turnonoutput = 0;
    turnonplot = 0;
end

% constants
samp_ratio = 3; %upsample ratio for determining fractional pixels (must be an odd integer)
down = 1/4;

if (useresults == 1)
    morepts = 1;
end

warning('off', 'all')

if (j02 < j01)
    error('j02 must be larger than or equal to j01')
end

if (j01 ~= round(j01) || j02 ~= round(j02))
    error('j01 and j02 must be integers')
end

if (minpts < 8)
    error('minpts must be larger than 8')
end

if ((samp_ratio-1) ~= 2*round((samp_ratio-1)/2))
    error('samp_ratio must be an odd integer')
end

%% *********************************************
% IR
% read in relative block-to-block shifts in pixels (rel_reg)
tmp = regexp(ir_fn,'\\');
pc = 1;
if (isempty(tmp))
    pc = 0;
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
if (pc == 1)
    fn_path2 = [fn_path trial '\'];
elseif (pc == 0)
    fn_path2 = [fn_path trial '/'];
end
fn_full2 = [fn_path2 'offset_values.csv'];
rel_reg = csvread(fn_full2);

%% *********************************************
% get size, bit depth, and image lists for blocks 
%   (m0,n0,m00,n00,bd2,tifnumlist)
listOftiffs = dir(fullfile(fn_path2,'*.tif'));
p2 = numel(listOftiffs);
clear listOftiffs
inc = 0;
p0 = 0;
tifnumlist = [];
while (inc < p2)
    p0 = p0 + 1;
    if (pc == 1)
        fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',p0)) '.tif'];
    elseif (pc == 0)
        fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03d',p0)) '.tif'];
    end
    if (exist(fn_full2,'file') == 2)
        info = imfinfo(fn_full2);
        m0 = info.Height;
        n0 = info.Width;
        bd2 = info.BitsPerSample;
        bd2 = bd2(1);
        inc = inc + 1;
        tifnumlist = [tifnumlist; p0];
    end
end
img2 = single(zeros(m0,n0));
img2 = imresize(img2,down);
[m00,n00] = size(img2);
clear img2

%% *********************************************
% order IR images by their entropy (IX)
H = [];
img2 = single(zeros(m0,n0,p0));
i = 0;
for inc = 1:p0
    if (pc == 1)
        fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',inc)) '.tif'];
    elseif (pc == 0)
        fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03d',inc)) '.tif'];
    end
    if (exist(fn_full2,'file') == 2)
        i = i + 1;
        img21 = single(imread(fn_full2,'tif'));
        img2(:,:,i) = img21;
        if (bd2 == 16)
            H = [H; entropy(uint16(img21))];
        elseif (bd2 == 8)
            H = [H; entropy(uint8(img21))];
        end
        clear img21
        pause(0.1)
    end
end
[~,IX] = sort(H,'descend');
clear H1 H
clear img2

%% *********************************************
% RGB -> R
% read/create reference image (image0tif), scaled reference image
%   (image1tif), and phase image (image0bin)
info = imfinfo(rgb_fn);
bd1 = info.BitsPerSample;
bd1 = bd1(1);
rgb1 = imread(rgb_fn);
img1 = imresize(rgb1,sc,'bicubic');
clear rgb1
fn_out = [fn_path trial '_RGB.tif'];
if (exist(fn_out,'file') ~= 2)
    imwrite(img1,fn_out,'tif','Compression','None');
end
[m1 n1 p1] = size(img1);

if (bd1 == 8)
    image0tif = uint8(ones(m1+2*(m0-1),n1+2*(n0-1),p1));
    tmp = rand(m1+2*(m0-1),n1+2*(n0-1),p1);
elseif (bd1 == 16)
    image0tif = uint16(ones(m1+2*(m0-1),n1+2*(n0-1),p1));
    tmp = rand(m1+2*(m0-1),n1+2*(n0-1),p1);
end
for i = 1:p1
    if (bd1 == 8) 
        image0tif(:,:,i) = uint8((mean2(double(img1(:,:,i))))*tmp(:,:,i));
    elseif (bd1 == 16)
        image0tif(:,:,i) = uint16((mean2(double(img1(:,:,i))))*tmp(:,:,i));
    end
end
clear tmp
image0tif(m0:(m0+m1-1),n0:(n0+n1-1),:) = img1;
clear img1
[m1x n1x ~] = size(image0tif);

image1tif = imresize(image0tif,down,'bicubic');
[m10 n10 ~] = size(image1tif);
if (memory_limited == 1)
    if (exist('image1.tif','file') ~= 2)
        imwrite(image1tif,[fn_path 'image1.tif'],'tif','Compression','None')
    end
    clear image1tif
    image1tif = zeros(1,1,p1);
end

%msk = true(m1x,n1x);

if (exist('image0.bin','file') ~= 2)
    for i = 1:p1
        i0fn = [fn_path 'image0' num2str(i) '.tif'];
        imwrite(image0tif(:,:,i),i0fn,'tif','Compression','None')
    end
    clear image0tif
    image0bin = single(zeros(m1x,n1x,p1));

    % phase image
    for i = 1:p1
        i0fn = [fn_path 'image0' num2str(i) '.tif'];
        tmp = single(imread(i0fn,'tif'));
        delete(i0fn)
        tmp1 = fft2(tmp);
        clear tmp
        tmp = atan2(imag(tmp1),real(tmp1));
        clear tmp1
        tmp1 = exp(1i*tmp);
        clear tmp
        image0bin(:,:,i) = real(ifft2(tmp1));
        clear tmp1
    end
    if (memory_limited == 1)
        for i = 1:p1
            multibandwrite(single(image0bin(:,:,i)),'image0.bin','bsq',[1,1,i],[m1x,n1x,p1],'precision','single','machfmt','ieee-le')
        end
        clear image0bin
        image0bin = zeros(m1x,n1x,p1);
    end
else
    clear image0bin
    image0bin = zeros(m1x,n1x,p1);
end

%% *********************************************
% if script was stopped midway, read in previous state when it is run later
restart_test = 1;
fn_out = [fn_path trial '_IR.tif'];
fn_zd = [fn_path trial '_dist.tif'];
fn_blknum = [fn_path trial '_blk.tif'];
if (exist(fn_out,'file') == 2)
    tmp = single(imread(fn_out,'tif'));
    img2_reg = single(zeros(m1+2*(m0-1),n1+2*(n0-1)));
    img2_reg(m0:(m0+m1-1),n0:(n0+n1-1)) = tmp;
    clear tmp
    tmp = uint16(imread(fn_zd,'tif'));
    zd = uint16(zeros(m1+2*(m0-1),n1+2*(n0-1)));
    zd(m0:(m0+m1-1),n0:(n0+n1-1)) = tmp;
    clear tmp
    tmp = uint16(imread(fn_blknum,'tif'));
    blknum = uint16(zeros(m1+2*(m0-1),n1+2*(n0-1)));
    blknum(m0:(m0+m1-1),n0:(n0+n1-1)) = tmp;
    clear tmp
else
    if (bd1 == 8)
        tmp = uint8(zeros(m1,n1));
    elseif (bd1 == 16)
        tmp = uint16(zeros(m1,n1));
    end
    imwrite(tmp,fn_out,'tif','Compression','None')
    tmp = uint16(zeros(m1,n1));
    imwrite(tmp,fn_zd,'tif','Compression','None')
    imwrite(tmp,fn_blknum,'tif','Compression','None')
    clear tmp
    img2_reg = single(zeros(m1x,n1x));
    zd = uint16((2^16-1)*ones(m1x,n1x));
    blknum = uint16(zeros(m1x,n1x));
    restart_test = 0;
end

fn_full2 = [fn_path2 'registration_results.csv'];
if (exist(fn_full2,'file') == 2)
    reg_results = csvread(fn_full2);
else
    reg_results = zeros(p0,5);
end

fn_full2 = [fn_path2 'Xstar.csv'];
if (exist(fn_full2,'file') == 2)
    Xstar00 = csvread(fn_full2);
else
    if (Ntran == 0)
        Xstar00 = zeros(p0,8);
    else
        Xstar00 = zeros(p0,sum(2*(1:(Ntran+1))));
    end
    restart_test = 0;
end

fn_full2 = [fn_path2 'init_reg.csv'];
if (exist(fn_full2,'file') == 2)
    init_reg00 = csvread(fn_full2);
else
    init_reg00 = zeros(p0,2);
    restart_test = 0;
end

%% *********************************************
% find initial translation for one subimage
if (restart_test == 0)
    if (memory_limited == 1)
        image1tif = imread([fn_path 'image1.tif'],'tif');
    end
    pass = 0;
    i = 1;
    p = 1;
    while (pass == 0)
        if (i <= p0)
            % create cube of scaled-down overlapping blocks from the reference image
            img1 = single(zeros(2*m00,2*n00,round(m10/(m10-(m00-1)))*round(n10/(n10-(n00-1)))));
            img1_corner = zeros(round(m10/(m10-(m00-1)))*round(n10/(n10-(n00-1))),2);
            inc = 0;
            stopm = 0;
            m11 = 1;
            m12 = m11 + 2*m00 - 1;
            if (m12 > m10)
                m12 = m10;
                stopm = stopm + 1;
            end
            while (stopm < 2)
                stopn = 0;
                n11 = 1;
                n12 = n11 + 2*n00 - 1;
                if (n12 > n10)
                    n12 = n10;
                    stopn = stopn + 1;
                end
                while (stopn < 2)
                    inc = inc + 1;
                    tmp1 = image1tif(m11:m12,n11:n12,p);
                    img1(1:(m12-m11+1),1:(n12-n11+1),inc) = tmp1;
                    clear tmp1
                    img1_corner(inc,:) = [m11 n11];

                    n11 = n12 - (n00-1) + 1;
                    n12 = n11 + 2*n00 - 1;
                    if (n12 >= n10)
                        n12 = n10;
                        stopn = stopn + 1;
                    end
                end
                m11 = m12 - (m00-1) + 1;
                m12 = m11 + 2*m00 - 1;
                if (m12 >= m10)
                    m12 = m10;
                    stopm = stopm + 1;
                end                
            end
            img1(:,:,(inc+1):end) = [];
            img1_corner((inc+1):end,:) = [];

            %i = i + 1;
            if (pc == 1)
                fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',tifnumlist(IX(i)))) '.tif'];
            elseif (pc == 0)
                fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03d',tifnumlist(IX(i)))) '.tif'];
            end
            if (exist(fn_full2,'file') == 2)
                % register downsampled image
                if (pc == 1)
                    fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',tifnumlist(IX(i)))) '.tif'];
                elseif (pc == 0)
                    fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03d',tifnumlist(IX(i)))) '.tif'];
                end
                img21 = single(imread(fn_full2,'tif'));
                img2 = imresize(img21,down);

                if (addsharp == 1)
                    %LOW_HIGH = stretchlim(uint16(img2),[0.001 0.999]);
                    %img2 = single(imadjust(uint16(img2),LOW_HIGH,[]));
                    %h = fspecial('unsharp',0.2);
                    %img2 = imfilter(img2,h);
                    %img2(img2<0) = 0;
                    %img2 = img2/max(img2(:))*(2^16-1);
                end

                clear img21
                [cctable_out fail] = init_pass(img1,img2,img1_corner,addsharp,gpu);
                clear img2
                clear init_pass
                corner1 = [cctable_out(1,5) cctable_out(1,6)];
                init_reg = corner1/down - [1/(2*down)+1 1/(2*down)+1];
                init_reg = round(init_reg);
                if (init_reg(1)<1)
                    init_reg(1) = 1;
                end
                if (init_reg(2)<1)
                    init_reg(2) = 1;
                end

                if (fail == 0)
                    %registration
                    img21 = single(imread(fn_full2,'tif'));
                    %repcnt = 0;
                    [Xstar,pass,img2_reg,zd,blknum,mxd1,mf1] = incremental_registration(image0bin(:,:,p),img2_reg,zd,blknum,tifnumlist(IX(i)),img21,max_shift,samp_ratio,init_reg,mxd,j01,j02,turnonoutput,turnonplot,0,filteron,m1x,n1x,p1,minpts,feather_width,gpu,morepts,p,memory_limited,Npoly,Ntran,addsharp,debugmode,usegui);
                    %[Xstar,pass,img2_reg,zd,blknum,mxd1,mf1] = incremental_registration(image0bin(:,:,p),img2_reg,zd,blknum,tifnumlist(IX(i)),img21,max_shift,samp_ratio,init_reg,mxd,j01,j02,turnonoutput,turnonplot,0,filteron,m1x,n1x,p1,minpts,feather_width,gpu,0,p,memory_limited,Npoly,Ntran,addsharp,debugmode,usegui);
                    clear incremental_registration
                    clear img21
                end
            end
            if (p >= p1)
                p = 1;
                if (pass == 0)
                    i = i + 1;
                end
            else
                p = p + 1;
            end
        else
            return
        end
    end
    fn_out = [fn_path trial '_IR.tif'];
    fn_zd = [fn_path trial '_dist.tif'];
    fn_blknum = [fn_path trial '_blk.tif'];
    if (bd2 == 16)
        imwrite(uint16(img2_reg(m0:(m0+m1-1),n0:(n0+n1-1))),fn_out,'tif','Compression','None');
    elseif (bd2 == 8)
        imwrite(uint8(img2_reg(m0:(m0+m1-1),n0:(n0+n1-1))),fn_out,'tif','Compression','None');
    end
    imwrite(uint16(zd(m0:(m0+m1-1),n0:(n0+n1-1))),fn_zd,'tif','Compression','None');
    imwrite(uint16(blknum(m0:(m0+m1-1),n0:(n0+n1-1))),fn_blknum,'tif','Compression','None');
    clear img1 img1_corner
    inc = tifnumlist(IX(i));
    delete(fn_full2);
    
    Xstar00(inc,:) = Xstar;
    init_reg00(inc,:) = init_reg;
    fn_full2 = [fn_path2 'Xstar.csv'];
    csvwrite(fn_full2,Xstar00)
    fn_full2 = [fn_path2 'init_reg.csv'];
    csvwrite(fn_full2,init_reg00)
    reg_results(inc,:) = [inc mxd1 mf1 p-1 morepts];
    fn_full2 = [fn_path2 'registration_results.csv'];
    csvwrite(fn_full2,reg_results)
    pass_test = [];
    if (memory_limited == 1)
        clear image1tif
        image1tif = zeros(1,1,p1);
    end
else
    tmp = sum((Xstar00==zeros(size(Xstar00))),2);
    if (Ntran == 0)
        pass_test = find(tmp ~= 8);
    else
        pass_test = find(tmp ~= sum(2*(1:(Ntran+1))));
    end
    clear tmp
end

%% *********************************************
% increment backward through blocks
direction = 0;
mult = 1;
if (restart_test == 1)
    inc = min(pass_test);
end
previ = inc;

for i = (inc-1):-1:1
    if (length(find(pass_test == i)) == 1)
        previ = i;
        mult = 1;
    else
        [init_reg00,Xstar00,mult,previ,img2_reg,zd,blknum,mxd1,mf1,boost,color] = reg_loop(image0bin,image1tif,img2_reg,zd,blknum,fn_path,fn_path2,trial,rel_reg,repcnt_max,init_reg00,Xstar00,mult,max_shift,samp_ratio,down,mxd,j01,j02,i,previ,direction,turnonoutput,turnonplot,useresults,pc,filteron,m1x,n1x,minpts,bd2,feather_width,gpu,morepts,memory_limited,Npoly,Ntran,addsharp,debugmode,usegui);
        clear reg_loop
        reg_results(i,:) = [i mxd1 mf1 color boost];
        fn_full2 = [fn_path2 'registration_results.csv'];
        csvwrite(fn_full2,reg_results)
    end  
    fn_full2 = [fn_path2 'Xstar.csv'];
    csvwrite(fn_full2,Xstar00)
    fn_full2 = [fn_path2 'init_reg.csv'];
    csvwrite(fn_full2,init_reg00)
end

% increment forward through blocks
direction = 1;
mult = 1;
if (restart_test == 1)
    inc = max(pass_test);
end

previ = inc;
for i = (inc+1):p0
    [init_reg00,Xstar00,mult,previ,img2_reg,zd,blknum,mxd1,mf1,boost,color] = reg_loop(image0bin,image1tif,img2_reg,zd,blknum,fn_path,fn_path2,trial,rel_reg,repcnt_max,init_reg00,Xstar00,mult,max_shift,samp_ratio,down,mxd,j01,j02,i,previ,direction,turnonoutput,turnonplot,useresults,pc,filteron,m1x,n1x,minpts,bd2,feather_width,gpu,morepts,memory_limited,Npoly,Ntran,addsharp,debugmode,usegui);
    clear reg_loop
    fn_full2 = [fn_path2 'Xstar.csv'];
    csvwrite(fn_full2,Xstar00)
    fn_full2 = [fn_path2 'init_reg.csv'];
    csvwrite(fn_full2,init_reg00)
    reg_results(i,:) = [i mxd1 mf1 color boost];
    fn_full2 = [fn_path2 'registration_results.csv'];
    csvwrite(fn_full2,reg_results)
end

% go to lowest-numbered registered block and increment forward
tmp = sum((Xstar00==zeros(size(Xstar00))),2);
if (Ntran == 0)
    pass_test = find(tmp ~= 8);
else
    pass_test = find(tmp ~= sum(2*(1:(Ntran+1))));
end
clear tmp
direction = 1;
previ = min(pass_test);
for i = (min(pass_test)+1):(inc-1)
    if (length(find(pass_test == i)) == 1)
        previ = i;
        mult = 1;
    else
        [init_reg00,Xstar00,mult,previ,img2_reg,zd,blknum,mxd1,mf1,boost,color] = reg_loop(image0bin,image1tif,img2_reg,zd,blknum,fn_path,fn_path2,trial,rel_reg,repcnt_max,init_reg00,Xstar00,mult,max_shift,samp_ratio,down,mxd,j01,j02,i,previ,direction,turnonoutput,turnonplot,0,pc,filteron,m1x,n1x,minpts,bd2,feather_width,gpu,morepts,memory_limited,Npoly,Ntran,addsharp,debugmode,usegui);
        clear reg_loop
        reg_results(i,:) = [i mxd1 mf1 color boost];
        fn_full2 = [fn_path2 'registration_results.csv'];
        csvwrite(fn_full2,reg_results)
    end
    fn_full2 = [fn_path2 'Xstar.csv'];
    csvwrite(fn_full2,Xstar00)
    fn_full2 = [fn_path2 'init_reg.csv'];
    csvwrite(fn_full2,init_reg00)
end

% go to higest-numbered registered block and increment backward
direction = 0;
previ = max(pass_test);
for i = (max(pass_test)-1):-1:(inc+1)
    if (length(find(pass_test == i)) == 1)
        previ = i;
        mult = 1;
    else
        [init_reg00,Xstar00,mult,previ,img2_reg,zd,blknum,mxd1,mf1,boost,color] = reg_loop(image0bin,image1tif,img2_reg,zd,blknum,fn_path,fn_path2,trial,rel_reg,repcnt_max,init_reg00,Xstar00,mult,max_shift,samp_ratio,down,mxd,j01,j02,i,previ,direction,turnonoutput,turnonplot,0,pc,filteron,m1x,n1x,minpts,bd2,feather_width,gpu,morepts,memory_limited,Npoly,Ntran,addsharp,debugmode,usegui);
        clear reg_loop
        reg_results(i,:) = [i mxd1 mf1 color boost];
        fn_full2 = [fn_path2 'registration_results.csv'];
        csvwrite(fn_full2,reg_results)
    end
    fn_full2 = [fn_path2 'Xstar.csv'];
    csvwrite(fn_full2,Xstar00)
    fn_full2 = [fn_path2 'init_reg.csv'];
    csvwrite(fn_full2,init_reg00)
end
clear rel_reg
clear image0bin image1tif img2_reg zd blknum
if (memory_limited == 1)
    delete('image1.tif')
    delete('image0.bin')
end
end

function [init_reg00,Xstar00,mult,previ,img2_reg,zd,blknum,mxd1,mf1,boost,color] = reg_loop(image0bin,image1tif,img2_reg,zd,blknum,fn_path,fn_path2,trial,rel_reg,repcnt_max,init_reg00,Xstar00,mult,max_shift0,samp_ratio,down,mxd,j01,j02,i,previ,direction,turnonoutput,turnonplot,useresults,pc,filteron,m2,n2,minpts,bd2,feather_width,gpu,morepts,memory_limited,Npoly,Ntran,addsharp,debugmode,usegui)
%% registration

if (usegui ~= 0)
    if (usegui == 1)
        handles = guidata(ArtRegister);
    elseif (usegui == 2)
        handles = guidata(HSIRegister);
    elseif (usegui == 3)
        handles = guidata(XrayRegister);
    end
    stopcomm = get(handles.stop,'Value');
    if (stopcomm == 1)
        set(handles.stop,'Value',0);
        error('Stop command issued')
    end
end

mxd1 = 0;
mf1 = 0;
color = 0;
boost = 0;
if (memory_limited == 1)
    image1tif = imread([fn_path 'image1.tif'],'tif');
end
[m10 n10 p1] = size(image1tif);

fn_out = [fn_path trial '_RGB.tif'];
info = imfinfo(fn_out);
m1 = info.Height;
n1 = info.Width;
bd1 = info.BitsPerSample;
bd1 = bd1(1);
        
if (pc == 1)
    fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
elseif (pc == 0)
    fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03d',i)) '.tif'];
end
if (exist(fn_full2,'file') == 2)
    reli = find(rel_reg(:,1) == i);
    Xstar = Xstar00(previ,:);
    init_reg = init_reg00(previ,:);
    if (direction == 1)
        %forward
        rel_reg0 = rel_reg(reli,2:3);
    elseif (direction == 0)
        %backward
        rel_reg0 = -rel_reg(reli+1,2:3);
    end

    %registration
    img21 = single(imread(fn_full2,'tif'));
    [m0 n0] = size(img21);
    img2 = imresize(img21,down);
    
    if (addsharp == 1)
        %LOW_HIGH = stretchlim(uint16(img2),[0.001 0.999]);
        %img2 = single(imadjust(uint16(img2),LOW_HIGH,[]));
        %h = fspecial('unsharp',0.2);
        %img2 = imfilter(img2,h);
        %img2(img2<0) = 0;
        %img2 = img2/max(img2(:))*(2^16-1);
    end
    
    if (Ntran == 0)
        search_test = 0;
    else
        search_test = 1;
    end
    
    repcnt = 0;
    pass = 0;
    mult0 = mult;
    if (morepts == 0)
        boost = 0;
    elseif (morepts == 1)
        boost = 1;
    end
    color = 1;
    while ((pass == 0) && (repcnt <= repcnt_max))
        max_shift = max_shift0*(1+repcnt);
        
        x0 = init_reg(2) + rel_reg0(2);
        y0 = init_reg(1) + rel_reg0(1);
        x1 = x0 + n0 - 1;
        y1 = y0 + m0 - 1;
        reg_out = [0 0];
        reg_out2 = [0 0];
        if (Ntran == 0)
            reg_out = [round(Xstar(5)*x0 + Xstar(6)*y0 + Xstar(7)*x0*y0 + Xstar(8)) round(Xstar(1)*x0 + Xstar(2)*y0 + Xstar(3)*x0*y0 + Xstar(4))];
            reg_out2 = [round(Xstar(5)*x1 + Xstar(6)*y1 + Xstar(7)*x1*y1 + Xstar(8)) round(Xstar(1)*x1 + Xstar(2)*y1 + Xstar(3)*x1*y1 + Xstar(4))];
        else
            pnt = 1;
            for i0 = 0:Ntran
                for j = 0:(Ntran-i0)
                    reg_out(2) = reg_out(2) + Xstar(pnt)*((x0.^i0).*(y0.^j));
                    reg_out(1) = reg_out(1) + Xstar(pnt+sum((1:(Ntran+1))))*((x0.^i0).*(y0.^j));
                    reg_out2(2) = reg_out2(2) + Xstar(pnt)*((x1.^i0).*(y1.^j));
                    reg_out2(1) = reg_out2(1) + Xstar(pnt+sum((1:(Ntran+1))))*((x1.^i0).*(y1.^j));
                    pnt = pnt + 1;
                end
            end
        end
        clear x1 y1
        reg_out(1) = round(reg_out(1));
        reg_out(2) = round(reg_out(2));
        mtmp1 = floor((reg_out(1)-m0*(mult-1)-max_shift*2)*down);
        if (mtmp1 < 1)
            mtmp1 = 1;
        end
        ntmp1 = floor((reg_out(2)-n0*(mult-1)-max_shift*2)*down);
        if (ntmp1 < 1)
            ntmp1 = 1;
        end
        reg_out = reg_out2;
        reg_out(1) = round(reg_out(1));
        reg_out(2) = round(reg_out(2));
        mtmp2 = floor((reg_out(1)+m0*(mult-1)+max_shift*2)*down);
        if (mtmp2 > m10)
            mtmp2 = m10;
        end
        ntmp2 = floor((reg_out(2)+n0*(mult-1)+max_shift*2)*down);
        if (ntmp2 > n10)
            ntmp2 = n10;
        end
        if ((mtmp1 >= mtmp2) && (mtmp2 > 1))
            mtmp1 = mtmp2 - round(3*m0*down);
            if (mtmp1 < 1)
                mtmp1 = 1;
            end
        elseif (mtmp1 >= mtmp2) && (mtmp2 < 1)
            mtmp2 = mtmp1 + round(3*m0*down);
            if (mtmp2 > m10)
                mtmp2 = m10;
            end
        end
        if ((ntmp1 >= ntmp2) && (ntmp2 > 1))
            ntmp1 = ntmp2 - round(3*n0*down);
            if (ntmp1 < 1)
                ntmp1 = 1;
            end
        elseif (ntmp1 >= ntmp2) && (ntmp2 < 1)
            ntmp2 = ntmp1 + round(3*n0*down);
            if (ntmp2 > n10)
                ntmp2 = n10;
            end
        end
        if ((useresults == 1) || (addsharp == 1))
            overlap_msk = (zd < (2^16-1));
            msk_tmp = imresize(overlap_msk,down,'bicubic');    
            if (bd1 == 8)
                img2_tmp = uint8(imresize(img2_reg,down,'bicubic'));
                %img1 = image1tif(mtmp1:mtmp2,ntmp1:ntmp2,color).*uint8(~msk_tmp(mtmp1:mtmp2,ntmp1:ntmp2)) + img2_tmp(mtmp1:mtmp2,ntmp1:ntmp2).*uint8(msk_tmp(mtmp1:mtmp2,ntmp1:ntmp2));
                img1 = img2_tmp(mtmp1:mtmp2,ntmp1:ntmp2).*uint8(msk_tmp(mtmp1:mtmp2,ntmp1:ntmp2));
                clear img2_tmp
            elseif (bd1 == 16)
                img2_tmp = uint16(imresize(img2_reg,down,'bicubic'));
                %img1 = image1tif(mtmp1:mtmp2,ntmp1:ntmp2,color).*uint16(~msk_tmp(mtmp1:mtmp2,ntmp1:ntmp2)) + img2_tmp(mtmp1:mtmp2,ntmp1:ntmp2).*uint16(msk_tmp(mtmp1:mtmp2,ntmp1:ntmp2));
                img1 = img2_tmp(mtmp1:mtmp2,ntmp1:ntmp2).*uint16(msk_tmp(mtmp1:mtmp2,ntmp1:ntmp2));
                clear msk_tmp img2_tmp
            end
            clear overlap_msk msk_tmp
        else
            img1 = image1tif(mtmp1:mtmp2,ntmp1:ntmp2,color);
        end
        img1_corner = [1 1];

        % register downsampled image
        if (search_test == 0)
            if (Ntran == 0)
                reg_out = [round(Xstar(5)*x0 + Xstar(6)*y0 + Xstar(7)*x0*y0 + Xstar(8)) round(Xstar(1)*x0 + Xstar(2)*y0 + Xstar(3)*x0*y0 + Xstar(4))];
            end
            init_reg0 = reg_out;
            fail = 0;
        else
            [cctable_out fail] = init_pass(img1,img2,img1_corner,addsharp,gpu);
            clear init_pass img1
            corner1 = [cctable_out(1,5) cctable_out(1,6)];
            corner1 = corner1 + [(mtmp1-1) (ntmp1-1)];
            init_reg0 = corner1/down - [1/(2*down)+1 1/(2*down)+1];
            init_reg0 = round(init_reg0);
            if (init_reg0(1)<1)
                init_reg0(1) = 1;
            end
            if (init_reg0(2)<1)
                init_reg0(2) = 1;
            end
        end

        if (fail == 0)
            if (boost == 0)
                [Xstar0,pass,img2_reg,zd,blknum,mxd1,mf1] = incremental_registration(image0bin(:,:,color),img2_reg,zd,blknum,i,img21,max_shift,samp_ratio,init_reg0,mxd,j01,j02,turnonoutput,turnonplot,useresults,filteron,m2,n2,p1,minpts,feather_width,gpu,morepts,color,memory_limited,Npoly,Ntran,addsharp,debugmode,usegui);
            else
                [Xstar0,pass,img2_reg,zd,blknum,mxd1,mf1] = incremental_registration(image0bin(:,:,color),img2_reg,zd,blknum,i,img21,max_shift,samp_ratio,init_reg0,mxd,j01,j02,turnonoutput,turnonplot,useresults,filteron,m2,n2,p1,minpts,feather_width,gpu,1,color,memory_limited,Npoly,Ntran,addsharp,debugmode,usegui);
            end
            clear incremental_registration
        else
            mult = mult + 1;
            if (mult > (mult0+1))
                mult = mult0 + 1;
            end
        end

        if (pass == 0)
            if (search_test == 0)
                if ((p1 > 1) && (color < p1))
                    repcnt
                    boost
                    color = color + 1   %change color band
                 else
                    color = 1;
                    search_test = 1;
                end
            else
                if ((p1 > 1) && (color < p1))
                    repcnt
                    boost
                    color = color + 1   %change color band
                elseif (boost == 0)
                    repcnt
                    boost = 1
                    color = 1
               else
                    repcnt = repcnt + 1 %increase search area
                    boost
                    color = 1
                    %{
                    mult = mult + 1
                    if (mult > (mult0+1))
                        mult = mult0 + 1;
                    end
                    %}
                end
            end
        else
            if (Ntran == 0)
                search_test = 0;
            else
                search_test = 1;
            end
            init_reg00(i,:) = init_reg0;
            Xstar00(i,:) = Xstar0;    
        end
    end
    clear img21
    if (pass == 1)
        delete(fn_full2);
        fn_out = [fn_path trial '_IR.tif'];
        if (bd2 == 16)
            imwrite(uint16(img2_reg(m0:(m0+m1-1),n0:(n0+n1-1))),fn_out,'tif','Compression','None');
        elseif (bd2 == 8)
            imwrite(uint8(img2_reg(m0:(m0+m1-1),n0:(n0+n1-1))),fn_out,'tif','Compression','None');
        end
        fn_zd = [fn_path trial '_dist.tif'];
        fn_blknum = [fn_path trial '_blk.tif'];
        imwrite(uint16(zd(m0:(m0+m1-1),n0:(n0+n1-1))),fn_zd,'tif','Compression','None');
        imwrite(uint16(blknum(m0:(m0+m1-1),n0:(n0+n1-1))),fn_blknum,'tif','Compression','None');
        mult = 1;
        previ = i;
    else
        if ((fail == 1) && (turnonoutput==1) && (usegui == 0))
            disp('Unable to obtain initial registration.');
            disp('FAIL');
        end
        mult = mult0 + 1;
        init_reg00(i,:) = zeros(1,2);
        if (Ntran == 0)
            Xstar00(i,:) = zeros(1,8);
        else
            Xstar00(i,:) = zeros(1,sum(2*(1:(Ntran+1))));
        end
    end
end
clear img2
end

function [cctable_out,fail] = init_pass(img1,img211,img1_corner,addsharp,gpu)
%% *********************************************
[m1sub n1sub p1sub] = size(img1);
[m2sub n2sub] = size(img211);

% phase image
img21 = single(img211);
clear img211
%mu21 = mean2(img211);
%sigma21 = std2(img211);
if (gpu == 1)
    img211 = gpuArray(img21);
    clear img21
    img21 = fft2(img211);
    clear img211
    img211 = atan2(imag(img21),real(img21));
    clear img21
    img21 = exp(1i*img211);
    clear img211
    img211 = real(ifft2(img21));
    clear img21
    img21 = gather(img211);
    clear img211
else
    img211 = fft2(img21);
    clear img21
    img21 = atan2(imag(img211),real(img211));
    clear img211
    img211 = exp(1i*img21);
    clear img21
    img21 = real(ifft2(img211));
    clear img211
end

cctable = zeros(p1sub,6);
for p = 1:p1sub
    img11 = single(img1(:,:,p));
    %mu11 = mean2(img11);
    %sigma11 = std2(img11);
    %img11 = (img11 - mu11)/sigma11*sigma21 + mu21;
    % phase image
    if (gpu == 1)
        img111 = gpuArray(img11);
        clear img11
        img11 = fft2(img111);
        clear img111
        img111 = atan2(imag(img11),real(img11));
        clear img11
        img11 = exp(1i*img111);
        clear img111
        img111 = real(ifft2(img11));
        clear img11
        img11 = gather(img111);
        clear img111
    else
        img111 = fft2(img11);
        clear img11
        img11 = atan2(imag(img111),real(img111));
        clear img111
        img111 = exp(1i*img11);
        clear img11
        img11 = real(ifft2(img111));
        clear img111
    end

    if (max(max(img21(:))) ~= min(img21(:)) && (m1sub>=m2sub) && (n1sub>=n2sub))
        cc = normxcorr2(img21,img11);
        cc(1:(m2sub-1),:) = [];
        cc(:,1:(n2sub-1)) = [];
        cc((end-m2sub+1):end,:) = [];
        cc(:,(end-n2sub+1):end) = [];
    else
        %cc = zeros(m2sub+m1sub-1,n2sub+n1sub-1);
        cc = 0;
    end
    clear img11
    
    if (addsharp ~= 1)
        sigma_cc = std2(cc);
        noise_floor = 4*sigma_cc;   %snf>2
    else
        noise_floor = 0;
    end
    [p1, imax] = max(cc(:));
    if (~isempty(imax))
        [ycc, xcc] = ind2sub(size(cc),imax(1));
    else
        ycc = [];
        xcc = [];
    end
    if (isempty(ycc))
        p1 = 0;
    end
    clear cc

    if (p1 > noise_floor)
        cctable(p,1) = p;
        cctable(p,2) = p1;
        cctable(p,3) = ycc(1);
        cctable(p,4) = xcc(1);
        cctable(p,5) = img1_corner(p,1) + ycc(1) - 1;
        cctable(p,6) = img1_corner(p,2) + xcc(1) - 1;
    else
        cctable(p,1) = 0;
        %cctable(p,2) = p1;
        cctable(p,2) = 0;
        cctable(p,3) = 0;
        cctable(p,4) = 0;
        cctable(p,5) = img1_corner(p,1) - 1;
        cctable(p,6) = img1_corner(p,2) - 1;
    end
    pause(0.1)
end
[~, ind] = max(cctable(:,2));

cctable_out = zeros(1,6);
cctable_out(1) = cctable(ind(1),1);
cctable_out(2) = cctable(ind(1),2);
cctable_out(3) = cctable(ind(1),3);
cctable_out(4) = cctable(ind(1),4);
cctable_out(5) = cctable(ind(1),5);
cctable_out(6) = cctable(ind(1),6);

fail = 0;
%if (sum(cctable(:,1)) == 0)
if (isempty(ind) || (cctable_out(2) == 0))
    fail = 1;
end
clear img1 img21 img1_corner
end

function [X,pass,img2_reg,zd,blknum,mxd1,mf1] = incremental_registration(image0bin,img2_reg,zd,blknum,blknum0,img21,max_shift,samp_ratio,init_reg,mxd,j01,j02,turnonoutput,turnonplot,useresults,filteron,m2,n2,p1,minpts,feather_width,gpu,morepts,color,memory_limited,Npoly,Ntran,addsharp,debugmode,usegui)
%% 

pass = 1;

[X mxd1 mf1 F2] = art_regx(image0bin,img2_reg,img21,init_reg,max_shift,samp_ratio,mxd,j01,j02,useresults,filteron,m2,n2,p1,minpts,gpu,morepts,color,memory_limited,Npoly,Ntran,turnonplot,addsharp,debugmode,usegui);
clear art_regx
if (usegui ~= 0)
    strfd = ['Final disparity: ', num2str(mxd1)];
    strnum = ['Number of feature point pairs: ', num2str(mf1)];
    if (mxd1 == 0)
        strpf = 'PASS';
    else
        strpf = 'FAIL';
    end
    if (usegui == 1)
        handles = guidata(ArtRegister);
    elseif (usegui == 2)
        handles = guidata(HSIRegister);
        set(handles.running,'String','Running ...');
    elseif (usegui == 3)
        handles = guidata(XrayRegister);
    end
    set(handles.text17,'String',strfd);
    set(handles.text18,'String',strnum);
    set(handles.text19,'String',strpf);
elseif (turnonoutput==1)
    disp(['Final disparity: ', num2str(mxd1)]);
    disp(['Number of feature point pairs: ', num2str(mf1)]);
    if (mxd1 == 0)
        disp('PASS');
    else
        disp('FAIL');
    end
end

if (mxd1 == 0)
    [img20,m21,n21] = apply_mapping(img21,init_reg,X,Ntran,m2,n2);
    [mimg21 nimg21] = size(img21);
    [tmp,~,~] = apply_mapping(2^15*ones(mimg21,nimg21),init_reg,X,Ntran,m2,n2);
    mskimg21 = (round(tmp) == 2^15);
    clear tmp
    % removed 20140222
    %{
    se = strel('disk',3);
    mskimg21 = imclose(mskimg21,se);
    se = strel('disk',3);
    mskimg21 = imerode(mskimg21,se);
    mskimg21 = imfill(mskimg21,'holes');
    img20 = img20.*mskimg21;
    %}
    
    [v1,u1] = meshgrid(0:(nimg21-1),0:(mimg21-1));
    zd0 =  single(sqrt((u1-((mimg21-1)/2)).^2 + (v1-((nimg21-1)/2)).^2));
    clear u1 v1
    [zd1,~,~] = apply_mapping(zd0,init_reg,X,Ntran,m2,n2);
    zd1(mskimg21 == 0) = sqrt(mimg21^2 + nimg21^2);
    zd1 = uint16(zd1);
    clear mskimg21 zd0
    %figure(101),imagesc(zd1)
    %pause
    clear apply_mapping
    [m0 n0] = size(img20);
    
    % mosaicing
    %removed20120222
    [img2_reg,zd,blknum] = mosaic1(img2_reg,zd,zd1,blknum,blknum0,img20,m21,n21);
    clear mosaic1
    clear img20 zd1

    if (usegui ~= 0)
        if (usegui == 1)
            handles = guidata(ArtRegister);
        elseif (usegui == 2)
            handles = guidata(HSIRegister);
            set(handles.running,'String','Running ...');
        elseif (usegui == 3)
            handles = guidata(XrayRegister);
        end
        tmp2 = single(img2_reg);
        tmp2 = tmp2/max(tmp2(:));
        tmp2 = 255*tmp2;
        %{
        [mtmp2 ntmp2] = size(tmp2);
        minx = floor(min(F2(:,1)))-256;
        if (minx < 1)
            minx = 1;
        end
        maxx = ceil(max(F2(:,1)))+256;
        if (maxx > ntmp2)
            maxx = ntmp2;
        end
        miny = floor(min(F2(:,2)))-256;
        if (miny < 1)
            miny = 1;
        end
        maxy = ceil(max(F2(:,2)))+256;
        if (maxy > mtmp2)
            maxy = mtmp2;
        end
        %}
        if (debugmode == 1)
            set(handles.axes5,'Visible','on');
            axes(handles.axes5)
            %imshow(uint8(tmp2(miny:maxy,minx:maxx)));
            imshow(uint8(tmp2(m21:(m21+m0-1),n21:(n21+n0-1))));
            hold on
            %plot(round(F2(:,1))-minx+1,round(F2(:,2)-miny+1),'r.')
            plot(round(F2(:,1))-n21+1,round(F2(:,2)-m21+1),'r.')
            hold off
        else
            set(handles.axes1,'Visible','on');
            axes(handles.axes1)
            imshow(uint8(tmp2(1:4:end,1:4:end)));
        end
        set(handles.text17,'String',strfd);
        set(handles.text18,'String',strnum);
        set(handles.text19,'String',strpf);
        clear tmp2
    elseif (turnonplot==1)
        tmp2 = single(img2_reg);
        tmp2 = tmp2/max(tmp2(:));
        tmp2 = 255*tmp2;
        figure(11)
        imshow(uint8(tmp2))
        clear tmp2
    end
else
    pass = 0;
end
clear img1 img21
pause(0.1)
end

function [X,mxd1,mf1,F2] = art_regx(image0bin,img2_reg,img1,init_reg,max_shift,samp_ratio,mxd,j01,j02,useresults,filteron,m2,n2,p1,minpts,gpu,morepts,color,memory_limited,Npoly,Ntran,turnonplot,addsharp,debugmode,usegui)

% constants
alloc = 1024;   %pre-allocation increment

% IR image
img11 = single(img1);
if (filteron == 1)
    img11 = medfilt2(img11,[3 3]);
end
[m1 n1] = size(img11);

% extract features from subregions for the frames of the first set of images
F1 = zeros(alloc,3);
pnt = 1;
for j = j01:j02
    [F] = total_features(single(img11),j,samp_ratio,morepts);
    clear total_features
    mt = size(F,1);
    mf1 = size(F1,1);    
    while ((mt+pnt-1) > mf1)
        F1 = cat(1,F1,zeros(alloc,3));
        mf1 = size(F1,1);
    end
    F1(pnt:(mt+pnt-1),:) = F;
    pnt = mt + pnt;
    clear F
end
F1 = F1(((F1(:,1)~=0) & (F1(:,2)~=0) & (F1(:,3)~=0)),:);
mf1 = size(F1,1);
clear img11

% phase image
img111 = single(img1);

if (addsharp == 1)
    %LOW_HIGH = stretchlim(uint16(img111),[0.001 0.999]);
    %img111 = single(imadjust(uint16(img111),LOW_HIGH,[]));
    %h = fspecial('unsharp',0.8);
    %img111 = imfilter(img111,h);
    %img111(img111<0) = 0;
    %img111 = img111/max(img111(:))*(2^16-1);
end
clear img1
%{
if (filteron == 1)
    img111 = medfilt2(img111,[3 3]);
end
%}
if (gpu == 1)
    img11 = gpuArray(img111);
    clear img111
    img1 = fft2(img11);
    clear img11
    img11 = atan2(imag(img1),real(img1));
    clear img1
    img1 = exp(1i*img11);
    clear img11
    img11 = real(ifft2(img1));
    clear img1
    img1 = gather(img11);
    clear img11
else
    img11 = fft2(img111);
    clear img111
    img1 = atan2(imag(img11),real(img11));
    clear img11
    img11 = exp(1i*img1);
    clear img1
    img1 = real(ifft2(img11));
    clear img11
end

if (useresults == 1)
    if (gpu == 1)
        img2_reg1 = gpuArray(img2_reg);
        clear img2_reg
        img2_reg = fft2(img2_reg1);
        clear img2_reg1
        img2_reg1 = atan2(imag(img2_reg),real(img2_reg));
        clear img2_reg
        img2_reg = exp(1i*img2_reg1);
        clear img2_reg1
        img2_reg1 = real(ifft2(img2_reg));
        clear img2_reg
        img2_reg = gather(img2_reg1);
        clear img2_reg1
    else
        img2_reg1 = fft2(img2_reg);
        clear img2_reg
        img2_reg = atan2(imag(img2_reg1),real(img2_reg1));
        clear img2_reg1
        img2_reg1 = exp(1i*img2_reg);
        clear img2_reg
        img2_reg = real(ifft2(img2_reg1));
        clear img2_reg1
    end
end

% pair up feature points
F2 = zeros(mf1,3);
img11 = single(zeros(m2,n2));
img11(init_reg(1):(init_reg(1)+m1-1),init_reg(2):(init_reg(2)+n1-1)) = img1;
F1(:,1:2) = F1(:,1:2) + repmat(([init_reg(2) init_reg(1)]-[1 1]),[mf1 1]);
clear img1

for i = mf1:-1:1
    %if ((F1(i,1) < 1) || (F1(i,2) < 1) || (F1(i,1) > n2) || (F1(i,2) > m2) || (msk(round(F1(i,2)),round(F1(i,1))) == 0))
    if ((F1(i,1) < 1) || (F1(i,2) < 1) || (F1(i,1) > n2) || (F1(i,2) > m2))
        F1(i,:) = [];
        F2(i,:) = [];
    else
        j = F1(i,3);
        lt = floor(F1(i,1))-((2^j)-1);
        rt = ceil(F1(i,1))+((2^j)-1);
        tp = floor(F1(i,2))-((2^j)-1);
        bt = ceil(F1(i,2))+((2^j)-1);
        tmp = img11(tp:bt,lt:rt);
        tmp1 = imresize(tmp,samp_ratio,'bicubic');
        clear tmp
        d = floor(samp_ratio/2)/samp_ratio;
        tmp1n = (lt-d):1/samp_ratio:(rt+d);
        [~, x1] = min(abs(tmp1n-F1(i,1)));
        clear tmp1n
        tmp1m = (tp-d):1/samp_ratio:(bt+d);
        [~, y1] = min(abs(tmp1m-F1(i,2)));
        clear tmp1m
        d = (((2^(j+1))-2+1)*samp_ratio-1)/2;
        tmp1 = tmp1((y1-d):(y1+d),(x1-d):(x1+d));

        lt = lt - max_shift;
        if (lt <= 1)
            lt = 1;
        elseif (lt >= (n2-2^(j+1)+1))
            lt = n2-2^(j+1)+1;
        end
        rt = rt + max_shift;
        if (rt >= n2)
            rt = n2;
        end
        tp = tp - max_shift;
        if (tp <= 1)
            tp = 1;
        elseif (tp >= (m2-2^(j+1)+1))
            tp = m2-2^(j+1)+1;
        end
        bt = bt + max_shift;
        if (bt >= m2)
            bt = m2;
        end

        % read data from disk
        if (memory_limited == 1)
            tmp = multibandread('image0.bin',[m2,n2,p1],'single',0,'bsq','ieee-le',{'Row','Range',[tp bt]},{'Column','Range',[lt rt]},{'Band','Direct',color});
        else
            tmp = image0bin(tp:bt,lt:rt);
        end
        % need to update to work with new code 20140306
        %{
        if (useresults == 1)
            img2_reg1 = single(img2_reg(tp:bt,lt:rt));
            overlap_msk = (zd < (2^16-1));
            overlap_msk1 = single(overlap_msk(tp:bt,lt:rt));
            clear overlap_msk
            %tmp = (img2_reg1.*overlap_msk1) + (tmp.*(~overlap_msk1));
            tmp = img2_reg1.*overlap_msk1;
            clear overlap_msk1 img2_reg1 
        end
        %}
        tmp2 = imresize(tmp,samp_ratio,'bicubic');
        clear tmp
        d = floor(samp_ratio/2)/samp_ratio;
        tmp2n = (lt-d):1/samp_ratio:(rt+d);
        tmp2m = (tp-d):1/samp_ratio:(bt+d);

        [mt1 nt1] = size(tmp1);
        [mt2 nt2] = size(tmp2);
        cc = abs(normxcorr2(tmp1,tmp2));
        if (addsharp ~= 1)
            sigma_cc = std2(cc);
            noise_floor = 4*sigma_cc;   %snf>2
        else
            noise_floor = 0;
        end
        mskcc = imregionalmax(cc,8);
        cc = cc.*single(mskcc);
        cc = cc((floor(mt1/2)+1):end,(floor(nt1/2)+1):end);
        [max_cc, imax] = max(cc(:));
        [ypeak, xpeak] = ind2sub(size(cc),imax(1));
        if (addsharp ~= 1)
            cnt = sum(sum(cc>(max_cc-sigma_cc)));
        else
            cnt = 1;
        end
        clear cc       
        if ((xpeak <= nt2) && (ypeak <= mt2) && (cnt == 1) && (max_cc > noise_floor))
            F2(i,1) = tmp2n(xpeak);
            F2(i,2) = tmp2m(ypeak);
            F2(i,3) = max_cc;
        else
            F1(i,:) = [];
            F2(i,:) = [];
        end
        clear tmp1 tmp2
        clear tmp2n tmp2m
    end
end

% remove bad pairs
[F1 F2 mxd1 mf1] = xy_disp(F1(:,1:2),F2,mxd,minpts,Npoly,turnonplot,debugmode,usegui);
clear xy_disp

X = [];
%if (mxd1 <= mxd)
if (mxd1 == 0)
    % transform image 1 to match image 2
    [X] = mapping(F2,F1,gpu,Ntran);
    clear mapping
end
clear F1
end

function [F] = total_features(img,j,samp_ratio,morepts)
%% total_features

[m0 n0] = size(img);

h = fspecial('unsharp',0.8);
img = imfilter(img,h);
img(img<0) = 0;

% Wavelet transform
h = fspecial('average',[1 2^j]);    %hsise = 1 x 2^j
LPx = imfilter(img,h);

h = [];
for q = 1:j
    h = [-1 h];
    h = [h 1];
end
LH = imfilter(LPx,h');
clear LPx
HPy = imfilter(img,h);
h = fspecial('average',[2^j 1]);    %hsize = 2^j x 1
HL = imfilter(HPy,h);
clear HPy

% Modulus of the wavelet transform
M = sqrt(LH.^2+HL.^2);
clear LH;
clear HL;
F0 = imregionalmax(M,8);

% Remove border pixels
hw = (2^j) + 1;
F0(1:hw,:) = 0;
F0(:,1:hw) = 0;
F0((m0-hw+1):m0,:) = 0;
F0(:,(n0-hw+1):n0) = 0;

[y1,x1] = find(F0==1);
yd = length(y1);
for i = 1:yd
    if (morepts == 0)
        tmp = F0((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)).*M((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j));
    elseif (morepts == 1)
        tmp = F0((y1(i)-2^(j-1)):(y1(i)+2^(j-1)),(x1(i)-2^(j-1)):(x1(i)+2^(j-1))).*M((y1(i)-2^(j-1)):(y1(i)+2^(j-1)),(x1(i)-2^(j-1)):(x1(i)+2^(j-1)));
    end
    
    [mx11,ind11] = max(tmp(:));
    if (mx11 > 0)
        xpeak = floor((ind11(1)-1)/size(tmp,1))+1;
        ypeak = ind11(1) - (xpeak-1)*size(tmp,1);
        msk11 = single(zeros(size(tmp)));
        msk11(ypeak,xpeak) = 1;
        if (morepts == 0)
            F0((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)) = F0((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)).*msk11;
        elseif (morepts == 1)
            F0((y1(i)-2^(j-1)):(y1(i)+2^(j-1)),(x1(i)-2^(j-1)):(x1(i)+2^(j-1))) = F0((y1(i)-2^(j-1)):(y1(i)+2^(j-1)),(x1(i)-2^(j-1)):(x1(i)+2^(j-1))).*msk11;
        end
        
    end
    clear tmp
end
clear x1 y1

[y1,x1] = find(F0==1);
clear F0
yd = length(y1);
F = zeros(yd,3);

pnt = 1;
% subsample and find subpixel maxima
for i = 1:yd
    lt = x1(i)-1;
    rt = x1(i)+1;
    tp = y1(i)-1;
    bt = y1(i)+1;
    tmp = M(tp:bt,lt:rt);
    tmp1 = imresize(tmp,samp_ratio,'bicubic');
    clear tmp
    d = floor(samp_ratio/2)/samp_ratio;
    tmp1n = (lt-d):1/samp_ratio:(rt+d);
    tmp1m = (tp-d):1/samp_ratio:(bt+d);
    max_tmp1 = max(tmp1(:));

    % if there are multiple maximum, pick the closest to the original
    [ymx,xmx] = find(tmp1==max_tmp1);
    clear tmp1
    ind_dist = 1;
    if (length(ymx) > 1)
        dist = zeros(length(ymx));
        for k = 1:length(ymx)
            dist(k) = (xmx(k)-x1(i))^2 + (ymx(k)-y1(i))^2;
        end
        [~, ind_dist] = max(dist(:));
    end

    F(pnt,1) = tmp1n(xmx(ind_dist));
    F(pnt,2) = tmp1m(ymx(ind_dist));
    F(pnt,3) = j;
    pnt = pnt + 1;
    clear tmp1n tmp1m
end
clear M
F = F(((F(:,1)~=0) & (F(:,2)~=0) & (F(:,3)~=0)),:);
clear img
end

function [F1 F2 mxd1 mf1] = xy_disp(F1,F2,mxd,minpts,Npoly,turnonplot,debugmode,usegui)
%% X/Y Disparity Minimumizing algorithm
% Usage: [F1 F2] = xy_disp(F1,F2,mxd,minpts)
%   inputs: F1,F2 -> two pairs of feature points
%           mxd -> maximum allowable disparity (pixels)
%           minpts -> mimumum number of final feature pairs
%   outputs: F1,F2 -> two pairs of optimal features points taken from the
%               original pairs
%

%csvwrite('f1_before54.csv',F1)
%csvwrite('f2_before54.csv',F2)

mxd_max = 1;
mf1 = size(F1,1);
minx = floor(min(F1(:,1)));
maxx = ceil(max(F1(:,1)));
miny = floor(min(F1(:,2)));
maxy = ceil(max(F1(:,2)));
midx = (maxx-minx)/2 + minx;
midy = (maxy-miny)/2 + miny;
ax = (mxd_max-mxd)/((maxx-midx)^Npoly);
ay = (mxd_max-mxd)/((maxy-midy)^Npoly);

xdist = F2(:,1) - F1(:,1);  %x-disperity
ydist = F2(:,2) - F1(:,2);  %y-disperity
mxx = Inf;  %maximum disperity
mxy = Inf;
test = 1;
%while ((mxx>=mxd || mxy>=mxd) && mf1>=minpts)
while ((mxx>0 || mxy>0) && mf1>=minpts)
    % fit best planes to x and y disparities
    % compute error between the displarites and the planes
    if (Npoly == 0)
        A = [F1(:,1) F1(:,2) F1(:,1).*F1(:,2) ones(mf1,1)];
        pA = pinv(A);
        clear A
        xx = pA*xdist;
        xy = pA*ydist;
        clear pA
    else
        %Aineq = diag(-ones(1,sum((1:(Npoly+1)))),0);
        %Aineq(1,:) = [];
        %bineq = zeros(sum((1:(Npoly+1))),1);
        %bineq(1) = [];
        Aineq = [];
        bineq = [];
        A = zeros(mf1,sum((1:(Npoly+1))));
        pnt = 1;
        for i = 0:Npoly
            for j = 0:(Npoly-i)
                A(:,pnt) = ((F1(:,1)).^i).*((F1(:,2)).^j);
                pnt = pnt + 1;
            end
        end
        opts = optimset('Display','off');
        xx = lsqlin(double(A),double(xdist),double(Aineq),double(bineq),[],[],[],[],[],opts);
        xy = lsqlin(double(A),double(ydist),double(Aineq),double(bineq),[],[],[],[],[],opts);
    end
    clear A
    errorx = xdist;
    errory = ydist;
    if (Npoly == 0)
        errorx = errorx - xx(1)*F1(:,1) - xx(2)*F1(:,2) - xx(3)*F1(:,1).*F1(:,2) - xx(4)*ones(mf1,1);
        errory = errory - xy(1)*F1(:,1) - xy(2)*F1(:,2) - xy(3)*F1(:,1).*F1(:,2) - xy(4)*ones(mf1,1);
    else
        pnt = 1;
        for i = 0:Npoly
            for j = 0:(Npoly-i)
                errorx = errorx - xx(pnt)*((F1(:,1)).^i).*((F1(:,2)).^j);
                errory = errory - xy(pnt)*((F1(:,1)).^i).*((F1(:,2)).^j);
                pnt = pnt + 1;
            end
        end
    end
    if (test == 0)
        [X,Y] = meshgrid(minx:20:maxx,miny:20:maxy);
        X = X(:);
        Y = Y(:);
        
        testx = zeros(length(X),1);
        testy = zeros(length(X),1);
        if (Npoly == 0)
            testx = testx + xx(1)*X + xx(2)*Y + xx(3)*X.*Y + xx(4)*ones(length(X),1);
            testy = testy + xy(1)*X + xy(2)*Y + xy(3)*X.*Y + xy(4)*ones(length(X),1);
        else
            pnt = 1;
            for i = 0:Npoly
                for j = 0:(Npoly-i)
                    testx = testx + xx(pnt)*(X.^i).*(Y.^j);
                    testy = testy + xy(pnt)*(X.^i).*(Y.^j);
                    pnt = pnt + 1;
                end
            end
        end
        
        figure(100)
        plot3(X,Y,testx,'r.')
        hold on
        plot3(F1(:,1),F1(:,2),xdist,'b*')
        grid on
        axis([minx maxx miny maxy])
        hold off
        figure(101)
        plot3(X,Y,testy,'r.')
        hold on
        plot3(F1(:,1),F1(:,2),ydist,'b*')
        grid on
        axis([minx maxx miny maxy])
        hold off
        test = 1;
    end
    
    % find the pair that contributes the largest error and remove it
    %if (Npoly == 0)
        errorx = abs(errorx) - mxd;
        errorx(errorx<0) = 0;
        errory = abs(errory) - mxd;
        errory(errory<0) = 0;
    %{
    else
        errorx = abs(errorx) - (abs(ax*((F1(:,1)-midx).^Npoly)) + mxd);   %threshold is a function of the distance to the center (fit is less accurate at the edges)
        errorx(errorx<0) = 0;
        errory = abs(errory) - (abs(ay*((F1(:,2)-midy).^Npoly)) + mxd);
        errory(errory<0) = 0;
    end
    %}
    [mxx indx] = max(abs(errorx));
    [mxy indy] = max(abs(errory));
    if (mxx(1)>=mxy(1))
        ind = indx(1);
    else
        ind = indy(1);
    end
    F1(ind,:) = [];
    F2(ind,:) = [];
    xdist(ind) = [];
    ydist(ind) = [];
    clear errorx errory
    mf1 = size(F1,1);
end
mxd1 = max(mxx,mxy);

test = 1;
if (mf1 >= minpts)
    if (debugmode == 1)
        %minx = floor(min(F1(:,1)));
        %maxx = ceil(max(F1(:,1)));
        %miny = floor(min(F1(:,2)));
        %maxy = ceil(max(F1(:,2)));
        [X,Y] = meshgrid(minx:20:maxx,miny:20:maxy);
        X = X(:);
        Y = Y(:);
        
        testx = zeros(length(X),1);
        testy = zeros(length(X),1);
        if (Npoly == 0)
            testx = testx + xx(1)*X + xx(2)*Y + xx(3)*X.*Y + xx(4)*ones(length(X),1);
            testy = testy + xy(1)*X + xy(2)*Y + xy(3)*X.*Y + xy(4)*ones(length(X),1);
        else
            pnt = 1;
            for i = 0:Npoly
                for j = 0:(Npoly-i)
                    testx = testx + xx(pnt)*(X.^i).*(Y.^j);
                    testy = testy + xy(pnt)*(X.^i).*(Y.^j);
                    pnt = pnt + 1;
                end
            end
        end
        
        if (usegui ~= 0)
            if (usegui == 1)
                handles = guidata(ArtRegister);
            elseif (usegui == 2)
                handles = guidata(HSIRegister);
            elseif (usegui == 3)
                handles = guidata(XrayRegister);
            end
            set(handles.axes2,'Visible','on');
            axes(handles.axes2)
            plot3(X,Y,testx,'r.')
            hold on
            plot3(F1(:,1),F1(:,2),xdist,'b*')
            grid on
            axis([minx maxx miny maxy])
            hold off
            set(handles.axes4,'Visible','on');
            axes(handles.axes4)
            plot3(X,Y,testy,'r.')
            hold on
            plot3(F1(:,1),F1(:,2),ydist,'b*')
            grid on
            axis([minx maxx miny maxy])
            hold off
        elseif (turnonplot==1)
            figure(11)
            subplot(121),plot3(X,Y,testx,'r.')
            hold on
            plot3(F1(:,1),F1(:,2),xdist,'b*')
            grid on
            axis([minx maxx miny maxy])
            hold off
            subplot(122),plot3(X,Y,testy,'r.')
            hold on
            plot3(F1(:,1),F1(:,2),ydist,'b*')
            grid on
            axis([minx maxx miny maxy])
            hold off
        end
        clear X Y testx testy
        
        if (test == 0)
            %minx = floor(min(F1(:,1)));
            %maxx = ceil(max(F1(:,1)));
            %miny = floor(min(F1(:,2)));
            %maxy = ceil(max(F1(:,2)));
            [X,Y] = meshgrid(minx:20:maxx,miny:20:maxy);
            X = X(:);
            Y = Y(:);

            testx = zeros(length(X),1);
            testy = zeros(length(X),1);
            if (Npoly == 0)
                testx = testx + xx(1)*X + xx(2)*Y + xx(3)*X.*Y + xx(4)*ones(length(X),1);
                testy = testy + xy(1)*X + xy(2)*Y + xy(3)*X.*Y + xy(4)*ones(length(X),1);
            else
                pnt = 1;
                for i = 0:Npoly
                    for j = 0:(Npoly-i)
                        testx = testx + xx(pnt)*(X.^i).*(Y.^j);
                        testy = testy + xy(pnt)*(X.^i).*(Y.^j);
                        pnt = pnt + 1;
                    end
                end
            end

            figure(102)
            plot3(X,Y,testx,'r.')
            hold on
            plot3(F1(:,1),F1(:,2),xdist,'b*')
            grid on
            axis([minx maxx miny maxy])
            hold off
            figure(103)
            plot3(X,Y,testy,'r.')
            hold on
            plot3(F1(:,1),F1(:,2),ydist,'b*')
            grid on
            axis([minx maxx miny maxy])
            hold off
            test = 1;
        end  
    end
end
clear xdist ydist

%csvwrite('f1_after54.csv',F1)
%csvwrite('f2_after54.csv',F2)
end

function [X] = mapping(F11,F21,gpu,Ntran)
%% mapping
% Transform second image so that it is registered with the first
%   Input: F11 -> control points for image 1
%          F21 -> control points for image 2
%   Outputs: z -> registered image

m2 = size(F21,1);
%{
if (gpu == 1)
    F11 = gpuArray(F11);
    F21 = gpuArray(F21);
    if (Ntran == 0)
        A = parallel.gpu.GPUArray.zeros(2*m2,8,'double');
    else
        A = parallel.gpu.GPUArray.zeros(2*m2,sum(2*(1:(Ntran+1))),'double');
    end
else
    if (Ntran == 0)
        A = zeros(2*m2,8);
    else
        A = zeros(2*m2,sum(2*(1:(Ntran+1))));
    end
end
%}
if (Ntran == 0)
    A = zeros(2*m2,8);
else
    A = zeros(2*m2,sum(2*(1:(Ntran+1))));
    %Aineq = diag(-ones(1,sum(2*(1:(Ntran+1)))),0);
    %Aineq(sum((1:(Ntran+1)))+1,:) = [];
    %Aineq(1,:) = [];
    %bineq = zeros(sum(2*(1:(Ntran+1))),1);
    %bineq(sum((1:(Ntran+1)))+1) = [];
    %bineq(1) = [];
    Aineq = [];
    bineq = [];
end
% bilinear transformation
Y = [F11(:,1); F11(:,2)];
x = F21(:,1);
y = F21(:,2);
if (Ntran == 0)
    A(1:m2,:) = [x y x.*y ones(m2,1) zeros(m2,4)];
    A((m2+1):2*m2,:) = [zeros(m2,4) x y x.*y ones(m2,1)];
    % pseudoinverse solution
    X = pinv(double(A))*double(Y);
else
    pnt = 1;
    for i = 0:Ntran
        for j = 0:(Ntran-i)
            A(1:m2,pnt) = (x.^i).*(y.^j);
            A((m2+1):(2*m2),(pnt+sum((1:(Ntran+1))))) = (x.^i).*(y.^j);
            pnt = pnt + 1;
        end
    end
    options = optimset('Display','off');
    X = lsqlin(double(A),double(Y),double(Aineq),double(bineq),[],[],[],[],[],options);
end
clear x y

% pseudoinverse solution
%{
if (gpu == 1)
    X = (A'*A)\(A')*Y;
    X = gather(X);
else
    X = pinv(double(A))*double(Y);
end
%}
clear A Y F11 F21
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

function [img2_reg,zd,blknum] = mosaic1(img2_reg,zd,zd1,blknum,blknum0,img20,m21,n21)
%% mosaicing

[m0,n0] = size(img20);
m22 = m21 + m0 - 1;
n22 = n21 + n0 - 1;
zd0 = zd(m21:m22,n21:n22);
msk = (zd1<zd0);
img2_reg(m21:m22,n21:n22) = img2_reg(m21:m22,n21:n22).*single(~msk) + img20.*single(msk);
zd(m21:m22,n21:n22) = zd(m21:m22,n21:n22).*uint16(~msk) + zd1.*uint16(msk);
blknum(m21:m22,n21:n22) = blknum(m21:m22,n21:n22).*uint16(~msk) + blknum0*uint16(msk);
clear msk img20
end
