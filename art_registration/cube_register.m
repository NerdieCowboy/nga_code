function [] = cube_register(cube_fn,j01,j02,repcnt_max,minpts,max_shift,filteron,morepts,usegui)
%% Multspecral Registration algorithm
% cube_register.m
% Author: Damon Conover
% Email: dconover@gwmail.gwu.edu
% Latest Revision: 15 December 2012
%
% Usage: cube_register(cube_fn,j01,j02,repcnt_max,minpts,max_shift,filteron,morepts,usegui)
%   inputs: cube_fn -> filename of multispecral cube
%           j01,j02 -> min/max scales for feature selection
%           repcnt_max -> repeat search region (note:  repcnt_max=2 means that 
%               max_shift*(1+0), max_shift*(1+1), max_shift*(1+2) will be used 
%               before moving on to the next template image)
%           max_shift -> maximum number of pixels that images can shift between blocks
%           minpts -> minimum number of final pairs
%           more_pts -> force maximum number of feature points to be extracted -> 1, otherwise -> 0
%           filteron -> 1 (median filter on), 0 (median filter off)
%           usegui -> use compiled code and ArtRegister GUI -> 1, otherwise -> 0
%   outputs: Two registered images will be written to the working directory
%

warning('off', 'all')
handles = guidata(MSIRegister);
set(handles.running,'Visible','on');
set(handles.running,'String','Running ...');
drawnow;
pause(0.1)
isOpen = matlabpool('size') > 0;
if (isOpen == 1)
    matlabpool close
end

% settings
samp_ratio = 3; %upsample ratio for determining fractional pixels (must be an odd integer)

% constants
mxd = 1/samp_ratio;  %maximum allowable disparity (pixels)
turnonoutput = 1;
if (usegui == 1)
    turnonoutput = 0;
end

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

tmp = regexp(cube_fn,'\\');
if (isempty(tmp))
    tmp = regexp(cube_fn,'\/');
end
tmp1 = tmp(end);
tmp2 = length(cube_fn);
fn_path = cube_fn(1:tmp1);
trial = cube_fn((tmp1+1):tmp2);
clear tmp tmp1 tmp2

% cube
fnh = [cube_fn '.hdr'];
fid = fopen(fnh);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
    [~,~,e] = regexp(line,'samples=','match','start','end');
    if (e>0)
        n1 = str2double(line((e+1):end));
    end
    [~,~,e] = regexp(line,'lines=','match','start','end');
    if (e>0)
        m1 = str2double(line((e+1):end));
    end
    [~,~,e] = regexp(line,'^bands=','match','start','end');
    if (e>0)
        bands1 = str2double(line((e+1):end));
    end
    [~,~,e] = regexp(line,'datatype=','match','start','end');
    if (e>0)
        precision1 = str2double(line((e+1):end));
        if (precision1 == 1)
            precision1 = 'uint8';
            elseif (precision1 == 4)
                precision1 = 'single';
            elseif (precision1 == 5)
                precision1 = 'double';
            elseif (precision1 == 12)
                precision1 = 'uint16';
        end
    end
    [~,~,e] = regexp(line,'interleave=','match','start','end');
    if (e>0)
        interleave1 = line((e+1):end);
    end
    [~,~,e] = regexp(line,'byteorder=','match','start','end');
    if (e>0)
        byteorder1 = str2double(line((e+1):end));
        if (byteorder1 == 0)
            byteorder1 = 'ieee-le';
        elseif (byteorder1 == 1)
            byteorder1 = 'ieee-be';
        end
    end
    [~,~,e] = regexp(line,'headeroffset=','match','start','end');
    if (e>0)
        offset1 = str2double(line((e+1):end));
    end
end
fclose(fid);
clear msk

% feature mask
fn_out = [fn_path trial '_mask.tif'];
if (exist(fn_out,'file') == 2)
    msk = logical(imread(fn_out)>0);
else
    msk = true(m1,n1);
end

cube1 = multibandread(cube_fn,[m1 n1 bands1],precision1,offset1,interleave1,byteorder1);
cube1(cube1==Inf) = 1;
cube1(cube1==-Inf) = -1;
cube1(isnan(cube1)) = 0;

ncores = str2double(get(handles.ncores,'String'));
matlabpool(ncores)
mxd0 = zeros(bands1-1);
mftot = zeros(bands1-1);
parfor i = 2:bands1
    repcnt = 0;
    while (repcnt <= repcnt_max)
        tmp1 = cube1(:,:,i);
        tmp2 = cube1(:,:,i-1);
        if (filteron == 1)
            tmp1 = medfilt2(tmp1,[3 3]);
            tmp2 = medfilt2(tmp2,[3 3]);
        end
        [F11 F21 mxd1] = art_regx(tmp1,tmp2,msk,max_shift,samp_ratio,mxd,j01,j02,minpts,morepts);
        mxd0(i-1) = mxd1;
        mftot(i-1) = size(F11,1);
        
        if (mxd1 <= mxd)
            F11(:,1) = F11(:,1) - round(n1/2);
            F21(:,1) = F21(:,1) - round(n1/2);
            F11(:,2) = F11(:,2) - round(m1/2);
            F21(:,2) = F21(:,2) - round(m1/2);
            F1{i-1} = F11;
            F2{i-1} = F21;
            repcnt = repcnt_max + 1;            
        else
            repcnt = repcnt + 1;
        end
    end
end
clear F11 F21 tmp1 tmp2

midband = floor(bands1/2)+1;
Yref = zeros(2*sum(mftot),1);
Atemp = zeros(2*sum(mftot),8*(bands1-1));
pnt = 1;
for i = 2:bands1
    j = i-1;    %j->ref, i->temp
    F11 = F1{i-1};
    F21 = F2{i-1};
    % x
    mf1 = size(F21,1);
    if (i == midband)
        Yref(pnt:(mf1+pnt-1),1) = -F11(:,1);
    elseif (i < midband)
        Atemp(pnt:(mf1+pnt-1),((i-1)*8+1):((i-1)*8+4)) = [F11(:,1) F11(:,2) F11(:,1).*F11(:,2) ones(mf1,1)];
    elseif (i > midband)
        Atemp(pnt:(mf1+pnt-1),((i-2)*8+1):((i-2)*8+4)) = [F11(:,1) F11(:,2) F11(:,1).*F11(:,2) ones(mf1,1)];
    end
    if (j == midband)
        Yref(pnt:(mf1+pnt-1),1) = F21(:,1);
    elseif (j < midband)
        Atemp(pnt:(mf1+pnt-1),((j-1)*8+1):((j-1)*8+4)) = [-F21(:,1) -F21(:,2) -F21(:,1).*F21(:,2) -ones(mf1,1)];
    elseif (j > midband)
        Atemp(pnt:(mf1+pnt-1),((j-2)*8+1):((j-2)*8+4)) = [-F21(:,1) -F21(:,2) -F21(:,1).*F21(:,2) -ones(mf1,1)];
    end
    pnt = mf1 + pnt;

    % y
    if (i == midband)
        Yref(pnt:(mf1+pnt-1),1) = -F11(:,2);
    elseif (i < midband)
        Atemp(pnt:(mf1+pnt-1),((i-1)*8+5):(i*8)) = [F11(:,1) F11(:,2) F11(:,1).*F11(:,2) ones(mf1,1)];
    elseif (i > midband)
        Atemp(pnt:(mf1+pnt-1),((i-2)*8+5):((i-1)*8)) = [F11(:,1) F11(:,2) F11(:,1).*F11(:,2) ones(mf1,1)];
    end
    if (j == midband)
        Yref(pnt:(mf1+pnt-1),1) = F21(:,2);
    elseif (j < midband)
        Atemp(pnt:(mf1+pnt-1),((j-1)*8+5):(j*8)) = [-F21(:,1) -F21(:,2) -F21(:,1).*F21(:,2) -ones(mf1,1)];
    elseif (j > midband)
        Atemp(pnt:(mf1+pnt-1),((j-2)*8+5):((j-1)*8)) = [-F21(:,1) -F21(:,2) -F21(:,1).*F21(:,2) -ones(mf1,1)];
    end
    pnt = mf1 + pnt;
    
    if (usegui == 1)
        strband = ['Band: ', num2str(i)];
        strdisp = ['Final disparity: ', num2str(mxd0(i-1))];
        strnum = ['Number of feature point pairs: ', num2str(mf1)];
        if (mxd0(i-1) <= mxd)
            strpf = 'PASS';
        else
            strpf = 'FAIL';
        end
        handles = guidata(MSIRegister);
        set(handles.text17,'String',strband);
        set(handles.text20,'String',strdisp);
        set(handles.text18,'String',strnum);
        set(handles.text19,'String',strpf);
        set(handles.running,'String','Running ...');
        drawnow;
    elseif (turnonoutput==1)
        disp(['Band: ', num2str(i)]);
        disp(['Final disparity: ', num2str(mxd0(i-1))]);
        disp(['Number of feature point pairs: ', num2str(mf1)]);
        if (mxd0(i-1) <= mxd)
            disp('PASS');
        else
            disp('FAIL');
        end
    end
    pause(1)
end
clear F1 F2 F11 F21

X = pinv(double(Atemp))*double(Yref);

cube2 = zeros(m1,n1,bands1);
cube2(:,:,midband) = cube1(:,:,midband);
parfor i = 1:bands1
    if (i < midband)
        Xstar = X(((i-1)*8+1):(i*8));
        Xstar = Xstar';
        cube2(:,:,i) = apply_mapping(cube1(:,:,i),[1 1],Xstar,1024,m1,n1);
    elseif (i > midband)
        Xstar = X(((i-2)*8+1):((i-1)*8));
        Xstar = Xstar';
        cube2(:,:,i) = apply_mapping(cube1(:,:,i),[1 1],Xstar,1024,m1,n1);
    end
end
matlabpool close

cube_fn = [fn_path trial '_registered'];
multibandwrite(cube2,cube_fn,interleave1,'precision',precision1,'machfmt',byteorder1);
fn_out = [fn_path trial '.hdr'];
fn_out2 = [fn_path trial '_registered.hdr'];
copyfile(fn_out,fn_out2);

clear cube1 cube2
handles = guidata(MSIRegister);
set(handles.text17,'String','');
set(handles.text20,'String','');
set(handles.text18,'String','');
set(handles.text19,'String','');
set(handles.running,'Visible','off');
set(handles.running,'String','');
drawnow;
pause(0.1)
end

function [F1,F2,mxd1] = art_regx(img1,img2,msk,max_shift,samp_ratio,mxd,j01,j02,minpts,morepts)

% constants
alloc = 1024;   %pre-allocation increment

% template image
img1 = double(img1);
%img1 = img1/max(img1(:));
%low_high = stretchlim(img1,[0.01 0.99]);
%img1 = imadjust(img1,low_high,[0 1],1);

% extract features from subregions for the frames of the first set of images
F1 = zeros(alloc,3);
pnt = 1;
for j = j01:j02
    [F] = total_features(double(img1),j,samp_ratio,morepts);
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

% phase image
FFT21 = fft2(img1);
clear img1
p1 = atan2(imag(FFT21),real(FFT21));
clear FFT21
FFT21 = exp(1i*p1);
clear p1
img1 = real(ifft2(FFT21));
clear FFT21

%reference image
img2 = double(img2);
img2 = img2/max(img2(:));
low_high = stretchlim(img2,[0.01 0.99]);
img2 = imadjust(img2,low_high,[0 1],1);
h = fspecial('unsharp',0.8);
img2 = imfilter(img2,h);

% phase image
FFT21 = fft2(img2);
clear img2
p1 = atan2(imag(FFT21),real(FFT21));
clear FFT21
FFT21 = exp(1i*p1);
clear p1
img2 = real(ifft2(FFT21));
clear FFT21
[m2 n2] = size(img2);

% pair up feature points
F2 = zeros(mf1,3);
for i = mf1:-1:1
    if (msk(round(F1(i,2)),round(F1(i,1))) == 0)
        F1(i,:) = [];
        F2(i,:) = [];
    else
        j = F1(i,3);
        lt = floor(F1(i,1))-((2^j)-1);
        rt = ceil(F1(i,1))+((2^j)-1);
        tp = floor(F1(i,2))-((2^j)-1);
        bt = ceil(F1(i,2))+((2^j)-1);
        tmp = img1(tp:bt,lt:rt);
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

        tmp = img2(tp:bt,lt:rt);
        tmp2 = imresize(tmp,samp_ratio,'bicubic');
        clear tmp
        d = floor(samp_ratio/2)/samp_ratio;
        tmp2n = (lt-d):1/samp_ratio:(rt+d);
        tmp2m = (tp-d):1/samp_ratio:(bt+d);

        [mt1 nt1] = size(tmp1);
        [mt2 nt2] = size(tmp2);
        cc = abs(normxcorr2(tmp1,tmp2));
        sigma_cc = std2(cc);
        noise_floor = 4*sigma_cc;   %snf>2
        mskcc = imregionalmax(cc,8);
        cc = cc.*double(mskcc);
        cc = cc((floor(mt1/2)+1):end,(floor(nt1/2)+1):end);
        [max_cc, imax] = max(cc(:));
        [ypeak, xpeak] = ind2sub(size(cc),imax(1));
        cnt = sum(sum(cc>(max_cc-sigma_cc)));
        clear cc
        if ((xpeak <= nt2) && (ypeak <= mt2) && (cnt == 1) && (max_cc > noise_floor))
            F2(i,1) = tmp2n(xpeak);
            F2(i,2) = tmp2m(ypeak);
            F2(i,3) = max_cc;
        else
            F1(i,:) = [];
            F2(i,:) = [];
        end
        clear tmp1 tmp2 tmp2n tmp2m
    end
end
clear img1 img2 msk

% remove bad pairs
mxd1 = Inf;
[F1 F2 mxd1] = xy_disp(F1(:,1:2),F2,mxd,minpts);
end

function [F] = total_features(img,j,samp_ratio,morepts)
%% total_features

%[mfull nfull] = size(img);
%img = imresize(img,1/sc2,'bicubic');
[m0 n0] = size(img);

h = fspecial('unsharp',0.8);
img = imfilter(img,h);

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
hw = (2^j)+2;
F0(1:hw,:) = 0;
F0(:,1:hw) = 0;
F0((m0-hw+1):m0,:) = 0;
F0(:,(n0-hw+1):n0) = 0;

[y1,x1] = find(F0==1);
yd = length(y1);
for i = 1:yd
    if (morepts == 0)
        tmp = F0((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)).*M((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j));
        %tmp = F0((y1(i)-2^(j+1)):(y1(i)+2^(j+1)),(x1(i)-2^(j+1)):(x1(i)+2^(j+1))).*M((y1(i)-2^(j+1)):(y1(i)+2^(j+1)),(x1(i)-2^(j+1)):(x1(i)+2^(j+1)));
    elseif (morepts == 1)
        tmp = F0((y1(i)-2^(j-1)):(y1(i)+2^(j-1)),(x1(i)-2^(j-1)):(x1(i)+2^(j-1))).*M((y1(i)-2^(j-1)):(y1(i)+2^(j-1)),(x1(i)-2^(j-1)):(x1(i)+2^(j-1)));
    end
    
    [mx11,ind11] = max(tmp(:));
    if (mx11 > 0)
        xpeak = floor((ind11(1)-1)/size(tmp,1))+1;
        ypeak = ind11(1) - (xpeak-1)*size(tmp,1);
        msk11 = zeros(size(tmp));
        msk11(ypeak,xpeak) = 1;
        if (morepts == 0)
            F0((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)) = F0((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)).*msk11;
            %F0((y1(i)-2^(j+1)):(y1(i)+2^(j+1)),(x1(i)-2^(j+1)):(x1(i)+2^(j+1))) = F0((y1(i)-2^(j+1)):(y1(i)+2^(j+1)),(x1(i)-2^(j+1)):(x1(i)+2^(j+1))).*msk11;
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

function [F1,F2,mxd1,mf1] = xy_disp(F1,F2,mxd,minpts)
%% X/Y Disparity Minimumizing algorithm
% Usage: [F1 F2] = xy_disp(F1,F2,mxd,minpts)
%   inputs: F1,F2 -> two pairs of feature points
%           mxd -> maximum allowable disparity (pixels)
%           minpts -> mimumum number of final feature pairs
%   outputs: F1,F2 -> two pairs of optimal features points taken from the
%               original pairs
%

mf1 = size(F1,1);

xdist = F2(:,1) - F1(:,1);  %x-disperity
ydist = F2(:,2) - F1(:,2);  %y-disperity

mxx = Inf;  %maximum disperity
mxy = Inf;
while ((mxx>=mxd || mxy>=mxd) && mf1>minpts)
    % fit best planes to x and y disparities
    % compute error between the displarites and the planes
    A = [F1(:,1) F1(:,2) ones(mf1,1)];
    pA = pinv(A);
    clear A
    x = pA*xdist;
    errorx = xdist - (x(1)*F1(:,1) + x(2)*F1(:,2) + x(3));
    A = [F1(:,1) F1(:,2) ones(mf1,1)];
    pA = pinv(A);
    x = pA*ydist;
    clear pA
    errory = ydist - (x(1)*F1(:,1) + x(2)*F1(:,2) + x(3));
    
    % find the pair that contributes the largest error and remove it
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
    clear errorx errory mdist xdistribution ydistribution
    mf1 = size(F1,1);
end
mxd1 = max(mxx,mxy);
clear xdist ydist
end

function [z] = apply_mapping(I2,init_reg,X,gridsz,mfull,nfull)
%% apply_mapping
% Transform second image so that it is registered with the first
%   Input: I2 -> image
%          X -> transform coefficients
%          gridsz -> 
%   Outputs: z -> registered image

z = zeros(mfull,nfull);
[m0 n0] = size(I2);
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
        x0 = x0 - round(n0/2);
        y0 = y0 - round(m0/2);
        u2 = X(1)*x0 + X(2)*y0 + X(3)*x0.*y0 + X(4);
        v2 = X(5)*x0 + X(6)*y0 + X(7)*x0.*y0 + X(8);
        clear x0 y0
        u2 = u2(:) + round(n0/2);
        v2 = v2(:) + round(m0/2);
        msk = (u2>=1).*(u2<=nfull).*(v2>=1).*(v2<=mfull);      
        u2(msk==0) = [];
        v2(msk==0) = [];
        tmp = I2((m11-init_reg(1)+1):(m12-init_reg(1)+1),(n11-init_reg(2)+1):(n12-init_reg(2)+1));
        tmp(msk==0) = [];
        clear msk
        F = TriScatteredInterp(u2(:),v2(:),tmp(:),'natural');
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
        if (m11~=init_reg(1))
            msk(1:4,:) = 0;
        end
        if (m12~=mend)
            msk((end-4+1):end,:) = 0;
        end
        if (n11~=init_reg(2))
            msk(:,1:4) = 0;
        end
        if (n12~=nend)
            msk(:,(end-4+1):end) = 0;
        end
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
clear I2
end
