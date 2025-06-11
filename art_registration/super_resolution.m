function [] = super_resolution(cube_fn,m1,n1,bands1,precision1,offset1,interleave1,byteorder1,sc,samp_ratio,mxd,j01,j02,ppthresh,repcnt_max)
%% 
% super_resolution.m
% Author: Damon Conover
% Email: dconover@gwmail.gwu.edu
% Latest Revision: 3 July 2012
%
% Usage: 
%   inputs: 
%   outputs: 
%

%settings
max_shift = 12;  %maximum number of pixels that images can shift between band
minpts = 128;   %minimum number of final pairs (must be >= 8)
gpu = 0;    %no gpu->0, gpu->1

% constants
turnonoutput = 1;
filteron = 0;
morepts = 0;    %force maximum number of feature points to be extracted -> 1, otherwise -> 0

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
pc = 1;
if (isempty(tmp))
    pc = 0;
    tmp = regexp(cube_fn,'\/');
end
tmp1 = tmp(end);
tmp2 = length(cube_fn);
fn_path = cube_fn(1:tmp1);
trial = cube_fn((tmp1+1):tmp2);
clear tmp tmp1 tmp2

% feature mask
fn_out = [fn_path trial '_mask.tif'];
if (exist(fn_out,'file') == 2)
    msk = logical(imread(fn_out)>0);
else
    msk = true(m1,n1);
end

% cube
cube1 = multibandread(cube_fn,[m1 n1 bands1],precision1,offset1,interleave1,byteorder1,{'Band','Range',[1 1 bands1]});
cube1(cube1==Inf) = 1;
cube1(cube1==-Inf) = -1;
cube1(isnan(cube1)) = 0;

Xstar = zeros(bands1,8);
Xstar(1,:) = [1 0 0 0 1 0 0 0];
for i = 2:bands1
    repcnt = 0;
    while (repcnt <= repcnt_max)
        %ppthresh0 = ppthresh*(2^repcnt);
        tmp1 = cube1(:,:,i);
        tmp2 = cube1(:,:,1);
        if (filteron == 1)
            tmp1 = medfilt2(tmp1,[3 3]);
            tmp2 = medfilt2(tmp2,[3 3]);
        end
        [Xstar(i,:) mxd1 mf1] = art_regx(tmp1,tmp2,msk,max_shift,samp_ratio,mxd,ppthresh,j01,j02,minpts,gpu,morepts);
        clear tmp1 tmp2
        if (turnonoutput==1)
            disp(['Band: ', num2str(i)]);
            disp(['Final disparity: ', num2str(mxd1)]);
            disp(['Number of feature point pairs: ', num2str(mf1)]);
            if (mxd1 <= mxd)
                disp('PASS');
            else
                disp('FAIL');
            end
        end
        if (mxd1 <= mxd)
            repcnt = repcnt_max + 1;
        else
            repcnt = repcnt + 1;
        end
    end
end

[x0,y0] = meshgrid(1:n1,1:m1);
F20 = zeros(m1*n1*bands1,4);
F20(:,3) = cube1(:);
clear cube1
F20(1:(m1*n1),1) = x0(:);
F20(1:(m1*n1),2) = y0(:);
F20(1:(m1*n1),4) = ones(m1*n1,1);
for i = 2:bands1
    u2 = Xstar(i,1)*x0 + Xstar(i,2)*y0 + Xstar(i,3)*x0.*y0 + Xstar(i,4);
    v2 = Xstar(i,5)*x0 + Xstar(i,6)*y0 + Xstar(i,7)*x0.*y0 + Xstar(i,8);
    F20(((i-1)*m1*n1+1):(i*m1*n1),1) = u2(:);
    F20(((i-1)*m1*n1+1):(i*m1*n1),2) = v2(:);
    F20(((i-1)*m1*n1+1):(i*m1*n1),4) = i*ones(m1*n1,1);
    clear u2 v2
end
clear x0 y0
%{
s=RandStream('mt19937ar');
RandStream.setDefaultStream(s);
rnum = rand(s,1,2) - 0.5;
F20(:,1) = F20(:,1) + rnum(1);
F20(:,2) = F20(:,2) + rnum(2);
%}
delta0 = zeros(1,2);
mx = 0;
for deltax = (-1/samp_ratio+0.01):0.01:(1/samp_ratio-0.01)
    for deltay = (-1/samp_ratio+0.01):0.01:(1/samp_ratio-0.01)
        tmp1 = F20(:,1) - deltax;
        tmp2 = F20(:,2) - deltay;
        test = sum(abs(sqrt((tmp1 - round(tmp1*samp_ratio)/samp_ratio).^2 + (tmp2 - round(tmp2*samp_ratio)/samp_ratio).^2)));
        clear tmp1 tmp2
        if (test > mx)
            delta0(1) = deltax;
            delta0(2) = deltay;
        end
    end
end
F20(:,1) = F20(:,1) - delta0(1);
F20(:,2) = F20(:,2) - delta0(2);

regsz = round(256/samp_ratio);
z = zeros(samp_ratio*m1,samp_ratio*n1);
mend = m1;
nend = n1;
stopm = 0;
m11 = 1;
m12 = m11 + regsz - 1;
if (m12 >= mend)
    m12 = mend;
    stopm = stopm + 1;
end
while (stopm < 2)
    stopn = 0;
    n11 = 1;
    n12 = n11 + regsz - 1;
    if (n12 >= nend)
        n12 = nend;
        stopn = stopn + 1;
    end
    while (stopn < 2)
        Ftmp = F20;
        Ftmp(Ftmp(:,1)<(n11-8),:) = [];
        Ftmp(Ftmp(:,2)<(m11-8),:) = [];
        Ftmp(Ftmp(:,1)>(n12+8),:) = [];
        Ftmp(Ftmp(:,2)>(m12+8),:) = [];
        mug1 = mean(Ftmp(Ftmp(:,4)==1,3));
        for i = 2:bands1
            if (sum(Xstar(i,:)==zeros(1,8)) == 8)
                Ftmp(Ftmp(:,4)==i,:) = [];
            else
                mugi = mean(Ftmp(Ftmp(:,4)==i,3));
                Ftmp(Ftmp(:,4)==i,3) = Ftmp(Ftmp(:,4)==i,3)*mug1/mugi;
            end
        end
        F = TriScatteredInterp(Ftmp(:,1),Ftmp(:,2),Ftmp(:,3),'natural');
        clear Ftmp
        [u1,v1] = meshgrid((n11-floor(samp_ratio/2)/samp_ratio):1/samp_ratio:(n12+floor(samp_ratio/2)/samp_ratio),(m11-floor(samp_ratio/2)/samp_ratio):1/samp_ratio:(m12+floor(samp_ratio/2)/samp_ratio));
        z0 = F(u1,v1);
        z0(isnan(z0)) = 0;
        clear F
        
        msk = (u1>=(1-floor(samp_ratio/2)/samp_ratio)).*(u1<=(n1+floor(samp_ratio/2)/samp_ratio)).*(v1>=(1-floor(samp_ratio/2)/samp_ratio)).*(v1<=(m1+floor(samp_ratio/2)/samp_ratio));      
        if (m11~=1)
            msk(1:2*samp_ratio,:) = 0;
        end
        if (m12~=mend)
            msk(((end-2*samp_ratio)+1):end,:) = 0;
        end
        if (n11~=1)
            msk(:,1:2*samp_ratio) = 0;
        end
        if (n12~=nend)
            msk(:,((end-2*samp_ratio)+1):end) = 0;
        end
        z(((m11-1)*samp_ratio+1):m12*samp_ratio,((n11-1)*samp_ratio+1):n12*samp_ratio) = z(((m11-1)*samp_ratio+1):m12*samp_ratio,((n11-1)*samp_ratio+1):n12*samp_ratio).*(~msk) + z0.*msk;
        clear z0 u1 v1 msk

        n11 = n12 - round(regsz/8) + 1;
        n12 = n11 + regsz - 1;
        if (n12 >= nend)
            n12 = nend;
            stopn = stopn + 1;
        end
    end
    m11 = m12 - round(regsz/8) + 1;
    m12 = m11 + regsz - 1;
    if (m12 >= mend)
        m12 = mend;
        stopm = stopm + 1;
    end
end
clear F20 Xstar
z(isnan(z)) = 0;

z = double(z);
z = imresize(z,sc/samp_ratio,'bicubic');
z = (z-min(z(:)))/(max(z(:))-min(z(:)));
low_high = stretchlim(z,[0.01 0.99]);
z = imadjust(z,low_high,[0 1],1);
fn_out = [trial '_super.tif'];
imwrite(uint16((2^16-1)*z),fn_out,'Compression','None')
clear z
end

function [X,mxd1,mf1] = art_regx(img1,img2,msk,max_shift,samp_ratio,mxd,pptthresh,j01,j02,minpts,gpu,morepts)

% constants
alloc = 1024;   %pre-allocation increment

% template image
img1 = double(img1);
img1 = img1/max(img1(:));
low_high = stretchlim(img1,[0.01 0.99]);
img1 = imadjust(img1,low_high,[0 1],1);
h = fspecial('unsharp',0.8);
img1 = imfilter(img1,h);

% extract features from subregions for the frames of the first set of images
F1 = zeros(alloc,3);
pnt = 1;
for j = j01:j02
    [F] = total_features(double(img1),j,pptthresh,samp_ratio,morepts);
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
if (gpu == 1)
    img1 = gpuArray(double(img1));
    FFT21 = fft2(img1);
    clear img1
    p1 = atan2(imag(FFT21),real(FFT21));
    clear FFT21
    FFT21 = exp(1i*p1);
    clear p1
    img1 = real(ifft2(FFT21));
    clear FFT21
    img1 = gather(img1);
else
    FFT21 = fft2(img1);
    clear img1
    p1 = atan2(imag(FFT21),real(FFT21));
    clear FFT21
    FFT21 = exp(1i*p1);
    clear p1
    img1 = real(ifft2(FFT21));
    clear FFT21
end

%reference image
img2 = double(img2);
img2 = img2/max(img2(:));
low_high = stretchlim(img2,[0.01 0.99]);
img2 = imadjust(img2,low_high,[0 1],1);
h = fspecial('unsharp',0.8);
img2 = imfilter(img2,h);
% phase image
if (gpu == 1)
    img2 = gpuArray(double(img2));
    FFT21 = fft2(img2);
    clear img2
    p1 = atan2(imag(FFT21),real(FFT21));
    clear FFT21
    FFT21 = exp(1i*p1);
    clear p1
    img2 = real(ifft2(FFT21));
    clear FFT21
    img2 = gather(img2);
else
    FFT21 = fft2(img2);
    clear img2
    p1 = atan2(imag(FFT21),real(FFT21));
    clear FFT21
    FFT21 = exp(1i*p1);
    clear p1
    img2 = real(ifft2(FFT21));
    clear FFT21
end
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
mf1 = size(F1,1);
clear img1 img2 msk

% remove bad pairs
mxd1 = Inf;
[F1 F2 mxd1 mf1] = xy_disp(F1(:,1:2),F2,mxd,minpts);

X = zeros(1,8);
if (mxd1 <= mxd)
    % transform image 1 to match image 2
    [X] = mapping(F2,F1,gpu);
end
end

function [F] = total_features(img,j,pptthresh,samp_ratio,morepts)
%% total_features

[m0 n0] = size(img);

% Wavelet transform
h = fspecial('average',[1 2^j]);    %hsise = 1 x 2^j
LPx = imfilter(img,h);

h = [];
for q = 1:j
    h = [-1 h];
    h = [h 1];
end
LH = imfilter(LPx,h');
clear LPx;
HPy = imfilter(img,h);
h = fspecial('average',[2^j 1]);    %hsize = 2^j x 1
HL = imfilter(HPy,h);
clear HPy;

% Modulus of the wavelet transform
M = sqrt(LH.^2+HL.^2);
clear LH;
clear HL;

F0 = imregionalmax(M,8);

% Remove border pixels
hw = (2^j);
F0(1:hw,:) = 0;
F0(:,1:hw) = 0;
F0((m0-hw+1):m0,:) = 0;
F0(:,(n0-hw+1):n0) = 0;
M(1:hw,:) = 0;
M(:,1:hw) = 0;
M((m0-hw+1):m0,:) = 0;
M(:,(n0-hw+1):n0) = 0;

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

% subsample and find subpixel maxima
pnt = 1;
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
clear M0 M
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

function [X] = mapping(F11,F21,gpu)
%% mapping
% Transform second image so that it is registered with the first
%   Input: F11 -> control points for image 1 from art_controlpts.m
%          F21 -> control points for image 2 from art_controlpts.m
%   Outputs: z -> registered image

m2 = size(F21,1);
if (gpu == 1)
    F11 = gpuArray(F11);
    F21 = gpuArray(F21);
    A = parallel.gpu.GPUArray.zeros(2*m2,8,'double');
else
    A = zeros(2*m2,8);
end
% bilinear transformation
Y = [F11(:,1); F11(:,2)];
x = F21(:,1);
y = F21(:,2);
A(1:m2,:) = [x y x.*y ones(m2,1) zeros(m2,4)];
A((m2+1):2*m2,:) = [zeros(m2,4) x y x.*y ones(m2,1)];
% pseudoinverse solution
if (gpu == 1)
    X = (A'*A)\(A')*Y;
    X = gather(X);
else
    X = pinv(double(A))*double(Y);
end
clear A Y F11 F21
end
