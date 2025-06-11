% Painting stitching tool
% mosaic_full.m
% Author: Damon Conover
% Email: dconover@gwu.edu
% Latest Revision: 25 February 2013
%
% Usage: [mosaic] = mosaic_full(fn_path,fn_out,mmpixel,feather_width,ncores,turnonoutput,quick)
%   inputs: fn_path -> path to easel CSV file
%           fn_out -> filename of resulting image (without extension)
%           mmpixel -> mm/pixel
%           feather_width -> how many pixels across the stitch to feather
%           ncores -> number of CPU cores
%           turnonoutput -> write image to disk (=1), otherwise(=0)
%           quick -> quick method (=1), otherwise (=0)
%   outputs: mosaic -> stitche image

function [] = mosaic_full(fn_path,fn_out,mmpixel,feather_width,ncores,turnonoutput,quick,usegui)

if (quick == 1)
    mosiac_quick(fn_path,fn_out,mmpixel,turnonoutput,usegui);
else

isOpen = matlabpool('size') > 0;
if (isOpen == 1)
    matlabpool close
end
ncores0 = feature('NumCores');
if (ncores > ncores0)
    ncores = ncores0;
end
matlabpool(ncores)
    
% settings
max_shift = 24;  %maximum number of pixels that images can shift between blocks
samp_ratio = 3; %upsample ratio for determining fractional pixels (must be an odd integer)
j0 = 4;
bord = 256;     %pixels
Ntran = 0;

if (usegui == 1)
    handles = guidata(PaintingStitch);
    set(handles.regtime,'String','0');
    set(handles.trantime,'String','0');
    set(handles.feathertime,'String','0');
    set(handles.ncores,'String',num2str(ncores));
    drawnow;
end
        
% Get all file names
everyImage = dir(fullfile(fn_path, '*.tif'));
offsetFile = dir(fullfile(fn_path, '*.csv'));

offset_val = csvread(fullfile(fn_path,offsetFile(1).name));
offset_val = offset_val(:,4:5);     %% x-shift(mm), y-shift(mm)
offset_val = offset_val/mmpixel;    %% x-shift(pixel), y-shift(pixel)
[num ~] = size(offset_val);

% Get image dimensions, icc profile, and bit depth
info = imfinfo(fullfile(fn_path, everyImage(1).name));
m0 = info.Height;
n0 = info.Width;
bd = info.BitsPerSample(1);
iccProfile = iccread(fullfile(fn_path, everyImage(1).name));
embeddedICC = iccwrite(iccProfile, fullfile(fn_path,'embeddedICCprofile.icc'));

offset_val = offset_val - repmat(min(offset_val)-1,num,1);
[~,colind] = max(abs(offset_val(:,1)));
colmx = round(offset_val(colind,1));
[~,rowind] = max(abs(offset_val(:,2)));
rowmx = round(offset_val(rowind,2));
offset_val(:,1) = -(offset_val(:,1) - 1) + colmx;
mtot = abs(rowmx)+m0-1;
ntot = abs(colmx)+n0-1;

tic
% identify and match points
if (bd == 8)
    imgAll = uint8(zeros(m0,n0,3,num));
elseif (bd == 16)
    imgAll = uint16(zeros(m0,n0,3,num));
end
pt_pair_set = zeros(1024,4);
ov_tot = zeros(num,1);
cnt = 1;
pnt1 = 1;
fn11 = [];
fn21 = [];
imgAll(:,:,:,1) = imread((fullfile(fn_path,everyImage(1).name)),'tif');
for i = 1:(num-1)
    msk1 = false(mtot+2*bord,ntot+2*bord);
    img1 = single(rgb2gray(imgAll(:,:,:,i)));
    [m0 n0] = size(img1);
    m11 = round(offset_val(pnt1,2)+bord);
    m21 = round(m11 + m0 - 1);
    n11 = round(offset_val(pnt1,1)+bord);
    n21 = round(n11 + n0 - 1);
    msk1(m11:m21,n11:n21) = 1;
    pnt2 = pnt1 + 1;
    for j = (i+1):num
        msk2 = false(mtot+2*bord,ntot+2*bord);
        if (i == 1)
            imgAll(:,:,:,j) = imread((fullfile(fn_path,everyImage(j).name)),'tif');
        end
        img2 = single(rgb2gray(imgAll(:,:,:,j)));
        [m0 n0] = size(img2);
        m12 = round(offset_val(pnt2,2)+bord);
        m22 = round(m12 + m0 - 1);
        n12 = round(offset_val(pnt2,1)+bord);
        n22 = round(n12 + n0 - 1);
        msk2(m12:m22,n12:n22) = 1;
        msk = and(msk1,msk2);
        tmp = msk(m12:m22,n12:n22);
        tmp(1,:) = 0;
        tmp(m0,:) = 0;
        tmp(:,1) = 0;
        tmp(:,n0) = 0;
        se = strel('square',(2*(2^j0)+1));
        msk21 = imerode(tmp,se);
        clear tmp
        ov = sum(msk21(:));
        if (ov > 0)
            fn11{cnt} = everyImage(i).name;
            fn21{cnt} = everyImage(j).name;
            ov_tot(pnt1) = ov_tot(pnt1) + ov;
            ov_tot(pnt2) = ov_tot(pnt2) + ov;
            pt_pair_set(cnt,1) = i;
            pt_pair_set(cnt,2) = j;
            pt_pair_set(cnt,3) = pnt1;
            pt_pair_set(cnt,4) = pnt2;
            cnt = cnt + 1;
        end
        pnt2 = pnt2 + 1;
        clear img2 msk msk2 msk21
    end
    pnt1 = pnt1 + 1;
    clear img1 msk1
end
pt_pair_set(cnt:end,:) = [];
[cnt,~] = size(pt_pair_set);
[~,ind] = max(ov_tot);

parfor k = 1:cnt
    %i = pt_pair_set(k,1);
    %j = pt_pair_set(k,2);
    pnt1 = pt_pair_set(k,3);
    pnt2 = pt_pair_set(k,4);
    msk1 = false(mtot+2*bord,ntot+2*bord);
    img1 = imread((fullfile(fn_path,fn11{k})),'tif');
    %img1 = imgAll(:,:,:,i);
    img1 = single(rgb2gray(img1));
    [m0 n0] = size(img1);
    m11 = round(offset_val(pnt1,2)+bord);
    m21 = round(m11 + m0 - 1);
    n11 = round(offset_val(pnt1,1)+bord);
    n21 = round(n11 + n0 - 1);
    msk1(m11:m21,n11:n21) = 1;
    
    msk2 = false(mtot+2*bord,ntot+2*bord);
    img2 = imread((fullfile(fn_path,fn21{k})),'tif');
    %img2 = imgAll(:,:,:,j);
    img2 = single(rgb2gray(img2));
    [m0 n0] = size(img2);
    m12 = round(offset_val(pnt2,2)+bord);
    m22 = round(m12 + m0 - 1);
    n12 = round(offset_val(pnt2,1)+bord);
    n22 = round(n12 + n0 - 1);
    msk2(m12:m22,n12:n22) = 1;
    msk = and(msk1,msk2);
    tmp = msk(m12:m22,n12:n22);
    tmp(1,:) = 0;
    tmp(m0,:) = 0;
    tmp(:,1) = 0;
    tmp(:,n0) = 0;
    se = strel('square',(2*(2^j0)+1));
    msk21 = imerode(tmp,se);
                 
    F2 = total_features(img2,msk21,j0,samp_ratio);
    mt = size(F2,1);
    F1 = F2 + repmat([(n12-n11),(m12-m11),0],mt,1);
    [F1,F2] = art_regx(img1,img2,F1,F2,samp_ratio,max_shift);
    F1(:,1) = F1(:,1) - round(n0/2);
    F2(:,1) = F2(:,1) - round(n0/2);
    F1(:,2) = F1(:,2) - round(m0/2);
    F2(:,2) = F2(:,2) - round(m0/2);
    mt = size(F2,1);
    
    pt_pairs0 = zeros(mt,7);
    pt_pairs0(:,1) = pnt1;
    pt_pairs0(:,2) = pnt2;
    pt_pairs0(:,3:4) = F1(:,1:2);
    pt_pairs0(:,5:7) = F2;
    pt_pairs1{k} = pt_pairs0;
end
clear tmp F2 F1 img2 msk msk2 msk21 img1 msk1 pt_pairs0

mf10 = zeros(cnt,1);
parfor k = 1:cnt
    mf10(k) = size(pt_pairs1{k},1);
end
mf1 = sum(mf10);
pt_pairs = zeros(mf1,7);
pnt = 1;
for k = 1:cnt
    pt_pairs(pnt:(pnt+mf10(k)-1),:) = pt_pairs1{k};
    pnt = pnt + mf10(k);
end
clear pt_pairs1

% compute coefficients
if (Ntran == 0)
    num_coeff = 4;
else
    num_coeff = sum((1:(Ntran+1)));
end
A = zeros(2*mf1,2*num_coeff*(num-1));
Y = zeros(2*mf1,1);
pnt1 = 1;
for i = 1:(num-1)
    for j = (i+1):num
        Rx = pt_pairs((pt_pairs(:,1)==i) & (pt_pairs(:,2)==j),3);
        Ry = pt_pairs((pt_pairs(:,1)==i) & (pt_pairs(:,2)==j),4);
        Tx = pt_pairs((pt_pairs(:,1)==i) & (pt_pairs(:,2)==j),5);
        Ty = pt_pairs((pt_pairs(:,1)==i) & (pt_pairs(:,2)==j),6);
        npt = size(Rx,1);
        if (Ntran == 0)
            if (i == ind)
                Y(pnt1:(pnt1+npt-1),1) = Rx;
                Y((mf1+pnt1):(mf1+pnt1+npt-1),1) = Ry;
            elseif (i < ind)
                A(pnt1:(pnt1+npt-1),(i-1)*num_coeff+1) = -Rx;
                A(pnt1:(pnt1+npt-1),(i-1)*num_coeff+2) = -Ry;
                A(pnt1:(pnt1+npt-1),(i-1)*num_coeff+3) = -Rx.*Ry;
                A(pnt1:(pnt1+npt-1),(i-1)*num_coeff+4) = -ones(npt,1);
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(i-1)*num_coeff+1) = -Rx;
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(i-1)*num_coeff+2) = -Ry;
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(i-1)*num_coeff+3) = -Rx.*Ry;
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(i-1)*num_coeff+4) = -ones(npt,1);
            elseif (i > ind)
                A(pnt1:(pnt1+npt-1),(i-2)*num_coeff+1) = -Rx;
                A(pnt1:(pnt1+npt-1),(i-2)*num_coeff+2) = -Ry;
                A(pnt1:(pnt1+npt-1),(i-2)*num_coeff+3) = -Rx.*Ry;
                A(pnt1:(pnt1+npt-1),(i-2)*num_coeff+4) = -ones(npt,1);
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(i-2)*num_coeff+1) = -Rx;
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(i-2)*num_coeff+2) = -Ry;
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(i-2)*num_coeff+3) = -Rx.*Ry;
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(i-2)*num_coeff+4) = -ones(npt,1);
            end
            if (j == ind)
                Y(pnt1:(pnt1+npt-1),1) = -Tx;
                Y((mf1+pnt1):(mf1+pnt1+npt-1),1) = -Ty;
            elseif (j < ind)
                A(pnt1:(pnt1+npt-1),(j-1)*num_coeff+1) = Tx;
                A(pnt1:(pnt1+npt-1),(j-1)*num_coeff+2) = Ty;
                A(pnt1:(pnt1+npt-1),(j-1)*num_coeff+3) = Tx.*Ty;
                A(pnt1:(pnt1+npt-1),(j-1)*num_coeff+4) = ones(npt,1);
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(j-1)*num_coeff+1) = Tx;
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(j-1)*num_coeff+2) = Ty;
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(j-1)*num_coeff+3) = Tx.*Ty;
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(j-1)*num_coeff+4) = ones(npt,1);
            elseif (j > ind)
                A(pnt1:(pnt1+npt-1),(j-2)*num_coeff+1) = Tx;
                A(pnt1:(pnt1+npt-1),(j-2)*num_coeff+2) = Ty;
                A(pnt1:(pnt1+npt-1),(j-2)*num_coeff+3) = Tx.*Ty;
                A(pnt1:(pnt1+npt-1),(j-2)*num_coeff+4) = ones(npt,1);
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(j-2)*num_coeff+1) = Tx;
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(j-2)*num_coeff+2) = Ty;
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(j-2)*num_coeff+3) = Tx.*Ty;
                A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(j-2)*num_coeff+4) = ones(npt,1);
            end
        else
            pnt = 1;
            for k = 0:Ntran
                for q = 0:(Ntran-k)
                    if (i == ind)
                        Y(pnt1:(pnt1+npt-1),1) = Rx;
                        Y((mf1+pnt1):(mf1+pnt1+npt-1),1) = Ry;
                    elseif (i < ind)
                        A(pnt1:(pnt1+npt-1),(i-1)*num_coeff+pnt) = -(Rx.^k).*(Ry.^q);
                        A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(i-1)*num_coeff+pnt) = -(Rx.^k).*(Ry.^q);
                    elseif (i > ind)
                        A(pnt1:(pnt1+npt-1),(i-2)*num_coeff+pnt) = -(Rx.^k).*(Ry.^q);
                        A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(i-2)*num_coeff+pnt) = -(Rx.^k).*(Ry.^q);
                    end
                    if (j == ind)
                        Y(pnt1:(pnt1+npt-1),1) = -Tx;
                        Y((mf1+pnt1):(mf1+pnt1+npt-1),1) = -Ty;
                    elseif (j < ind)
                        A(pnt1:(pnt1+npt-1),(j-1)*num_coeff+pnt) = (Tx.^k).*(Ty.^q);
                        A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(j-1)*num_coeff+pnt) = (Tx.^k).*(Ty.^q);
                    elseif (j > ind)
                        A(pnt1:(pnt1+npt-1),(j-2)*num_coeff+pnt) = (Tx.^k).*(Ty.^q);
                        A((mf1+pnt1):(mf1+pnt1+npt-1),num_coeff*(num-1)+(j-2)*num_coeff+pnt) = (Tx.^k).*(Ty.^q);
                    end
                    pnt = pnt + 1;
                end
            end
        end
        pnt1 = pnt1 + npt;
    end 
    clear Rx Ry Tx Ty
end
X = pinv(double(A))*double(Y);
X = reshape(X,num_coeff,2*(num-1))';
clear A Y pt_pairs
regtime = toc;
if (usegui == 1)
    handles = guidata(PaintingStitch);
    set(handles.regtime,'String',num2str(regtime));
    drawnow;
end

trantime = 0;
feathertime = 0;
% apply transformation
if (bd == 8)
    mosaic = uint8(zeros(mtot+2*bord,ntot+2*bord,3));
elseif (bd == 16)
    mosaic = uint16(zeros(mtot+2*bord,ntot+2*bord,3));
end
%img1 = imread((fullfile(fn_path,everyImage(ind).name)),'tif');
img1 = imgAll(:,:,:,ind);
[m0,n0,~] = size(img1);
m1 = round(offset_val(ind,2)+bord);
m2 = round(m1 + m0 - 1);
n1 = round(offset_val(ind,1)+bord);
n2 = round(n1 + n0 - 1);
[x0,y0] = meshgrid(1:n0,1:m0);
x0 = x0 - round(n0/2);
y0 = y0 - round(m0/2);
mosaic(m1:m2,n1:n2,:) = img1;
if (usegui == 1)
    handles = guidata(PaintingStitch);
    set(handles.axes1,'Visible','on');
    axes(handles.axes1)
    imshow(mosaic(1:4:end,1:4:end,:));
    drawnow;
end
        
x02 = n0 - round(n0/2);
y02 = m0 - round(m0/2);
tmp = single(sqrt(x0.^2 + y0.^2));
dist_from_center = single((sqrt(x02^2 + y02^2)+1000)*ones(mtot+2*bord,ntot+2*bord));
dist_from_center(m1:m2,n1:n2) = tmp;
clear x0 y0 tmp
clear img1
for i = 1:num
    if (i ~= ind)
        %img1 = imread((fullfile(fn_path,everyImage(i).name)),'tif');
        img1 = imgAll(:,:,:,i);
        if (i < ind)
            X1 = X(i,:);
            X2 = X((num-1+i),:);
        elseif (i > ind)
            X1 = X((i-1),:);
            X2 = X((num-1+i-1),:);
        end
        X0 = [X1 X2];
        
        vall = 1:n0;
        vall = vall - round(n0/2);
        uall = 1:m0;
        uall = uall - round(m0/2);
        x01 = 1 - round(n0/2);
        x02 = n0 - round(n0/2);
        y01 = 1 - round(m0/2);
        y02 = m0 - round(m0/2);
        bc = zeros(2*m0+2*n0,2);
        bc(1:m0,:) = [x01*ones(m0,1) uall'];
        bc((m0+1):(2*m0),:) = [x02*ones(m0,1) uall'];
        bc((2*m0+1):(2*m0+n0),:) = [vall' y01*ones(n0,1)];
        bc((2*m0+n0+1):(2*m0+2*n0),:) = [vall' y02*ones(n0,1)];
        if (Ntran == 0)
            v0 = X0(1)*bc(:,1) + X0(2)*bc(:,2) + X0(3)*bc(:,1).*bc(:,2) + X0(4)*ones(2*m0+2*n0,1);
            u0 = X0(5)*bc(:,1) + X0(6)*bc(:,2) + X0(7)*bc(:,1).*bc(:,2) + X0(8)*ones(2*m0+2*n0,1);
        else
            u0 = zeros(2*m0+2*n0,1);
            v0 = zeros(2*m0+2*n0,1);
            pnt = 1;
            for k = 0:Ntran
                for j = 0:(Ntran-k)
                    v0 = v0 + X0(pnt)*((bc(:,1).^k).*(bc(:,2).^j));
                    u0 = u0 + X0(pnt+sum((1:(Ntran+1))))*((bc(:,1).^k).*(bc(:,2).^j));
                    pnt = pnt + 1;
                end
            end
        end
        clear bc uall vall
        u0 = u0 + round(m0/2);
        v0 = v0 + round(n0/2);
        n211 = ceil(min(v0));
        n221 = floor(max(v0));
        m211 = ceil(min(u0));
        m221 = floor(max(u0));
        %mosaic0 = single(zeros((m221-m211+1),(n221-n211+1),3));
        %dist_from_center0 = single(zeros((m221-m211+1),(n221-n211+1),3));
        mosaic00 = single(zeros((m221-m211+1),(n221-n211+1),12));
        dist_from_center00 = single(zeros((m221-m211+1),(n221-n211+1),12));
        clear u0 v0
        tmp = [1 1 round((m0-bord)/2) round((m0-bord)/2)];
        m01 = [tmp tmp tmp];
        tmp = [round((m0+bord)/2) round((m0+bord)/2) m0 m0];
        m02 = [tmp tmp tmp];
        tmp = [1 round((n0-bord)/2) 1 round((n0-bord)/2)];
        n01 = [tmp tmp tmp];
        tmp = [round((n0+bord)/2) n0 round((n0+bord)/2) n0];
        n02 = [tmp tmp tmp];
        parfor k = 1:12
            [mosaic00(:,:,k),dist_from_center00(:,:,k)] = apply_mapping(img1(:,:,ceil(k/4)),m01(k),m02(k),n01(k),n02(k),X0,Ntran);
        end
        clear img1
        mosaic0 = single(zeros((m221-m211+1),(n221-n211+1),3));
        dist_from_center0 = single(zeros((m221-m211+1),(n221-n211+1),3));
        pnt = 0;
        for k = 1:3
            for q = 1:4
                pnt = pnt + 1;
                msk = (mosaic00(:,:,pnt) > 0);
                mosaic0(:,:,k) = mosaic0(:,:,k).*(~msk) + mosaic00(:,:,pnt).*msk;
                dist_from_center0(:,:,k) = dist_from_center0(:,:,k).*(~msk) + dist_from_center00(:,:,q).*msk;
            end
        end
        clear mosaic00 dist_from_center00
        m21 = m211 + m1 - 1;
        m22 = m221 + m1 - 1;
        n21 = n211 + n1 - 1;
        n22 = n221 + n1 - 1;
        
        % compare distance from center with mosaic distances
        msk1 = (mosaic0(:,:,1) ~= 0);
        msk2 = (mosaic(m21:m22,n21:n22,1) ~= 0);
        mskd = and((dist_from_center0(:,:,1) < dist_from_center(m21:m22,n21:n22)),msk1);
        dist_from_center(m21:m22,n21:n22) = dist_from_center(m21:m22,n21:n22).*(~mskd) + dist_from_center0(:,:,1).*mskd;
        %{
        imagesc(dist_from_center)
        axis image
        pause(5)
        %}
        clear dist_from_center0
        trantime1 = toc;
        trantime = trantime + trantime1 - regtime;
        if (usegui == 1)
            handles = guidata(PaintingStitch);
            set(handles.trantime,'String',num2str(trantime));
            drawnow;
        end
        %{
        % mosaic
        se = strel('disk',round(feather_width/2));
        new_msk1 = and(and((imdilate(mskd,se) - mskd),msk1),msk2);   %feather region in overall mosaic
        new_msk2 = and(and((imdilate(~mskd,se) - (~mskd)),msk1),msk2);  %feather region in new data
        clear msk1 msk2
        new_msk2(and(new_msk1,new_msk2)) = 0;   %if the masks ovelap, use mosaic info
        new_msk = or(new_msk1,new_msk2);
        
        se = strel('square',3);
        C1 = imdilate(new_msk,se) - new_msk;
        A1 = and(C1,imdilate(new_msk1,se));
        [A1y,A1x] = find(A1==1);    %feather boundary in overall mosaic
        B1 = and(C1,imdilate(new_msk2,se));
        [B1y,B1x] = find(B1==1);    %feather boundary in new data
        clear C1
        clear new_msk1 new_msk2
        [Cy,Cx] = find(new_msk==1);
        %dA = dsearchn1([A1y A1x],[Cy Cx]);
        %dB = dsearchn1([B1y B1x],[Cy Cx]);
        dA = dsearchn1([A1y A1x],[Cy Cx]);
        dB = dsearchn1([B1y B1x],[Cy Cx]);
        test = ((dA > 2*feather_width) | (dB > 2*feather_width));
        Cy(test) = [];
        Cx(test) = [];
        dA(test) = [];
        dB(test) = [];
        clear test
        feather_regionA = single(zeros((m221-m211+1),(n221-n211+1)));
        feather_regionB = single(zeros((m221-m211+1),(n221-n211+1)));
        if (~isempty([A1y A1x]) && ~isempty([B1y B1x]) && ~isempty([Cy Cx]))
            tmp1 = dA./(dA + dB);
            tmp2 = dB./(dA + dB);
            for p = 1:length(Cy)
                feather_regionA(Cy(p),Cx(p)) = tmp1(p);
                feather_regionB(Cy(p),Cx(p)) = tmp2(p);
            end
            clear tmp
        end
        clear kA kB dA dB
        clear Cx Cy new_msk A1y A1x B1y B1x
        new_msk = or((feather_regionA > 0),(feather_regionB > 0));
        %}
        mosaic_full0 = single(mosaic(m21:m22,n21:n22,:));
        %{
        parfor k = 1:3
            mosaic_full0(:,:,k) = (mosaic_full0(:,:,k).*(~mskd) + mosaic0(:,:,k).*mskd).*(~new_msk) + mosaic_full0(:,:,k).*feather_regionB + mosaic0(:,:,k).*feather_regionA;
        end
        %}
        
        
        parfor k = 1:3
            mosaic_full0(:,:,k) = (mosaic_full0(:,:,k).*(~mskd) + mosaic0(:,:,k).*mskd);
        end
        
        
        
        
        clear feather_regionA feather_regionB mskd new_msk mosaic0
        if (bd == 8)
            mosaic(m21:m22,n21:n22,:) = uint8(mosaic_full0);
        elseif (bd == 16)
            mosaic(m21:m22,n21:n22,:) = uint16(mosaic_full0);
        end
        clear mosaic_full0
        
        if (usegui == 1)
            handles = guidata(PaintingStitch);
            set(handles.axes1,'Visible','on');
            axes(handles.axes1)
            imshow(mosaic(1:4:end,1:4:end,:));
            drawnow;
        end
        feathertime1 = toc;
        feathertime = feathertime + feathertime1 - trantime1;
        regtime = feathertime1;
        %{
        if (usegui == 1)
            handles = guidata(PaintingStitch);
            set(handles.feathertime,'String',num2str(feathertime));
            drawnow;
        end
        %}
    end
end
matlabpool close
mosaic(1:bord,:,:) = [];
mosaic((mtot+1):end,:,:) = [];
mosaic(:,1:bord,:) = [];
mosaic(:,(ntot+1):end,:) = [];
if (usegui == 1)
    handles = guidata(PaintingStitch);
    set(handles.text15,'String','');
    set(handles.axes1,'Visible','on');
    axes(handles.axes1)
    imshow(mosaic(1:4:end,1:4:end,:));
    drawnow;
end
if (turnonoutput == 1)
    imwrite(mosaic,fullfile(fn_path,[fn_out '.tif']),'tif','Compression','None');
    ICCfid = fopen(fullfile(fn_path,'embeddedICCprofile.icc'));
    rawProfileBytes = fread(ICCfid,Inf,'uint8=>uint8');
    fclose(ICCfid);
    img = Tiff(fullfile(fn_path,[fn_out '.tif']), 'r+');
    img.setTag('ICCProfile', rawProfileBytes);
    img.rewriteDirectory();
    img.close();
    delete(fullfile(fn_path, 'embeddedICCprofile.icc'));
end
clear mosaic dist_from_center imgAll
end
end

function [F] = total_features(img,msk,j,samp_ratio)
%% total_features

img = single(img);

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
M = sqrt(LH.^2 + HL.^2);
clear LH;
clear HL;
F0 = imregionalmax(M,8);
F0 = F0.*msk;

[y1,x1] = find(F0==1);
yd = length(y1);
for i = 1:yd
    tmp = F0((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)).*M((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j));
    [mx11,ind11] = max(tmp(:));
    if (mx11 > 0)
        xpeak = floor((ind11(1)-1)/size(tmp,1))+1;
        ypeak = ind11(1) - (xpeak-1)*size(tmp,1);
        msk11 = single(zeros(size(tmp)));
        msk11(ypeak,xpeak) = 1;
        F0((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)) = F0((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)).*msk11;
    end
end
clear x1 y1 tmp

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
end
clear M tmp tmp1 tmp1n tmp1m
F = F(((F(:,1)~=0) & (F(:,2)~=0) & (F(:,3)~=0)),:);
clear img
end

function [F1,F2] = art_regx(img1,img2,F1,F2,samp_ratio,max_shift)

img1 = single(img1);
[m1 n1] = size(img1);
img2 = single(img2);
mf1 = size(F2,1);

for i = mf1:-1:1
    j = F2(i,3);
    lt = floor(F2(i,1))-((2^j)-1);
    rt = ceil(F2(i,1))+((2^j)-1);
    tp = floor(F2(i,2))-((2^j)-1);
    bt = ceil(F2(i,2))+((2^j)-1);
    tmp = img2(tp:bt,lt:rt);
    tmp1 = imresize(tmp,samp_ratio,'bicubic');
    clear tmp
    d = floor(samp_ratio/2)/samp_ratio;
    tmp1n = (lt-d):1/samp_ratio:(rt+d);
    [~, x1] = min(abs(tmp1n-F2(i,1)));
    clear tmp1n
    tmp1m = (tp-d):1/samp_ratio:(bt+d);
    [~, y1] = min(abs(tmp1m-F2(i,2)));
    clear tmp1m
    d = (((2^(j+1))-2+1)*samp_ratio-1)/2;
    tmp1 = tmp1((y1-d):(y1+d),(x1-d):(x1+d));

    lt = floor(F1(i,1))-((2^j)-1)-max_shift;
    if (lt < 1)
        lt = 1;
    end
    rt = ceil(F1(i,1))+((2^j)-1)+max_shift;
    if (rt > n1)
        rt = n1;
    end
    tp = floor(F1(i,2))-((2^j)-1)-max_shift;
    if (tp < 1)
        tp = 1;
    end
    bt = ceil(F1(i,2))+((2^j)-1)+max_shift;
    if (bt > m1)
        bt = m1;
    end
    tmp = img1(tp:bt,lt:rt);
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
    cc = cc.*single(mskcc);
    cc = cc((floor(mt1/2)+1):end,(floor(nt1/2)+1):end);
    [max_cc, imax] = max(cc(:));
    [ypeak, xpeak] = ind2sub(size(cc),imax(1));
    cnt = sum(sum(cc>(max_cc-sigma_cc)));
    clear cc       
    if ((xpeak <= nt2) && (ypeak <= mt2) && (cnt == 1) && (max_cc > noise_floor))
        F1(i,1) = tmp2n(xpeak);
        F1(i,2) = tmp2m(ypeak);
        F1(i,3) = max_cc;
    else
        F1(i,:) = [];
        F2(i,:) = [];
    end
    clear tmp1 tmp2
    clear tmp2n tmp2m
end

end

function [z,zd] = apply_mapping(I2,m01,m02,n01,n02,X,Ntran)
%% apply_mapping
% Transform second image so that it is registered with the first
%   Input: I2 -> image
%          X -> transform coefficients
%   Outputs: z -> registered image

[m0,n0] = size(I2);
mend = m02;
nend = n02;

vall = 1:n0;
vall = vall - round(n0/2);
uall = 1:m0;
uall = uall - round(m0/2);
x01 = 1 - round(n0/2);
x02 = n0 - round(n0/2);
y01 = 1 - round(m0/2);
y02 = m0 - round(m0/2);
bc = zeros(2*m0+2*n0,2);
bc(1:m0,:) = [x01*ones(m0,1) uall'];
bc((m0+1):(2*m0),:) = [x02*ones(m0,1) uall'];
bc((2*m0+1):(2*m0+n0),:) = [vall' y01*ones(n0,1)];
bc((2*m0+n0+1):(2*m0+2*n0),:) = [vall' y02*ones(n0,1)];
if (Ntran == 0)
    v0 = X(1)*bc(:,1) + X(2)*bc(:,2) + X(3)*bc(:,1).*bc(:,2) + X(4)*ones(2*m0+2*n0,1);
    u0 = X(5)*bc(:,1) + X(6)*bc(:,2) + X(7)*bc(:,1).*bc(:,2) + X(8)*ones(2*m0+2*n0,1);
else
    u0 = zeros(2*m0+2*n0,1);
    v0 = zeros(2*m0+2*n0,1);
    pnt = 1;
    for k = 0:Ntran
        for j = 0:(Ntran-k)
            v0 = v0 + X(pnt)*((bc(:,1).^k).*(bc(:,2).^j));
            u0 = u0 + X(pnt+sum((1:(Ntran+1))))*((bc(:,1).^k).*(bc(:,2).^j));
            pnt = pnt + 1;
        end
    end
end
clear bc uall vall
u0 = u0 + round(m0/2);
v0 = v0 + round(n0/2);
n211 = ceil(min(v0));
n221 = floor(max(v0));
m211 = ceil(min(u0));
m221 = floor(max(u0));
z = single(zeros((m221-m211+1),(n221-n211+1)));
clear u0 v0
zd = single((sqrt(x02^2 + y02^2)+1000)*ones((m221-m211+1),(n221-n211+1)));

m11 = m01;
m12 = mend;
n11 = n01;
n12 = nend;
[x0,y0] = meshgrid(n11:n12,m11:m12);
x0 = x0 - round(n0/2);
y0 = y0 - round(m0/2);
if (Ntran == 0)
    v2 = X(1)*x0 + X(2)*y0 + X(3)*x0.*y0 + X(4);
    u2 = X(5)*x0 + X(6)*y0 + X(7)*x0.*y0 + X(8);
else
    u2 = zeros(size(y0));
    v2 = zeros(size(x0));
    pnt = 1;
    for i = 0:Ntran
        for j = 0:(Ntran-i)
            v2 = v2 + X(pnt)*((x0.^i).*(y0.^j));
            u2 = u2 + X(pnt+sum((1:(Ntran+1))))*((x0.^i).*(y0.^j));
            pnt = pnt + 1;
        end
    end
end
%tmpd0 = sqrt(x0.^2 + y0.^2);
clear x0 y0
u2 = u2(:) + round(m0/2);
v2 = v2(:) + round(n0/2);
n21 = ceil(min(v2));
n22 = floor(max(v2));
m21 = ceil(min(u2));
m22 = floor(max(u2));
[v1,u1] = meshgrid(n21:n22,m21:m22);
tmp = I2(m11:m12,n11:n12);
F = TriScatteredInterp(double(v2),double(u2),double(tmp(:)),'nearest');
%performance tradeoff (use original, not the transformed, pixel distances)
%Fd = TriScatteredInterp(double(v2),double(u2),double(tmpd0(:)),'natural');
clear tmp u2 v2
%clear tmpd0
z0 = single(F(v1,u1));
%zd0 = single(Fd(v1,u1));

clear u1 v1 F
%clear Fd
[v1,u1] = meshgrid((n21-n211+1):(n22-n211+1),(m21-m211+1):(m22-m211+1));
zd0 = single(sqrt((u1-round(m0/2)).^2 + (v1-round(n0/2)).^2));
clear u1 v1
z0(isnan(z0)) = 0;
z0(isinf(z0)) = 0;
msk = (z0~=0);
zd0(isnan(zd0)) = sqrt(x02^2 + y02^2)+1000;
zd0(isinf(zd0)) = sqrt(x02^2 + y02^2)+1000;
msk2 = (zd0 < sqrt(x02^2 + y02^2)+1000);
msk2 = and(msk2,msk);
z((m21-m211+1):(m22-m211+1),(n21-n211+1):(n22-n211+1)) = z((m21-m211+1):(m22-m211+1),(n21-n211+1):(n22-n211+1)).*(~msk) + z0.*msk;
zd((m21-m211+1):(m22-m211+1),(n21-n211+1):(n22-n211+1)) = zd((m21-m211+1):(m22-m211+1),(n21-n211+1):(n22-n211+1)).*(~msk2) + zd0.*msk2;
clear z0 zd0 msk msk2 I2
end

function [d] = dsearchn1(x,xi)
%DSEARCHN N-D nearest point search (modified version of MathWorks function).
%   K = DSEARCHN(X,XI) performs the search without using a triangulation.
%   With large X and small XI, this approach is faster and uses much
%   less memory. 
%
%   [K,D] = DSEARCHN(X,...) also returns the distances D to the closest
%   points. D is a column vector of length p. 

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.7.4.10 $  $Date: 2010/11/22 02:46:36 $

clear A1y A1x
d = zeros(size(xi,1),1);
parfor i = 1:size(xi,1)
    yi = repmat(xi(i,:),size(x,1),1);
    d(i) = min(sum((x-yi).^2,2));
end
d = sqrt(d);
end

function [] = mosiac_quick(fn_path,fn_out,mmpixel,turnonoutput,usegui)

bord = 256;     %pixels
bord = 2*floor(bord/2);

% Get all file names
everyImage = dir(fullfile(fn_path, '*.tif'));
offsetFile = dir(fullfile(fn_path, '*.csv'));

offset_val = csvread(fullfile(fn_path,offsetFile(1).name));
offset_val = offset_val(:,4:5);     %% x-shift(mm), y-shift(mm)
offset_val = offset_val/mmpixel;    %% x-shift(pixel), y-shift(pixel)
[num ~] = size(offset_val);

% Get image dimensions, icc profile, and bit depth
info = imfinfo(fullfile(fn_path,everyImage(1).name));
m0 = info.Height;
n0 = info.Width;
bd = info.BitsPerSample(1);
iccProfile = iccread(fullfile(fn_path,everyImage(1).name));
embeddedICC = iccwrite(iccProfile, fullfile(fn_path, 'embeddedICCprofile.icc'));

offset_val = offset_val - repmat(min(offset_val)-1,num,1);
[~,colind] = max(abs(offset_val(:,1)));
colmx = round(offset_val(colind,1));
[~,rowind] = max(abs(offset_val(:,2)));
rowmx = round(offset_val(rowind,2));
offset_val(:,1) = -(offset_val(:,1) - 1) + colmx;
mtot = abs(rowmx)+m0-1;
ntot = abs(colmx)+n0-1;

if (bd == 8)
    mosaic = uint8(zeros(mtot+2*bord,ntot+2*bord,3));
elseif (bd == 16)
    mosaic = uint16(zeros(mtot+2*bord,ntot+2*bord,3));
end

pnt1 = 1;
for i = 1:num
    img1 = imread((fullfile(fn_path,everyImage(i).name)),'tif');
    [m0,n0,~] = size(img1);
    m1 = round(offset_val(pnt1,2)+bord);
    m2 = round(m1 + m0 - 1);
    n1 = round(offset_val(pnt1,1)+bord);
    n2 = round(n1 + n0 - 1);
    mosaic(m1:m2,n1:n2,:) = img1;
    pnt1 = pnt1 + 1;
    clear img1
end
mosaic(1:bord,:,:) = [];
mosaic((mtot+1):end,:,:) = [];
mosaic(:,1:bord,:) = [];
mosaic(:,(ntot+1):end,:) = [];

if (usegui == 1)
    handles = guidata(PaintingStitch);
    set(handles.text15,'String','');
    set(handles.axes1,'Visible','on');
    axes(handles.axes1)
    imshow(mosaic(1:4:end,1:4:end,:));
    drawnow;
end
if (turnonoutput == 1)
    imwrite(mosaic,fullfile(fn_path,[fn_out '.tif']),'tif','Compression','None');
    ICCfid = fopen(fullfile(fn_path,'embeddedICCprofile.icc'));
    rawProfileBytes = fread(ICCfid,Inf, 'uint8=>uint8');
    fclose(ICCfid);
    img = Tiff(fullfile(fn_path,[fn_out '.tif']),'r+');
    img.setTag('ICCProfile',rawProfileBytes);
    img.rewriteDirectory();
    img.close();
    delete(fullfile(fn_path, 'embeddedICCprofile.icc'));
end
clear mosaic

end
