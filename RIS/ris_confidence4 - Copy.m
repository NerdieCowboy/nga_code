clear all
close all
clc

fn_path = 'C:\data\Pacino\RIS\';
features_table_fn = 'C:\damon\dissertation\features_pigments_throughRenasaince.csv';
cube_fn = 'sub_sc_crop_cal_Pacino_ChristMajesty_vnir_xnir';
derivative_fn = [fn_path 'slope_' cube_fn];
derivative2_fn = [fn_path 'slope_slope_' cube_fn];
derivative3_fn = [fn_path 'slope_slope_slope_' cube_fn];
%abs_derivative_fn = [fn_path 'absorbance_slope_' cube_fn];
%abs_derivative2_fn = [fn_path 'absorbance_slope_slope_' cube_fn];
%abs_derivative3_fn = [fn_path 'absorbance_slope_slope_slope_' cube_fn];
spectra_lib = 'resampled_EarlyRennaisance3';
lib_derivative_fn = [fn_path 'slope_' spectra_lib];
lib_derivative2_fn = [fn_path 'slope_slope_' spectra_lib];
lib_derivative3_fn = [fn_path 'slope_slope_slope_' spectra_lib];
%abs_lib_derivative_fn = [fn_path 'absorbance_slope_' spectra_lib];
%abs_lib_derivative2_fn = [fn_path 'absorbance_slope_slope_' spectra_lib];
%abs_lib_derivative3_fn = [fn_path 'absorbance_slope_slope_slope_' spectra_lib];
spectra_names = {'lamp black', 'bone black', 'burnt umber', ...
    'red ochre', 'Yellow Ochre', 'gypsum', 'chalk', 'Madder lake', ...
    'Realgar', 'Malachite', 'Orpiment', 'Indigo', 'Azurite', ...
    'Red lead', 'vermilion', 'Green earth', 'verdigris', 'lead white', ...
    'ultramarine', 'Naples Yellow', 'Smalt', 'indian yellow', ...
    'Copper resinate', 'lead tin yellow', 'Van Dyke brown', ...
    'Carmine lake', 'cobalt blue', 'raw umber', 'Cd Red', ...
    'titanium white', 'Cd yellow', 'zinc white'};
pk_dist_range = 80;    %nm
bw_dist_range = 100;    %nm
abs_dist_range = 150;   %nm
pkdelta_thresh1 = 0.05;
pkdelta_thresh2 = 0.02;
blksz = 8;
%w = 35;

lib_test = 0;
lambdaon = 0;
lambda = [];
fnh = [derivative_fn '.hdr'];
fid = fopen(fnh);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
    if (lambdaon == 0)
        [~,~,e] = regexp(line,'^samples=','match','start','end');
        if (e>0)
            n0 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^lines=','match','start','end');
        if (e>0)
            m0 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^bands=','match','start','end');
        if (e>0)
            p0 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^datatype=','match','start','end');
        if (e>0)
            datatype = str2double(line((e+1):end));
            if (datatype == 1)
                datatype = 'uint8';
                elseif (datatype == 4)
                    datatype = 'single';
                elseif (datatype == 5)
                    datatype = 'double';
                elseif (datatype == 12)
                    datatype = 'uint16';
            end
        end
        [~,~,e] = regexp(line,'^interleave=','match','start','end');
        if (e>0)
            interleave = line((e+1):end);
        end
        [~,~,e] = regexp(line,'^byteorder=','match','start','end');
        if (e>0)
            byteorder = str2double(line((e+1):end));
            if (byteorder == 0)
                byteorder = 'ieee-le';
            elseif (byteorder == 1)
                byteorder = 'ieee-be';
            end
        end
        [~,~,e] = regexp(line,'^headeroffset=','match','start','end');
        if (e>0)
            headeroffset = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^wavelength={','match','start','end');
        if (e>0)
            lambdaon = 1;
            tmp1 = line((e+1):end);
            tmp2 = regexp(tmp1,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda = [lambda str2double(tmp2(i))];
                end
                [~,~,e1] = regexp(tmp2(i),'}$','match','start','end');
                if (e1{1} > 0)
                    tmp3 = cell2mat(tmp2(i));
                    tmp3 = tmp3(1:(e1{1}-1));
                    lambda = [lambda str2double(tmp3)];
                end
            end
            [~,~,e] = regexp(line,'}$','match','start','end');
            if (e>0)
                lambdaon = 0;
            end
        end
    end
    if (lambdaon == 1)
        [~,~,e] = regexp(line,'}$','match','start','end');
        if (e>0)
            lambdaon = 0;
            tmp1 = line(1:(e-1));
            tmp2 = regexp(tmp1,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda = [lambda str2double(tmp2(i))];
                end
            end
        else
            tmp2 = regexp(line,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda = [lambda str2double(tmp2(i))];
                end
            end
        end
    end
    clear msk tmp1 tmp2 line tmp3
end
fclose(fid);
cube = multibandread([fn_path cube_fn],[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
if (lib_test == 1)
    cube = reshape(cube,[m0 1 n0]);
end
der1 = multibandread(derivative_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
if (lib_test == 1)
    der1 = reshape(der1,[m0 1 n0]);
end
der2 = multibandread(derivative2_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
if (lib_test == 1)
    der2 = reshape(der2,[m0 1 n0]);
end
der3 = multibandread(derivative3_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
if (lib_test == 1)
    der3 = reshape(der3,[m0 1 n0]);
end
%{
% data cube with local mean subtracted off
data_muremov = single(zeros(m0,n0,p0));
for p = 1:p0
    p011 = p - floor(w/2);
    if (p011 < 1)
        p011 = 1;
    end
    p012 = p + floor(w/2);
    if (p012 > p0)
        p012 = p0;
    end
    tmp = cube(:,:,p011:p012);
    data_muremov(:,:,p) = single(cube(:,:,p) - mean(tmp,3));
    clear tmp
end
%}
% compute black/white features
databw = cube;
lambdabw = lambda;
databw(:,:,(lambda<360)) = [];
databw(:,:,(lambda>830)) = [];
lambdabw(lambdabw<360) = [];
lambdabw(lambdabw>830) = [];
p0bw = length(lambdabw);

[cie] = getCieStruct(lambdabw);
d50 = cieIllD( [5000], cie );
d50(isnan(d50)) = 0;

databw = shiftdim(databw,2);
databw = databw(:);
databw = reshape(databw,[p0bw m0*n0]);
XYZ = ref2XYZ(databw,cie.cmf10deg,d50);
clear databw
imgLab = XYZ2Lab(XYZ,whitepoint('d50')*100);
clear XYZ
imgLab = shiftdim(imgLab,1);
imgLab = reshape(imgLab,[m0 n0 3]);

ab = sqrt(imgLab(:,:,2).^2 + imgLab(:,:,3).^2);
brit = imgLab(:,:,1);
mxL = max(brit(:));
low_high =  stretchlim(brit/mxL,[0.005 0.995]);
brit1 = imadjust(brit/mxL,low_high,[]);
clear brit
mxab = max(ab(:));
low_high =  stretchlim(ab/mxab,[0.005 0.995]);
ab1 = imadjust(ab/mxab,low_high,[]);
clear ab
ab1 = abs(1 - ab1);
bwf = brit1.*ab1;
clear brit1 ab1

pnt = 1;
features = zeros(m0*n0*20,6);
[fsz,~] = size(features);
for m = 1:blksz:m0
    mend = m + blksz - 1;
    if (mend > m0)
        mend = m0;
    end
    for n = 1:blksz:n0
        [m n]
        nend = n + blksz - 1;
        if (nend > n0)
            nend = n0;
        end
        [features0,tot_ris_table,Lp,Lf] = match_features(cube(m:mend,n:nend,:),der1(m:mend,n:nend,:),der2(m:mend,n:nend,:),der3(m:mend,n:nend,:),fn_path,features_table_fn,spectra_lib,spectra_names,pk_dist_range,bw_dist_range,abs_dist_range,pkdelta_thresh1,pkdelta_thresh2);
        features0(:,1) = features0(:,1) + m - 1;
        features0(:,2) = features0(:,2) + n - 1;
        [fsz0,~] = size(features0);
        while ((pnt + fsz0 - 1) > fsz)
            features = [features;zeros(m0*n0*20,6)];
            [fsz,~] = size(features);
        end
        features(pnt:(pnt + fsz0 - 1),:) = features0;
        pnt = pnt + fsz0;
        clear features0
    end
end
features(features(:,1)==0,:) = [];

lib_test = 1;
lambdaon = 0;
lambda = [];
fnh = [lib_derivative_fn '.hdr'];
fid = fopen(fnh);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
    if (lambdaon == 0)
        [~,~,e] = regexp(line,'^samples=','match','start','end');
        if (e>0)
            n01 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^lines=','match','start','end');
        if (e>0)
            m01 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^bands=','match','start','end');
        if (e>0)
            p01 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^datatype=','match','start','end');
        if (e>0)
            datatype = str2double(line((e+1):end));
            if (datatype == 1)
                datatype = 'uint8';
                elseif (datatype == 4)
                    datatype = 'single';
                elseif (datatype == 5)
                    datatype = 'double';
                elseif (datatype == 12)
                    datatype = 'uint16';
            end
        end
        [~,~,e] = regexp(line,'^interleave=','match','start','end');
        if (e>0)
            interleave = line((e+1):end);
        end
        [~,~,e] = regexp(line,'^byteorder=','match','start','end');
        if (e>0)
            byteorder = str2double(line((e+1):end));
            if (byteorder == 0)
                byteorder = 'ieee-le';
            elseif (byteorder == 1)
                byteorder = 'ieee-be';
            end
        end
        [~,~,e] = regexp(line,'^headeroffset=','match','start','end');
        if (e>0)
            headeroffset = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^wavelength={','match','start','end');
        if (e>0)
            lambdaon = 1;
            tmp1 = line((e+1):end);
            tmp2 = regexp(tmp1,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda = [lambda str2double(tmp2(i))];
                end
                [~,~,e1] = regexp(tmp2(i),'}$','match','start','end');
                if (e1{1} > 0)
                    tmp3 = cell2mat(tmp2(i));
                    tmp3 = tmp3(1:(e1{1}-1));
                    lambda = [lambda str2double(tmp3)];
                end
            end
            [~,~,e] = regexp(line,'}$','match','start','end');
            if (e>0)
                lambdaon = 0;
            end
        end
    end
    if (lambdaon == 1)
        [~,~,e] = regexp(line,'}$','match','start','end');
        if (e>0)
            lambdaon = 0;
            tmp1 = line(1:(e-1));
            tmp2 = regexp(tmp1,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda = [lambda str2double(tmp2(i))];
                end
            end
        else
            tmp2 = regexp(line,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda = [lambda str2double(tmp2(i))];
                end
            end
        end
    end
    clear msk tmp1 tmp2 line tmp3
end
fclose(fid);
cube = multibandread([fn_path spectra_lib],[m01,n01,p01],datatype,headeroffset,interleave,byteorder);
if (lib_test == 1)
    cube = reshape(cube,[m01 1 n01]);
end
der1 = multibandread(lib_derivative_fn,[m01,n01,p01],datatype,headeroffset,interleave,byteorder);
if (lib_test == 1)
    der1 = reshape(der1,[m01 1 n01]);
end
der2 = multibandread(lib_derivative2_fn,[m01,n01,p01],datatype,headeroffset,interleave,byteorder);
if (lib_test == 1)
    der2 = reshape(der2,[m01 1 n01]);
end
der3 = multibandread(lib_derivative3_fn,[m01,n01,p01],datatype,headeroffset,interleave,byteorder);
if (lib_test == 1)
    der3 = reshape(der3,[m01 1 n01]);
end
[featuresL,tot_ris_tableL,LpL,LfL] = match_features(cube,der1,der2,der3,fn_path,features_table_fn,spectra_lib,spectra_names,pk_dist_range,bw_dist_range,abs_dist_range,pkdelta_thresh1,pkdelta_thresh2);

clear cube lib der1 lib_der1 der2 lib_der2 der3 lib_der3
% normalize
for j = 1:Lf
    ind = ((features(:,4) == j)&(features(:,6) ~= Inf));
    if (sum(ind) > 0)
        tmp = features(ind,6);
        mx = max(tmp);
        mn = min(tmp);
        tmp = (tmp-mn)/(mx-mn);
        sigma = std(tmp);
        mu = mean(tmp);
        tmp = (tmp - mu)/sigma;
        features(ind,6) = tmp;
        clear tmp ind
        
        ind = ((featuresL(:,4) == j)&(featuresL(:,6) ~= Inf));
        if (sum(ind) > 0)
            tmp = featuresL(ind,6);
            tmp = (tmp-mn)/(mx-mn);
            tmp = (tmp - mu)/sigma;
            featuresL(ind,6) = tmp;
        end
        clear tmp ind
    end
end
tmp = features(:,6);
tmp(isinf(tmp)) = [];
mx1 = max(tmp);
mn1 = min(tmp);
clear tmp
tmp = featuresL(:,6);
tmp(isinf(tmp)) = [];
mx2 = max(tmp);
mn2 = min(tmp);
mx = max(mx1,mx2);
mn = min(mn1,mn2);
clear tmp
features(:,6) = (features(:,6)-mn)/(mx-mn);
featuresL(:,6) = (featuresL(:,6)-mn)/(mx-mn);

tmp = features(:,6);
tmp(isinf(tmp)) = [];
low_high1 = stretchlim(tmp,[0.01 0.99]);
clear tmp
tmp = featuresL(:,6);
tmp(isinf(tmp)) = [];
low_high2 = stretchlim(tmp,[0.01 0.99]);
mx = max(low_high1(2),low_high2(2));
mn = min(low_high1(1),low_high2(1));
clear tmp
features(:,6) = (features(:,6)-mn)/(mx-mn);
featuresL(:,6) = (featuresL(:,6)-mn)/(mx-mn);
ind = (features(:,6)<0) & (features(:,6)~=Inf);
features(ind,6) = 0;
ind = (features(:,6)>1) & (features(:,6)~=Inf);
features(ind,6) = 1;

ncores = feature('numCores');
if (ncores > 1)
    isOpen = matlabpool('size') > 0;
    if (isOpen == 1)
        matlabpool close
    end
    matlabpool(ncores)
end
%{
score = single(zeros(m0,n0,Lp));
tot = sum(tot_ris_table,2);
parfor m = 1:m0
    ind = ((features(:,1) == m)&(features(:,6) ~= Inf));
    tmp1 = features(ind,:);
    for n = 1:n0
        ind = (tmp1(:,2) == n);
        tmp2 = tmp1(ind,:);
        for i = 1:Lp
            %ind = (features(:,1) == m)&(features(:,2) == n)&(features(:,3) == i)&(features(:,6) ~= Inf);
            ind = (tmp2(:,3) == i);
            %tmp = features(ind,6);
            tmp = tmp2(ind,6);
            [Ltmp,~] = size(tmp);
            tmp = [tmp;ones((tot(i)-Ltmp),1)];
            score(m,n,i) = 1-(sum(tmp)/tot(i));
        end
    end
end

scoreL = single(zeros(Lp,Lp));
parfor m = 1:Lp
    if (tot(m) > 0)
        ind = ((featuresL(:,1) == m)&(featuresL(:,6) ~= Inf));
        tmp1 = featuresL(ind,:);
        for n = 1:Lp
            if (tot(n) > 0)
                ind = (tmp1(:,3) == n);
                tmp2 = tmp1(ind,6);
                [Ltmp,~] = size(tmp2);
                tmp2 = [tmp2;ones((tot(n)-Ltmp),1)];
                scoreL(m,n) = 1-(sum(tmp2)/tot(n));
            end
        end
    end
end
%}
score = single(zeros(m0,n0,Lp));
tot = sum(tot_ris_table,2);
parfor m = 1:m0
    ind = ((features(:,1) == m)&(features(:,6) ~= Inf));
    tmp1 = features(ind,:);
    for n = 1:n0
        ind = (tmp1(:,2) == n);
        tmp2 = tmp1(ind,:);
        for i = 1:Lp
            %ind = (features(:,1) == m)&(features(:,2) == n)&(features(:,3) == i)&(features(:,6) ~= Inf);
            ind = (tmp2(:,3) == i);
            %tmp = features(ind,6);
            tmp = tmp2(ind,6);
            [Ltmp,~] = size(tmp);
            tmp = [tmp;ones((tot(i)-Ltmp),1)];
            score(m,n,i) = 1-(sum(tmp)/tot(i));
        end
    end
end

scoreL = single(zeros(Lp,Lp));
parfor m = 1:Lp
    if (tot(m) > 0)
        ind = ((featuresL(:,1) == m)&(featuresL(:,6) ~= Inf));
        tmp1 = featuresL(ind,:);
        for n = 1:Lp
            if (tot(n) > 0)
                ind = (tmp1(:,3) == n);
                tmp2 = tmp1(ind,6);
                [Ltmp,~] = size(tmp2);
                tmp2 = [tmp2;ones((tot(n)-Ltmp),1)];
                scoreL(m,n) = 1-(sum(tmp2)/tot(n));
            end
        end
    end
end

if (ncores > 1)
    isOpen = matlabpool('size') > 0;
    if (isOpen == 1)
        matlabpool close
    end
end

for i = 1:Lp
    figure,imagesc(score(:,:,i))
    axis image
    axis off
    title(spectra_names{i})
    caxis([0 1])
    colormap(gray)
    fn = [fn_path cube_fn '_' spectra_names{i} '.tif'];
    imwrite(uint16((2^16-1)*score(:,:,i)),fn,'tif','Compression','None')
end
