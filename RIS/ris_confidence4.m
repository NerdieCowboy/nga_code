clear all
close all
clc

fn_path = 'C:\data\Pacino\RIS\';
features_table_fn = 'C:\damon\dissertation\features_pigments_throughRenasaince.csv';
%cube_fn = 'sub_sc_crop_cal_Pacino_ChristMajesty_vnir_xnir';
%cube_fn = 'sc_crop_cal_Pacino_ChristMajesty_vnir_xnir';
cube_fn = 'spatialavg_sub_sc_crop_cal_Pacino_ChristMajesty_vnir_xnir';
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
pk_dist_range = 50;    %nm
bw_dist_range = 100;    %nm
abs_dist_range = 150;   %nm
pkdelta_thresh1 = 0.02;
pkdelta_thresh2 = 0.02;
blksz = 8;
low_penalty = 0.5;
%w = 35;

%features
pigment_family_col = 1;
pigment_col = 2;
color_col = 3;
reflectance_col = 4;
transition_edge_col = 5;
transition_edge_fwhm_col = 6;
transition_edge_asym_col = 7;
absorption_col = 8;
strong_absorption_center_col = 9;
strong_absorption_width_col = 10;
gradual_slope_col = 11;
interesting_col = 12;
organic_col = 13;
xrf_cols = [14 113];

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
%ab1 = abs(1 - ab1);
%bwf = brit1.*ab1;
blkf = sqrt(brit1.^2 + ab1.^2);
whtf = sqrt(brit1.^2 + (1-ab1).^2);
clear brit1 ab1

% replace missing features with a penalty proportational to the
% "difficulty" in indentifying it
pnt = 1;
features = zeros(m0*n0*20,7);
%{
missing_penalty = single(zeros(1024,5));
[Lmp,~] = size(missing_penalty);
pnt_mp = 1;
%}
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
        [features0,tot_ris_table,Lp,Lf,pigment_ris_table,~] = match_features(cube(m:mend,n:nend,:),der1(m:mend,n:nend,:),der2(m:mend,n:nend,:),der3(m:mend,n:nend,:),fn_path,features_table_fn,spectra_lib,spectra_names,pk_dist_range,bw_dist_range,abs_dist_range,pkdelta_thresh1,pkdelta_thresh2,0);
        %{
        [Lmp1,~] = size(missing_penalty1);
        while ((pnt_mp + Lmp1 - 1) > Lmp)
           missing_penalty = [missing_penalty;single(Inf(1024,5))];
           [Lmp,~] = size(missing_penalty);
        end
        missing_penalty(pnt:(pnt+Lmp1-1),:) = missing_penalty1;
        pnt_mp = pnt_mp + Lmp1;
        %}
        features0(:,1) = features0(:,1) + m - 1;
        features0(:,2) = features0(:,2) + n - 1;
        [fsz0,~] = size(features0);
        while ((pnt + fsz0 - 1) > fsz)
            features = [features;zeros(m0*n0*20,7)];
            [fsz,~] = size(features);
        end
        features(pnt:(pnt + fsz0 - 1),:) = features0;
        pnt = pnt + fsz0;
        clear features0
    end
end
features(features(:,1)==0,:) = [];
%missing_penalty(missing_penalty(:,1)==0,:) = [];
%clear missing_penalty1
%{
missing_penalty_whole = single(zeros(Lp*Lf*100,4));
cnt = 1;
for i = 1:Lp
    for j = 1:Lf
        tmp = pigment_ris_table(i,j);
        tmp = tmp{1};
        if (~strcmp(tmp,'0'))
            tmp1 = regexp(tmp, '\;', 'split');
            Ltmp1 = length(tmp1);
            for k = 1:Ltmp1
                tmp_mp_ind = (missing_penalty(:,1)==i) & (missing_penalty(:,2)==j) & (missing_penalty(:,3)==k);
                tmp_mp = missing_penalty(tmp_mp_ind,:);
                pkdelta1 = sum(tmp_mp(:,4).*tmp_mp(:,5))/sum(tmp_mp(:,5));
                if (j == (reflectance_col-3))
                    missing_penalty_whole(cnt,:) = [i j k pkdelta1];
                    cnt = cnt + 1;
                elseif (j == (transition_edge_col-3))
                    missing_penalty_whole(cnt,:) = [i j k pkdelta1];
                    cnt = cnt + 1;
                    missing_penalty_whole(cnt,:) = [i j+1 k pkdelta1];
                    cnt = cnt + 1;
                    missing_penalty_whole(cnt,:) = [i j+2 k pkdelta1];
                    cnt = cnt + 1;
                elseif (j == (absorption_col-3))
                    missing_penalty_whole(cnt,:) = [i j k pkdelta1];
                    cnt = cnt + 1;
                elseif (j == (strong_absorption_center_col-3))
                    missing_penalty_whole(cnt,:) = [i j k pkdelta1];
                    cnt = cnt + 1;
                    missing_penalty_whole(cnt,:) = [i j+1 k pkdelta1];
                    cnt = cnt + 1;
                elseif (j == (gradual_slope_col-3))
                    missing_penalty_whole(cnt,:) = [i j k 1];
                    cnt = cnt + 1;
                elseif (j == (interesting_col-3))
                    missing_penalty_whole(cnt,:) = [i j k 1];
                    cnt = cnt + 1;
                end
            end
        end
    end
end
missing_penalty_whole(cnt:end,:) = [];

for j = 1:Lf
    tmp_mp_ind = (missing_penalty_whole(:,2)==j);
    if (sum(tmp_mp_ind) > 0)
        tmp_mp = missing_penalty_whole(tmp_mp_ind,:);
        tmp_mp = tmp_mp/max(tmp_mp);
        low_high = stretchlim(tmp_mp,[0.01 0.99]);
        tmp_mp = (tmp_mp - low_high(1))/(low_high(2) - low_high(1));
        tmp_mp = tmp_mp/2 + 0.5;
        if ((j ~= (gradual_slope_col-3)) && (j ~= (interesting_col-3)))
            missing_penalty_whole(tmp_mp_ind,4) = abs(tmp_mp);
        end
    end
end
%}
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
[featuresL,tot_ris_tableL,LpL,LfL,~,missing_penalty] = match_features(cube,der1,der2,der3,fn_path,features_table_fn,spectra_lib,spectra_names,pk_dist_range,bw_dist_range,abs_dist_range,pkdelta_thresh1,pkdelta_thresh2,1);

cnt = 1;
for i = 1:Lp
    for j = 1:Lf
        tmp = pigment_ris_table(i,j);
        tmp = tmp{1};
        if (~strcmp(tmp,'0'))
            tmp_mp_ind1 = ((missing_penalty(:,2)==j) & (missing_penalty(:,4)~=Inf));
            tmp_mp1 = missing_penalty(tmp_mp_ind1,4);
            if (~isempty(tmp_mp1))
                mx_mp1 = max(tmp_mp1);
            else
                mx_mp1 = 1;
            end
            
            tmp1 = regexp(tmp, '\;', 'split');
            Ltmp1 = length(tmp1);
            for k = 1:Ltmp1
                tmp_mp_ind = (missing_penalty(:,1)==i) & (missing_penalty(:,2)==j) & (missing_penalty(:,3)==k);
                tmp_mp = missing_penalty(tmp_mp_ind,4)/mx_mp1;
                low_high = stretchlim(tmp_mp,[0.01 0.99]);
                tmp_mp = (tmp_mp - low_high(1))/(low_high(2) - low_high(1));
                missing_penalty(tmp_mp_ind,4) = tmp_mp*(1-low_penalty) + low_penalty;
            end
        end
    end
end
clear der1 lib_der1 der2 lib_der2 der3 lib_der3

% normalize
for j = 1:Lf
    ind1 = ((features(:,4) == j)&(features(:,6) ~= Inf)&~isnan(features(:,6)));
    ind2 = ((featuresL(:,4) == j)&(featuresL(:,6) ~= Inf)&~isnan(featuresL(:,6)));
    if (sum(ind1) > 0)
        tmp1 = features(ind1,6);
        tmp2 = featuresL(ind2,6);
        mx1 = max(tmp1);
        mn1 = min(tmp1);
        mx2 = max(tmp2);
        mn2 = min(tmp2);
        mx = max(mx1,mx2);
        mn = min(mn1,mn2);
        tmp1 = (tmp1-mn)/(mx-mn);
        tmp2 = (tmp2-mn)/(mx-mn);
        sigma = std(tmp1);
        mu = mean(tmp1);
        tmp1 = (tmp1 - mu)/sigma + mu/sigma;
        tmp2 = (tmp2 - mu)/sigma + mu/sigma;
        features(ind1,6) = tmp1;
        featuresL(ind2,6) = tmp2;
        clear tmp1 tmp2 ind1 ind2
    end
end
tmp = features(:,6);
tmp(isinf(tmp)) = [];
tmp(isnan(tmp)) = [];
mx1 = max(tmp);
clear tmp
tmp = featuresL(:,6);
tmp(isinf(tmp)) = [];
tmp(isnan(tmp)) = [];
mx2 = max(tmp);
mx = max(mx1,mx2);
clear tmp
features(:,6) = features(:,6)/mx;
featuresL(:,6) = featuresL(:,6)/mx;

tmp = features(:,6);
tmp(isinf(tmp)) = [];
tmp(isnan(tmp)) = [];
low_high1 = stretchlim(tmp,[0.01 0.99]);
clear tmp
tmp = featuresL(:,6);
tmp(isinf(tmp)) = [];
tmp(isnan(tmp)) = [];
low_high2 = stretchlim(tmp,[0.01 0.99]);
mx = max(low_high1(2),low_high2(2));
clear tmp
features(:,6) = features(:,6)/mx;
featuresL(:,6) = featuresL(:,6)/mx;
ind = (features(:,6)<0) & (features(:,6)~=Inf);
features(ind,6) = 0;
ind = (features(:,6)>1) & (features(:,6)~=Inf);
features(ind,6) = 1;
ind = isnan(features(:,6));
features(ind,6) = 1;
ind = (featuresL(:,6)<0) & (featuresL(:,6)~=Inf);
featuresL(ind,6) = 0;
ind = (featuresL(:,6)>1) & (featuresL(:,6)~=Inf);
featuresL(ind,6) = 1;
ind = isnan(featuresL(:,6));
featuresL(ind,6) = 1;
clear ind

ncores = feature('numCores');
if (ncores > 1)
    isOpen = matlabpool('size') > 0;
    if (isOpen == 1)
        matlabpool close
    end
    matlabpool(ncores)
end

score = single(zeros(m0,n0,Lp));
tot = sum(tot_ris_table,2);
parfor m = 1:m0
    ind = (features(:,1) == m);
    tmp1 = features(ind,:);
    for n = 1:n0
        ind = (tmp1(:,2) == n);
        tmp2 = tmp1(ind,:);
        for i = 1:Lp
            ind = (tmp2(:,3) == i);
            tmp3 = tmp2(ind,:);
            tmp_score = 0;
            for j = 1:Lf
                tmp_prt = pigment_ris_table(i,j);
                tmp_prt = tmp_prt{1};
                if (~strcmp(tmp_prt,'0'))
                    tmp1_prt = regexp(tmp_prt, '\;', 'split');
                    Ltmp1 = length(tmp1_prt);
                    for k = 1:Ltmp1
                        ind = (tmp3(:,4) == j) & (tmp3(:,5) == k);
                        tmp4 = tmp3(ind,:);
                        if (~isempty(tmp4))
                            if (isinf(tmp4(1,6)))
                                ind_mp = ((missing_penalty(:,1)==i) & (missing_penalty(:,2)==j) & (missing_penalty(:,3)==k));
                                if (sum(ind_mp)>0)
                                    tmp_score = tmp_score + missing_penalty(ind_mp,4);
                                end
                            else
                                tmp_score = tmp_score + tmp4(1,6);
                            end
                        else
                            ind_mp = ((missing_penalty(:,1)==i) & (missing_penalty(:,2)==j) & (missing_penalty(:,3)==k));
                            tmp_score = tmp_score + missing_penalty(ind_mp,4);
                        end
                    end         
                end
            end
            if (tot(i) > 0)
                score(m,n,i) = 1-(tmp_score/tot(i));
            else
                score(m,n,i) = 0;
            end
        end
    end
end

scoreL = single(zeros(Lp,Lp));
parfor m = 1:Lp
    ind = (featuresL(:,1) == m);
    tmp1 = featuresL(ind,:);
    for n = 1:1
        ind = (tmp1(:,2) == n);
        tmp2 = tmp1(ind,:);
        for i = 1:Lp
            ind = (tmp2(:,3) == i);
            tmp3 = tmp2(ind,:);
            tmp_score = 0;
            for j = 1:Lf
                tmp_prt = pigment_ris_table(i,j);
                tmp_prt = tmp_prt{1};
                if (~strcmp(tmp_prt,'0'))
                    tmp1_prt = regexp(tmp_prt, '\;', 'split');
                    Ltmp1 = length(tmp1_prt);
                    for k = 1:Ltmp1
                        ind = (tmp3(:,4) == j) & (tmp3(:,5) == k);
                        tmp4 = tmp3(ind,:);
                        if (~isempty(tmp4))
                            if (isinf(tmp4(1,6)))
                                ind_mp = ((missing_penalty(:,1)==i) & (missing_penalty(:,2)==j) & (missing_penalty(:,3)==k));
                                if (sum(ind_mp)>0)
                                    tmp_score = tmp_score + missing_penalty(ind_mp,4);
                                end
                            else
                                tmp_score = tmp_score + tmp4(1,6);
                            end
                        else
                            ind_mp = ((missing_penalty(:,1)==i) & (missing_penalty(:,2)==j) & (missing_penalty(:,3)==k));
                            tmp_score = tmp_score + missing_penalty(ind_mp,4);
                        end
                    end         
                end
            end
            if (tot(i) > 0)
                scoreL(m,i) = 1-(tmp_score/tot(i));
            else
                scoreL(m,i) = 0;
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

fid = fopen(features_table_fn);
features_table = textscan(fid,'%s','delimiter',',','EmptyValue',NaN);
fclose(fid);
tmp = reshape(features_table{1},[xrf_cols(2),33]);
tmp = tmp';
pigment_family = tmp(2:33,pigment_family_col);
clear tmp

sthresh = single(zeros(32,1));
for i = 1:32
    pigment_family0 = pigment_family(i);
    family_msk = strcmp(pigment_family,pigment_family0);
    tmp = scoreL(i,:)';
    tmp2 = scoreL(:,i);
    tmp3 = tmp.*tmp2;
    
    tmp3(family_msk==1) = 0;
    [sthresh(i),~] = max(tmp3);
    if (sthresh(i) == 0)
        sthresh(i) = Inf;
    end
end

%{
for i = 1:Lp
    tmp = score(:,:,i);
    msk = (tmp>sthresh(i));
    tmp = tmp.*single(msk);
    figure,imagesc(tmp)
    axis image
    axis off
    title(spectra_names{i})
    %caxis([0 1])
    colormap(gray)
    fn = [fn_path 'thesh_' cube_fn '_' spectra_names{i} '.tif'];
    imwrite(uint16((2^16-1)*tmp),fn,'tif','Compression','None')
end
%}
% confidence
conf = single(zeros(m0,n0,Lp));
class_family_name = false(m0,n0,Lp);
for m = 1:m0
    for n = 1:n0
        tmp = score(m,n,:);
        tmp = tmp(:);
        for i = 1:Lp
            pigment_family0 = pigment_family(i);
            family_msk = strcmp(pigment_family,pigment_family0);
            
            tmp2 = scoreL(:,i);
            tmp3 = tmp.*tmp2;

            [mx,ind] = max(tmp3);
            if (ind == i)
                tmp3(i) = 0;
                [mx0,~] = max(tmp3);
                
                tmp3(family_msk==1) = 0;
                [mx1,~] = max(tmp3);
                if (mx0 ~= mx1)
                    class_family_name(m,n,i) = 1;
                end
                
                conf(m,n,i) = abs(mx - mx1);
            else
                conf(m,n,i) = 0;
            end
        end
    end
end

% normalize confidence
tmp = conf;
tmp = tmp(:);
tmp(tmp==0) = [];
mx = max(tmp);
tmp = tmp/mx;
low_high = stretchlim(tmp,[0.01 0.99]);
clear tmp
conf = conf/mx;
conf = conf/low_high(2);
conf(conf>1) = 1;

for i = 1:Lp
    %{
    tmp = score(:,:,i);
    tmp2 = scoreL(i,i);
    tmp3 = tmp*tmp2;
    msk = (tmp3>sthresh(i));
    %}
    tmp = score(:,:,i);
    %tmp = tmp.*single(msk);
    figure,subplot(121),imagesc(tmp)
    axis image
    axis off
    caxis([0 1])
    title(spectra_names{i})
    colormap(gray)
    
    fn = [fn_path 'score_' cube_fn '_' spectra_names{i} '.tif'];
    imwrite(uint16((2^16-1)*tmp),fn,'tif','Compression','None')
    
    tmp = conf(:,:,i);
    %tmp = tmp.*single(msk);
    subplot(122),imagesc(tmp)
    axis image
    axis off
    caxis([0 1])
    colormap(gray)
    
    fn = [fn_path 'confidence_' cube_fn '_' spectra_names{i} '.tif'];
    imwrite(uint16((2^16-1)*tmp),fn,'tif','Compression','None')
end
