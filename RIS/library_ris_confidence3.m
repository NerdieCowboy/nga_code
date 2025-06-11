clear all
close all
clc

%load('C:\data\Pacino\RIS\contremov_crop_Pacino_ChristMajesty_vnir_xnir.mat')
features_table_fn = 'C:\damon\dissertation\features_pigments_throughRenasaince.csv';
cube_fn = 'C:\data\Pacino\RIS\resampled_EarlyRennaisance2';
derivative_fn = 'C:\data\Pacino\RIS\slope_resampled_EarlyRennaisance2';
derivative2_fn = 'C:\data\Pacino\RIS\slope_slope_resampled_EarlyRennaisance2';
derivative3_fn = 'C:\data\Pacino\RIS\slope_slope_slope_resampled_EarlyRennaisance2';
hsz = 3;
dthresh = 0.01;
pk_dist_range = 40;    %bands
bw_dist_range = 50;    %nm

lambdaon = 0;
lambda = [];
fnh = [cube_fn '.hdr'];
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
lib = multibandread(cube_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
lib = reshape(lib,[m0 1 n0]);
der1 = multibandread(derivative_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
der1 = reshape(der1,[m0 1 n0]);
der2 = multibandread(derivative2_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
der2 = reshape(der2,[m0 1 n0]);
der3 = multibandread(derivative3_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
der3 = reshape(der3,[m0 1 n0]);
[m0,n0,p0] = size(lib);

% list of element in each pigment
%[~,~,N] = size(hsi_map);
%hsi_map = (hsi_map+1)/2;
pigment_element_table = csvread(features_table_fn,1,10);

% spectral features for each pigment
fid = fopen(features_table_fn);
features_table = textscan(fid,'%s','delimiter',',','EmptyValue',NaN);
fclose(fid);
tmp = reshape(features_table{1},[110,33]);
tmp = tmp';
pigments = tmp(2:33,1);
colors = tmp(2:33,2);
pigment_ris_table = tmp(2:33,3:9);
clear tmp
[Lp,Lf] = size(pigment_ris_table);

% derivative features
%   reflectance peaks (zero-crossing), absorption peaks (zero-crossing), 
%   transition edges (peak), transition edge fwhm
% reflectance features
%   reflectance peak fwhm
rpeak_list = single(zeros(m0*n0,4)); %[m,n,lambda,band]
apeak_list = single(zeros(m0*n0,3)); %[m,n,lambda]
edge_list = single(zeros(m0*n0,4)); %[m,n,lambda,band]
[Ltotr,~] = size(rpeak_list);
[Ltota,~] = size(apeak_list);
[Ltote,~] = size(edge_list);
ptotr = 1;
ptota = 1;
ptote = 1;
for p = (hsz+ceil(hsz/2)):(p0-hsz-floor(hsz/2))
    %tmp = multibandread(derivative_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Range',[p-1 2 p+1]});
    %tmp1 = multibandread(derivative2_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Range',[p-1 p+1]});
    tmp = der1(:,:,(p-1):2:(p+1));
    tmp1 = der2(:,:,(p-1):(p+1));
    
    snr_pass = (abs(tmp1(:,:,2)) > dthresh);
    mskr = (tmp(:,:,1) >= 0) & (tmp(:,:,2) <= 0) & (snr_pass == 1);
    mska = (tmp(:,:,1) <= 0) & (tmp(:,:,2) >= 0) & (snr_pass == 1);
    clear tmp snr_pass
    [r,c] = find(mskr==1);
    L = length(r);
    while ((ptotr + L - 1) > Ltotr)
        rpeak_list = [rpeak_list;single(zeros(round(sqrt(m0*n0)),4))];
        [Ltotr,~] = size(rpeak_list);
    end
    for i = 1:L
        rpeak_list(ptotr,:) = [r(i) c(i) lambda(p) p];
        ptotr = ptotr + 1;
    end
    clear mskr
    [r,c] = find(mska==1);
    L = length(r);
    while ((ptota + L - 1) > Ltota)
        apeak_list = [apeak_list;single(zeros(round(sqrt(m0*n0)),3))];
        [Ltota,~] = size(apeak_list);
    end
    for i = 1:L
        apeak_list(ptota,:) = [r(i) c(i) lambda(p)];
        ptota = ptota + 1;
    end
    clear mska
    
    %tmp2 = multibandread(derivative3_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Direct',p});
    tmp2 = der3(:,:,p);
    snr_pass = (abs(tmp2) > dthresh);
    clear tmp2
    mske = (tmp1(:,:,1) >= 0) & (tmp1(:,:,3) <= 0) & (snr_pass == 1);
    clear tmp1 snr_pass
    [r,c] = find(mske==1);
    L = length(r);
    while ((ptote + L - 1) > Ltote)
        edge_list = [edge_list;single(zeros(round(sqrt(m0*n0)),4))];
        [Ltote,~] = size(edge_list);
    end
    for i = 1:L
        edge_list(ptote,:) = [r(i) c(i) lambda(p) p];
        ptote = ptote + 1;
    end
    clear mske
end
clear r c
rpeak_list(rpeak_list(:,1)==0,:) = [];
apeak_list(apeak_list(:,1)==0,:) = [];
edge_list(edge_list(:,1)==0,:) = [];

% match features
%clear features
tot_ris_table = zeros(Lp,Lf);
for i = 1:Lp
    for j = 1:Lf
        tmp = pigment_ris_table(i,j);
        tmp = tmp{1};
        if (~strcmp(tmp,'0'))
            tmp1 = regexp(tmp, '\;', 'split');
            Ltmp1 = length(tmp1);
            tot_ris_table(i,j) = tot_ris_table(i,j) + Ltmp1;
            %features(1:m0,1:n0,i,j) = {Inf(Ltmp1,1)};
        else
            %features(1:m0,1:n0,i,j) = {Inf(1,1)};
        end
    end
end

% match spectral features
%cube = multibandread(cube_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
%der1 = multibandread(derivative_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
features = Inf(1024,6);
[fsz,~] = size(features);
pklist = zeros(1024,1);
pnt = 1;
Lf_list = [1 3 5];
Lf_list2 = [6 7];
for m = 1:m0
    indj1m = (rpeak_list(:,1) == m);
    indj3m = (edge_list(:,1) == m);
    indj5m = (apeak_list(:,1) == m);
    rtmpj1m = rpeak_list(indj1m,1:4);
    rtmpj3m = edge_list(indj3m,1:4);
    rtmpj5m = apeak_list(indj5m,1:3);
    tmpm = lib(m,:,:);
    tmpm = reshape(tmpm,[n0,p0]);
    for n = 1:n0
        [m n]
        indj1 = (rtmpj1m(:,2) == n);
        indj3 = (rtmpj3m(:,2) == n);
        indj5 = (rtmpj5m(:,2) == n);
        rtmpj1 = rtmpj1m(indj1,1:4);
        rtmpj3 = rtmpj3m(indj3,1:4);
        rtmpj5 = rtmpj5m(indj5,1:3);
        if ((sum(indj1) > 0)||(sum(indj3) > 0)||(sum(indj5) > 0))
            for i = 1:Lp
                for j1 = 1:length(Lf_list)
                    j = Lf_list(j1);
                    tmp = pigment_ris_table(i,j);
                    tmp = tmp{1};
                    if (~strcmp(tmp,'0'))
                        tmp1 = regexp(tmp, '\;', 'split');
                        Ltmp1 = length(tmp1);
                        for k = 1:Ltmp1
                            md = single(str2double(tmp1{k}));
                            %lb = md - pk_dist_range;
                            %ub = md + pk_dist_range;

                            if (j == 1)
                                if (sum(indj1) > 0)
                                    rdist = abs(rtmpj1(:,3) - md);
                                    %[i j k]
                                    [mnq,q] = min(rdist);
                                    if (mnq(1) > pk_dist_range)
                                        mnq(1) = Inf;
                                    end
                                    features(pnt,:) = [m n i j k mnq(1)];
                                    pnt = pnt + 1;
                                    if (pnt > fsz)
                                        features = [features;Inf(1024,6)];
                                        pklist = [pklist;zeros(1024,1)];
                                    end
                                    [fsz,~] = size(features);
                                    features(pnt,:) = [m n i j+1 k Inf];
                                    pklist(pnt) = rtmpj1(q(1),4);
                                    pnt = pnt + 1;
                                    if (pnt > fsz)
                                        features = [features;Inf(1024,6)];
                                        pklist = [pklist;zeros(1024,1)];
                                    end
                                    [fsz,~] = size(features);
                                end
                            elseif (j == 3)
                                if (sum(indj3) > 0)
                                    rdist = abs(rtmpj3(:,3) - md);
                                    %[i j k]
                                    [mnq,q] = min(rdist);
                                    if (mnq(1) > pk_dist_range)
                                        mnq(1) = Inf;
                                    end
                                    features(pnt,:) = [m n i j k mnq(1)];
                                    pnt = pnt + 1;
                                    if (pnt > fsz)
                                        features = [features;Inf(1024,6)];
                                        pklist = [pklist;zeros(1024,1)];
                                    end
                                    [fsz,~] = size(features);
                                    features(pnt,:) = [m n i j+1 k Inf];
                                    pklist(pnt) = rtmpj3(q(1),4);
                                    pnt = pnt + 1;
                                    if (pnt > fsz)
                                        features = [features;Inf(1024,6)];
                                        pklist = [pklist;zeros(1024,1)];
                                    end
                                    [fsz,~] = size(features);
                                end
                            elseif (j == 5)
                                if (sum(indj5) > 0)
                                    rdist = abs(rtmpj5(:,3) - md);
                                    %[i j k]
                                    [mnq,q] = min(rdist);
                                    if (mnq(1) > pk_dist_range)
                                        mnq(1) = Inf;
                                    end
                                    features(pnt,:) = [m n i j k mnq(1)];
                                    pnt = pnt + 1;
                                    if (pnt > fsz)
                                        features = [features;Inf(1024,6)];
                                        pklist = [pklist;zeros(1024,1)];
                                    end
                                    [fsz,~] = size(features);
                                end 
                            end
                        end
                    end
                end
            end
        end
        for i = 1:Lp
            for j1 = 1:length(Lf_list2)
                j = Lf_list2(j1);
                tmp = pigment_ris_table(i,j);
                tmp = tmp{1};
                if (~strcmp(tmp,'0'))
                    tmp1 = regexp(tmp, '\;', 'split');
                    Ltmp1 = length(tmp1);
                    for k = 1:Ltmp1
                        tmp2 = regexp(tmp1{k}, '\:', 'split');
                        lb = str2double(tmp2{1});
                        if (lb < lambda(1))
                            lb = lambda(1);
                        end
                        ub = str2double(tmp2{2});
                        if (ub > lambda(end))
                            ub = lambda(end);
                        end
                        tmp = abs(lambda - lb);
                        [~,ind] = min(tmp);
                        lb1 = ind;
                        tmp = abs(lambda - ub);
                        [~,ind] = min(tmp);
                        ub1 = ind;

                        G1 = lib(i,lb1:ub1);
                        G1 = reshape(G1,[1 (ub1-lb1+1)]);
                        F1 = tmpm(n,lb1:ub1);
                        tmpcc = sum((G1-mean(G1)).*(F1-mean(F1)))/sqrt(sum((G1-mean(G1)).^2).*sum((F1-mean(F1)).^2));
                        tmpcc = 1 - (tmpcc + 1)/2;

                        features(pnt,:) = [m n i j k tmpcc];
                        pnt = pnt + 1;
                        if (pnt > fsz)
                            features = [features;Inf(1024,6)];
                            pklist = [pklist;zeros(1024,1)];
                        end
                        [fsz,~] = size(features);
                    end
                end
            end
        end
        clear F G F1 G1
    end
    clear tmpm tmp
end
clear tmp ind rtmp
features(pnt:end,:) = [];
[szf,~] = size(features);
pklist(pnt:end) = [];

%bw from table
Lf_list = [2 4];
for j1 = 1:length(Lf_list)
    j = Lf_list(j1);
    ind = (features(:,4)==j);
    ind1 = (1:szf)';
    tmp = features(ind,:);
    ind1 = ind1(ind);
    for ii = 1:length(tmp)
        m = tmp(ii,1);
        n = tmp(ii,2);
        i = tmp(ii,3);
        k = tmp(ii,5);
        tmpBW = pigment_ris_table(i,j);
        tmpBW = tmpBW{1};
        if (~strcmp(tmpBW,'0'))
            tmpBW1 = regexp(tmpBW, '\;', 'split');
            if (~strcmp(tmpBW1{k},'0'))
                mdBW = single(str2double(tmpBW1{k}));

                %measured bw
                pk_p = pklist(ind1(ii));
                pk1 = pk_p - bw_dist_range;
                pk2 = pk_p + bw_dist_range;
                if (pk1 < 1)
                    pk1 = 1;
                    pk0 = pk_p;
                elseif (pk2 > p0)
                    pk2 = p0;
                    pk0 = ceil((pk2 - pk1 + 1)/2);
                else
                    pk0 = ceil((pk2 - pk1 + 1)/2);
                end
                pk_m = features(ind1(ii),1);
                pk_n = features(ind1(ii),2);
                if (j == 2)
                    cube_tmp = lib(pk_m,pk_n,pk1:pk2);
                elseif (j == 4)
                    cube_tmp = der1(pk_m,pk_n,pk1:pk2);
                end
                cube_tmp = cube_tmp(:);
                pk_val = cube_tmp(pk0);
                cube_tmp = cube_tmp - pk_val/2;
                lambda_tmp = lambda(pk1:pk2);
                lambda_tmp = lambda_tmp(:);
                lambda_tmp = (lambda_tmp(2:(pk0-1))-lambda_tmp(1:(pk0-2)))/2+lambda_tmp(1:(pk0-2));
                zc1 = (cube_tmp(1:(pk0-2))<=0) & (cube_tmp(2:(pk0-1))>=0);
                zc1 = lambda_tmp.*zc1;
                zc1(zc1==0)=[];
                bw1 = Inf;
                if (length(zc1)>0)
                    lambda_tmp = lambda(pk1:pk2);
                    zc1 = zc1(end);
                    bw1 = 2*abs(zc1-lambda_tmp(pk0));
                end
                lambda_tmp = lambda(pk1:pk2);
                lambda_tmp = lambda_tmp(:);
                lambda_tmp = (lambda_tmp((pk0+2):end)-lambda_tmp((pk0+1):(end-1)))/2+lambda_tmp((pk0+1):(end-1));
                zc2 = (cube_tmp((pk0+1):(end-1))>=0) & (cube_tmp((pk0+2):end)<=0);
                zc2 = lambda_tmp.*zc2;
                zc2(zc2==0)=[];
                bw2 = Inf;
                if (length(zc2)>0)
                    lambda_tmp = lambda(pk1:pk2);
                    zc2 = zc2(1);
                    bw2 = 2*abs(zc2-lambda_tmp(pk0));
                end
                bw = abs(min(bw1,bw2) - mdBW);
                if (bw > (2*bw_dist_range))
                    bw = Inf;
                end
                features(ind1(ii),:) = [m n i j k bw];
            end
        end
    end
end
clear ind ind1 tmp

