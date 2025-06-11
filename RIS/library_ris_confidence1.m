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
pk_dist_range = 100;    %bands
bw_dist_range = 100;    %nm

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
features = Inf(1024,6);
[fsz,~] = size(features);
pklist = zeros(1024,1);
pnt = 1;
Lf_list = [1 3 5];
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
                lb = md - pk_dist_range;
                ub = md + pk_dist_range;
                
                if (j == 1)
                    ind = ((rpeak_list(:,3) >= lb) & (rpeak_list(:,3) <= ub));
                    if (sum(ind) > 0)
                        rtmp = rpeak_list(ind,1:4);
                        [Lpr,~] = size(rtmp);
                        rdist = abs(rtmp(:,3) - md);
                        for q = 1:Lpr
                            %tmp = features(rtmp(q,1),rtmp(q,2),i,j);
                            %tmp = tmp{1};
                            indf = (features(:,1)==rtmp(q,1))&(features(:,2)==rtmp(q,2))&(features(:,3)==i)&(features(:,4)==j)&(features(:,5)==k);
                            if (sum(indf) == 0)
                                features(pnt,:) = [rtmp(q,1) rtmp(q,2) i j k rdist(q)];
                                [fsz,~] = size(features);
                                pnt = pnt + 1;
                                if (pnt > fsz)
                                    features = [features;Inf(1024,6)];
                                    pklist = [pklist;zeros(1024,1)];
                                end
                                features(pnt,:) = [rtmp(q,1) rtmp(q,2) i j+1 k Inf];
                                pklist(pnt) = rtmp(q,4);
                                [fsz,~] = size(features);
                                pnt = pnt + 1;
                                if (pnt > fsz)
                                    features = [features;Inf(1024,6)];
                                    pklist = [pklist;zeros(1024,1)];
                                end
                            else
                                tmp = features(indf,6);
                                if (rdist(q) < tmp)
                                    features(indf,6) = rdist(q);
                                    pklist(indf+1) = rtmp(q,4);
                                end
                            end
                        end
                        clear rtmp
                    end
                    clear ind
                elseif (j == 5)
                    ind = ((apeak_list(:,3) >= lb) & (apeak_list(:,3) <= ub));
                    if (sum(ind) > 0)
                        rtmp = apeak_list(ind,1:3);
                        [Lpr,~] = size(rtmp);
                        rdist = abs(rtmp(:,3) - md);
                        for q = 1:Lpr
                            indf = (features(:,1)==rtmp(q,1))&(features(:,2)==rtmp(q,2))&(features(:,3)==i)&(features(:,4)==j)&(features(:,5)==k);
                            if (sum(indf) == 0)
                                features(pnt,:) = [rtmp(q,1) rtmp(q,2) i j k rdist(q)];
                                [fsz,~] = size(features);
                                pnt = pnt + 1;
                                if (pnt > fsz)
                                    features = [features;Inf(1024,6)];
                                    pklist = [pklist;zeros(1024,1)];
                                end
                            else
                                tmp = features(indf,6);
                                if (rdist(q) < tmp)
                                    features(indf,6) = rdist(q);
                                end
                            end
                        end
                        clear rtmp
                    end
                    clear ind
                elseif (j == 3)
                    ind = ((edge_list(:,3) >= lb) & (edge_list(:,3) <= ub));
                    if (sum(ind) > 0)
                        rtmp = edge_list(ind,1:4);
                        [Lpr,~] = size(rtmp);
                        rdist = abs(rtmp(:,3) - md);
                        for q = 1:Lpr
                            indf = (features(:,1)==rtmp(q,1))&(features(:,2)==rtmp(q,2))&(features(:,3)==i)&(features(:,4)==j)&(features(:,5)==k);
                            if (sum(indf) == 0)
                                features(pnt,:) = [rtmp(q,1) rtmp(q,2) i j k rdist(q)];
                                [fsz,~] = size(features);
                                pnt = pnt + 1;
                                if (pnt > fsz)
                                    features = [features;Inf(1024,6)];
                                    pklist = [pklist;zeros(1024,1)];
                                end
                                features(pnt,:) = [rtmp(q,1) rtmp(q,2) i j+1 k Inf];
                                pklist(pnt) = rtmp(q,4);
                                [fsz,~] = size(features);
                                pnt = pnt + 1;
                                if (pnt > fsz)
                                    features = [features;Inf(1024,6)];
                                    pklist = [pklist;zeros(1024,1)];
                                end
                            else
                                tmp = features(indf,6);
                                if (rdist(q) < tmp)
                                    features(indf,6) = rdist(q);
                                    pklist(indf+1) = rtmp(q,4);
                                end
                            end
                        end
                        clear rtmp
                    end
                    clear ind
                end
            end
        end
    end
end

% compute correlation with library spectra
Lf_list = [6 7];
for m = 1:m0
    tmpm = lib(m,:,:);
    tmpm = reshape(tmpm,[n0,p0]);
    for n = 1:n0
        for i = 1:Lp
            for j1 = 1:length(Lf_list)
                j = Lf_list(j1);
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
                        
                        G1 = lib(i,:,lb1:ub1);
                        G1 = reshape(G1,[1 (ub1-lb1+1)]);
                        F1 = tmpm(n,lb1:ub1);
                        tmpcc = sum((G1-mean(G1)).*(F1-mean(F1)))/sqrt(sum((G1-mean(G1)).^2).*sum((F1-mean(F1)).^2));
                        tmpcc = 1 - (tmpcc + 1)/2;
                        
                        features(pnt,:) = [m n i j k tmpcc];
                        [fsz,~] = size(features);
                        pnt = pnt + 1;
                        if (pnt > fsz)
                            features = [features;Inf(1024,6)];
                            pklist = [pklist;zeros(1024,1)];
                        end
                    end
                end
            end
        end
    end
    clear tmpm F G F1 G1 tmp
end
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
                %bw = abs(min(bw1,bw2));
                features(ind1(ii),:) = [m n i j k bw];
            end
        end
    end
end
clear ind ind1 tmp

sigma = Inf(1024,4);
mu = Inf(1024,4);
[sz,~] = size(sigma);
pnt = 1;
for i = 1:Lp
    for j = 1:Lf
        for k = 1:tot_ris_table(i,j)
            ind = (features(:,3) == i)&(features(:,4) == j)&(features(:,5) == k)&(features(:,6)~=Inf);
            tmp = features(ind,6);
            sigma(pnt,:) = [i j k std(tmp)];
            mu(pnt,:) = [i j k mean(tmp)];
            tmp1 = (features(ind,6) - mean(tmp))/std(tmp);
            
            [h,~] = hist(tmp1,1001);
            nout = min(tmp1):0.001:max(tmp1);
            cdf = zeros(1,length(h));
            lim = [0 0];
            lim1 = 0;
            lim2 = 0;
            for ii = 1:length(h)
                cdf(ii) = sum(h(1:ii));
                if ((cdf(ii) >= 0.01*sum(h)) && (lim1 == 0))
                    lim(1) = nout(i);
                    lim1 = 1;
                end
                if ((cdf(ii) >= 0.99*sum(h)) && (lim2 == 0))
                    lim(2) = nout(ii);
                    lim2 = 1;
                end
            end
            tmp1 = (tmp1 - lim(1))/(lim(2) - lim(1))/2;
            tmp1(tmp1<0) = 0;
            tmp1(tmp1>1) = 1;
            features(ind,6) = tmp1;
            pnt = pnt + 1;
            if (pnt > sz)
                sigma = [sigma;single(Inf(1024,4))];
                mu = [mu;single(Inf(1024,4))];
                [sz,~] = size(sigma);
            end
        end
    end
end
sigma(pnt:end,:) = [];
mu(pnt:end,:) = [];

score = ones(m0,n0,Lp);
for m = 1:m0
    for n = 1:n0
        for i = 1:Lp
            ind = (features(:,1) == m)&(features(:,2) == n)&(features(:,3) == i)&(features(:,6) ~= Inf);
            tmp = features(ind,6);
            [Ltmp,~] = size(tmp);
            tot = sum(tot_ris_table(i,:));
            tmp = [tmp;ones(tot-Ltmp,1)];
            score(m,n,i) = sum(tmp)/tot;
        end
    end
end

score = reshape(score,[m0,Lp]);
imagesc(score)
axis image

