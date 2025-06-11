clear all
close all
clc

load('C:\data\Pacino\RIS\contremov_crop_Pacino_ChristMajesty_vnir_xnir.mat')
features_table_fn = 'C:\damon\dissertation\features_pigments_throughRenasaince.csv';
%cube_fn = 'C:\data\Pacino\crop_Pacino_ChristMajesty_vnir_xnir';
derivative_fn = 'C:\data\Pacino\sc_slope_crop_Pacino_ChristMajesty_vnir_xnir';
derivative2_fn = 'C:\data\Pacino\sc_slope_slope_crop_Pacino_ChristMajesty_vnir_xnir';
derivative3_fn = 'C:\data\Pacino\sc_slope_slope_slope_crop_Pacino_ChristMajesty_vnir_xnir';
hsz = 3;
dthresh = 0.01;

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

% list of element in each pigment
[~,~,N] = size(hsi_map);
hsi_map = (hsi_map+1)/2;
pigment_element_table = csvread(features_table_fn,1,10);

% spectral features for each pigment
fid = fopen(features_table_fn);
features_table = textscan(fid,'%s','delimiter',',','EmptyValue',NaN);
fclose(fid);
tmp = reshape(features_table{1},[110,35]);
tmp = tmp';
pigments = tmp(2:35,1);
colors = tmp(2:35,2);
pigment_ris_table = tmp(2:35,3:9);
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
    tmp = multibandread(derivative_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Range',[p-1 2 p+1]});
    tmp1 = multibandread(derivative2_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Range',[p-1 p+1]});
    
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
    
    tmp2 = multibandread(derivative3_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Direct',p});
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
tot_ris_table = zeros(Lp,Lf);
obs_ris_table = zeros(m0,n0,Lp);
% match spectral features
for i = 1:Lp
    for j = 1:Lf
        tmp = pigment_ris_table(i,j);
        tmp = tmp{1};
        if (~strcmp(tmp,'0'))
            tmp1 = regexp(tmp, '\;', 'split');
            Ltmp1 = length(tmp1);
            tot_ris_table(i,j) = tot_ris_table(i,j) + Ltmp1;
            for k = 1:Ltmp1
                tmp2 = regexp(tmp1{k}, '\:', 'split');
                lb = str2double(tmp2{1});
                ub = str2double(tmp2{2});
                
                if (j == 1)
                    ind = ((rpeak_list(:,3) >= lb) & (rpeak_list(:,3) <= ub));
                    rtmp = rpeak_list(ind,1:2);
                    [Lpr,~] = size(rtmp);
                    prev_match = zeros(m0,n0);
                    if (sum(ind) > 0)
                        for q = 1:Lpr
                            if (prev_match(rtmp(q,1),rtmp(q,2)) == 0)
                                obs_ris_table(rtmp(q,1),rtmp(q,2),i) = obs_ris_table(rtmp(q,1),rtmp(q,2),i) + 1;
                                prev_match(rtmp(q,1),rtmp(q,2)) = 1;
                            end
                        end
                    end
                    clear rtmp ind prev_match
                elseif (j==3)
                    ind = ((edge_list(:,3) >= lb) & (edge_list(:,3) <= ub));
                    etmp = edge_list(ind,1:2);
                    [Lpe,~] = size(etmp);
                    prev_match = zeros(m0,n0);
                    if (sum(ind) > 0)
                        for q = 1:Lpe
                            if (prev_match(etmp(q,1),etmp(q,2)) == 0)
                                obs_ris_table(etmp(q,1),etmp(q,2),i) = obs_ris_table(etmp(q,1),etmp(q,2),i) + 1;
                                prev_match(etmp(q,1),etmp(q,2)) = 1;
                            end
                        end
                    end
                    clear etmp ind prev_match
                elseif (j==5)
                    ind = ((apeak_list(:,3) >= lb) & (apeak_list(:,3) <= ub));
                    atmp = apeak_list(ind,1:2);
                    [Lpa,~] = size(atmp);
                    prev_match = zeros(m0,n0);
                    if (sum(ind) > 0)
                        for q = 1:Lpa
                            if (prev_match(atmp(q,1),atmp(q,2)) == 0)
                                obs_ris_table(atmp(q,1),atmp(q,2),i) = obs_ris_table(atmp(q,1),atmp(q,2),i) + 1;
                                prev_match(atmp(q,1),atmp(q,2)) = 1;
                            end
                        end
                    end
                    clear atmp ind prev_match
                end
            end
        end
    end
end
tot_ris_table(:,2) = 0;
tot_ris_table(:,4) = 0;
tot_ris_table(:,6) = 0;
tot_ris_table(:,7) = 0;
tot_ris_table = sum(tot_ris_table,2);
tot_ris_table = reshape(tot_ris_table,[1 1 Lp]);
tot_ris_table = repmat(tot_ris_table,[m0 n0 1]);
%tot_score = (pigment_matches + obs_ris_table)./(pigment_total + tot_ris_table);
tot_score = (pigment_matches + obs_ris_table);
obs_ris_table = obs_ris_table./tot_ris_table;
clear tot_ris_table
obs_ris_table(isnan(obs_ris_table)) = 0;

for i = 1:Lp
    %figure,imagesc(tot_score(:,:,i))
    figure,imagesc(obs_ris_table(:,:,i))
    %figure,imshow(obs_ris_table(:,:,i)==1)
    title(pigments{i})
    axis image
    axis off
    colormap gray
    %caxis([0 max(max(tot_score(:,:,i)))])
    caxis([0 max(max(obs_ris_table(:,:,i)))])
    colorbar
    pause
end

%{
for i = 1:N
    tmp = hsi_map(:,:,i);
    tmp = tmp(:);
    tmp(isnan(tmp)) = [];
    stretchlim(tmp,[0.01 0.99])
    max(tmp)
    disp(spectra_names{i})
    pause
end

T = sum(abs(hsi_map),3);
T = repmat(T,[1 1 N]);
%Q = (T - abs(hsi_map))./T;
Q = hsi_map./T;
mx = max(Q(:));
mn = min(Q(:));
Q = (Q-mn)/(mx-mn);

for i = 1:N
    tmp = Q(:,:,i);
    imagesc(tmp,[0 1])
    axis image
    title(spectra_names{i})
    colorbar
    pause
end
%}
