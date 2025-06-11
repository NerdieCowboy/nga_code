%function [pigment_ris_table,tot_ris_table] = define_features(cube_fn,cube_path,spectra_names)

%{
clear all
close all
clc
cube_path = 'C:\data\Pacino\RIS\';
cube_fn = 'resampled_EarlyRennaisance2';
spectra_names = {'lamp black', 'bone black', 'burnt umber', ...
    'red ochre', 'Yellow Ochre', 'gypsum', 'chalk', 'Madder lake', ...
    'Realgar', 'Malachite', 'Orpiment', 'Indigo', 'Azurite', ...
    'Red lead', 'vermilion', 'Green earth', 'verdigris', 'lead white', ...
    'ultramarine', 'Naples Yellow', 'Smalt', 'indian yellow', ...
    'Copper resinate', 'lead tin yellow', 'Van Dyke brown', ...
    'Carmine lake', 'cobalt blue', 'raw umber', 'Cd Red', ...
    'titanium white', 'Cd yellow', 'zinc white'};
pk_dist_range = 20;    %nm
bw_dist_range = 50;    %bands
%}
derivative_fn = [cube_path 'slope_' cube_fn];
derivative2_fn = [cube_path 'slope_slope_' cube_fn];
derivative3_fn = [cube_path 'slope_slope_slope_' cube_fn];
hsz = 3;
%dthresh = 0.01;
dark_thresh = 0.05;
bw_dist_range = 50;    %bands
Lf = 5;

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
cube = multibandread(cube_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
cube = reshape(cube,[m0 1 n0]);
der1 = multibandread(derivative_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
der1 = reshape(der1,[m0 1 n0]);
der2 = multibandread(derivative2_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
der2 = reshape(der2,[m0 1 n0]);
der3 = multibandread(derivative3_fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
der3 = reshape(der3,[m0 1 n0]);
[m0,n0,p0] = size(cube);
Lp = m0;

% derivative features
%   reflectance peaks (zero-crossing), absorption peaks (zero-crossing), 
%   transition edges (peak), transition edge fwhm
% reflectance features
%   reflectance peak fwhm
rpeak_list = single(zeros(m0*n0,4)); %[m,n,lambda,band]
apeak_list = single(zeros(m0*n0,4)); %[m,n,lambda]
edge_list = single(zeros(m0*n0,4)); %[m,n,lambda,band]
[Ltotr,~] = size(rpeak_list);
[Ltota,~] = size(apeak_list);
[Ltote,~] = size(edge_list);
ptotr = 1;
ptota = 1;
ptote = 1;
%mx_cube = max(cube(:));
%mx_der1 = max(der1(:));
mu_cube = mean(cube,3);
mu_der1 = mean(der1,3);
mu_der2 = mean(der2,3);
mu_der3 = mean(der3,3);
sigma_cube = std(cube,0,3);
sigma_der1 = std(der1,0,3);
sigma_der2 = std(der2,0,3);
sigma_der3 = std(der3,0,3);
prevr = zeros(m0,n0);
preva = zeros(m0,n0);
preve = zeros(m0,n0);
for p = (hsz+ceil(hsz/2)):(p0-hsz-floor(hsz/2))
    tmp0 = cube(:,:,(p-1):(p+1));
    tmp = der1(:,:,(p-1):(p+1));
    tmp1 = der2(:,:,(p-1):(p+1));
    %dthresh1 = dthresh*mean(tmp0,3);
    %dthresh21 = dthresh*mean(tmp,3);
    dthresh21 = 0;
    dthresh22 = mu_der1;
    
    %snr_pass = ((abs(tmp1(:,:,1)) > dthresh1) | (abs(tmp1(:,:,3)) > dthresh1)) & (mu_cube > dark_thresh);
    snr_pass = (abs(tmp1(:,:,2)) > 0.5*sigma_der2) & (mu_cube > dark_thresh);
    mskr = (tmp(:,:,1) >= 0) & (tmp(:,:,3) <= 0) & (snr_pass == 1);
    mska = (tmp(:,:,1) <= 0) & (tmp(:,:,3) >= 0) & (snr_pass == 1);
    clear snr_pass
    [r,c] = find(mskr==1);
    L = length(r);
    while ((ptotr + L - 1) > Ltotr)
        rpeak_list = [rpeak_list;single(zeros(round(sqrt(m0*n0)),4))];
        [Ltotr,~] = size(rpeak_list);
    end
    for i = 1:L
        %if ((tmp0(r(i),c(i),1) <= tmp0(r(i),c(i),2)) && (tmp0(r(i),c(i),3) <= tmp0(r(i),c(i),2)))
        if (prevr(r(i),c(i)) == 0)
            rpeak_list(ptotr,:) = [r(i) c(i) lambda(p) p];
            prevr(r(i),c(i)) = 1;
            ptotr = ptotr + 1;
        else
            prevr(r(i),c(i)) = 0;
        end
    end
    clear mskr
    [r,c] = find(mska==1);
    L = length(r);
    while ((ptota + L - 1) > Ltota)
        apeak_list = [apeak_list;single(zeros(round(sqrt(m0*n0)),4))];
        [Ltota,~] = size(apeak_list);
    end
    for i = 1:L
        %if ((tmp0(r(i),c(i),1) >= tmp0(r(i),c(i),2)) && (tmp0(r(i),c(i),3) >= tmp0(r(i),c(i),2)))
        if (preva(r(i),c(i)) == 0)
            apeak_list(ptota,:) = [r(i) c(i) lambda(p) p];
            preva(r(i),c(i)) = 1;
            ptota = ptota + 1;
        else
            preva(r(i),c(i)) = 0;
        end
    end
    clear mska
    
    tmp2 = der3(:,:,p);
    %snr_pass = ((abs(tmp2) > dthresh21) & (tmp(:,:,2) > dthresh22)) & (mu_cube > dark_thresh);
    snr_pass = (abs(tmp2) > 0.25*sigma_der3) & (mu_cube > dark_thresh);
    clear tmp2
    mske = (tmp(:,:,2) > 0) &(tmp1(:,:,1) >= 0) & (tmp1(:,:,3) <= 0) & (snr_pass == 1);
    %clear tmp1 snr_pass
    [r,c] = find(mske==1);
    L = length(r);
    while ((ptote + L - 1) > Ltote)
        edge_list = [edge_list;single(zeros(round(sqrt(m0*n0)),4))];
        [Ltote,~] = size(edge_list);
    end
    for i = 1:L
        %if ((tmp(r(i),c(i),1) <= tmp(r(i),c(i),2)) && (tmp(r(i),c(i),3) <= tmp(r(i),c(i),2)))
        if (preve(r(i),c(i)) == 0)
            edge_list(ptote,:) = [r(i) c(i) lambda(p) p];
            preve(r(i),c(i)) = 1;
            ptote = ptote + 1;
        else
            preve(r(i),c(i)) = 0;
        end
    end
    mske = (tmp(:,:,2) <= 0) &(tmp1(:,:,1) >= 0) & (tmp1(:,:,3) <= 0) & (snr_pass == 1);
    [r,c] = find(mske==1);
    L = length(r);
    while ((ptota + L - 1) > Ltota)
        apeak_list = [apeak_list;single(zeros(round(sqrt(m0*n0)),4))];
        [Ltota,~] = size(apeak_list);
    end
    for i = 1:L
        if (preva(r(i),c(i)) == 0)
            apeak_list(ptota,:) = [r(i) c(i) lambda(p) p];
            preva(r(i),c(i)) = 1;
            ptota = ptota + 1;
        else
            preva(r(i),c(i)) = 0;
        end
    end
    mske = (tmp(:,:,2) >= 0) &(tmp1(:,:,1) <= 0) & (tmp1(:,:,3) >= 0) & (snr_pass == 1);
    clear tmp1 snr_pass
    [r,c] = find(mske==1);
    L = length(r);
    while ((ptotr + L - 1) > Ltotr)
        rpeak_list = [rpeak_list;single(zeros(round(sqrt(m0*n0)),4))];
        [Ltotr,~] = size(rpeak_list);
    end
    for i = 1:L
        if (prevr(r(i),c(i)) == 0)
            rpeak_list(ptotr,:) = [r(i) c(i) lambda(p) p];
            prevr(r(i),c(i)) = 1;
            ptotr = ptotr + 1;
        else
            prevr(r(i),c(i)) = 0;
        end
    end
    clear mske
end
clear r c
rpeak_list(rpeak_list(:,1)==0,:) = [];
apeak_list(apeak_list(:,1)==0,:) = [];
edge_list(edge_list(:,1)==0,:) = [];

% define features
tot_ris_table = zeros(Lp,Lf);
for i = 1:Lp
    j = 1;
    ind = (rpeak_list(:,1)==i);
    tmp = rpeak_list(ind,:);
    tmp = tmp(:,3);
    Ltmp = length(tmp);
    tot_ris_table(i,j) = Ltmp;
    strtmp = '';
    if (Ltmp > 0)
        strtmp = num2str(tmp(1));
        for k = 2:Ltmp
            strtmp = [strtmp ';' num2str(tmp(k))];
        end
    end
    pigment_ris_table{i,j} = strtmp;
    clear strtmp
    tmp = rpeak_list(ind,:);
    tmp = tmp(:,4);
    strtmp = '';
    if (Ltmp > 0)
        strtmp = num2str(tmp(1));
        for k = 2:Ltmp
            strtmp = [strtmp ';' num2str(tmp(k))];
        end
    end
    pigment_ris_table_band{i,j} = strtmp;
    clear strtmp
    
    j = 3;
    ind = (edge_list(:,1)==i);
    tmp = edge_list(ind,:);
    tmp = tmp(:,3);
    Ltmp = length(tmp);
    tot_ris_table(i,j) = Ltmp;
    strtmp = '';
    if (Ltmp > 0)
        strtmp = num2str(tmp(1));
        for k = 2:Ltmp
            strtmp = [strtmp ';' num2str(tmp(k))];
        end
    end
    pigment_ris_table{i,j} = strtmp;
    clear strtmp
    tmp = edge_list(ind,:);
    tmp = tmp(:,4);
    strtmp = '';
    if (Ltmp > 0)
        strtmp = num2str(tmp(1));
        for k = 2:Ltmp
            strtmp = [strtmp ';' num2str(tmp(k))];
        end
    end
    pigment_ris_table_band{i,j} = strtmp;
    clear strtmp
    
    j = 5;
    ind = (apeak_list(:,1)==i);
    tmp = apeak_list(ind,:);
    tmp = tmp(:,3);
    Ltmp = length(tmp);
    tot_ris_table(i,j) = Ltmp;
    strtmp = '';
    if (Ltmp > 0)
        strtmp = num2str(tmp(1));
        for k = 2:Ltmp
            strtmp = [strtmp ';' num2str(tmp(k))];
        end
    end
    pigment_ris_table{i,j} = strtmp;
    clear strtmp
end

%bw from table
Lf_list = [2 4];
for i = 1:Lp
    for j1 = 1:length(Lf_list)
        j = Lf_list(j1);
        tmpBW = pigment_ris_table_band(i,j-1);
        if (~isempty(tmpBW))
            tmpBW1 = regexp(tmpBW, '\;', 'split');
            tmpBW1 = tmpBW1{1};
            strtmp = '';
            for k = 1:tot_ris_table(i,j-1)
                pk_p = single(str2double(tmpBW1{k}));
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
                if (j == 2)
                    cube_tmp = cube(i,1,pk1:pk2);
                elseif (j == 4)
                    cube_tmp = der1(i,1,pk1:pk2);
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
                bw = min(bw1,bw2);
                %{
                if (bw > (2*bw_dist_range))
                    bw = Inf;
                end
                %}
                if (isinf(bw))
                    bw = 0;
                end
                if (isempty(strtmp))
                    strtmp = num2str(bw);
                else
                    strtmp = [strtmp ';' num2str(bw)];
                end
                pigment_ris_table{i,j} = strtmp;
            end
        end
    end
end
clear tmp

fid = fopen('define_features.csv','w');
%tmp = 'Pigment	Color,Reflectance peak (nm),Reflectance peak - FWHM (nm),Transition edge (nm),Transition edge - FWHM (nm),Absorption bands (nm),"interesting" spectral region (nm),"interesting" 1st derivative region (nm)';
tmp1 = 'Pigment	Color,Reflectance peak (nm),Reflectance peak - FWHM (nm),Transition edge (nm),Transition edge - FWHM (nm),Absorption bands (nm)';
fprintf(fid,'%s\n',char(tmp1));
for i = 1:Lp
    tmp1 = char(spectra_names{i});
    for j = 1:Lf
        tmp2 = pigment_ris_table(i,j);
        tmp2 = tmp2{1};
        if (isempty(tmp2))
            tmp1 = [tmp1 ',0'];
        else
            tmp1 = [tmp1 ',' char(tmp2)];
        end
    end 
    fprintf(fid,'%s\n',char(tmp1));
end
clear tmp1
fclose(fid);

for i = 1:Lp
    lib1 = cube(i,:,:);
    lib1 = lib1(:);
    figure
    subplot(311),plot(lambda,lib1)
    hold on
    title(spectra_names{i})
    axis([425 2500 0 max(lib1)])
    
    j = 1;
    tmp = pigment_ris_table(i,j);
    if (~isempty(tmp))
        tmp1 = regexp(tmp, '\;', 'split');
        tmp1 = tmp1{1};
        for k = 1:tot_ris_table(i,j)
            [~,ind] = min(abs(lambda - str2double(tmp1(k))));
            plot(str2double(tmp1(k)),lib1(ind(1)),'rx')
            hold on
        end
    end
    j = 5;
    tmp = pigment_ris_table(i,j);
    if (~isempty(tmp))
        tmp1 = regexp(tmp, '\;', 'split');
        tmp1 = tmp1{1};
        for k = 1:tot_ris_table(i,j)
            [~,ind] = min(abs(lambda - str2double(tmp1(k))));
            plot(str2double(tmp1(k)),lib1(ind(1)),'gx')
            hold on
        end
    end
    hold off

    lib_der11 = der1(i,:,:);
    lib_der11 = lib_der11(:);
    lib_der21 = der2(i,:,:);
    lib_der21 = lib_der21(:);
    subplot(312),plot(lambda,lib_der11)
    hold on
    axis([425 2500 min(lib_der11) max(lib_der11)])
    j = 3;
    tmp = pigment_ris_table(i,j);
    if (~isempty(tmp))
        tmp1 = regexp(tmp, '\;', 'split');
        tmp1 = tmp1{1};
        for k = 1:tot_ris_table(i,j)
            [~,ind] = min(abs(lambda - str2double(tmp1(k))));
            plot(str2double(tmp1(k)),lib_der11(ind(1)),'rx')
            hold on
        end
    end
    hold off
    subplot(313),plot(lambda,lib_der21)
    axis([425 2500 min(lib_der21) max(lib_der21)])
    %pause
end
