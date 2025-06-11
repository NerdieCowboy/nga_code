%function [pigment_ris_table,tot_ris_table] = define_features(cube_fn,cube_path,spectra_names)

clear all
close all
clc
tic

cube_path = 'C:\data\Pacino\RIS\';
%cube_fn = 'resampled_EarlyRennaisance2';
cube_fn = 'resampled_EarlyRennaisance3';
spectra_names = {'lamp black', 'bone black', 'burnt umber', ...
    'red ochre', 'Yellow Ochre', 'gypsum', 'chalk', 'Madder lake', ...
    'Realgar', 'Malachite', 'Orpiment', 'Indigo', 'Azurite', ...
    'Red lead', 'vermilion', 'Green earth', 'verdigris', 'lead white', ...
    'ultramarine', 'Naples Yellow', 'Smalt', 'indian yellow', ...
    'Copper resinate', 'lead tin yellow', 'Van Dyke brown', ...
    'Carmine lake', 'cobalt blue', 'raw umber', 'Cd Red', ...
    'titanium white', 'Cd yellow', 'zinc white'};
pk_dist_range = 50;    %nm
bw_dist_range = 100;    %bands
pkdelta_thresh1 = 0.02;
pkdelta_thresh2 = 0.02;
w = 35;

derivative_fn = [cube_path 'slope_' cube_fn];
derivative2_fn = [cube_path 'slope_slope_' cube_fn];
derivative3_fn = [cube_path 'slope_slope_slope_' cube_fn];
hsz = 3;
dthresh = 0.001;
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
    if (lambdaon == 0)
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
[N,M] = meshgrid(1:n0,1:m0);
Lp = m0;

% derivative features
%   reflectance peaks (zero-crossing), absorption peaks (zero-crossing), 
%   transition edges (peak), transition edge fwhm
% reflectance features
%   reflectance peak fwhm
rpeak_listt = single(zeros(m0*n0,5)); %[m,n,lambda,band]
rpeak_list2t = single(zeros(m0*n0,5)); %[m,n,lambda,band]
apeak_listt = single(zeros(m0*n0,5)); %[m,n,lambda]
apeak_list2t = single(zeros(m0*n0,5)); %[m,n,lambda]
[Ltotr,~] = size(rpeak_listt);
[Ltotr2,~] = size(rpeak_list2t);
[Ltota,~] = size(apeak_listt);
[Ltota2,~] = size(apeak_list2t);
ptotr = 1;
ptotr2 = 1;
ptota = 1;
ptota2 = 1;
for p = (hsz+ceil(hsz/2)):(p0-hsz-floor(hsz/2))
    tmp = der1(:,:,(p-1):(p+1));
    tmp1 = der2(:,:,(p-1):(p+1));
    tmp2 = der3(:,:,(p-1):(p+1));
    
    % find peaks in 1st derivative
    msk2p = (tmp1(:,:,1) >= 0) & (tmp1(:,:,3) <= 0) & (abs(tmp2(:,:,2)) > dthresh); %zero-crossing in 2nd derivative(+ -> -)
    %[r,c] = find(msk2p==1);
    r = M(msk2p);
    c = N(msk2p);
    L = length(r);
    while ((ptotr2 + L - 1) > Ltotr2)
        rpeak_list2t = [rpeak_list2t;single(zeros(L,5))];
        [Ltotr2,~] = size(rpeak_list2t);
    end
    for i = 1:L
        rpeak_list2t(ptotr2,:) = [r(i) c(i) lambda(p) p der1(r(i),c(i),p)];
        ptotr2 = ptotr2 + 1;
    end
    
    msk2n = (tmp1(:,:,1) <= 0) & (tmp1(:,:,3) >= 0) & (abs(tmp2(:,:,2)) > dthresh); %zero-crossing in 2nd derivative(- -> +)
    %[r,c] = find(msk2n==1);
    r = M(msk2n);
    c = N(msk2n);
    L = length(r);
    while ((ptota2 + L - 1) > Ltota2)
        apeak_list2t = [apeak_list2t;single(zeros(L,5))];
        [Ltota2,~] = size(apeak_list2t);
    end
    for i = 1:L
        apeak_list2t(ptota2,:) = [r(i) c(i) lambda(p) p der1(r(i),c(i),p)];
        ptota2 = ptota2 + 1;
    end
    
    % find peaks in data
    msk1p = (tmp(:,:,1) >= 0) & (tmp(:,:,3) <= 0) & (abs(tmp1(:,:,2)) > dthresh); %zero-crossing in 1st derivative(+ -> -)
    %[r,c] = find(msk1p==1);
    r = M(msk1p);
    c = N(msk1p);
    L = length(r);
    while ((ptotr + L - 1) > Ltotr)
        rpeak_listt = [rpeak_listt;single(zeros(L,5))];
        [Ltotr,~] = size(rpeak_listt);
    end
    for i = 1:L
        rpeak_listt(ptotr,:) = [r(i) c(i) lambda(p) p cube(r(i),c(i),p)];
        ptotr = ptotr + 1;
    end
    
    msk1n = (tmp(:,:,1) <= 0) & (tmp(:,:,3) >= 0) & (abs(tmp1(:,:,2)) > dthresh); %zero-crossing in 1st derivative(- -> +)
    %[r,c] = find(msk1n==1);
    r = M(msk1n);
    c = N(msk1n);
    L = length(r);
    while ((ptota + L - 1) > Ltota)
        apeak_listt = [apeak_listt;single(zeros(L,5))];
        [Ltota,~] = size(apeak_listt);
    end
    for i = 1:L
        apeak_listt(ptota,:) = [r(i) c(i) lambda(p) p cube(r(i),c(i),p)];
        ptota = ptota + 1;
    end
end
clear r c
rpeak_listt(rpeak_listt(:,1)==0,:) = [];
rpeak_list2t(rpeak_list2t(:,1)==0,:) = [];
apeak_listt(apeak_listt(:,1)==0,:) = [];
apeak_list2t(apeak_list2t(:,1)==0,:) = [];
[Ltotr1,~] = size(rpeak_listt);
[Ltotr2,~] = size(rpeak_list2t);
[Ltota1,~] = size(apeak_listt);
[Ltota2,~] = size(apeak_list2t);

rpeak_list = single(zeros(m0*n0,4)); %[m,n,lambda,band]
apeak_list = single(zeros(m0*n0,4)); %[m,n,lambda]
edge_list = single(zeros(m0*n0,4)); %[m,n,lambda,band]
[Ltotr,~] = size(rpeak_list);
[Ltota,~] = size(apeak_list);
[Ltote,~] = size(edge_list);
ptotr = 1;
ptota = 1;
ptote = 1;
for p = (hsz+ceil(hsz/2)):(p0-hsz-floor(hsz/2))
    p
    % find peaks in 1st derivative
    Lv = 1:Ltotr2;
    ind1 = (rpeak_list2t(:,4)==p);
    num = Lv(ind1);
    rp2 = rpeak_list2t(num,:);
    ind2 = (rpeak_list2t(:,4) == (p+1));
    num2 = Lv(ind2);
    rp2p = rpeak_list2t(num2,:);
    ind3 = (rpeak_list2t(:,4) == (p-1));
    num3 = Lv(ind3);
    rp2n = rpeak_list2t(num3,:);
    [L,~] = size(rp2);
    while ((ptote + L - 1) > Ltote)
        edge_list = [edge_list;single(zeros(L,4))];
        [Ltote,~] = size(edge_list);
    end
    while ((ptota + L - 1) > Ltota)
        apeak_list = [apeak_list;single(zeros(L,4))];
        [Ltota,~] = size(apeak_list);
    end
    for i = 1:L
        r = rp2(i,1);
        c = rp2(i,2);
        val = rp2(i,5);
        
        ind2 = (rp2p(:,1) == r) & (rp2p(:,2) == c);
        val2 = rp2p(ind2,5);
        ind3 = (rp2n(:,1) == r) & (rp2n(:,2) == c);
        val3 = rp2n(ind3,5);
        pass = 1;
        if (~isempty(val2))
            if (val2 < val)
                pass = 1;
            else
                pass = 0;
            end
        end
        if (~isempty(val3))
            if (val3 < val)
                pass = 1;
            else
                pass = 0;
            end
        end
        if (pass == 1)
            Lv = 1:Ltota2;
            ind0 = Lv((apeak_list2t(:,1) == r) & (apeak_list2t(:,2) == c));
            tmpLv = apeak_list2t(ind0,:);
            Lv = 1:size(tmpLv,1);
            ind = Lv(tmpLv(:,4) < p);
            [~,ind1] = max(tmpLv(ind,4) - p);
            if (~isempty(ind1))
                %v1 = der1(r,c,apeak_list2t(ind0(ind(ind1)),4));
                v1 = apeak_list2t(ind0(ind(ind1)),5);
            else
                v1 = NaN;
            end
            ind = Lv(tmpLv(:,4) > p);
            [~,ind1] = min(tmpLv(ind,4) - p);
            if (~isempty(ind1))
                %v2 = der1(r,c,apeak_list2t(ind0(ind(ind1)),4));
                v2 = apeak_list2t(ind0(ind(ind1)),5);
            else
                v2 = NaN;
            end
            %pkdelta1 = abs(v1 - der1(r,c,p));
            %pkdelta2 = abs(v2 - der1(r,c,p));
            pkdelta1 = abs(v1 - val);
            pkdelta2 = abs(v2 - val);
            pkdelta = max(pkdelta1,pkdelta2);
            
            if (val > 0)
                if (pkdelta >= pkdelta_thresh2)
                    % if the peak is greater than 0
                    edge_list(ptote,:) = [r c lambda(p) p];
                    ptote = ptote + 1;
                end
            elseif (val <= 0)
                if (pkdelta >= (pkdelta_thresh2/2))
                    apeak_list(ptota,:) = [r c lambda(p) p];
                    ptota = ptota + 1;
                end
            end
        end
    end
    clear Lv ind1 ind2 ind3 num num2 num3
    
    % find valleys in 1st derivative
    Lv = 1:Ltota2;
    ind1 = (apeak_list2t(:,4)==p);
    num = Lv(ind1);
    rp2 = apeak_list2t(num,:);
    ind2 = (apeak_list2t(:,4) == (p+1));
    num2 = Lv(ind2);
    rp2p = apeak_list2t(num2,:);
    ind3 = (apeak_list2t(:,4) == (p-1));
    num3 = Lv(ind3);
    rp2n = apeak_list2t(num3,:);
    [L,~] = size(rp2);
    while ((ptota + L - 1) > Ltota)
        apeak_list = [apeak_list;single(zeros(L,4))];
        [Ltota,~] = size(apeak_list);
    end
    for i = 1:L
        r = rp2(i,1);
        c = rp2(i,2);
        val = rp2(i,5);
        
        ind2 = (rp2p(:,1) == r) & (rp2p(:,2) == c);
        val2 = rp2p(ind2,5);
        ind3 = (rp2n(:,1) == r) & (rp2n(:,2) == c);
        val3 = rp2n(ind3,5);
        pass = 1;
        if (~isempty(val2))
            if (val2 > val)
                pass = 1;
            else
                pass = 0;
            end
        end
        if (~isempty(val3))
            if (val3 > val)
                pass = 1;
            else
                pass = 0;
            end
        end
        if (pass == 1)
            Lv = 1:Ltotr2;
            ind0 = Lv((rpeak_list2t(:,1) == r) & (rpeak_list2t(:,2) == c));
            tmpLv = rpeak_list2t(ind0,:);
            Lv = 1:size(tmpLv,1);
            ind = Lv(tmpLv(:,4) < p);
            [~,ind1] = max(tmpLv(ind,4) - p);
            if (~isempty(ind1))
                %p1 = der1(r,c,rpeak_list2t(ind0(ind(ind1)),4));
                p1 = rpeak_list2t(ind0(ind(ind1)),5);
            else
                p1 = NaN;
            end
            ind = Lv(tmpLv(:,4) > p);
            [~,ind1] = min(tmpLv(ind,4) - p);
            if (~isempty(ind1))
                %p2 = der1(r,c,rpeak_list2t(ind0(ind(ind1)),4));
                p2 = rpeak_list2t(ind0(ind(ind1)),5);
            else
                p2 = NaN;
            end
            %pkdelta1 = abs(p1 - der1(r,c,p));
            %pkdelta2 = abs(p2 - der1(r,c,p));
            pkdelta1 = abs(p1 - val);
            pkdelta2 = abs(p2 - val);
            pkdelta = max(pkdelta1,pkdelta2);
            
            if (val > 0)
                if (pkdelta >= (pkdelta_thresh2/2))
                    % if the peak is greater than 0
                    apeak_list(ptota,:) = [r c lambda(p) p];
                    ptota = ptota + 1;
                end
            end
        end
    end
    clear Lv ind1 ind2 ind3 num num2 num3
    
    % find peaks in data
    Lv = 1:Ltotr1;
    ind1 = (rpeak_listt(:,4)==p);
    num = Lv(ind1);
    rp2 = rpeak_listt(num,:);
    ind2 = (rpeak_listt(:,4) == (p+1));
    num2 = Lv(ind2);
    rp2p = rpeak_listt(num2,:);
    ind3 = (rpeak_listt(:,4) == (p-1));
    num3 = Lv(ind3);
    rp2n = rpeak_listt(num3,:);
    [L,~] = size(rp2);
    while ((ptotr + L - 1) > Ltotr)
        rpeak_list = [rpeak_list;single(zeros(L,4))];
        [Ltotr,~] = size(rpeak_list);
    end
    for i = 1:L
        r = rp2(i,1);
        c = rp2(i,2);
        val = rp2(i,5);
        
        ind2 = (rp2p(:,1) == r) & (rp2p(:,2) == c);
        val2 = rp2p(ind2,5);
        ind3 = (rp2n(:,1) == r) & (rp2n(:,2) == c);
        val3 = rp2n(ind3,5);
        pass = 1;
        if (~isempty(val2))
            if (val2 < val)
                pass = 1;
            else
                pass = 0;
            end
        end
        if (~isempty(val3))
            if (val3 < val)
                pass = 1;
            else
                pass = 0;
            end
        end
        if (pass == 1)
            Lv = 1:Ltota1;
            ind0 = Lv((apeak_listt(:,1) == r) & (apeak_listt(:,2) == c));
            tmpLv = apeak_listt(ind0,:);
            Lv = 1:size(tmpLv,1);
            ind = Lv(tmpLv(:,4) < p);
            [~,ind1] = max(tmpLv(ind,4) - p);
            if (~isempty(ind1))
                %v1 = der1(r,c,apeak_listt(ind0(ind(ind1)),4));
                v1 = apeak_listt(ind0(ind(ind1)),5);
            else
                v1 = NaN;
            end
            ind = Lv(tmpLv(:,4) > p);
            [~,ind1] = min(tmpLv(ind,4) - p);
            if (~isempty(ind1))
                %v2 = der1(r,c,apeak_listt(ind0(ind(ind1)),4));
                v2 = apeak_listt(ind0(ind(ind1)),5);
            else
                v2 = NaN;
            end
            %pkdelta1 = abs(v1 - cube(r,c,p));
            %pkdelta2 = abs(v2 - cube(r,c,p));
            pkdelta1 = abs(v1 - val);
            pkdelta2 = abs(v2 - val);
            pkdelta = max(pkdelta1,pkdelta2);
            
            if (pkdelta >= pkdelta_thresh1)
                % if the peak is greater than 0
                rpeak_list(ptotr,:) = [r c lambda(p) p];
                ptotr = ptotr + 1;
            end
        end
    end
    clear Lv ind1 ind2 ind3 num num2 num3
    
    % find valleys in data
    Lv = 1:Ltota1;
    ind1 = (apeak_listt(:,4)==p);
    num = Lv(ind1);
    rp2 = apeak_listt(num,:);
    ind2 = (apeak_listt(:,4) == (p+1));
    num2 = Lv(ind2);
    rp2p = apeak_listt(num2,:);
    ind3 = (apeak_listt(:,4) == (p-1));
    num3 = Lv(ind3);
    rp2n = apeak_listt(num3,:);
    [L,~] = size(rp2);
    while ((ptota + L - 1) > Ltota)
        apeak_list = [apeak_list;single(zeros(L,4))];
        [Ltota,~] = size(apeak_list);
    end
    for i = 1:L
        r = rp2(i,1);
        c = rp2(i,2);
        val = rp2(i,5);
        
        ind2 = (rp2p(:,1) == r) & (rp2p(:,2) == c);
        val2 = rp2p(ind2,5);
        ind3 = (rp2n(:,1) == r) & (rp2n(:,2) == c);
        val3 = rp2n(ind3,5);
        pass = 1;
        if (~isempty(val2))
            if (val2 > val)
                pass = 1;
            else
                pass = 0;
            end
        end
        if (~isempty(val3))
            if (val3 > val)
                pass = 1;
            else
                pass = 0;
            end
        end
        if (pass == 1)
            Lv = 1:Ltotr1;
            ind0 = Lv((rpeak_listt(:,1) == r) & (rpeak_listt(:,2) == c));
            tmpLv = rpeak_listt(ind0,:);
            Lv = 1:size(tmpLv,1);
            ind = Lv(tmpLv(:,4) < p);
            [~,ind1] = max(tmpLv(ind,4) - p);
            if (~isempty(ind1))
                %p1 = der1(r,c,rpeak_listt(ind0(ind(ind1)),4));
                p1 = rpeak_listt(ind0(ind(ind1)),5);
            else
                p1 = NaN;
            end
            ind = Lv(tmpLv(:,4) > p);
            [~,ind1] = min(tmpLv(ind,4) - p);
            if (~isempty(ind1))
                %p2 = der1(r,c,rpeak_listt(ind0(ind(ind1)),4));
                p2 = rpeak_listt(ind0(ind(ind1)),5);
            else
                p2 = NaN;
            end
            %pkdelta1 = abs(p1 - cube(r,c,p));
            %pkdelta2 = abs(p2 - cube(r,c,p));
            pkdelta1 = abs(p1 - val);
            pkdelta2 = abs(p2 - val);
            pkdelta = max(pkdelta1,pkdelta2);
            
            if (pkdelta >= pkdelta_thresh1/2)
                % if the peak is greater than 0
                apeak_list(ptota,:) = [r c lambda(p) p];
                ptota = ptota + 1;
            end
        end
    end
    clear Lv ind1 ind2 ind3 num num2 num3
end
clear r c val
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
    
    j = 2;
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
for i = 1:Lp
    j = 3;
    tmpBW = pigment_ris_table_band(i,j-1);
    if (~isempty(tmpBW))
        tmpBW1 = regexp(tmpBW, '\;', 'split');
        tmpBW1 = tmpBW1{1};
        strtmp = '';
        strtmp2 = '';
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
            cube_tmp = der1(i,1,pk1:pk2);
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
                %bw1 = 2*abs(zc1-lambda_tmp(pk0));
                bw1 = abs(zc1-lambda_tmp(pk0));
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
                %bw2 = 2*abs(zc2-lambda_tmp(pk0));
                bw2 = abs(zc2-lambda_tmp(pk0));
            end
            %bw = min(bw1,bw2);
            bw = bw1 + bw2;
            asym = bw2 - bw1;
            %{
            if (bw > (2*bw_dist_range))
                bw = Inf;
            end
            %}
            if (isinf(bw))
                bw = 0;
                asym = 0;
            end
            if (isempty(strtmp))
                strtmp = num2str(bw);
                strtmp2 = num2str(asym);
            else
                strtmp = [strtmp ';' num2str(bw)];
                strtmp2 = [strtmp2 ';' num2str(asym)];
            end
            pigment_ris_table{i,j} = strtmp;
            pigment_ris_table{i,j+1} = strtmp2;
        end
    end
end
clear tmp

% strong absorption features
tmp = repmat(mean(cube,3),[1 1 p0]);
msk = (abs(der1) < 0.01) & (cube < tmp);
clear tmp
msk1 = uint16(zeros(m0,n0,p0));
[X,Y] = meshgrid(1:n0,1:m0);
abs_features = zeros(m0*n0*10,4);
[Labsf,~] = size(abs_features);
pnt = 1;
for p = 2:(p0-1)
    tmp = (msk1(:,:,(p-1)) + uint16(msk(:,:,p))).*uint16(msk(:,:,p));
    msk1(:,:,p) = tmp;
    tmp1 = single(tmp).*single((tmp > 0) & (msk(:,:,(p+1)) == 0));
    tmp2 = (tmp1 == 0);
    tmp1(tmp2) = 1;
    abs_dist0 = lambda(p) - lambda(p - tmp1 + 1);
    abs_dist0(tmp2) = 0;
    abs_pos0 = lambda(round(p - tmp1/2));
    abs_pos0(tmp2) = 0;
    tmp3 = [X(:) Y(:) abs_dist0(:) abs_pos0(:)];
    clear tmp tmp1 tmp2 abs_dist0 abs_pos0
    tmp3(tmp3(:,3) == 0,:) = [];
    [L3,~] = size(tmp3);
    while ((pnt + L3 - 1) > Labsf)
        abs_features = [abs_features;zeros(m0*n0*10,4)];
        [Labsf,~] = size(abs_features);
    end
    abs_features(pnt:(pnt+L3-1),:) = tmp3;
    pnt = pnt + L3;
    clear tmp3
end
clear msk msk1
abs_features(abs_features(:,1) == 0,:) = [];
%remove features of regions < 50 nm wide
abs_features(abs_features(:,3) < 50,:) = [];

% gradual (unspecified value) slope feature
%{
slope_diff_thresh = 0.005;
msk1 = uint16(zeros(m0,n0,p0));
slope_features = zeros(m0*n0*10,4);
[Lslopef,~] = size(slope_features);
pnt = 1;
for p = 2:(p0-1)
    p
    tmp = abs(der1(:,:,p) - der1(:,:,(p-1)));
    tmp = (tmp < slope_diff_thresh);
    tmpe = abs(der1(:,:,(p+1)) - der1(:,:,p));
    tmpe = (tmpe < slope_diff_thresh);
    msk1(:,:,p) = (msk1(:,:,(p-1)) + uint16(tmp)).*uint16(tmp);
    
    tmp1 = single(msk1(:,:,p)).*single((msk1(:,:,p) > 0) & (tmpe == 0));
    tmp2 = (tmp1 == 0);
    tmp1(tmp2) = 1;
    slope_dist0 = lambda(p) - lambda(p - tmp1 + 1);
    slope_dist0(tmp2) = 0;
    slope_pos0 = lambda(round(p - tmp1/2));
    slope_pos0(tmp2) = 0;
    tmp3 = [X(:) Y(:) slope_dist0(:) slope_pos0(:)];
    clear tmp tmpe tmp1 tmp2 slope_dist0 slope_pos0
    tmp3(tmp3(:,3) == 0,:) = [];
    [L3,~] = size(tmp3);
    while ((pnt + L3 - 1) > Lslopef)
        slope_features = [slope_features;zeros(m0*n0*10,4)];
        [Lsloef,~] = size(slope_features);
    end
    slope_features(pnt:(pnt+L3-1),:) = tmp3;
    pnt = pnt + L3;
    clear tmp3
end
clear msk1 X Y
slope_features(slope_features(:,1) == 0,:) = [];
%remove features of regions < 200 nm wide
slope_features(slope_features(:,3) < 200,:) = [];
%}
%{
msk = ((der1 > 0) & (abs(der2) < 0.01));
msk1 = uint16(zeros(m0,n0,p0));
[X,Y] = meshgrid(1:n0,1:m0);
slope_features = zeros(m0*n0*10,4);
[Labsf,~] = size(slope_features);
pnt = 1;
for p = 2:(p0-1)
    msk01 = (msk(:,:,p) & (sign(der1(:,:,(p-1))) == sign(der1(:,:,p))));
    msk02 = (msk(:,:,p) & (sign(der1(:,:,(p+1))) == sign(der1(:,:,p))));
    tmp = (msk1(:,:,(p-1)) + uint16(msk01)).*uint16(msk01);
    clear msk01
    msk1(:,:,p) = tmp;
    tmp1 = single(tmp).*single((tmp > 0) & (msk02 == 0));
    clear msk02
    tmp2 = (tmp1 == 0);
    tmp1(tmp2) = 1;
    abs_dist0 = lambda(p) - lambda(p - tmp1 + 1);
    abs_dist0(tmp2) = 0;
    abs_pos0 = lambda(round(p - tmp1/2));
    abs_pos0(tmp2) = 0;
    tmp3 = [X(:) Y(:) abs_dist0(:) abs_pos0(:)];
    clear tmp tmp1 tmp2 abs_dist0 abs_pos0
    tmp3(tmp3(:,3) == 0,:) = [];
    [L3,~] = size(tmp3);
    while ((pnt + L3 - 1) > Labsf)
        slope_features = [slope_features;zeros(m0*n0*10,4)];
        [Labsf,~] = size(slope_features);
    end
    slope_features(pnt:(pnt+L3-1),:) = tmp3;
    pnt = pnt + L3;
    clear tmp3
end
clear msk msk1
slope_features(slope_features(:,1) == 0,:) = [];
%remove features of regions < 200 nm wide
slope_features(slope_features(:,3) < 200,:) = [];
%}

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
toc

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
            plot(str2double(tmp1(k)),lib1(ind(1)),'kx')
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
            plot(str2double(tmp1(k)),lib1(ind(1)),'ko')
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
    j = 2;
    tmp = pigment_ris_table(i,j);
    if (~isempty(tmp))
        tmp1 = regexp(tmp, '\;', 'split');
        tmp1 = tmp1{1};
        for k = 1:tot_ris_table(i,j)
            [~,ind] = min(abs(lambda - str2double(tmp1(k))));
            plot(str2double(tmp1(k)),lib_der11(ind(1)),'kx')
            hold on
        end
    end
    hold off

    subplot(313),plot(lambda,lib_der21)
    axis([425 2500 min(lib_der21) max(lib_der21)])

    %{
    muremov = data_muremov(i,:,:);
    muremov1 = muremov(:);
    subplot(414),plot(lambda,muremov1)
    axis([425 2500 min(muremov1) max(muremov1)])
    %}
    %{
    dataA1 = dataA(i,:,:);
    dataA1 = dataA1(:);
    subplot(413),plot(lambda,dataA1)
    axis([425 2500 min(dataA1) max(dataA1)])
    dataA1 = der1A(i,:,:);
    dataA1 = dataA1(:);
    subplot(414),plot(lambda,dataA1)
    axis([425 2500 min(dataA1) max(dataA1)])
    %}
    %pause
end
