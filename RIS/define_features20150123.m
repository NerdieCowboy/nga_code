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
%dark_thresh = 0.05;
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
    
    % find peaks in 1st derivative
    msk2p = (tmp1(:,:,1) >= 0) & (tmp1(:,:,3) <= 0); %zero-crossing in 2nd derivative(+ -> -)
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
    
    msk2n = (tmp1(:,:,1) <= 0) & (tmp1(:,:,3) >= 0); %zero-crossing in 2nd derivative(- -> +)
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
    msk1p = (tmp(:,:,1) >= 0) & (tmp(:,:,3) <= 0); %zero-crossing in 1st derivative(+ -> -)
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
    
    msk1n = (tmp(:,:,1) <= 0) & (tmp(:,:,3) >= 0); %zero-crossing in 1st derivative(- -> +)
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
    pkdelta_thresh1 = 0.05;
    pkdelta_thresh2 = 0.02;
    
    ind1 = (rpeak_listt(:,4) == p);
    r1 = rpeak_listt(ind1,1);
    c1 = rpeak_listt(ind1,2);
    v1 = rpeak_listt(ind1,5);
    L = length(r1);
    Lv = 1:Ltotr1;
    num = Lv(ind1);
    for i = 1:L
        ind2 = (rpeak_listt(:,4) == (p+1)) & (rpeak_listt(:,1) == r1(i)) & (rpeak_listt(:,2) == c1(i));
        v2 = rpeak_listt(ind2,5);
        ind3 = (rpeak_listt(:,4) == (p-1)) & (rpeak_listt(:,1) == r1(i)) & (rpeak_listt(:,2) == c1(i));
        v3 = rpeak_listt(ind3,5);
        if (~isempty(v2))
            if (v2 >= v1(i))
                rpeak_listt(num(i),:) = zeros(1,5);
            end
        end
        if (~isempty(v3))
            if (v3 >= v1(i))
                rpeak_listt(num(i),:) = zeros(1,5);
            end
        end
    end
    clear Lv
    
    ind1 = (apeak_listt(:,4) == p);
    r1 = apeak_listt(ind1,1);
    c1 = apeak_listt(ind1,2);
    v1 = apeak_listt(ind1,5);
    L = length(r1);
    Lv = 1:Ltota1;
    num = Lv(ind1);
    for i = 1:L
        ind2 = (apeak_listt(:,4) == (p+1)) & (apeak_listt(:,1) == r1(i)) & (apeak_listt(:,2) == c1(i));
        v2 = apeak_listt(ind2,5);
        ind3 = (apeak_listt(:,4) == (p-1)) & (apeak_listt(:,1) == r1(i)) & (apeak_listt(:,2) == c1(i));
        v3 = apeak_listt(ind3,5);
        if (~isempty(v2))
            if (v2 <= v1(i))
                apeak_listt(num(i),:) = zeros(1,5);
            end
        end
        if (~isempty(v3))
            if (v3 <= v1(i))
                apeak_listt(num(i),:) = zeros(1,5);
            end
        end
    end
    clear Lv
    
    ind1 = (rpeak_list2t(:,4) == p);
    r1 = rpeak_list2t(ind1,1);
    c1 = rpeak_list2t(ind1,2);
    v1 = rpeak_list2t(ind1,5);
    L = length(r1);
    Lv = 1:Ltotr2;
    num = Lv(ind1);
    for i = 1:L
        ind2 = (rpeak_list2t(:,4) == (p+1)) & (rpeak_list2t(:,1) == r1(i)) & (rpeak_list2t(:,2) == c1(i));
        v2 = rpeak_list2t(ind2,5);
        ind3 = (rpeak_list2t(:,4) == (p-1)) & (rpeak_list2t(:,1) == r1(i)) & (rpeak_list2t(:,2) == c1(i));
        v3 = rpeak_list2t(ind3,5);
        if (~isempty(v2))
            if (v2 >= v1(i))
                rpeak_list2t(num(i),:) = zeros(1,5);
            end
        end
        if (~isempty(v3))
            if (v3 >= v1(i))
                rpeak_list2t(num(i),:) = zeros(1,5);
            end
        end
    end
    clear Lv
    
    ind1 = (apeak_list2t(:,4) == p);
    r1 = apeak_list2t(ind1,1);
    c1 = apeak_list2t(ind1,2);
    v1 = apeak_list2t(ind1,5);
    L = length(r1);
    Lv = 1:Ltota2;
    num = Lv(ind1);
    for i = 1:L
        ind2 = (apeak_list2t(:,4) == (p+1)) & (apeak_list2t(:,1) == r1(i)) & (apeak_list2t(:,2) == c1(i));
        v2 = apeak_list2t(ind2,5);
        ind3 = (apeak_list2t(:,4) == (p-1)) & (apeak_list2t(:,1) == r1(i)) & (apeak_list2t(:,2) == c1(i));
        v3 = apeak_list2t(ind3,5);
        if (~isempty(v2))
            if (v2 <= v1(i))
                apeak_list2t(num(i),:) = zeros(1,5);
            end
        end
        if (~isempty(v3))
            if (v3 <= v1(i))
                apeak_list2t(num(i),:) = zeros(1,5);
            end
        end
    end
    clear Lv
    
    % find peaks in 1st derivative
    %num = find(rpeak_list2t(:,4)==p);
    Lv = 1:Ltotr2;
    num = Lv((rpeak_list2t(:,4)==p) & (rpeak_list2t(:,5) > 0));
    clear Lv
    L = length(num);
    while ((ptote + L - 1) > Ltote)
        edge_list = [edge_list;single(zeros(L,4))];
        [Ltote,~] = size(edge_list);
    end
    for i = 1:L
        r = rpeak_list2t(num(i),1);
        c = rpeak_list2t(num(i),2);
        
        %ind = find((apeak_list2t(:,1) == r) & (apeak_list2t(:,2) == c) & (apeak_list2t(:,4) < p));
        Lv = 1:Ltota2;
        ind = Lv((apeak_list2t(:,1) == r) & (apeak_list2t(:,2) == c) & (apeak_list2t(:,4) < p));
        clear Lv
        [~,ind1] = max(apeak_list2t(ind,4) - p);
        if (~isempty(ind1))
            v1 = der1(r,c,apeak_list2t(ind(ind1),4));
        else
            v1 = NaN;
        end
        %ind = find((apeak_list2t(:,1) == r) & (apeak_list2t(:,2) == c) & (apeak_list2t(:,4) > p));
        [L1,~] = size(apeak_list2t);
        Lv = 1:L1;
        ind = Lv((apeak_list2t(:,1) == r) & (apeak_list2t(:,2) == c) & (apeak_list2t(:,4) > p));
        clear Lv
        [~,ind1] = min(apeak_list2t(ind,4) - p);
        if (~isempty(ind1))
            v2 = der1(r,c,apeak_list2t(ind(ind1),4));
        else
            v2 = NaN;
        end
        
        pkdelta1 = abs(v1 - der1(r,c,p));
        pkdelta2 = abs(v2 - der1(r,c,p));
        pkdelta = max(pkdelta1,pkdelta2);
        if (pkdelta >= pkdelta_thresh2)
            % if the peak is greater than 0
            edge_list(ptote,:) = [r c lambda(p) p];
            ptote = ptote + 1;
        end
    end
    
    Lv = 1:Ltotr2;
    num = Lv((rpeak_list2t(:,4)==p) & (rpeak_list2t(:,5) <= 0));
    clear Lv
    L = length(num);
    while ((ptota + L - 1) > Ltota)
        apeak_list = [apeak_list;single(zeros(L,4))];
        [Ltota,~] = size(apeak_list);
    end
    for i = 1:L
        r = rpeak_list2t(num(i),1);
        c = rpeak_list2t(num(i),2);
        
        %ind = find((apeak_list2t(:,1) == r) & (apeak_list2t(:,2) == c) & (apeak_list2t(:,4) < p));
        Lv = 1:Ltota2;
        ind = Lv((apeak_list2t(:,1) == r) & (apeak_list2t(:,2) == c) & (apeak_list2t(:,4) < p));
        clear Lv
        [~,ind1] = max(apeak_list2t(ind,4) - p);
        if (~isempty(ind1))
            v1 = der1(r,c,apeak_list2t(ind(ind1),4));
        else
            v1 = NaN;
        end
        %ind = find((apeak_list2t(:,1) == r) & (apeak_list2t(:,2) == c) & (apeak_list2t(:,4) > p));
        [L1,~] = size(apeak_list2t);
        Lv = 1:L1;
        ind = Lv((apeak_list2t(:,1) == r) & (apeak_list2t(:,2) == c) & (apeak_list2t(:,4) > p));
        clear Lv
        [~,ind1] = min(apeak_list2t(ind,4) - p);
        if (~isempty(ind1))
            v2 = der1(r,c,apeak_list2t(ind(ind1),4));
        else
            v2 = NaN;
        end
        
        pkdelta1 = abs(v1 - der1(r,c,p));
        pkdelta2 = abs(v2 - der1(r,c,p));
        pkdelta = max(pkdelta1,pkdelta2);
        if (pkdelta >= (pkdelta_thresh2/2))
            apeak_list(ptota,:) = [r c lambda(p) p];
            ptota = ptota + 1;
        end
    end
    
    %num = find(apeak_list2t(:,4)==p);
    Lv = 1:Ltota2;
    num = Lv((apeak_list2t(:,4)==p) & (apeak_list2t(:,5)>0));
    clear Lv
    L = length(num);
    while ((ptota + L - 1) > Ltota)
        apeak_list = [apeak_list;single(zeros(L,4))];
        [Ltota,~] = size(apeak_list);
    end
    for i = 1:L
        r = apeak_list2t(num(i),1);
        c = apeak_list2t(num(i),2);
        
        %ind = find((rpeak_list2t(:,1) == r) & (rpeak_list2t(:,2) == c) & (rpeak_list2t(:,4) < p));
        [L1,~] = size(rpeak_list2t);
        Lv = 1:L1;
        ind = Lv((rpeak_list2t(:,1) == r) & (rpeak_list2t(:,2) == c) & (rpeak_list2t(:,4) < p));
        clear Lv
        [~,ind1] = max(rpeak_list2t(ind,4) - p);
        if (~isempty(ind1))
            p1 = der1(r,c,rpeak_list2t(ind(ind1),4));
        else
            p1 = NaN;
        end
        %ind = find((rpeak_list2t(:,1) == r) & (rpeak_list2t(:,2) == c) & (rpeak_list2t(:,4) > p));
        [L1,~] = size(rpeak_list2t);
        Lv = 1:L1;
        ind = Lv((rpeak_list2t(:,1) == r) & (rpeak_list2t(:,2) == c) & (rpeak_list2t(:,4) > p));
        clear Lv
        [~,ind1] = min(rpeak_list2t(ind,4) - p);
        if (~isempty(ind1))
            p2 = der1(r,c,rpeak_list2t(ind(ind1),4));
        else
            p2 = NaN;
        end
        
        pkdelta1 = abs(p1 - der1(r,c,p));
        pkdelta2 = abs(p2 - der1(r,c,p));
        pkdelta = max(pkdelta1,pkdelta2);
        if (pkdelta >= (pkdelta_thresh2/2))
            % if the valley is greater than 0
            apeak_list(ptota,:) = [r c lambda(p) p];
            ptota = ptota + 1;
        end
    end
    
    % find peaks in data
    %num = find(rpeak_listt(:,4)==p);
    Lv = 1:Ltotr1;
    num = Lv(rpeak_listt(:,4)==p);
    clear Lv
    L = length(num);
    while ((ptotr + L - 1) > Ltotr)
        rpeak_list = [rpeak_list;single(zeros(L,4))];
        [Ltotr,~] = size(rpeak_list);
    end
    for i = 1:L
        r = rpeak_listt(num(i),1);
        c = rpeak_listt(num(i),2);
        
        %ind = find((apeak_listt(:,1) == r) & (apeak_listt(:,2) == c) & (apeak_listt(:,4) < p));
        [L1,~] = size(apeak_listt);
        Lv = 1:L1;
        ind = Lv((apeak_listt(:,1) == r) & (apeak_listt(:,2) == c) & (apeak_listt(:,4) < p));
        clear Lv
        [~,ind1] = max(apeak_listt(ind,4) - p);
        if (~isempty(ind1))
            v1 = cube(r,c,apeak_listt(ind(ind1),4));
        else
            v1 = NaN;
        end
        %ind = find((apeak_listt(:,1) == r) & (apeak_listt(:,2) == c) & (apeak_listt(:,4) > p));
        [L1,~] = size(apeak_listt);
        Lv = 1:L1;
        ind = Lv((apeak_listt(:,1) == r) & (apeak_listt(:,2) == c) & (apeak_listt(:,4) > p));
        clear Lv
        [~,ind1] = min(apeak_listt(ind,4) - p);
        if (~isempty(ind1))
            v2 = cube(r,c,apeak_listt(ind(ind1),4));
        else
            v2 = NaN;
        end
        
        pkdelta1 = abs(v1 - cube(r,c,p));
        pkdelta2 = abs(v2 - cube(r,c,p));
        pkdelta = max(pkdelta1,pkdelta2);
        if (pkdelta >= pkdelta_thresh1)
            rpeak_list(ptotr,:) = [r c lambda(p) p];
            ptotr = ptotr + 1;
        end
    end
    
    %num = find(apeak_listt(:,4)==p);
    Lv = 1:Ltota1;
    num = Lv(apeak_listt(:,4)==p);
    clear Lv
    L = length(num);
    while ((ptota + L - 1) > Ltota)
        apeak_list = [apeak_list;single(zeros(L,4))];
        [Ltota,~] = size(apeak_list);
    end
    for i = 1:L
        r = apeak_listt(num(i),1);
        c = apeak_listt(num(i),2);
        
        %ind = find((rpeak_listt(:,1) == r) & (rpeak_listt(:,2) == c) & (rpeak_listt(:,4) < p));
        [L1,~] = size(rpeak_listt);
        Lv = 1:L1;
        ind = Lv((rpeak_listt(:,1) == r) & (rpeak_listt(:,2) == c) & (rpeak_listt(:,4) < p));
        clear Lv
        [~,ind1] = max(rpeak_listt(ind,4) - p);
        if (~isempty(ind1))
            p1 = cube(r,c,rpeak_listt(ind(ind1),4));
        else
            p1 = NaN;
        end
        %ind = find((rpeak_listt(:,1) == r) & (rpeak_listt(:,2) == c) & (rpeak_listt(:,4) > p));
        [L1,~] = size(rpeak_listt);
        Lv = 1:L1;
        ind = Lv((rpeak_listt(:,1) == r) & (rpeak_listt(:,2) == c) & (rpeak_listt(:,4) > p));
        clear Lv
        [~,ind1] = min(rpeak_listt(ind,4) - p);
        if (~isempty(ind1))
            p2 = cube(r,c,rpeak_listt(ind(ind1),4));
        else
            p2 = NaN;
        end
        
        pkdelta1 = abs(p1 - cube(r,c,p));
        pkdelta2 = abs(p2 - cube(r,c,p));
        pkdelta = max(pkdelta1,pkdelta2);
        if (pkdelta >= pkdelta_thresh1)
            apeak_list(ptota,:) = [r c lambda(p) p];
            ptota = ptota + 1;
        end
    end
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
