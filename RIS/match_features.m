function [features,tot_ris_table,Lp,Lf,pigment_ris_table,missing_penalty] = match_features(cube,der1,der2,der3,cube_path,features_table_fn,spectra_lib,spectra_names,pk_dist_range,bw_dist_range,abs_dist_range,pkdelta_thresh1,pkdelta_thresh2,lib_test)

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

lib_derivative_fn = [cube_path 'slope_' spectra_lib];
lib_derivative2_fn = [cube_path 'slope_slope_' spectra_lib];
lib_derivative3_fn = [cube_path 'slope_slope_slope_' spectra_lib];
hsz = 3;
dthresh = 0.001;
dark_thresh = 0.05;
wiggle = 5; %bands

% spectral library
fn = [cube_path spectra_lib];
lambdaon = 0;
lambda = [];
fnh = [fn '.hdr'];
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
lib = multibandread(fn,[m01,n01,p01],datatype,headeroffset,interleave,byteorder);
lib_der1 = multibandread(lib_derivative_fn,[m01,n01,p01],datatype,headeroffset,interleave,byteorder);
lib_der2 = multibandread(lib_derivative2_fn,[m01,n01,p01],datatype,headeroffset,interleave,byteorder);
lib_der3 = multibandread(lib_derivative3_fn,[m01,n01,p01],datatype,headeroffset,interleave,byteorder);
[m0,n0,p0] = size(cube);
[N,M] = meshgrid(1:n0,1:m0);

pigment_element_table = csvread(features_table_fn,1,xrf_cols(1)-1);

% spectral features for each pigment
fid = fopen(features_table_fn);
features_table = textscan(fid,'%s','delimiter',',','EmptyValue',NaN);
fclose(fid);
tmp = reshape(features_table{1},[xrf_cols(2),33]);
tmp = tmp';
pigment_family = tmp(2:33,pigment_family_col);
pigments = tmp(2:33,pigment_col);
colors = tmp(2:33,color_col);
pigment_ris_table = tmp(2:33,reflectance_col:interesting_col);
clear tmp
[Lp,Lf] = size(pigment_ris_table);

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
    %msk2p = (tmp1(:,:,1) >= 0) & (tmp1(:,:,3) <= 0) & (abs(tmp2(:,:,2)) > (pkdelta_thresh2/8)); %zero-crossing in 2nd derivative(+ -> -)
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
    
    %msk2n = (tmp1(:,:,1) <= 0) & (tmp1(:,:,3) >= 0) & (abs(tmp2(:,:,2)) > (pkdelta_thresh2/8)); %zero-crossing in 2nd derivative(- -> +)
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
    %msk1p = (tmp(:,:,1) >= 0) & (tmp(:,:,3) <= 0) & (abs(tmp1(:,:,2)) > (pkdelta_thresh1/8)); %zero-crossing in 1st derivative(+ -> -)
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
    
    %msk1n = (tmp(:,:,1) <= 0) & (tmp(:,:,3) >= 0) & (abs(tmp1(:,:,2)) > (pkdelta_thresh1/8)); %zero-crossing in 1st derivative(- -> +)
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

rpeak_list = single(zeros(m0*n0,5)); %[m,n,lambda,band,depth]
apeak_list = single(zeros(m0*n0,5)); %[m,n,lambda,band,depth]
edge_list = single(zeros(m0*n0,5)); %[m,n,lambda,band,depth]
[Ltotr,~] = size(rpeak_list);
[Ltota,~] = size(apeak_list);
[Ltote,~] = size(edge_list);
ptotr = 1;
ptota = 1;
ptote = 1;
for p = (hsz+ceil(hsz/2)):(p0-hsz-floor(hsz/2))
    % find peaks in 1st derivative
    %Lv = 1:Ltotr2;
    ind1 = (rpeak_list2t(:,4)==p);
    %num = Lv(ind1);
    %rp2 = rpeak_list2t(num,:);
    rp2 = rpeak_list2t(ind1,:);
    ind2 = (rpeak_list2t(:,4) == (p+1));
    %num2 = Lv(ind2);
    %rp2p = rpeak_list2t(num2,:);
    rp2p = rpeak_list2t(ind2,:);
    ind3 = (rpeak_list2t(:,4) == (p-1));
    %num3 = Lv(ind3);
    %rp2n = rpeak_list2t(num3,:);
    rp2n = rpeak_list2t(ind3,:);
    [L,~] = size(rp2);
    while ((ptote + L - 1) > Ltote)
        edge_list = [edge_list;single(zeros(L,5))];
        [Ltote,~] = size(edge_list);
    end
    while ((ptota + L - 1) > Ltota)
        apeak_list = [apeak_list;single(zeros(L,5))];
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
                    edge_list(ptote,:) = [r c lambda(p) p pkdelta];
                    ptote = ptote + 1;
                end
            elseif (val <= 0)
                if (pkdelta >= (pkdelta_thresh2/2))
                    apeak_list(ptota,:) = [r c lambda(p) p pkdelta];
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
        apeak_list = [apeak_list;single(zeros(L,5))];
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
                    apeak_list(ptota,:) = [r c lambda(p) p pkdelta];
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
        rpeak_list = [rpeak_list;single(zeros(L,5))];
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
                rpeak_list(ptotr,:) = [r c lambda(p) p pkdelta];
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
        apeak_list = [apeak_list;single(zeros(L,5))];
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
                apeak_list(ptota,:) = [r c lambda(p) p pkdelta];
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
    tmp3 = [Y(:) X(:) abs_dist0(:) abs_pos0(:)];
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
features = Inf(1024,7);
[fsz,~] = size(features);
pklist = zeros(1024,1);
pnt = 1;
Lf_list = [reflectance_col-3 transition_edge_col-3 absorption_col-3 strong_absorption_center_col-3];
Lf_list2 = [gradual_slope_col-3 interesting_col-3];
for m = 1:m0
    indj1m = (rpeak_list(:,1) == m);
    indj3m = (edge_list(:,1) == m);
    indj5m = (apeak_list(:,1) == m);
    ind_absm = (abs_features(:,1) == m);
    rtmpj1m = rpeak_list(indj1m,1:5);
    rtmpj3m = edge_list(indj3m,1:5);
    rtmpj5m = apeak_list(indj5m,1:5);
    rtmp_absm = abs_features(ind_absm,:);
    tmpm = cube(m,:,:);
    tmpm = reshape(tmpm,[n0,p0]);
    tmpm1 = der1(m,:,:);
    tmpm1 = reshape(tmpm1,[n0,p0]);
    for n = 1:n0
        indj1 = (rtmpj1m(:,2) == n);
        indj3 = (rtmpj3m(:,2) == n);
        indj5 = (rtmpj5m(:,2) == n);
        ind_abs = (rtmp_absm(:,2) == n);
        rtmpj1 = rtmpj1m(indj1,1:5);
        rtmpj3 = rtmpj3m(indj3,1:5);
        rtmpj5 = rtmpj5m(indj5,1:5);
        rtmp_abs = rtmp_absm(ind_abs,:);
        if ((sum(indj1) > 0)||(sum(indj3) > 0)||(sum(indj5) > 0)||(sum(ind_abs) > 0))
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
                            
                            if (j == (strong_absorption_center_col-3))
                                if (sum(ind_abs) > 0)
                                    tmpWidth = pigment_ris_table(i,j+1);
                                    tmpWidth = tmpWidth{1};
                                    tmp2 = regexp(tmpWidth, '\;', 'split');
                                    mdw = single(str2double(tmp2{k}));
                                    rdist = abs(rtmp_abs(:,4) - md);
                                    wdist = abs(rtmp_abs(:,3) - mdw);
                                    
                                    [mnq,q] = min(rdist);
                                    if (mnq(1) > abs_dist_range)
                                        mnq(1) = Inf;
                                        wdist(q(1)) = Inf;
                                    end
                                    features(pnt,:) = [m n i j k mnq(1) rtmp_abs(q(1),3)];
                                    pnt = pnt + 1;
                                    if (pnt > fsz)
                                        features = [features;Inf(1024,7)];
                                        pklist = [pklist;zeros(1024,1)];
                                    end
                                    [fsz,~] = size(features);
                                    features(pnt,:) = [m n i j+1 k wdist(q(1)) rtmp_abs(q(1),3)];
                                    pnt = pnt + 1;
                                    if (pnt > fsz)
                                        features = [features;Inf(1024,7)];
                                        pklist = [pklist;zeros(1024,1)];
                                    end
                                    [fsz,~] = size(features);
                                end
                            elseif (j == (reflectance_col-3))
                                if (sum(indj1) > 0)
                                    rdist = abs(rtmpj1(:,3) - md);
                                    %[i j k]
                                    [mnq,q] = min(rdist);
                                    if (mnq(1) > pk_dist_range)
                                        mnq(1) = Inf;
                                    end
                                    features(pnt,:) = [m n i j k mnq(1) rtmpj1(q(1),5)];
                                    pnt = pnt + 1;
                                    if (pnt > fsz)
                                        features = [features;Inf(1024,7)];
                                        pklist = [pklist;zeros(1024,1)];
                                    end
                                    [fsz,~] = size(features);
                                end
                            elseif (j == (transition_edge_col-3))
                                if (sum(indj3) > 0)
                                    rdist = abs(rtmpj3(:,3) - md);
                                    %[i j k]
                                    [mnq,q] = min(rdist);
                                    if (mnq(1) > pk_dist_range)
                                        mnq(1) = Inf;
                                    end
                                    features(pnt,:) = [m n i j k mnq(1) rtmpj3(q(1),5)];
                                    pnt = pnt + 1;
                                    if (pnt > fsz)
                                        features = [features;Inf(1024,7)];
                                        pklist = [pklist;zeros(1024,1)];
                                    end
                                    [fsz,~] = size(features);
                                    features(pnt,:) = [m n i j+1 k Inf rtmpj3(q(1),5)];
                                    pklist(pnt) = rtmpj3(q(1),4);
                                    pnt = pnt + 1;
                                    if (pnt > fsz)
                                        features = [features;Inf(1024,7)];
                                        pklist = [pklist;zeros(1024,1)];
                                    end
                                    [fsz,~] = size(features);
                                    features(pnt,:) = [m n i j+2 k Inf rtmpj3(q(1),5)];
                                    pnt = pnt + 1;
                                    if (pnt > fsz)
                                        features = [features;Inf(1024,7)];
                                        pklist = [pklist;zeros(1024,1)];
                                    end
                                end
                            elseif (j == (absorption_col-3))
                                if (sum(indj5) > 0)
                                    rdist = abs(rtmpj5(:,3) - md);
                                    %[i j k]
                                    [mnq,q] = min(rdist);
                                    if (mnq(1) > pk_dist_range)
                                        mnq(1) = Inf;
                                    end
                                    features(pnt,:) = [m n i j k mnq(1) rtmpj5(q(1),5)];
                                    pnt = pnt + 1;
                                    if (pnt > fsz)
                                        features = [features;Inf(1024,7)];
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
                        
                        if (j == 8)
                            %gradual slope regions
                            F1 = tmpm(n,lb1:ub1);
                            if (mean(F1) > dark_thresh)
                                A = ones((ub1-lb1+1),2);
                                A(:,1) = (lb1:ub1)';
                                pA = pinv(A);
                                coeff = pA*(F1');
                                fi = A*coeff;
                                if (coeff(1) > (0.2/(ub1 - lb1 + 1)))
                                    %1-r^2
                                    tmpcc1 = sum((F1' - fi).^2)/sum((F1' - mean(F1)).^2);
                                    clear fi F1

                                    G1 = lib(i,lb1:ub1);
                                    coeff = pA*(G1');
                                    clear pA
                                    gi = A*coeff;
                                    clear A
                                    %1-r^2
                                    tmpcc2 = sum((G1' - gi).^2)/sum((G1' - mean(G1)).^2);
                                    clear gi G1

                                    tmpcc = tmpcc1 - tmpcc2;
                                    if (tmpcc < 0)
                                        tmpcc = 0;
                                    end
                                else
                                    tmpcc = Inf;
                                end
                            else
                                tmpcc = Inf;
                            end
                        elseif (j == 9)
                            %"interesting" regions
                            tmpcc1 = zeros(2*wiggle+1,1);
                            inc = 1;
                            for w = -wiggle:wiggle
                                lb1 = lb1 + w;
                                ub1 = ub1 + w;
                                if ((lb1>=1) && (ub1<=p0))
                                    G1 = lib_der1(i,lb1:ub1);
                                    F1 = tmpm1(n,lb1:ub1);
                                    G1 = reshape(G1,[1 (ub1-lb1+1)]);
                                    tmpcc1(inc) = sum((G1-mean(G1)).*(F1-mean(F1)))/sqrt(sum((G1-mean(G1)).^2).*sum((F1-mean(F1)).^2));
                                    clear G1 F1
                                end
                                inc = inc + 1;
                            end
                            tmpcc = max(tmpcc1);
                            tmpcc = 1 - (tmpcc + 1)/2;
                        end
                        features(pnt,:) = [m n i j k tmpcc 1];
                        pnt = pnt + 1;
                        if (pnt > fsz)
                            features = [features;Inf(1024,7)];
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
Lf_list = [transition_edge_fwhm_col-3];
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
        tmpAsym = pigment_ris_table(i,j+1);
        tmpAsym = tmpAsym{1};
        if (~strcmp(tmpBW,'0'))
            tmpBW1 = regexp(tmpBW, '\;', 'split');
            tmpAsym1 = regexp(tmpAsym, '\;', 'split');
            if (~strcmp(tmpBW1{k},'0'))
                mdBW = single(str2double(tmpBW1{k}));
                mdAsym = single(str2double(tmpAsym1{k}));

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
                cube_tmp = der1(pk_m,pk_n,pk1:pk2);
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
                bw = abs(bw1 + bw2 - mdBW);
                %bw = abs(min(bw1,bw2) - mdBW);
                asym = abs(bw2 - bw1 - mdAsym);
                if (bw > (2*bw_dist_range))
                    bw = Inf;
                    asym = Inf;
                end
                features(ind1(ii),1:6) = [m n i j k bw];
                features(ind1(ii)+1,1:6) = [m n i j+1 k asym];
            end
        end
    end
end
clear ind ind1 tmp

% compute missing data penalty
if (lib_test == 1)
missing_penalty = single(Inf(Lp*Lf*100,4));
cnt = 1;
for i = 1:Lp
    for j = 1:Lf
        tmp = pigment_ris_table(i,j);
        tmp = tmp{1};
        if (~strcmp(tmp,'0'))
            tmp1 = regexp(tmp, '\;', 'split');
            Ltmp1 = length(tmp1);
            for k = 1:Ltmp1
                ind_pk1 = ((features(:,1)==i) & (features(:,3)==i) & (features(:,4)==j) & (features(:,5)==k));
                missing_penalty(cnt,:) = [i j k features(ind_pk1,7)];
                cnt = cnt + 1;
            end
        end
    end
end
missing_penalty(cnt:end,:) = [];
else
    missing_penalty = [];
end
