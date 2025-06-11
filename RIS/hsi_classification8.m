clear all
close all
clc

opt1 = 2; %0->none, 1->mean removed, 2->derivative
opt2 = 1; %1->correlation, 2->spectral angle
fn_path = 'C:\data\Pacino\';
trial = 'contremov_crop_Pacino_ChristMajesty_vnir_xnir';
spectra_lib = 'contremov_resampled_EarlyRennaisance1';
spectra_names = {'bone_black', 'cobalt_blue', 'natural_raw_umber', ...
    'verdigris', 'burnt_umber', 'Madder_lake', 'Cd_Red', ...
    'vermilion', 'titanium_white', 'Yellow_Ochre', 'lead_tin_yellow', ...
    'Naples_Yellow', 'Cd_yellow', 'Azurite', 'ultramarine', ...
    'indian_yellow', 'Malachite', 'Indigo', 'Precipitated_chalk', ...
    'lead_white', 'lamp_black', 'Lapis_ultramarine', 'Smalt', ...
    'VanDyke_brown', 'Green_earth', 'Realgar', 'Red_lead', 'Orpiment', ...
    'gypsum','Copper resinate', 'Carmine lake', 'red ochre'};

roi = [1 430 650; ...
    1 650 1010; ...
    1 1010 1600; ...
    2 430 630; ...
    2 630 750; ...
    2 950 1010; ...
    2 1010 1250; ...
    3 430 630; ...
    3 630 1010; ...
    3 1010 1600; ...
    4 430 650; ...
    4 1010 1350; ...
    5 430 650; ...
    5 650 1010; ...
    5 1010 1600; ...
    6 430 650; ...
    6 650 1010; ...
    6 1010 1125; ...
    7 430 650; ...
    7 650 1010; ...
    7 1010 1125; ...
    8 430 650; ...
    8 650 1010; ...
    8 1010 1600; ...
    9 430 650; ...
    9 650 1010; ...
    9 1010 1120; ...
    9 1250 1600; ...
    10 430 650; ...
    10 650 1010; ...
    11 430 650; ...
    11 650 1010; ...
    11 1010 1100; ...
    12 430 650; ...
    12 650 1010; ...
    12 1010 1100; ...
    13 430 650; ...
    13 650 1010; ...
    13 1010 1100; ...
    14 430 650; ...
    14 650 875; ...
    14 2200 2400; ...
    15 430 650; ...
    15 650 1010; ...
    15 1010 1150; ...
    16 430 650; ...
    16 650 1010; ...
    17 430 650; ...
    17 650 800; ...
    17 800 1010; ...
    17 1010 1600; ...
    18 430 650; ...
    18 650 1010; ...
    18 1010 1250; ...
    19 430 650; ...
    19 650 1010; ...
    19 1010 1120; ...
    19 1250 1600; ...
    20 430 650; ...
    20 650 1010; ...
    20 1010 1120; ...
    20 1250 1600; ...
    21 430 650; ...
    21 650 1010; ...
    21 1010 1600; ...
    22 430 650; ...
    22 650 1010; ...
    22 1010 1150; ...
    23 430 630; ...
    23 630 750; ...
    23 950 1010; ...
    23 1010 1300; ...
    24 430 650; ...
    24 650 1010; ...
    24 1010 1300; ...
    25 430 650; ...
    25 650 1010; ...
    25 1010 1300; ...
    26 430 650; ...
    26 650 800; ...
    27 430 650; ...
    27 650 800; ...
    28 430 650; ...
    28 650 800; ...
    29 430 650; ...
    29 1800 2100; ...
    30 430 650; ...
    30 650 1130; ...
    31 430 650; ...
    32 430 650; ...
    32 625 810; ...
    32 810 1010];
Lsn =  length(spectra_names);
%w = 27;
w = 69;
corr_thresh = 0.6;

% spectral library
fn = [fn_path spectra_lib];
lambdaon = 0;
lambda1 = [];
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
        [~,~,e] = regexp(line,'^wavelength={','match','start','end');
        if (e>0)
            lambdaon = 1;
            tmp1 = line((e+1):end);
            tmp2 = regexp(tmp1,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda1 = [lambda1 str2double(tmp2(i))];
                end
                [~,~,e1] = regexp(tmp2(i),'}$','match','start','end');
                if (e1{1} > 0)
                    tmp3 = cell2mat(tmp2(i));
                    tmp3 = tmp3(1:(e1{1}-1));
                    lambda1 = [lambda1 str2double(tmp3)];
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
                    lambda1 = [lambda1 str2double(tmp2(i))];
                end
            end
        else
            tmp2 = regexp(line,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda1 = [lambda1 str2double(tmp2(i))];
                end
            end
        end
    end
    clear msk tmp1 tmp2 line tmp3
end
fclose(fid);
lambda1 = lambda1';
lib = multibandread(fn,[m01,n01,p01],datatype,headeroffset,interleave,byteorder);
[mr,nr] = size(roi);
roi1 = zeros(mr,nr);
roisz = zeros(Lsn,1);
if (mr > 0)
    roi1(:,1) = roi(:,1);
    roi1(:,2) = ones(1,mr);
    roi1(:,3) = length(lambda1)*ones(1,mr);
    for i = 1:mr
        tmp = abs(lambda1 - roi(i,2));
        [~,ind] = min(tmp);
        roi1(i,2) = ind;
        tmp = abs(lambda1 - roi(i,3));
        [~,ind] = min(tmp);
        clear tmp
        roi1(i,3) = ind;
        roisz(roi1(i,1)) = roisz(roi1(i,1)) + roi1(i,3) - roi1(i,2) + 1;
    end
end

if (opt1 == 1)
% spectra library with local mean subtracted off
lib1 = zeros(m01,n01);
%tmp = zeros(m01,w);
%cnt = 1;
for n = 1:n01
%for n = 1:n01
    %tmp(:,cnt) = lib(:,n);
    n011 = n - floor(w/2);
    if (n011 < 1)
        n011 = 1;
    end
    n012 = n + floor(w/2);
    if (n012 > n01)
        n012 = n01;
    end
    %if (n >= w)
        lib1(:,n) = lib(:,n) - mean(lib(:,n011:n012),2);
        %lib1(:,(n-floor(w/2))) = lib(:,(n-floor(w/2))) - mean(tmp,2);
    %end
    %{
    cnt = cnt + 1;
    if (cnt > w)
        cnt = 1;
    end
    %}
end
%clear tmp
multibandwrite(single(lib1),[fn '_muremoved'],'bsq',[1,1,1],[m01,n01,p01]);
copyfile(fnh,[fn '_muremoved.hdr']);
end

if (opt1 == 2)
% library slope
lib2 = zeros(m01,n01);
hsz_lp = 3;
for n = (hsz_lp+ceil(hsz_lp/2)):(n01-hsz_lp-floor(hsz_lp/2))
    tmp = lib(:,(n-hsz_lp-floor(hsz_lp/2)):(n+hsz_lp+floor(hsz_lp/2)));
    Llr = 3*hsz_lp;
    mpt = ceil(Llr/2);
    lt = mean(tmp(:,1:hsz_lp),2);
    rt = mean(tmp(:,(Llr-hsz_lp+1):Llr),2);
    lib2(:,n) = rt - lt;
end
clear tmp lt rt
clear lib
multibandwrite(single(lib2),[fn '_slope'],'bsq',[1,1,1],[m01,n01,p01]);
copyfile(fnh,[fn '_slope.hdr']);
end

% data cube
fn = [fn_path trial];
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
lambda = lambda';
data = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);

% data cube with local mean subtracted off
if (opt1 == 1)
fn2 = [fn_path 'muremoved_avg' num2str(w) '_' trial];
fn2h = [fn2 '.hdr'];
data2 = single(zeros(m0,n0,p0));
end
bnd_msk = true(m0,n0);
%for p = (floor(w/2)+1):(p0-floor(w/2))
for p = 1:p0
    p011 = p - floor(w/2);
    if (p011 < 1)
        p011 = 1;
    end
    p012 = p + floor(w/2);
    if (p012 > p0)
        p012 = p0;
    end
    %tmp = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Range',[(p - floor(w/2)) (p + floor(w/2))]});
    %tmp = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Range',[p011 p012]});
    tmp = data(:,:,p011:p012);
    %msk = (tmp(:,:,(floor(w/2)+1)) == 0);
    p011 = p - floor(w/2);
    if (p011 < 1)
        msk = (tmp(:,:,1) == 0);
    else
        msk = (tmp(:,:,(floor(w/2)+1)) == 0);
    end
    bnd_msk = and(bnd_msk,msk);
    clear msk
    %data2(:,:,p) = single(tmp(:,:,(floor(w/2)+1)) - mean(tmp,3));
    if (opt1 == 1)
    if (p011 < 1)
        data2(:,:,p) = single(tmp(:,:,1) - mean(tmp,3));
    else
        data2(:,:,p) = single(tmp(:,:,(floor(w/2)+1)) - mean(tmp,3));
    end
    end
end
clear tmp
bnd_msk = (bnd_msk == 0);
if (opt1 == 1)
multibandwrite(single(data2),fn2,'bsq',[1,1,1],[m0,n0,p0]);
fid1 = fopen(fn2h,'w');
fprintf(fid1,'ENVI\n');
fprintf(fid1,'description = {}\n');
fprintf(fid1,'samples = %u\n',n0);
fprintf(fid1,'lines = %u\n',m0);
fprintf(fid1,'bands = %u\n',p0);
fprintf(fid1,'header offset = 0\n');
fprintf(fid1,'file type = ENVI Standard\n');
fprintf(fid1,'data type = 4\n');
fprintf(fid1,'interleave = bsq\n');
fprintf(fid1,'byte order = 0\n');
fprintf(fid1,'Wavelength = {');
for i = 1:(p0-1)
    fprintf(fid1,'%f, ',lambda(i));
end
fprintf(fid1,'%f}\n',lambda(p0));
fclose(fid1);
end

%{
cc = single(zeros(p0,p0));
for p1 = 1:(p0-1)
    for p2 = (p1+1):p0
        tmp1 = data2(:,:,p1);
        tmp1 = tmp1(:);
        tmp1(bnd_msk(:)==0) = [];
        tmp2 = data2(:,:,p2);
        tmp2 = tmp2(:);
        tmp2(bnd_msk(:)==0) = [];
        cc(p1,p2) = sum(tmp1.*tmp2)/sqrt(sum(tmp1.^2)*sum(tmp2.^2));
        clear tmp1 tmp2
    end
end
cc1 = cc';
cc_tot = cc + cc1;
%clear cc cc1
for i = 1:p0
    cc_tot(i,i) = 1;
end
cc_tot(isnan(cc_tot)) = 0;

for i = 1:p0
    figure(100),plot(cc_tot(i,:))
    pause(0.1)
end
%}

if (opt1 == 2)
% data cube derivative
fn2 = [fn_path 'slope_' trial];
fn2h = [fn2 '.hdr'];
data3 = single(zeros(m0,n0,p0));
for p = (hsz_lp+ceil(hsz_lp/2)):(p0-hsz_lp-floor(hsz_lp/2))
    %tmp = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Range',[(p-hsz_lp-floor(hsz_lp/2)) (p+hsz_lp+floor(hsz_lp/2))]});
    tmp = data(:,:,(p-hsz_lp-floor(hsz_lp/2)):(p+hsz_lp+floor(hsz_lp/2)));
    Llr = 3*hsz_lp;
    mpt = ceil(Llr/2);
    lt = mean(tmp(:,:,1:hsz_lp),3);
    rt = mean(tmp(:,:,(Llr-hsz_lp+1):Llr),3);
    data3(:,:,p) = single(rt - lt);
end
clear tmp lt rt
multibandwrite(single(data3),fn2,'bsq',[1,1,1],[m0,n0,p0]);
fid1 = fopen(fn2h,'w');
fprintf(fid1,'ENVI\n');
fprintf(fid1,'description = {}\n');
fprintf(fid1,'samples = %u\n',n0);
fprintf(fid1,'lines = %u\n',m0);
fprintf(fid1,'bands = %u\n',p0);
fprintf(fid1,'header offset = 0\n');
fprintf(fid1,'file type = ENVI Standard\n');
fprintf(fid1,'data type = 4\n');
fprintf(fid1,'interleave = bsq\n');
fprintf(fid1,'byte order = 0\n');
fprintf(fid1,'Wavelength = {');
for i = 1:(p0-1)
    fprintf(fid1,'%f, ',lambda(i));
end
fprintf(fid1,'%f}\n',lambda(p0));
fclose(fid1);
end

% compute similarity between each library spectra
ambiguity_map = single(zeros(m01,m01));
[mr1,~] = size(roi1);
for i1 = 1:(m01-1)
    if (mr1 > 0)
        roi11 = roi1((roi1(:,1) == i1),:);
    else
        roi11 = [];
    end
    [mr11,~] = size(roi11);
    if (mr11 == 0)
        roi11 = [i1 1 length(lambda1)];
        [mr11,~] = size(roi11);
    end
    G = [];
    for j = 1:mr11
        if (opt1 == 0)
            G1 = lib(i1,roi11(j,2):roi11(j,3));
        elseif (opt1 == 1)
            G1 = lib1(i1,roi11(j,2):roi11(j,3));
        elseif (opt1 == 2)
            G1 = lib2(i1,roi11(j,2):roi11(j,3));
        else
            G1 = lib1(i1,roi11(j,2):roi11(j,3));
        end
        %G = [G (G1-mean(G1))/std(G1)*((roi1(j,3)-roi1(j,2)+1)/roisz(roi1(j,1)))];
        G = [G G1];
        %G = (G1-mean(G1))/std(G1);
    end
    for i2 = (i1+1):m01
        %F = [];
        tmpcc = 0;
        pnt = 1;
        for j = 1:mr11
            if (opt1 == 0)
                F1 = lib(i2,roi11(j,2):roi11(j,3));
            elseif (opt1 == 1)
                F1 = lib1(i2,roi11(j,2):roi11(j,3));
            elseif (opt1 == 2)
                F1 = lib2(i2,roi11(j,2):roi11(j,3));
            else
                F1 = lib1(i2,roi11(j,2):roi11(j,3));
            end
            %F = [F (F1-mean(F1))/std(F1)*((roi1(j,3)-roi1(j,2)+1)/roisz(roi1(j,1)))];
            %F = (F1-mean(F1))/std(F1);
            LF1 = length(F1);
            G11 = G(pnt:(pnt+LF1-1));
            pnt = pnt + LF1;
            FG = single(sum((G11-mean(G11)).*(F1-mean(F1)))/(std(G11)*std(F1)));
            %tmpcc = tmpcc + single(sum(F.*G)*((roi1(j,3)-roi1(j,2)+1)/roisz(roi1(j,1))));
            %tmpcc = tmpcc + single(FG*((roi11(j,3)-roi11(j,2)+1)/roisz(roi11(j,1))));
            tmpcc = tmpcc + single(FG/roisz(roi11(j,1)));
        end
        %ambiguity_map(i1,i2) = single(sum(F.*G)/sqrt(sum(F.^2)*sum(G.^2)));
        ambiguity_map(i1,i2) = single(tmpcc);
        %ambiguity_map(i1,i2) = (pi - abs(single(acos(F*G'/norm(F)/norm(G)))))/pi;
    end
end
figure,imagesc(ambiguity_map)

if (opt1 == 2)
ambiguity_map2 = single(zeros(m01,m01));
[mr1,~] = size(roi1);
for i1 = 1:(m01-1)
    if (mr1 > 0)
        roi11 = roi1((roi1(:,1) == i1),:);
    else
        roi11 = [];
    end
    [mr11,~] = size(roi11);
    if (mr11 == 0)
        roi11 = [i1 1 length(lambda1)];
        [mr11,~] = size(roi11);
    end
    G = [];
    for j = 1:mr11
        if (opt1 == 0)
            G1 = lib(i,roi11(j,2):roi11(j,3));
        elseif (opt1 == 1)
            G1 = lib1(i,roi11(j,2):roi11(j,3));
        elseif (opt1 == 2)
            G1 = lib2(i,roi11(j,2):roi11(j,3));
        else
            G1 = lib1(i,roi11(j,2):roi11(j,3));
        end

        G = [G (G1-mean(G1))/std(G1)*((roi1(j,3)-roi1(j,2)+1)/roisz(roi1(j,1)))];
    end
    for i2 = (i1+1):m01
        F = [];
        for j = 1:mr11
            if (opt1 == 0)
                F1 = lib(i2,roi11(j,2):roi11(j,3));
            elseif (opt1 == 1)
                F1 = lib1(i2,roi11(j,2):roi11(j,3));
            elseif (opt1 == 2)
                F1 = lib2(i2,roi11(j,2):roi11(j,3));
            else
                F1 = lib1(i2,roi11(j,2):roi11(j,3));
            end
            F = [F (F1-mean(F1))/std(F1)*((roi1(j,3)-roi1(j,2)+1)/roisz(roi1(j,1)))];
        end
        %ambiguity_map2(i1,i2) = single(sum(F.*G)/sqrt(sum(F.^2)*sum(G.^2)));
        ambiguity_map2(i1,i2) = (pi - abs(single(acos(F*G'/norm(F)/norm(G)))))/pi;
    end
end
figure,imagesc(ambiguity_map2)
end

ncores = feature('numCores');
if (ncores > 1)
    isOpen = matlabpool('size') > 0;
    if (isOpen == 1)
        matlabpool close
    end
    matlabpool(ncores)
end

% matched filter(?)
hsi_map = single(zeros(m0,n0,m01));
for m = 1:m0
    if (opt1 == 0)
        tmpm = data(m,:,:);
    elseif (opt1 == 1)
        tmpm = data2(m,:,:);
    elseif (opt1 == 2)
        tmpm = data3(m,:,:);
    end
    
    tmpm = reshape(tmpm,[n0,p0]);
    hsi_map0 = single(zeros(1,n0,m01));
    parfor n = 1:n0
        for i = 1:m01
            [mr1,~] = size(roi1);
            if (mr1 > 0)
                roi11 = roi1((roi1(:,1) == i),:);
            else
                roi11 = [];
            end
            [mr11,nr11] = size(roi11);
            if (mr11 == 0)
                roi11 = [i 1 length(lambda1)];
                [mr11,nr11] = size(roi11);
            end
            G = [];
            F = [];
            tmpcc = 0;
            for j = 1:mr11
                if (opt1 == 0)
                    G1 = lib(i,roi11(j,2):roi11(j,3));
                elseif (opt1 == 1)
                    G1 = lib1(i,roi11(j,2):roi11(j,3));
                elseif (opt1 == 2)
                    G1 = lib2(i,roi11(j,2):roi11(j,3));
                else
                    G1 = lib1(i,roi11(j,2):roi11(j,3));
                end
                
                %G = [G (G1-mean(G1))/std(G1)*((roi1(j,3)-roi1(j,2)+1)/roisz(roi1(j,1)))];
                %G = (G1-mean(G1))/std(G1);
                F1 = tmpm(n,roi11(j,2):roi11(j,3));
                %muG = sum(abs(G1(:)));
                %muF = sum(abs(F1(:)));
                %F = [F F1*muG/muF];
                %F = [F (F1-mean(F1))/std(F1)*((roi1(j,3)-roi1(j,2)+1)/roisz(roi1(j,1)))];
                %F = (F1-mean(F1))/std(F1);
                FG = single(sum((G1-mean(G1)).*(F1-mean(F1)))/(std(G1)*std(F1)));
                %tmpcc = tmpcc + single(sum(F.*G)*((roi1(j,3)-roi1(j,2)+1)/roisz(roi1(j,1))));
                tmpcc = tmpcc + single(FG/roisz(roi11(j,1)));
            end
            if (opt2 == 1)
                %hsi_map0(1,n,i) = single(sum(F.*G)/sqrt(sum(F.^2)*sum(G.^2)));
                %hsi_map0(1,n,i) = single(sum(F.*G));
                hsi_map0(1,n,i) = single(tmpcc);
            elseif (opt2 == 2)
                hsi_map0(1,n,i) = (pi - abs(single(acos(F*G'/norm(F)/norm(G)))))/pi;
            end
        end
    end
    hsi_map(m,:,:) = single(hsi_map0);
    clear hsi_map0 tmpm xc F G F1 G1
end

if (ncores > 1)
    isOpen = matlabpool('size') > 0;
    if (isOpen == 1)
        matlabpool close
    end
end

% save/display results
for i = 1:m01
    tmp = hsi_map(:,:,i).*single(bnd_msk);
    
    figure,imagesc(tmp)
    axis image
    colormap gray
    title(spectra_names{i})
    
    tmp1 = tmp;
    tmp1(tmp1 == 0) = [];
    %mx = max(tmp1(:));
    %mn = min(tmp1(:));
    mx = 1;
    mn = -1;
    clear tmp1
    tmp = (tmp-mn)/(mx-mn)*(2^16-1);
    fn = [fn_path trial '_' spectra_names{i} '_' num2str(w) '.tif'];
    imwrite(uint16(tmp),fn,'tif','Compression','None')
    clear tmp
end

pass_thresh = zeros(m01,1);
for i = 1:m01
    hsz_lp = 1;
    tmp0 = hsi_map(:,:,i);
    tmp0 = tmp0.*single(bnd_msk);
    tmp = tmp0(:);
    tmp(bnd_msk(:)==0) = [];
    if (opt2 == 1)
        [h1,x1] = hist(tmp,-(1+hsz_lp*0.01):0.02:(1+hsz_lp*0.01));
    elseif (opt2 == 2)
        [h1,x1] = hist(tmp,-(hsz_lp*0.01):0.02:(1+hsz_lp*0.01));
    end
    h1 = h1/sum(h1);
    
    peak_list = single(zeros(length(x1),2)); %[m,n,lambda]
    valley_list = single(zeros(length(x1),2)); %[m,n,lambda]
    ptot = 1;
    ptot2 = 1;
    change_thresh = [];
    for p = 2:(length(x1)-1)
        tmp = h1((p-1):(p+1));
        tmp1 = tmp(2) - tmp(1);
        tmp2 = tmp(3) - tmp(2);
        msk = (tmp1 >= 0) & (tmp2 <= 0) & (tmp1 ~= tmp2);
        if (msk == 1)
            peak_list(ptot,:) = [x1(p) tmp(2)];
            ptot = ptot + 1;
        end
        clear msk
        msk = (tmp1 <= 0) & (tmp2 >= 0) & (tmp1 ~= tmp2);
        if (msk == 1)
            valley_list(ptot2,:) = [x1(p) tmp(2)];
            ptot2 = ptot2 + 1;
        end
        clear msk
        clear tmp1 tmp2 tmp
    end
    peak_list(peak_list(:,1)==0,:) = [];
    valley_list(valley_list(:,1)==0,:) = [];
    clear msk
        
    % find closest peak closest for each valley
    [mxp,mpind] = max(peak_list(:,1));
    peak_list(peak_list(:,2)<0.01,:) = [];
    [Lv,~] = size(valley_list);
    for j = Lv:-1:1
        [~,vpind] = min(abs(peak_list(:,1) - valley_list(j,1)));
        if (valley_list(j,2) > 0.5*peak_list(vpind(1),2))
            valley_list(j,:) = [];
        end
    end
    valley_list = valley_list((valley_list(:,1)<mxp),:);
    plot(x1,h1,'r');
    hold on
    plot(peak_list(:,1),peak_list(:,2),'bx')
    hold on
    plot(valley_list(:,1),valley_list(:,2),'bo')
    hold off
    
    if (mxp > corr_thresh)
        valley_list0 = valley_list(:,1);
        pass_thresh(i) = max(valley_list0);
    else
        pass_thresh(i) = Inf;
    end
    
    spectra_names{i}
    pause
end

% create maps
class_map = uint8(zeros(m0,n0,m01));
for i = 1:m01
    figure
    subplot(121),imagesc(hsi_map(:,:,i))
    axis image
    colormap gray
    title(spectra_names{i})
    
    class_map(:,:,i) = 255*uint8(hsi_map(:,:,i) >= pass_thresh(i));
    subplot(122),imshow(class_map(:,:,i))
    
    fn = [fn_path trial '_map_' spectra_names{i} '_' num2str(w) '_' 'opt' num2str(opt2) '.tif'];
    imwrite(uint8(class_map(:,:,i)),fn,'tif','Compression','None')
    clear tmp
end
save([fn_path trial '.mat'],'hsi_map','spectra_names','bnd_msk','class_map')
