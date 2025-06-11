function [map1,mapp,A1,fwhm1,c1,r1,element_list,LN,c,xrf_table,xray_label] = xrf_classification(fnpath,trial,xrf_table_fn,m1,m2,n1,n2,w,hsz_lp,N,cal_energy,cal_band,minE,maxE,ncores)

warning('off','MATLAB:rankDeficientMatrix');

lp_mn = 5;
%lp_p = 5;
dcnt0 = 3;
dcnt1 = lp_mn^2;

% read in previous calibration coefficients
fnc = 'XRFMapping.txt';
fid = fopen(fnc);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
    [~,~,e] = regexp(line,'^a1=','match','start','end');
    if (e>0)
        a1 = str2double(line((e+1):end));
    end
    [~,~,e] = regexp(line,'^a0=','match','start','end');
    if (e>0)
        a0 = str2double(line((e+1):end));
    end
end
fclose(fid);

b1 = 0.0118;    %fwhm fit coeffs
b0 = 0.0947;
fn_fit = [fnpath 'fit_' trial];
element_list = {'H'; 'He'; 'Li'; 'Be'; 'B'; 'C'; 'N'; 'O'; 'F'; ...
    'Ne'; 'Na'; 'Mg'; 'Al'; 'Si'; 'P'; 'S'; 'Cl'; 'Ar'; 'K'; ...
    'Ca'; 'Sc'; 'Ti'; 'V'; 'Cr'; 'Mn'; 'Fe'; 'Co'; 'Ni'; 'Cu'; ...
    'Zn'; 'Ga'; 'Ge'; 'As'; 'Se'; 'Br'; 'Kr'; 'Rb'; 'Sr'; 'Y'; ...
    'Zr'; 'Nb'; 'Mo'; 'Tc'; 'Ru'; 'Rh'; 'Pd'; 'Ag'; 'Cd'; 'In'; ...
    'Sn'; 'Sb'; 'Te'; 'I'; 'Xe'; 'Cs'; 'Ba'; 'La'; 'Ce'; 'Pr'; ...
    'Nd'; 'Pm'; 'Sm'; 'Eu'; 'Gd'; 'Tb'; 'Dy'; 'Ho'; 'Er'; 'Tm'; ...
    'Yb'; 'Lu'; 'Hf'; 'Ta'; 'W'; 'Re'; 'Os'; 'Ir'; 'Pt'; 'Au'; ...
    'Hg'; 'Tl'; 'Pb'; 'Bi'; 'Po'; 'At'; 'Rn'; 'Fr'; 'Ra'; 'Ac'; ...
    'Th'; 'Pa'; 'U'; 'Np'; 'Pu'; 'Am'; 'Cm'; 'Bk'; 'Cf'; 'Es'; ...
    'Fm'};
xray_label = {'Kalpha';'Kbeta';'Malpha';'Mbeta';'Lalpha';'Lbeta'};

xrf_table = csvread(xrf_table_fn,1,2);
fn = [fnpath trial];
fnh = [fnpath trial '.hdr'];
fid = fopen(fnh);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
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
end
if (m2 == 0)
    m2 = m0;
end
if (n2 == 0)
    n2 = n0;
end
fclose(fid);
lambda = single(1:p0);
lambda = lambda(:);

% find peaks (1st derivatives = 0, 2nd derivative = -)
peak_list = single(zeros(m0*n0,3)); %[m,n,lambda]
snr_pass = false(m0,n0,p0);
[Ltot,~] = size(peak_list);
ptot = 1;
for p = (hsz_lp+ceil(hsz_lp/2)):(p0-hsz_lp-floor(hsz_lp/2))
    tmp = single(multibandread(fn,[m0,n0,p0],'single',0,'bsq','ieee-le',{'Band','Range',[p-hsz_lp-floor(hsz_lp/2) p+hsz_lp+floor(hsz_lp/2)]}));
    Llr = 3*hsz_lp;
    mpt = ceil(Llr/2);
    lt = mean(tmp(:,:,1:hsz_lp),3);
    md = mean(tmp(:,:,(mpt-floor(hsz_lp/2)):(mpt+floor(hsz_lp))),3);
    rt = mean(tmp(:,:,(Llr-hsz_lp+1):Llr),3);
    tmp1 = md - lt;
    tmp2 = rt - md;
    snr_pass(:,:,p) = (abs(tmp1) > dcnt0) | (abs(tmp2) > dcnt0);
    msk = (tmp1 > 0) & (tmp2 < 0) & snr_pass(:,:,p);
    lambda_interp = lambda(p-hsz_lp) + (lambda(p+hsz_lp) - lambda(p-hsz_lp))*(abs(tmp1)./abs(tmp1-tmp2));
    clear tmp tmp1 tmp2 delta1 delta2
    [r,c] = find(msk==1);
    L = length(r);
    while ((ptot + L - 1) > Ltot)
        peak_list = [peak_list;single(zeros(m0*n0,3))];
        [Ltot,~] = size(peak_list);
    end
    for i = 1:L
        peak_list(ptot,:) = [r(i) c(i) lambda_interp(r(i),c(i))];
        ptot = ptot + 1;
    end
end
peak_list(peak_list(:,1)==0,:) = [];
clear lambda_interp msk
clear snr_pass

% calibration coefficients
[n,xout] = hist(peak_list(:,3),21920);
clear peak_list
%figure,plot(xout,n)
[a1,a0,r2] = auto_cal(xout,n,a1,a0,xrf_table,minE,maxE,N);
clear n xout
A = [cal_band ones(length(cal_energy),1)];
coeff = pinv(A)*cal_energy;
a1 = coeff(1);
a0 = coeff(2);
figure,plot(cal_band,cal_energy,'k*')
hold on
plot(cal_band,a1*cal_band+a0,'r');
hold off
ybar = mean(cal_energy);
sigma2 = sum((cal_energy-ybar).^2);
sse = sum((a1*cal_band+a0-cal_energy).^2);
r2 = 1 - (sse/sigma2);
clear peak_list
lambda0 = single(1:p0)*a1 + a0;
lambda0 = lambda0(:);

% spatially-smoothed data cube
mu_mn = single(zeros(m0,n0,p0));
cube = multibandread(fn,[m0,n0,p0],'single',0,'bsq','ieee-le');
for m = 1:m0
    for n = 1:n0
        m11 = (m-floor(lp_mn/2));
        if (m11 < 1)
            m11 = 1;
        end
        m21 = (m+floor(lp_mn/2));
        if (m21 > m0)
            m21 = m0;
        end
        n11 = (n-floor(lp_mn/2));
        if (n11 < 1)
            n11 = 1;
        end
        n21 = (n+floor(lp_mn/2));
        if (n21 > n0)
            n21 = n0;
        end
        tmp = cube(m11:m21,n11:n21,:);
        tmp = reshape(tmp,[(m21-m11+1)*(n21-n11+1),p0]);
        %tmp21 = mean(tmp,1);
        tmp21 = sum(tmp,1);
        mu_mn(m,n,:) = reshape(tmp21,[1,1,p0]);
        clear tmp21
    end
end
clear cube
fn1 = [fnpath 'avg' num2str(lp_mn) '_' trial];
fn1h = [fn1 '.hdr'];
multibandwrite(single(mu_mn),fn1,'bsq',[1,1,1],[m0,n0,p0]);
fid1 = fopen(fn1h,'w');
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
    fprintf(fid1,'%f, ',lambda0(i));
end
fprintf(fid1,'%f}\n',lambda0(p0));
fclose(fid1);
%{
% spectral average
mu_p = single(zeros(m0,n0,p0));
tmp = single(zeros(m0,n0,lp_p));
cnt = 1;
cube = multibandread(fn,[m0,n0,p0],'single',0,'bsq','ieee-le');
for p = 1:p0
    tmp(:,:,cnt) = cube(:,:,p);
    if (p >= lp_p)
        mu_p(:,:,(p - floor(lp_p/2))) = mean(tmp,3);
    end
    
    cnt = cnt + 1;
    if (cnt > lp_p)
        cnt = 1;
    end
end
clear tmp cube
%}
% find peaks (1st derivatives = 0, 2nd derivative = -)
snr_pass = false(m0,n0,p0);
for p = (hsz_lp+ceil(hsz_lp/2)):(p0-hsz_lp-floor(hsz_lp/2))
    tmp = mu_mn(:,:,(p-hsz_lp-floor(hsz_lp/2)):(p+hsz_lp+floor(hsz_lp/2)));
    Llr = 3*hsz_lp;
    mpt = ceil(Llr/2);
    lt = mean(tmp(:,:,1:hsz_lp),3);
    md = mean(tmp(:,:,(mpt-floor(hsz_lp/2)):(mpt+floor(hsz_lp))),3);
    rt = mean(tmp(:,:,(Llr-hsz_lp+1):Llr),3);
    tmp1 = md - lt;
    tmp2 = rt - md;
    
    snr_pass(:,:,p) = (abs(tmp1) > dcnt1) | (abs(tmp2) > dcnt1);
end
lambda = single(1:p0)*a1 + a0;
lambda = lambda(:);
clear mu_mn

% element mapping
[L,c] = size(xrf_table);
element_maps = false(m0,n0,L);
for i = 1:c
    tmp = xrf_table(:,i);
    tmp(tmp <= minE) = 0;
    tmp(tmp > maxE) = 0;
    xrf_table(:,i) = tmp;
    clear tmp
end
for i = 1:length(N)
    if (sum(N == i) == 0)
        xrf_table(i,:) = zeros(1,c);
    end
end

for i = 1:L
    pass = uint8(zeros(m0,n0)); 
    tmp = xrf_table(i,:);   %test each row of xrf table
    validE = sum(tmp>0);
    if (validE > 0)
        for j = 1:c
            if (tmp(j) ~= 0)
                e1 = tmp(j) - a1*floor(w/2);
                e2 = tmp(j) + a1*floor(w/2);
                
                % test if peak is encompassed in a different peak
                p1 = floor((e1 - a0)/a1);
                p2 = ceil((e2 - a0)/a1);
                test = (sum(snr_pass(:,:,p1:p2),3)>0);
                pass = pass + uint8(test);
                clear test
            end
        end
    end
    clear tmp
    tmp = (pass > 0); %no peak requirement, only look for counts
    element_maps(:,:,i) = tmp;
    clear tmp
end
clear pass
LN = length(N);
clear snr_pass

for i = 1:LN
    figure,imshow(element_maps(:,:,N(i)))
    title(element_list{N(i)})
    imwrite(uint8((2^8-1)*element_maps(:,:,N(i))),[fnpath trial '_map_' element_list{N(i)} '.tif'],'tif','Compression','None')
end

xrf_table0 = zeros(L,c);
for i = 1:LN
    xrf_table0(N(i),:) = xrf_table(N(i),:);
end
[L,G] = size(xrf_table0);
xrf_table0 = xrf_table0';
xrf_table0 = xrf_table0(:);
ov = logical(eye(L*G,L*G));
for G1 = 1:(L*G-1)
    sigma0 = 0.4247*(xrf_table0(G1)*b1 + b0);
    for G2 = (G1+1):(L*G)
        if ((xrf_table0(G1) ~= 0) && (xrf_table0(G2) ~= 0))
            if (((xrf_table0(G2)+2*sigma0)>=(xrf_table0(G1)-2*sigma0)) && ((xrf_table0(G2)+2*sigma0)<=(xrf_table0(G1)+2*sigma0)))
                ov(G1,G2) = 1;
                ov(G2,G1) = 1;
            end
            if (((xrf_table0(G2)-2*sigma0)>=(xrf_table0(G1)-2*sigma0)) && ((xrf_table0(G2)-2*sigma0)<=(xrf_table0(G1)+2*sigma0)))
                ov(G1,G2) = 1;
                ov(G2,G1) = 1;
            end
        end
    end
end
cntG1 = 1;
lasttot = 1;
while(cntG1 <= (L*G))
    tmp = ov(cntG1,:);
    ind = find(tmp);
    tot = length(ind);
    if (tot > lasttot);
        for k = 1:tot
            if (ind(k) ~= cntG1)
                tmp2 = ov(ind(k),:);
                tmp = or(tmp,tmp2);
            end
        end
        ov(cntG1,:) = tmp;
        lasttot = tot;
    else
        lasttot = 1;
        cntG1 = cntG1 + 1;
    end
    clear tmp tmp2
end
ov = ov.*(~eye(L*G,L*G));
clear xrf_table0

% sum peaks
if (ncores > 1)
    isOpen = matlabpool('size') > 0;
    if (isOpen == 1)
        matlabpool close
    end
    matlabpool(ncores)
end

mapp = zeros(m0,n0,LN*c);
map1 = zeros(m0,n0,LN);
A1 = zeros(m0,n0,LN*c);
fwhm1 = zeros(m0,n0,LN*c);
c1 = zeros(m0,n0,LN*c);
r1 = ones(m0,n0,LN*c);
tmp = zeros(1,n0,p0);
multibandwrite(single(tmp),fn_fit,'bip',[1 1 1],[m0,n0,p0]);
clear tmp
parfor m = m1:m2
    spect1 = single(multibandread(fn,[m0,n0,p0],'single',0,'bsq','ieee-le',{'Row','Direct',m}));
    %spect1 = single(mu_p(m,:,:));
    spect2 = single(zeros(1,n0,p0));
    mapp1 = zeros(1,n0,LN*c);
    A11 = zeros(1,n0,LN*c);
    fwhm11 = zeros(1,n0,LN*c);
    c11 = zeros(1,n0,LN*c);
    r11 = ones(1,n0,LN*c);
    for n = n1:n2
        spect = spect1(1,n,:);
        spect = spect(:);
        N0 = [];
        Ni = [];
        Npeak = [];
        e0 = [];
        fwhm = [];
        for i = 1:LN
            if (element_maps(m,n,N(i)) == 1)
                tmp = xrf_table(N(i),:);
                for j = 1:c
                    if (tmp(j) ~= 0)
                        N0 = [N0 N(i)];
                        Ni = [Ni i];
                        Npeak =[Npeak j];
                        e0 = [e0 xrf_table(N(i),j)];     %KeV
                        fwhm = [fwhm (b1*e0(end)+b0)];     %KeV
                    end
                end
            end
        end
        G = length(e0);
        if (G > 0)
            [map0,A0,fwhm0,c0,r0] = xrf_map(spect,ov,c,a0,a1,N0,Npeak,e0,fwhm,p0,xrf_table);
            for i = 1:G
                mapp1(1,n,(Ni(i)-1)*c+Npeak(i)) = map0(i);
                A11(1,n,(Ni(i)-1)*c+Npeak(i)) = A0(i);
                fwhm11(1,n,(Ni(i)-1)*c+Npeak(i)) = fwhm0(i);
                c11(1,n,(Ni(i)-1)*c+Npeak(i)) = c0(i);
                r11(1,n,(Ni(i)-1)*c+Npeak(i)) = r0(i);
                                
                tmp = (A0(i))*exp(-((lambda-e0(i)).^2)/(2*(0.4247*fwhm0(i))^2));
                tmp = reshape(tmp,[1,1,p0]);
                spect2(1,n,:) = spect2(1,n,:) + tmp;
            end
        end
    end
    mapp(m,:,:) = mapp1;
    map11 = zeros(1,n0,LN);
    for i = 1:LN
        tmp = xrf_table(N(i),:);
        for j = 1:c
            if (tmp(j) ~= 0)
                msk = (mapp1(:,:,(i-1)*c+j)>0);
                map11(:,:,i) = map11(:,:,i) + mapp1(:,:,(i-1)*c+j).*single(msk);
            end
        end
    end
    map1(m,:,:) = map11;
    A1(m,:,:) = A11;
    fwhm1(m,:,:) = fwhm11;
    c1(m,:,:) = c11;
    r1(m,:,:) = r11;
    multibandwrite(single(spect2),fn_fit,'bip',[m 1 1],[m0,n0,p0]);
end
clear tmp spect1 spect2 mu_lambda

fid1 = fopen([fn_fit '.hdr'],'w');
fprintf(fid1,'ENVI\n');
fprintf(fid1,'description = {}\n');
fprintf(fid1,'samples = %u\n',n0);
fprintf(fid1,'lines = %u\n',m0);
fprintf(fid1,'bands = %u\n',p0);
fprintf(fid1,'header offset = 0\n');
fprintf(fid1,'file type = ENVI Standard\n');
fprintf(fid1,'data type = 4\n');
fprintf(fid1,'interleave = bip\n');
fprintf(fid1,'byte order = 0\n');
fprintf(fid1,'Wavelength = {');
for i = 1:(p0-1)
    fprintf(fid1,'%f, ',lambda(i));
end
fprintf(fid1,'%f}\n',lambda(p0));
fclose(fid1);

if (ncores > 1)
    isOpen = matlabpool('size') > 0;
    if (isOpen == 1)
        matlabpool close
    end
end

for i = 1:LN
    figure,imagesc(map1(:,:,i))
    axis image
    colormap gray
    title(element_list{N(i)})
    mx = max(max(map1(:,:,i)));
    imwrite(uint16(map1(:,:,i)/mx*(2^16-1)),[fnpath trial '_sum_' element_list{N(i)} '.tif'],'tif','Compression','None')
    %{
    msk = (map1(:,:,i) > 0);
    imwrite(uint8((2^8-1)*msk),[fnpath trial '_final_map_' element_list{N(i)} '.tif'],'tif','Compression','None')
    %}
end

for i = 1:LN
    tmp = xrf_table(N(i),:);
    for j = 1:c
        if (tmp(j) ~= 0)
            ij = (i-1)*c+j;
            figure,imagesc(mapp(:,:,ij))
            axis image
            colormap gray
            title([element_list{N(i)} ' - ' xray_label{j}])
            mx = max(max(mapp(:,:,ij)));
            imwrite(uint16(mapp(:,:,ij)/mx*(2^16-1)),[fnpath trial '_sum_' element_list{N(i)} '_' xray_label{j} '.tif'],'tif','Compression','None')
        end
    end
end
warning('on','MATLAB:rankDeficientMatrix');
end

function [map1,A1,fwhm1,c1,r1] = xrf_map(spect1,ov,c,a0,a1,N0,Npeak,e0,fwhm,p0,xrf_table)

%e0 -> peaks found in each spectrum
%fwhm -> estimated fwhm of each peak
%N0 -> element number corresponding to each peak

warning('off','MATLAB:rankDeficientMatrix');

sigma = 0.4247*fwhm;    %KeV
N = length(e0);     %number of peaks in the spectrum
lambda = single(1:p0)*a1 + a0;        %KeV
lambda = lambda(:);

map1 = zeros(1,N);
A1 = zeros(1,N);
fwhm1 = zeros(1,N);
c1 = zeros(1,N);
r1 = ones(1,N);
for i = 1:N
    tmp = ov((N0(i)-1)*c+Npeak(i),:);
    ind = find(tmp);
    e00 = [e0(i)];
    for ii = 1:length(ind)
        N01 = floor(ind(ii)/c)+1;
        if (sum(N0==N01) > 0)
            Npeak1 = ind(ii)-c*floor((ind(ii)-1)/c);
            e00 = [e00 xrf_table(N01,Npeak1)];
        end
    end
    Aest = spect1(floor((e00-a0)/a1))';

    sigma00 = sigma(i);
    e1 = min(e00) - 2*sigma00; %KeV
    e2 = max(e00) + 2*sigma00; %KeV
    lambda1 = floor((e1-a0)/a1);    %channel
    if (lambda1 < 1)
        lambda1 = 1;
    end
    lambda2 = ceil((e2-a0)/a1);     %channel
    if (lambda2 > p0)
        lambda2 = p0;
    end
    lambda11 = lambda(lambda1:lambda2); %KeV
    spect = spect1(lambda1:lambda2);
    %N00 = length(e00);
    
    %{
    cest = 0;
    L = length(lambda11);
    resid_min = Inf;
        
    for j = 1:(2^N00)
        ind = dec2bin((j-1),N00)-'0';
        e001 = e00;
        e001(ind==1) = [];
        Aest001 = Aest;
        Aest001(ind==1) = [];
        if (~isempty(e001))
            parguess = [fwhm(i),cest,Aest001];
            [pars0,resid0] = nlinfit([lambda11;e001';L],spect,@gauss_func,parguess);
            %[pars0,resid0] = nlinfit([lambda11;e001';L],spect,@lorentz_func,parguess);
            if (resid0 < resid_min)
                resid_min = resid0;
                pars = pars0;
                best_e00 = e001;
            end
        end
    end
    
    if (best_e00(1) == e0(i))
        tmp = pars(1);
        if ((tmp > 0) && (tmp < 0.5))
            fwhm1(i) = tmp;
        elseif (tmp < 0)
            fwhm1(i) = 0;
        else
            fwhm1(i) = fwhm(i);
        end 
        tmp = pars(2);
        if (tmp < 0)
            c1(i) = 0;
        else
            c1(i) = tmp;
        end
        tmp = pars(3);
        if (tmp < 0)
            A1(i) = 0;
        else
            A1(i) = tmp;
        end
        sigma00 = 0.4247*fwhm1(i);

        e1 = e0(i) - 2*sigma00; %KeV
        e2 = e0(i) + 2*sigma00; %KeV
        lambda1 = floor((e1-a0)/a1);    %channel
        if (lambda1 < 1)
            lambda1 = 1;
        end
        lambda2 = ceil((e2-a0)/a1);     %channel
        if (lambda2 > p0)
            lambda2 = p0;
        end
        lambda11 = lambda(lambda1:lambda2); %KeV
        spect = spect1(lambda1:lambda2);
        f = A1(i)*exp(-((lambda11-e0(i)).^2)/(2*sigma00^2)) + c1(i);
        r1(i) = sum(abs(spect - f))/sum(spect);
        f = f - c1(i);
        f(f<0) = 0;
        map1(i) = sum(f);
    else
        map1(i) = 0;
        A1(i) = 0;
        fwhm1(i) = 0;
        c1(i) = 0;
        r1(i) = Inf;
    end
    %}
    
    cest = 0;
    parguess = [fwhm(i),cest,Aest];
    L = length(lambda11);
    pars = nlinfit([lambda11;e00';L],spect,@gauss_func,parguess);
    %pars = nlinfit([lambda11;e00';L],spect,@lorentz_func,parguess);
    tmp = pars(1);
    bad = 0;
    if ((tmp > 0) && (tmp < 2*fwhm(i)))
        fwhm1(i) = tmp;
        %{
    elseif (tmp < 0)
        fwhm1(i) = 0;
    else
        fwhm1(i) = fwhm(i);
        %}
    else
        fwhm1(i) = 0;
        bad = 1;
    end 
    tmp = pars(2);
    if (tmp < 0)
        c1(i) = 0;
    else
        c1(i) = tmp;
    end
    tmp = pars(3);
    if (tmp < 0)
        A1(i) = 0;
        bad = 1;
    else
        A1(i) = tmp;
    end
    sigma00 = 0.4247*fwhm1(i);
    
    if (bad == 0)
        e1 = e0(i) - 2*sigma00; %KeV
        e2 = e0(i) + 2*sigma00; %KeV
        lambda1 = floor((e1-a0)/a1);    %channel
        if (lambda1 < 1)
            lambda1 = 1;
        end
        lambda2 = ceil((e2-a0)/a1);     %channel
        if (lambda2 > p0)
            lambda2 = p0;
        end
        lambda11 = lambda(lambda1:lambda2); %KeV
        spect = spect1(lambda1:lambda2);
        f = A1(i)*exp(-((lambda11-e0(i)).^2)/(2*sigma00^2)) + c1(i);
        r1(i) = sum(abs(spect - f))/sum(spect);
        f = f - c1(i);
        f(f<0) = 0;
        map1(i) = sum(f);
    elseif (bad == 1)
        r1(i) = 1;
        map1(i) = 0;
        A1(i) = 0;
        fwhm1(i) = 0;
        c1(i) = 0;
    end
end
end

function f = gauss_func(pars,lambda11)
    L = lambda11(end);
    lambda = lambda11(1:L);
    e0 = lambda11((L+1):(end-1));
    N = length(e0);
    fwhm = pars(1);
    c = abs(pars(2));
    A = abs(pars(3:end));
    sigma = abs(0.4247*fwhm);
    f = zeros(L,1);
    for i = 1:N
        f = f + abs(A(i))*exp(-((lambda-e0(i)).^2)/(2*sigma^2));
    end
    f = f + c;
end

function f = lorentz_func(pars,lambda11)
    L = lambda11(end);
    lambda = lambda11(1:L);
    e0 = lambda11((L+1):(end-1));
    N = length(e0);
    fwhm = pars(1);
    c = abs(pars(2));
    A = abs(pars(3:end));
    gamma = fwhm/2;
    f = zeros(L,1);
    for i = 1:N
        f = f + abs(A(i))*(gamma^2./((lambda-e0(i)).^2 + gamma^2));
    end
    f = f + c;
end
