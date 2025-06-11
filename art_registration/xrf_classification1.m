function [map1,mapp,A1,fwhm1,c1,r1,element_list,LN,c,xrf_table,xray_label,a1,a0,r2] = xrf_classification1(fnpath,trial,xrf_table_fn,m1,m2,n1,n2,N,minE,maxE,a1,a0,update_coeffs,ncores)

if (ncores > 1)
    isOpen = matlabpool('size') > 0;
    if (isOpen == 1)
        matlabpool close
    end
    matlabpool(ncores)
end

w = 7;  %threshold for peak finder (bands)
hsz_lp = 3;    %~0.07 KeV
lp_mn = 5;
dcnt1 = lp_mn^2;

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
fnh = [fn '.hdr'];
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

if (update_coeffs == 1)
    [a1,a0,r2] = auto_cal(a1,a0,xrf_table,minE,maxE,N,fn,hsz_lp);
else
    r2 = NaN;
end
lambda0 = single(1:p0)*a1 + a0;
lambda0 = lambda0(:);

% spatially-smoothed data cube
mu_mn = zeros(m0,n0,p0,'single');
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
        xrf_table(i,:) = zeros(1,c,'single');
    end
end

for i = 1:L
    pass = zeros(m0,n0,'uint8'); 
    tmp = xrf_table(i,:);   %test each row of xrf table
    validE = sum(tmp>0);
    if (validE > 0)
        % only check the alphas
        for j = 1:2:c
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
%{
for i = 1:LN
    imwrite(uint8((2^8-1)*element_maps(:,:,N(i))),[fnpath trial '_map_' element_list{N(i)} '.tif'],'tif','Compression','None')
end
%}
xrf_table0 = zeros(L,c,'single');
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
            if (((xrf_table0(G2)+3*sigma0)>=(xrf_table0(G1)-3*sigma0)) && ((xrf_table0(G2)+3*sigma0)<=(xrf_table0(G1)+3*sigma0)))
                ov(G1,G2) = 1;
                ov(G2,G1) = 1;
            end
            if (((xrf_table0(G2)-3*sigma0)>=(xrf_table0(G1)-3*sigma0)) && ((xrf_table0(G2)-3*sigma0)<=(xrf_table0(G1)+3*sigma0)))
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
mapp = zeros(m0,n0,LN*c,'single');
map1 = zeros(m0,n0,LN,'single');
A1 = zeros(m0,n0,LN*c,'single');
fwhm1 = zeros(m0,n0,LN*c,'single');
c1 = zeros(m0,n0,LN*c,'single');
r1 = ones(m0,n0,LN*c,'single');
tmp = zeros(1,n0,p0,'single');
multibandwrite(single(tmp),fn_fit,'bip',[1 1 1],[m0,n0,p0]);
clear tmp
parfor m = m1:m2
    spect1 = single(multibandread(fn,[m0,n0,p0],'single',0,'bsq','ieee-le',{'Row','Direct',m}));
    %spect1 = single(mu_p(m,:,:));
    spect2 = zeros(1,n0,p0,'single');
    mapp1 = zeros(1,n0,LN*c,'single');
    A11 = zeros(1,n0,LN*c,'single');
    fwhm11 = zeros(1,n0,LN*c,'single');
    c11 = zeros(1,n0,LN*c,'single');
    r11 = ones(1,n0,LN*c,'single');
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
            [map0,A0,fwhm0,c0,r0] = xrf_map(spect,ov,a0,a1,N0,Npeak,e0,fwhm,p0,xrf_table);
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
    map11 = zeros(1,n0,LN,'single');
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

for i = 1:LN
    mx = max(max(map1(:,:,i)));
    f = figure('visible','off');
    a = axes;
    image(map1(:,:,i),'cdatamapping','scaled');
    caxis([0 mx]);
    colormap gray
    colorbar('peer',a);
    cb=colorbar;
    ylabel(cb,'sum counts');
    axis image
    axis off
    saveas(f,[fnpath trial '_sum_' element_list{N(i)} '.png'])
    close(f)
    
    imwrite(uint16(map1(:,:,i)/mx*(2^16-1)),[fnpath trial '_sum_' element_list{N(i)} '.tif'],'tif','Compression','None')
end

for i = 1:LN
    tmp = xrf_table(N(i),:);
    for j = 1:c
        if (tmp(j) ~= 0)
            ij = (i-1)*c+j;
            
            mx = max(max(mapp(:,:,ij)));
            f = figure('visible','off');
            a = axes;
            image(mapp(:,:,ij),'cdatamapping','scaled');
            caxis([0 mx]);
            colormap gray
            colorbar('peer',a);
            cb=colorbar;
            ylabel(cb,'sum counts');
            axis image
            axis off
            saveas(f,[fnpath trial '_sum_' element_list{N(i)} '_' xray_label{j} '.png'])
            close(f)
            
            imwrite(uint16(mapp(:,:,ij)/mx*(2^16-1)),[fnpath trial '_sum_' element_list{N(i)} '_' xray_label{j} '.tif'],'tif','Compression','None')
        end
    end
end

% confidence
% define overlap tighter than above
xrf_table0 = zeros(L,c,'single');
for i = 1:LN
    xrf_table0(N(i),:) = xrf_table(N(i),:);
end
[L,G] = size(xrf_table0);
xrf_table0 = xrf_table0';
xrf_table0 = xrf_table0(:);
ov = logical(eye(L*G,L*G));
for G1 = 1:(L*G-1)
    fwhm0 = xrf_table0(G1)*b1 + b0;
    for G2 = (G1+1):(L*G)
        if ((xrf_table0(G1) ~= 0) && (xrf_table0(G2) ~= 0))
            if (((xrf_table0(G2)+fwhm0/2)>=(xrf_table0(G1)-fwhm0/2)) && ((xrf_table0(G2)+fwhm0/2)<=(xrf_table0(G1)+fwhm0/2)))
                ov(G1,G2) = 1;
                ov(G2,G1) = 1;
            end
            if (((xrf_table0(G2)-fwhm0/2)>=(xrf_table0(G1)-fwhm0/2)) && ((xrf_table0(G2)-fwhm0/2)<=(xrf_table0(G1)+fwhm0/2)))
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

A0 = 2;
hsz = 3;
element_amount_map = zeros(m0,n0,LN,'single');
conf = zeros(m0,n0,LN,'single');
for i = 1:LN
    % K-alpha
    e0 = xrf_table(N(i),1);
    if (e0 ~= 0)
        fwhm0 = b0 + b1*e0;
        sigma0 = 0.4247*fwhm0;
        lambda = (e0-3*sigma0):a1:(e0+3*sigma0);
        thresh1 = sum(A0*exp(-((lambda-e0).^2)/(2*sigma0^2)));
    else
        thresh1 = 0;
    end
    
    % K-beta
    e0 = xrf_table(N(i),2);
    if (e0 ~= 0)
        fwhm0 = b0 + b1*e0;
        sigma0 = 0.4247*fwhm0;
        lambda = (e0-3*sigma0):a1:(e0+3*sigma0);
        thresh2 = sum(A0*exp(-((lambda-e0).^2)/(2*sigma0^2)));
    else
        thresh2 = 0;
    end
    
    % L-alpha
    e0 = xrf_table(N(i),5);
    if (e0 ~= 0)
        fwhm0 = b0 + b1*e0;
        sigma0 = 0.4247*fwhm0;
        lambda = (e0-3*sigma0):a1:(e0+3*sigma0);
        thresh3 = sum(A0*exp(-((lambda-e0).^2)/(2*sigma0^2)));
    else
        thresh3 = 0;
    end
    
    % L-beta
    e0 = xrf_table(N(i),6);
    if (e0 ~= 0)
        fwhm0 = b0 + b1*e0;
        sigma0 = 0.4247*fwhm0;
        lambda = (e0-3*sigma0):a1:(e0+3*sigma0);
        thresh4 = sum(A0*exp(-((lambda-e0).^2)/(2*sigma0^2)));
    else
        thresh4 = 0;
    end

    if ((thresh1 > 0) && (thresh2 > 0))
        alpha = mapp(:,:,(i-1)*c+1);
        alpha(alpha<thresh1) = 0;
        beta = mapp(:,:,(i-1)*c+2);
        beta(beta<thresh2) = 0;
        
        tmp = r1(:,:,(i-1)*c+1);
        tmp(tmp==0) = 1;
        tmp(tmp<0) = 0;
        tmp(tmp>1) = 1;
        fit_alpha = 1 - tmp;
        clear tmp
        tmp = r1(:,:,(i-1)*c+2);
        tmp(tmp==0) = 1;
        tmp(tmp<0) = 0;
        tmp(tmp>1) = 1;
        fit_beta = 1 - tmp;
        clear tmp
    elseif ((thresh3 > 0) && (thresh4 > 0))
        alpha = mapp(:,:,(i-1)*c+5);
        alpha(alpha<thresh3) = 0;
        beta = mapp(:,:,(i-1)*c+6);
        beta(beta<thresh4) = 0;
        
        thresh1 = thresh3;
        thresh2 = thresh4;
        
        tmp = r1(:,:,(i-1)*c+5);
        tmp(tmp==0) = 1;
        tmp(tmp<0) = 0;
        tmp(tmp>1) = 1;
        fit_alpha = 1 - tmp;
        clear tmp
        tmp = r1(:,:,(i-1)*c+6);
        tmp(tmp==0) = 1;
        tmp(tmp<0) = 0;
        tmp(tmp>1) = 1;
        fit_beta = 1 - tmp;
        clear tmp
    end
    
    alpha1 = (alpha>=thresh1);
    L = bwlabel(alpha1,8);
    h = hist(L(:),(max(L(:))+1));
    for j = 0:max(L(:))
        if (h(j+1) == 1)
            L(L==j) = 0;
        end
    end
    msk = (L>0);
    alpha1 = alpha1.*msk;
    clear msk
    
    beta1 = (beta>=thresh2);
    L = bwlabel(beta1,8);
    h = hist(L(:),(max(L(:))+1));
    for j = 0:max(L(:))
        if (h(j+1) == 1)
            L(L==j) = 0;
        end
    end
    msk = (L>0);
    beta1 = beta1.*msk;
    clear msk

    if (sum(beta1(:)) > 0)
        tmp = and(alpha1,beta1) + and(~alpha1,~beta1);
        h = fspecial('average',[hsz hsz]);
        spatial_corr = imfilter(tmp,h);
        clear tmp
    else
        spatial_corr = zeros(m0,n0,'single');
    end

    if ((xrf_table(N(i),1) ~= 0) || (xrf_table(N(i),2) ~= 0))
        ov_alpha_tot = sum(ov((N(i)-1)*c+1,:));
        ov_beta_tot = sum(ov((N(i)-1)*c+2,:));
    elseif ((xrf_table(N(i),5) ~= 0) || (xrf_table(N(i),6) ~= 0))
        ov_alpha_tot = sum(ov((N(i)-1)*c+5,:));
        ov_beta_tot = sum(ov((N(i)-1)*c+6,:));
    end
    fit_alpha = fit_alpha.*double(alpha1);
    fit_beta = fit_beta.*double(beta1);
    mn = min(fit_alpha,fit_beta);
    mx = max(fit_alpha,fit_beta);
    
    if ((ov_alpha_tot == 0) && (ov_beta_tot == 0))
        conf(:,:,i) = single(fit_alpha).*single(alpha>thresh1);
        element_amount_map(:,:,i) = single(alpha + beta).*single(alpha>thresh1);
    elseif (ov_alpha_tot == 0)
        conf(:,:,i) = single(fit_alpha).*single(alpha>thresh1);
        element_amount_map(:,:,i) = single(alpha).*single(alpha>thresh1);
    elseif (ov_beta_tot == 0)
        conf(:,:,i) = single(fit_beta).*single(beta>thresh2);
        element_amount_map(:,:,i) = single(beta).*single(beta>thresh2);
    else
        conf(:,:,i) = single((mx - mn).*spatial_corr + mn).*single(alpha>thresh1);
        if (mean(fit_alpha(:)) >= mean(fit_beta(:)))
            element_amount_map(:,:,i) = single(alpha).*single(alpha>thresh1);
        else
            element_amount_map(:,:,i) = single(beta).*single(beta>thresh2);
        end
    end
    
    f = figure('visible','off');
    a = axes;
    image(element_amount_map(:,:,i),'cdatamapping','scaled');
    caxis([0 max(max(element_amount_map(:,:,i)))]);
    colormap gray
    axis image
    axis off
    colorbar('peer',a);
    cb=colorbar;
    ylabel(cb,'sum counts')
    saveas(f,[fnpath trial '_amountmap_' element_list{N(i)} '.png'])
    close(f)
    
    f = figure('visible','off');
    a = axes;
    image(conf(:,:,i),'cdatamapping','scaled');
    caxis([0 1]);
    colormap gray
    axis image
    axis off
    colorbar('peer',a);
    cb=colorbar;
    ylabel(cb,'confidence')
    saveas(f,[fnpath trial '_confidence_' element_list{N(i)} '.png'])
    close(f)
end

if (ncores > 1)
    isOpen = matlabpool('size') > 0;
    if (isOpen == 1)
        matlabpool close
    end
end

end

function [map1,A1,fwhm1,c1,r1] = xrf_map(spect1,ov,a0,a1,N0,Npeak,e0,fwhm,p0,xrf_table)

%e0 -> peaks found in each spectrum
%fwhm -> estimated fwhm of each peak
%N0 -> element number corresponding to each peak

[Lxt,c] = size(xrf_table);
sigma = 0.4247*fwhm;    %KeV
N = length(e0);     %number of peaks in the spectrum
lambda = single(1:p0)*a1 + a0;        %KeV
lambda = lambda(:);

test = 1:N;
map1 = zeros(1,N,'single');
A1 = zeros(1,N,'single');
fwhm1 = zeros(1,N,'single');
sigma1 = zeros(1,N,'single');
c1 = zeros(1,N,'single');
r1 = ones(1,N,'single');
while (~isempty(test))
    tmp = logical(ov((N0(1)-1)*c+Npeak(1),:));
    tmp2 = 1:(Lxt*c);
    ind = tmp2(tmp);
    clear tmp tmp2
    e00 = [e0(1)];
    fwhm0 = [fwhm(1)];
    ovE = [test(1)];   %index in N0
    testi = [1];
    for ii = 1:length(ind)
        N01 = floor(ind(ii)/c)+1;
        Npeak1 = ind(ii)-c*floor((ind(ii)-1)/c);
        tmp1 = (N0==N01);
        tmp2 = 1:(length(N0));
        Nind = tmp2(tmp1);
        tmp1 = (Npeak==Npeak1);
        Pind = tmp2(tmp1);
        NPind = intersect(Nind,Pind);
        if (length(NPind) > 0)
            NPind = NPind(1);
            e00 = [e00 e0(NPind)];
            fwhm0 = [fwhm0 fwhm(NPind)];
            ovE = [ovE test(NPind)];
            testi = [testi NPind];
            clear Nind Pind NPind tmp1 tmp2
        end
    end
    Aest = spect1(floor((e00-a0)/a1))';
    sigma0 = 0.4247*fwhm0;
    LovE = length(ovE);
    
    sigma00 = mean(sigma0);
    e1 = min(e00) - 3*sigma00; %KeV
    e2 = max(e00) + 3*sigma00; %KeV
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
    
    %cest = 0;
    cest = min(spect);
    parguess = [mean(fwhm0),cest,Aest];
    L = length(lambda11);
    pars = nlinfit([lambda11;e00';L],spect,@gauss_func,parguess);
    %pars = nlinfit([lambda11;e00';L],spect,@lorentz_func,parguess);
    
    for ii = 1:LovE
        tmp = pars(1);
        bad = 0;
        if ((tmp > 0) && (tmp < 2*mean(fwhm0)))
            fwhm1(ovE(ii)) = tmp;
            sigma1(ovE(ii)) = 0.4247*tmp;
        else
            fwhm1(ovE(ii)) = 0;
            sigma1(ovE(ii)) = 0;
            %bad = 1;
        end 
        tmp = pars(2);
        if (tmp < 0)
            c1(ovE(ii)) = 0;
        else
            c1(ovE(ii)) = tmp;
        end
        tmp = pars(ii+2);
        if (tmp < 0)
            A1(ovE(ii)) = 0;
        else
            A1(ovE(ii)) = tmp;
        end
    
        if (bad == 0)
            e1 = e00(ii) - 2*sigma00; %KeV
            e2 = e00(ii) + 2*sigma00; %KeV
            %e1 = e0(ii) - fwhm1(1)/2; %KeV
            %e2 = e0(ii) + fwhm1(1)/2; %KeV
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
            f = A1(ovE(ii))*exp(-((lambda11-e00(ii)).^2)/(2*(sigma1(ovE(ii)))^2)) + c1(ovE(ii));
            r1(ovE(ii)) = sum(abs(spect - f))/sum(spect);
            f = f - c1(ovE(ii));
            f(f<0) = 0;
            map1(ovE(ii)) = sum(f);
        elseif (bad == 1)
            r1(ovE(ii)) = 1;
            map1(ovE(ii)) = 0;
            A1(ovE(ii)) = 0;
            fwhm1(ovE(ii)) = 0;
            sigma1(ovE(ii)) = 0;
            c1(ovE(ii)) = 0;
        end
    end
    for ii = LovE:-1:1
        test(testi(ii)) = [];
        N0(testi(ii)) = [];
        Npeak(testi(ii)) = [];
        fwhm(testi(ii)) = [];
        sigma(testi(ii)) = [];
        e0(testi(ii)) = [];
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
    f = zeros(L,1,'single');
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
    f = zeros(L,1,'single');
    for i = 1:N
        f = f + abs(A(i))*(gamma^2./((lambda-e0(i)).^2 + gamma^2));
    end
    f = f + c;
end

function [a1,a0,r2] = auto_cal(a1,a0,xrf_table,minE,maxE,N,fn,hsz_lp)

dcnt0 = 3;
dcnt1 = 30;

fnh = [fn '.hdr'];
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
fclose(fid);
lambda = single(1:p0);
lambda = lambda(:);

% find peaks (1st derivatives = 0, 2nd derivative = -)
peak_list = zeros(m0*n0,3,'single'); %[m,n,lambda]
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
        peak_list = [peak_list;zeros(m0*n0,3,'single')];
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

[~,c] = size(xrf_table);
for i = 1:c
    tmp = xrf_table(:,i);
    tmp(tmp <= minE) = 0;
    tmp(tmp > maxE) = 0;
    xrf_table(:,i) = tmp;
    clear tmp
end
for i = 1:length(N)
    if (sum(N == i) == 0)
        xrf_table(i,:) = zeros(1,c,'single');
    end
end
xrf_table = xrf_table(:);

peak_list = zeros(1024*1024*2,2,'single');
[Ltot,~] = size(peak_list);
ptot = 1;
for p = 14:(length(n)-13)
    tmp = n((p-13):(p+13));
    mx = max(tmp);
    lt = mean(tmp(5));
    md = mean(tmp(14));
    rt = mean(tmp(23));
    tmp1 = md - lt;
    tmp2 = rt - md;
    msk = (tmp1 > 0) & (tmp2 < 0) & (tmp(14) > dcnt1) & (tmp(14) == mx);
    clear tmp tmp1 tmp2
    [r,~] = find(msk==1);
    L = length(r);
    while ((ptot + L - 1) > Ltot)
        peak_list = [peak_list;zeros(1024*1024*2,2,'single')];
        [Ltot,~] = size(peak_list);
    end
    for i = 1:L
        peak_list(ptot,:) = [p xout(p)];
        ptot = ptot + 1;
    end
end
peak_list(peak_list(:,1)==0,:) = [];
%{
figure,plot(xout,n)
hold on
plot(peak_list(:,2),n(peak_list(:,1)),'rx')
hold off
%}
cal_band0 = peak_list(:,2);
cal_energy0 = cal_band0*a1 + a0;
cal_band = [];
cal_energy = [];
for i = 1:length(cal_energy0)
    % find closest match for first peak
    [d1,ind1] = min(abs(cal_energy0(i) - xrf_table));
    if (d1 < 0.1)
        % is the closest match a better fit to a different peak
        [~,ind2] = min(abs(cal_energy0 - xrf_table(ind1(1))));
        if (i == ind2(1))
            cal_band = [cal_band; cal_band0(i)];
            cal_energy = [cal_energy; xrf_table(ind1(1))];
        end
    end
end
clear peak_list cal_band0 cal_energy0 xrf_table n xout

A = [cal_band ones(length(cal_energy),1,'single')];
coeff = pinv(A)*cal_energy;
a1 = coeff(1);
a0 = coeff(2);
%{
figure,plot(cal_band,cal_energy,'k*')
hold on
plot(cal_band,a1*cal_band+a0,'r');
hold off
%}
ybar = mean(cal_energy);
sigma2 = sum((cal_energy-ybar).^2);
sse = sum((a1*cal_band+a0-cal_energy).^2);
r2 = 1 - (sse/sigma2);
clear ybar sigma2 sse A coeff
end
