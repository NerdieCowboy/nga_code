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
        xrf_table(i,:) = zeros(1,c);
    end
end
xrf_table = xrf_table(:);

peak_list = single(zeros(1024*1024*2,2));
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
        peak_list = [peak_list;single(zeros(1024*1024*2),2)];
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

A = [cal_band ones(length(cal_energy),1)];
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
