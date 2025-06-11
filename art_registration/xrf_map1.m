function [map1,A1,fwhm1,c1,r1] = xrf_map1(spect1,a0,a1,e0,fwhm,p0)

sigma = 0.4247*fwhm;    %KeV
N = length(e0);
e1 = min(e0) - 2*sigma; %KeV
e2 = max(e0) + 2*sigma; %KeV
lambda1 = floor((e1-a0)/a1);    %channel
lambda2 = ceil((e2-a0)/a1);     %channel
lambda = (1:p0)*a1 + a0;        %KeV
lambda = lambda(:);
lambda11 = lambda(lambda1:lambda2); %KeV
spect = spect1(lambda1:lambda2);

cest = sum(spect)/(lambda2-lambda1+1);
Aest = (max(spect) - cest)*ones(1,N);
parguess = [fwhm,cest,Aest];
L = length(lambda11);
pars = nlinfit([lambda11;e0';L],spect,@gauss_func,parguess);
f = gauss_func(pars,[lambda11;e0';L]);
r1 = sum(abs(spect - f))/sum(spect);
c1 = pars(2);
f = f - c1;
f(f<0) = 0;
%map1 = sum(f);
fwhm1 = pars(1);
sigma1 = 0.4247*fwhm1;
A1 = pars(3:end);
map1 = zeros(1,N);
for i = 1:N
    map1(i) = sum(A1(i)*exp(-((lambda11-e0(i)).^2)/(2*sigma1^2)));
end
end

function f = gauss_func(pars,lambda11)
    L = lambda11(end);
    lambda = lambda11(1:L);
    e0 = lambda11((L+1):(end-1));
    N = length(e0);
    fwhm = pars(1);
    c = pars(2);
    A = pars(3:end);
    sigma = 0.4247*fwhm;
    f = zeros(L,1);
    for i = 1:N
        f = f + A(i)*exp(-((lambda-e0(i)).^2)/(2*sigma^2));
    end
    f = f + c;
end
