function [map1,A1,fwhm1,c1,r1] = xrf_map0(spect1,a0,a1,e0,fwhm,p0)

sigma = 0.4247*fwhm;    %KeV
e1 = e0 - 3*sigma;
e2 = e0 + 3*sigma;
lambda1 = floor((e1-a0)/a1);
lambda2 = ceil((e2-a0)/a1);    
lambda = (1:p0)*a1 + a0;
lambda = lambda(:);
lambda11 = lambda(lambda1:lambda2);
spect = spect1(lambda1:lambda2);

cest = sum(spect)/(lambda2-lambda1+1);
Aest = max(spect) - cest;
parguess = [Aest,fwhm,cest];
pars = nlinfit([lambda11;e0],spect,@gauss_func,parguess);
f = gauss_func(pars,[lambda11;e0]);
r1 = sum(abs(spect - f))/sum(spect);
f = f - pars(3);
f(f<0) = 0;
map1 = sum(f);
A1 = pars(1);
fwhm1 = pars(2);
c1 = pars(3);
end

function f = gauss_func(pars,lambda11)
    lambda = lambda11(1:(end-1));
    e0 = lambda11(end);
    A = pars(1);
    fwhm = pars(2);
    c = pars(3);
    sigma = 0.4247*fwhm;
    f = A*exp(-((lambda-e0).^2)/(2*sigma^2)) + c;
end
