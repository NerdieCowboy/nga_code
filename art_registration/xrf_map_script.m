function [map1,A1,fwhm1,c1,r1] = xrf_map(spect1,ov,c,a0,a1,N0,Npeak,e0,fwhm,p0,xrf_table)

%e0 -> peaks found in each spectrum
%fwhm -> estimated fwhm of each peak
%N0 -> element number corresponding to each peak

sigma = 0.4247*fwhm;    %KeV
N = length(e0);     %number of peaks in the spectrum
lambda = single(1:p0)*a1 + a0;        %KeV
lambda = lambda(:);

map1 = zeros(1,N);
A1 = zeros(1,N);
fwhm1 = zeros(1,N);
c1 = zeros(1,N);
r1 = zeros(1,N);
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
    N00 = length(e00);
    
    %cest = sum(spect)/(lambda2-lambda1+1);
    cest = 0;
    %Aest = (max(spect) - cest)*ones(1,N00);
    %Aest = max(spect)*ones(1,N00);
    parguess = [fwhm(i),cest,Aest];
    L = length(lambda11);
    pars = nlinfit([lambda11;e00';L],spect,@gauss_func,parguess);
    %pars = nlinfit([lambda11;e00';L],spect,@lorentz_func,parguess);
    tmp = pars(1);
    %if ((tmp > 0) && (tmp < 0.5))
    if (tmp > 0)
        fwhm1(i) = tmp;
    else
        fwhm1(i) = 0;
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
        f = f + A(i)*exp(-((lambda-e0(i)).^2)/(2*sigma^2));
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
        f = f + A(i)*(gamma^2./((lambda-e0(i)).^2 + gamma^2));
    end
    f = f + c;
end