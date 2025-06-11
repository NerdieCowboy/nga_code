function [map1,A1,fwhm1,c1,map21,A21,fwhm21,c2,map22,A22,fwhm22] = xrf_map(trial,fnpath,a0,a1,e0,e11,e12,fwhm,ncores)

fn = [fnpath trial];
%fn = 'ApostlePaulXRF_full';
%a0 = 1.345;
%a1 = 0.0137;

%e0 = 6.403;     %Fe K-alpha (KeV)
%fwhm = 0.2; %keV
sigma = 0.4247*fwhm;    %KeV
e1 = e0 - 3*sigma;
e2 = e0 + 3*sigma;
e3 = e11 - 3*sigma;
e4 = e12 + 3*sigma;
lambda1 = floor((e1-a0)/a1);
lambda2 = ceil((e2-a0)/a1);
lambda3 = floor((e3-a0)/a1);
lambda4 = ceil((e4-a0)/a1);
    
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
fclose(fid);
lambda = (1:p0)*a1 + a0;
lambda = lambda(:);
lambda11 = lambda(lambda1:lambda2);   %keV
lambda12 = lambda(lambda3:lambda4);   %keV

isOpen = matlabpool('size') > 0;
if (isOpen == 1)
    matlabpool close
end
%ncores = feature('numCores');
matlabpool(ncores)

map1 = zeros(m0,n0);
A1 = zeros(m0,n0);
fwhm1 = zeros(m0,n0);
c1 = zeros(m0,n0);
map21 = zeros(m0,n0);
A21 = zeros(m0,n0);
fwhm21 = zeros(m0,n0);
c2 = zeros(m0,n0);
map22 = zeros(m0,n0);
A22 = zeros(m0,n0);
fwhm22 = zeros(m0,n0);
parfor m = 1:m0
    cube = single(multibandread(fn,[m0,n0,p0],'single',0,'bsq','ieee-le',{'Row','Direct',m},{'Band','Range',[lambda1 lambda2]}));
    cube2 = single(multibandread(fn,[m0,n0,p0],'single',0,'bsq','ieee-le',{'Row','Direct',m},{'Band','Range',[lambda3 lambda4]}));
    last_pars = zeros(1,4);
    %last_pars2 = zeros(1,7);
    for n = 1:n0
        cube1 = cube(1,n,:);
        cube1 = cube1(:);
        cube21 = cube2(1,n,:);
        cube21 = cube21(:);
        if (sum(cube1) ~= 0)
            if (length(cube1) == length(lambda11))
                %figure(1)
                %plot(lambda12,cube21,'k');
                %hold on

                cest2 = sum(cube21)/(lambda4-lambda3+1);
                Aest2 = max(cube1) - cest2;
                parguess2 = [Aest2,e11,fwhm,cest2,Aest2,e12,fwhm];
                if (sum(last_pars) == 0)
                    cest = sum(cube1)/(lambda2-lambda1+1);
                    Aest = max(cube1) - cest;
                    parguess = [Aest,e0,fwhm,cest];
                else
                    parguess = last_pars;
                    %parguess2 = last_pars2;
                end
         
                %plot(lambda,gauss_func(parguess,lambda),'r')
                %hold on
                pars = nlinfit(lambda11,cube1,@gauss_func,parguess);
                pars2 = nlinfit(lambda12,cube21,@gauss_func2,parguess2);
                pars21 = pars2(1:4);
                pars22 = [pars2(5:7) pars2(4)];
                %plot(lambda12,gauss_func(pars21,lambda12),'b')
                %hold on
                %plot(lambda12,gauss_func(pars22,lambda12),'g')
                %hold off

                f = gauss_func(pars,lambda11);
                map1(m,n) = sum(f - pars(4));
                A1(m,n) = pars(1);
                fwhm1(m,n) = pars(3);
                c1(m,n) = pars(4);
                
                f21 = gauss_func(pars21,lambda12);
                f22 = gauss_func(pars22,lambda12);
                                
                map21(m,n) = sum(f21 - pars2(4));
                A21(m,n) = pars2(1);
                fwhm21(m,n) = pars2(3);
                c2(m,n) = pars2(4);
                map22(m,n) = sum(f22 - pars2(4));
                A22(m,n) = pars2(5);
                fwhm22(m,n) = pars2(7);
                
                last_pars = pars;
                %last_pars2 = pars2;
                %figure(2)
                %imagesc(map1)
            end
            %pause(0.1)
        end
    end
end
clear cube cube1
isOpen = matlabpool('size') > 0;
if (isOpen == 1)
    matlabpool close
end

end

function f = gauss_func(pars,lambda)
    A = pars(1);
    e0 = pars(2);
    fwhm = pars(3);
    c = pars(4);
    sigma = 0.4247*fwhm;
    f = A*exp(-((lambda-e0).^2)/(2*sigma^2)) + c;
end

function f = gauss_func2(pars,lambda)
    A1 = pars(1);
    e01 = pars(2);
    fwhm1 = pars(3);
    c = pars(4);
    A2 = pars(5);
    e02 = pars(6);
    fwhm2 = pars(7);
    sigma1 = 0.4247*fwhm1;
    sigma2 = 0.4247*fwhm2;
    f = A1*exp(-((lambda-e01).^2)/(2*sigma1^2)) + A2*exp(-((lambda-e02).^2)/(2*sigma2^2)) + c;
end
