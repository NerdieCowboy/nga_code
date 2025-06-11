function [A1,A2,fwhm1,fwhm2,c1,r1] = overall_xrf_map2(trial,fnpath,a0,a1,e1,e2,fwhm,ncores)

fn = [fnpath trial];
sigma = 0.4247*fwhm;
e11 = e1 - 3*sigma;
e21 = e2 + 3*sigma;
lambda1 = floor((e11-a0)./a1);
lambda2 = ceil((e21-a0)./a1);
    
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

if (ncores > 1)
    isOpen = matlabpool('size') > 0;
    if (isOpen == 1)
        matlabpool close
    end
    matlabpool(ncores)
end

A1 = zeros(m0,n0);
A2 = zeros(m0,n0);
fwhm1 = zeros(m0,n0);
fwhm2 = zeros(m0,n0);
c1 = zeros(m0,n0);
r1 = zeros(m0,n0);
if (ncores > 1)
    parfor m = 1:m0
        cube = single(multibandread(fn,[m0,n0,p0],'single',0,'bsq','ieee-le',{'Row','Direct',m}));
        for n = 1:n0
            cube1 = cube(1,n,:);
            cube1 = cube1(:);
            if (sum(cube1) ~= 0)
                if (length(cube1) == length(lambda))
                    cest = sum(cube1)/(lambda2-lambda1+1);
                    Aest = max(cube1) - cest;
                    parguess = [Aest,Aest,fwhm,fwhm,cest];
                    pars = nlinfit([lambda;e1;e2],cube1,@gauss_func,parguess);
                    f = gauss_func(pars,[lambda;e1;e2]);
                    r1(m,n) = sum(abs(cube1 - f))/sum(cube1);
                    A1(m,n) = pars(1);
                    A2(m,n) = pars(2);
                    fwhm1(m,n) = pars(3);
                    fwhm2(m,n) = pars(4);
                    c1(m,n) = pars(5);
                end
            end
        end
    end
else
    
end
clear cube cube1
if (ncores > 1)
    isOpen = matlabpool('size') > 0;
    if (isOpen == 1)
        matlabpool close
    end
end
end

function f = gauss_func(pars,lambda11)
    lambda = lambda11(1:(end-2));
    e1 = lambda11(end-1);
    e2 = lambda11(end);
    A1 = pars(1);
    A2 = pars(2); 
    fwhm1 = pars(3);
    fwhm2 = pars(4);
    c = pars(5);
    sigma1 = 0.4247*fwhm1;
    sigma2 = 0.4247*fwhm2;
    f = A1*exp(-((lambda-e1).^2)/(2*sigma1^2)) + A2*exp(-((lambda-e2).^2)/(2*sigma2^2)) + c;
end
