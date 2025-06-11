function [map1,A1,fwhm1,c1,r1] = overall_xrf_map(trial,fnpath,a0,a1,e0,fwhm,ncores)

fn = [fnpath trial];
sigma = 0.4247*fwhm;    %KeV
L = length(e0);
%e1 = e0 - 3*sigma;
%e2 = e0 + 3*sigma;
%lambda1 = floor((e1-a0)./a1);
%lambda2 = ceil((e2-a0)./a1);
    
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

map1 = zeros(m0,n0,L);
A1 = zeros(m0,n0,L);
fwhm1 = zeros(m0,n0,L);
c1 = zeros(m0,n0,L);
r1 = zeros(m0,n0,L);
if (ncores > 1)
    parfor m = 1:m0
        cube = single(multibandread(fn,[m0,n0,p0],'single',0,'bsq','ieee-le',{'Row','Direct',m}}));
        for n = 1:n0
            cube1 = cube(1,n,:);
            cube1 = cube1(:);
            if (sum(cube1) ~= 0)
                if (length(cube1) == length(lambda))
                    cest = sum(cube1)/(lambda2-lambda1+1);
                    Aest = max(cube1) - cest;
                    parguess = [Aest,fwhm,cest];
                    pars = nlinfit([lambda;e0],cube1,@gauss_func,parguess);
                    f = gauss_func(pars,[lambda;e0]);
                    r1(m,n) = sum(abs(cube1 - f))/sum(cube1);
                    f = f - pars(3);
                    f(f<0) = 0;
                    map1(m,n) = sum(f);
                    A1(m,n) = pars(1);
                    fwhm1(m,n) = pars(2);
                    c1(m,n) = pars(3);
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
    lambda = lambda11(1:(end-1));
    e0 = lambda11(end);
    A = pars(1);
    fwhm = pars(2);
    c = pars(3);
    sigma = 0.4247*fwhm;
    f = A*exp(-((lambda-e0).^2)/(2*sigma^2)) + c;
end
