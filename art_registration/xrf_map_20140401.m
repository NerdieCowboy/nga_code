function [map1,A1,fwhm1,c1] = xrf_map(trial,fnpath,a0,a1,e0,fwhm,ncores)

fn = [fnpath trial];
sigma = 0.4247*fwhm;    %KeV
e1 = e0 - 3*sigma;
e2 = e0 + 3*sigma;
lambda1 = floor((e1-a0)/a1);
lambda2 = ceil((e2-a0)/a1);
    
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
lambda11 = lambda(lambda1:lambda2);

isOpen = matlabpool('size') > 0;
if (isOpen == 1)
    matlabpool close
end
matlabpool(ncores)

map1 = zeros(m0,n0);
A1 = zeros(m0,n0);
fwhm1 = zeros(m0,n0);
c1 = zeros(m0,n0);
parfor m = 1:m0
    cube = single(multibandread(fn,[m0,n0,p0],'single',0,'bsq','ieee-le',{'Row','Direct',m},{'Band','Range',[lambda1 lambda2]}));
    %last_pars = zeros(1,4);
    for n = 1:n0
        cube1 = cube(1,n,:);
        cube1 = cube1(:);
        if (sum(cube1) ~= 0)
            if (length(cube1) == length(lambda11))
                %{
                figure(1)
                plot(lambda11,cube1,'k');
                hold on
                %}
                %if (sum(last_pars) == 0)
                    %cest = sum(cube1)/(lambda2-lambda1+1);
                    %Aest = max(cube1) - cest;
                    %parguess = [Aest,e0,fwhm,cest];
                    %parguess = [Aest,fwhm,cest];
                %else
                %    parguess = last_pars;
                %end
         
                %[pars,r,J,~,~] = nlinfit([lambda11 e0],cube1,@gauss_func,parguess);
                %ci = nlparci(pars,r,'Jacobian',J);
                %if ((pars(1) >= ci(1,1)) || (pars(1) <= ci(1,2)) || (pars(2) >= ci(2,1)) || (pars(2) <= ci(2,2)) ||  (pars(3) >= ci(3,1)) || (pars(3) <= ci(3,2)) || (pars(4) >= ci(4,1)) || (pars(4) <= ci(4,2)))
                    %last_pars = pars;
                %else
                    %last_pars = zeros(1,4);
                    cest = sum(cube1)/(lambda2-lambda1+1);
                    Aest = max(cube1) - cest;
                    %parguess = [Aest,e0,fwhm,cest];
                    parguess = [Aest,fwhm,cest];
                    pars = nlinfit(lambda11,cube1,@gauss_func,parguess);
                %end
                %{
                plot(lambda11,gauss_func(parguess,lambda11),'r')
                hold on
                plot(lambda11,gauss_func(pars,lambda11),'b')
                hold off
                %}
                f = gauss_func(pars,lambda11) - pars(4);
                f(f<0) = 0;
                map1(m,n) = sum(f);
                A1(m,n) = pars(1);
                fwhm1(m,n) = pars(3);
                c1(m,n) = pars(4);
                
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
    %e0 = pars(2);
    %fwhm = pars(3);
    fwhm = pars(2);
    %c = pars(4);
    c = pars(3);
    sigma = 0.4247*fwhm;
    f = A*exp(-((lambda-e0).^2)/(2*sigma^2)) + c;
end
%{
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
%}