clear all
close all
clc

fn_path = 'C:\data\Pacino\RIS\';
trial = 'sub_sc_crop_cal_Pacino_ChristMajesty_vnir_xnir';

fn = [fn_path trial];
lambdaon = 0;
lambda = [];
fnh = [fn '.hdr'];
fid = fopen(fnh);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
    if (lambdaon == 0)
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
        [~,~,e] = regexp(line,'^datatype=','match','start','end');
        if (e>0)
            datatype = str2double(line((e+1):end));
            if (datatype == 1)
                datatype = 'uint8';
                elseif (datatype == 4)
                    datatype = 'single';
                elseif (datatype == 5)
                    datatype = 'double';
                elseif (datatype == 12)
                    datatype = 'uint16';
            end
        end
        [~,~,e] = regexp(line,'^interleave=','match','start','end');
        if (e>0)
            interleave = line((e+1):end);
        end
        [~,~,e] = regexp(line,'^byteorder=','match','start','end');
        if (e>0)
            byteorder = str2double(line((e+1):end));
            if (byteorder == 0)
                byteorder = 'ieee-le';
            elseif (byteorder == 1)
                byteorder = 'ieee-be';
            end
        end
        [~,~,e] = regexp(line,'^headeroffset=','match','start','end');
        if (e>0)
            headeroffset = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^wavelength={','match','start','end');
        if (e>0)
            lambdaon = 1;
            tmp1 = line((e+1):end);
            tmp2 = regexp(tmp1,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda = [lambda str2double(tmp2(i))];
                end
                [~,~,e1] = regexp(tmp2(i),'}$','match','start','end');
                if (e1{1} > 0)
                    tmp3 = cell2mat(tmp2(i));
                    tmp3 = tmp3(1:(e1{1}-1));
                    lambda = [lambda str2double(tmp3)];
                end
            end
            [~,~,e] = regexp(line,'}$','match','start','end');
            if (e>0)
                lambdaon = 0;
            end
        end
    end
    if (lambdaon == 1)
        [~,~,e] = regexp(line,'}$','match','start','end');
        if (e>0)
            lambdaon = 0;
            tmp1 = line(1:(e-1));
            tmp2 = regexp(tmp1,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda = [lambda str2double(tmp2(i))];
                end
            end
        else
            tmp2 = regexp(line,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda = [lambda str2double(tmp2(i))];
                end
            end
        end
    end
    clear msk tmp1 tmp2 line tmp3
end
fclose(fid);
data = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
data(:,:,(lambda<360)) = [];
data(:,:,(lambda>830)) = [];
lambda(lambda<360) = [];
lambda(lambda>830) = [];
p0 = length(lambda);

[cie] = getCieStruct(lambda);
d50 = cieIllD( [5000], cie );
d50(isnan(d50)) = 0;

data = shiftdim(data,2);
data = data(:);
data = reshape(data,[p0 m0*n0]);

XYZ = ref2XYZ(data,cie.cmf10deg,d50);
%{
imgXYZ = shiftdim(XYZ,1);
imgXYZ = reshape(imgXYZ,[m0 n0 3]);
figure,imagesc(imgXYZ(:,:,1))
axis image
colormap gray
figure,imagesc(imgXYZ(:,:,2))
axis image
colormap gray
figure,imagesc(imgXYZ(:,:,3))
axis image
colormap gray
%}
clear data
imgLab = XYZ2Lab(XYZ,whitepoint('d50')*100);
clear XYZ

%convert 3 x m0n0 to m0 x n0 x 3
imgLab = shiftdim(imgLab,1);
imgLab = reshape(imgLab,[m0 n0 3]);
%{
figure,imagesc(imgLab(:,:,1))
axis image
colormap gray
figure,imagesc(imgLab(:,:,2))
axis image
colormap gray
figure,imagesc(imgLab(:,:,3))
axis image
colormap gray
ab = sqrt(imgLab(:,:,2).^2 + imgLab(:,:,3).^2);
figure,imagesc(ab)
axis image
colormap gray

brit = imgLab(:,:,1);
mxL = max(brit(:));
low_high =  stretchlim(brit/mxL,[0.005 0.995]);
brit1 = imadjust(brit/mxL,low_high,[]);

mxab = max(ab(:));
low_high =  stretchlim(ab/mxab,[0.005 0.995]);
ab1 = imadjust(ab/mxab,low_high,[]);
ab1 = abs(1 - ab1);

figure,imagesc(brit1.*ab1)
axis image
colormap gray
%}
[R,G,B] = Lab2RGB(imgLab(:,:,1),imgLab(:,:,2),imgLab(:,:,3));
%{
Bmat = zeros(3,m0*n0);
Bmat(1,:) = R(:);
Bmat(2,:) = G(:);
Bmat(3,:) = B(:);
a = [1;1;1];
P = a*a'/(a'*a);
tmp = (P*Bmat - Bmat);
dist = sqrt(tmp(1,:).^2 + tmp(2,:).^2 + tmp(3,:).^2);
clear tmp
dist = reshape(dist,m0,n0);
%}
clear imgLab
RGB = uint16(zeros(m0,n0,3));
RGB(:,:,1) = uint16((2^16-1)*R);
clear R
RGB(:,:,2) = uint16((2^16-1)*G);
clear G
RGB(:,:,3) = uint16((2^16-1)*B);
clear B

cmf = cie.cmf2deg;
figure,plot(lambda,cmf(:,1),'r')
hold on
plot(lambda,cmf(:,2),'g')
hold on
plot(lambda,cmf(:,3),'b')
hold off

figure,imshow(RGB)
imwrite(RGB,'pacino_RIS2RGB.tif','tif','Compression','None');

RGB = double(RGB);
I = mean(RGB,3);
I = I/(2^16-1);
Imx = max(I(:));
R0 = RGB(:,:,1)./I;
G0 = RGB(:,:,2)./I;
B0 = RGB(:,:,3)./I;
clear RGB
%R0(I <= 0.1*Imx) = 0;   %hue variations are not perceivable at very low luminance (Itti et al. 1998)
%G0(I <= 0.1*Imx) = 0;
%B0(I <= 0.1*Imx) = 0;
clear I
R1 = R0 - (G0 + B0)/2;
G1 = G0 - (R0 + B0)/2;
B1 = B0 - (R0 + G0)/2;
Y1 = -(abs(R0 - G0)/2 + B1);
clear R0 G0 B0
R1(R1<0) = 0; %negative values are set to 0 (Itti et al. 1998)
G1(G1<0) = 0;
B1(B1<0) = 0;
Y1(Y1<0) = 0;
figure,imagesc(R1)
axis image
colormap gray
figure,imagesc(G1)
axis image
colormap gray
figure,imagesc(B1)
axis image
colormap gray
figure,imagesc(Y1)
axis image
colormap gray
%{
RG = R1 - G1;
YB = Y1 - B1;
figure,imagesc(RG)
axis image
colormap gray
figure,imagesc(YB)
axis image
colormap gray
%}