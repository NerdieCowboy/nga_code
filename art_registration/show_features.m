function [img1] = show_features(img1,j)
%% total_features

[m0 n0 p0] = size(img1);
if (p0 == 3)
    img = double(rgb2gray(img1));
else
    img = double(img1);
end

h = fspecial('unsharp',0.8);
img = imfilter(img,h);

% Wavelet transform
h = fspecial('average',[1 2^j]);    %hsise = 1 x 2^j
LPx = imfilter(img,h);

h = [];
for q = 1:j
    h = [-1 h];
    h = [h 1];
end
LH = imfilter(LPx,h');
clear LPx;
HPy = imfilter(img,h);
clear img
h = fspecial('average',[2^j 1]);    %hsize = 2^j x 1
HL = imfilter(HPy,h);
clear HPy;

% Modulus of the wavelet transform
M = sqrt(LH.^2+HL.^2);
clear LH;
clear HL;
F0 = imregionalmax(M,8);

% Remove border pixels
hw = (2^j);
F0(1:hw,:) = 0;
F0(:,1:hw) = 0;
F0((m0-hw+1):m0,:) = 0;
F0(:,(n0-hw+1):n0) = 0;

[y1,x1] = find(F0==1);
yd = length(y1);
for i = 1:yd
    tmp = F0((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)).*M((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j));
    [mx11,ind11] = max(tmp(:));
    if (mx11 > 0)
        xpeak = floor((ind11(1)-1)/size(tmp,1))+1;
        ypeak = ind11(1) - (xpeak-1)*size(tmp,1);
        msk11 = zeros(size(tmp));
        msk11(ypeak,xpeak) = 1;
        F0((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)) = F0((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)).*msk11;
    end
    clear tmp
end
clear x1 y1
[y1,x1] = find(F0==1);
yd = length(y1);

if (p0 == 1)
    img1 = 0.5*double(img1);
    if (max(img1(:)) > 255)
        img1 = img1/(2^16-1);
    else
        img1 = img1/(2^8-1);
    end
    low_high = stretchlim(img1(:),[0.001 0.999]);
    img1 = 0.5*(2^16-1)*imadjust(img1,low_high,[0.25 1],1);
    for i = 1:yd
        img1((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j)) = 2*img1((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j));
        %img1((y1(i)-2):(y1(i)+2),(x1(i)-2):(x1(i)+2)) = 0;
    end
elseif (p0 == 3)
    hsv = rgb2hsv(double(img1));
    hsv(:,:,3) = 0.5*hsv(:,:,3);
    for i = 1:yd
        hsv((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j),3) = 2*hsv((y1(i)-2^j):(y1(i)+2^j),(x1(i)-2^j):(x1(i)+2^j),3);
    end
    img1 = hsv2rgb(hsv);
    clear hsv
end

imshow(uint16(img1))

end
