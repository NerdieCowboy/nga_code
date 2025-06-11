function [] = rough_mosaic(ir_fn,usegui)

if (usegui == 0)
    reply = input('Do the template images need to be aligned? Y/N [N]: ', 's');
    if isempty(reply)
        reply = 'N';
    end
else
    choice = questdlg('Do the template images need to be aligned?','Rough Mosaic','Yes','No','No');
    switch choice
        case 'Yes'
            reply = 'Y';
        case 'No'
            reply = 'N';
    end
end
if (strcmpi(reply, 'Y') == 1)
    if (usegui == 0)
        ov = input('What is the minimum overlap between images (in pixels)? [32]: ', 's');
        if isempty(ov)
            ov = '32';
        end
        ov = str2double(ov);    %overlap (pixels)
    else
        prompt = {'What is the minimum overlap between images (in pixels)?'};
        dlg_title = 'Rough Mosaic';
        num_lines = 1;
        def = {'32'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        ov = str2double(answer{1});
    end
    
warning('off','all')

numbest = 16;

tmp = regexp(ir_fn,'\\');
pc = 1;
if (isempty(tmp))
    pc = 0;
    tmp = regexp(ir_fn,'\/');
end
tmp1 = tmp(end);
tmp2 = length(ir_fn);
fn_path = ir_fn(1:tmp1);
fn = ir_fn((tmp1+1):tmp2);
clear tmp tmp1 tmp2
tmp = regexp(fn,'\.');
trial = fn(1:(tmp-1));
clear tmp

if (pc == 1)
    fn_path2 = [fn_path trial '\'];
elseif (pc == 0)
    fn_path2 = [fn_path trial '/'];
end
listOftiffs = dir(fullfile(fn_path2,'*.tif'));
p2 = numel(listOftiffs);
clear listOftiffs
inc = 0;
p0 = 0;
while (inc < p2)
    p0 = p0 + 1;
    if (pc == 1)
        fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',p0)) '.tif'];
    elseif (pc == 0)
        fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03d',p0)) '.tif'];
    end
    if (exist(fn_full2,'file') == 2)
        info = imfinfo(fn_full2);
        m0 = info.Height;
        n0 = info.Width;
        bd2 = info.BitsPerSample;
        bd2 = bd2(1);
        inc = inc + 1;
    end
end

nblk1 = floor(n0/2)-floor(2*ov)+1;
nblk2 = floor(n0/2)+ceil(2*ov);
mblk1 = floor(m0/2)-floor(2*ov)+1;
mblk2 = floor(m0/2)+ceil(2*ov);
min_overlap = min((nblk2-nblk1)*ov,(mblk2-mblk1)*ov);
nsearch1 = floor(n0/2)-floor(2*2*ov)+1;
nsearch2 = floor(n0/2)+ceil(2*2*ov);
msearch1 = floor(m0/2)-floor(2*2*ov)+1;
msearch2 = floor(m0/2)+ceil(2*2*ov);

%IR
rel_reg = zeros(p0,4);
cnt = 0;
first = 1;
for inc = 1:p0
    if (pc == 1)
        fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',inc)) '.tif'];
    elseif (pc == 0)
        fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03d',inc)) '.tif'];
    end
        
    if (exist(fn_full2,'file') == 2)
        cnt = cnt + 1;
        img2 = double(imread(fn_full2,'tif'));
        
        best = zeros(numbest*4,3);
        if (first == 0)
            %top    
            temp = img2(1:ov,nblk1:nblk2);
            [m2 n2] = size(temp);

            temp2 = img2_1(:,nsearch1:nsearch2);
            
            try 
                cc = normxcorr2(temp,temp2);
            catch
                cc = 0;
            end

            for a = 1:numbest
                [max_cc, imax] = max(cc(:));
                [ypeak, xpeak] = ind2sub(size(cc),imax(1));
                corr_offset = [ (ypeak-m2) (xpeak-n2) ];
                best(a,1) = cc(ypeak,xpeak);
                cc(ypeak,xpeak) = 0;
                best(a,2) = corr_offset(1) - 1;
                best(a,3) = nblk1 - (nsearch1 + corr_offset(2) - 1);
                if (best(a,2) < 0)
                    best(a,1) = 0;
                end
            end
            clear temp cc img1sub
            
            %bottom
            temp = img2((m0-ov+1):m0,nblk1:nblk2);
            [m2 n2] = size(temp);

            temp2 = img2_1(:,nsearch1:nsearch2);
            
            try 
                cc = normxcorr2(temp,temp2);
            catch
                cc = 0;
            end

            for a = (numbest+1):2*numbest
                [max_cc, imax] = max(cc(:));
                [ypeak, xpeak] = ind2sub(size(cc),imax(1));
                corr_offset = [ (ypeak-m2) (xpeak-n2) ];
                best(a,1) = cc(ypeak,xpeak);
                cc(ypeak,xpeak) = 0;
                best(a,2) = corr_offset(1) - (m0-ov+1);
                best(a,3) = nblk1 - (nsearch1 + corr_offset(2) - 1);
                if (best(a,2) > 0)
                    best(a,1) = 0;
                end
            end
            clear temp temp2 cc img1sub

            %left
            temp = img2(mblk1:mblk2,1:ov);
            [m2 n2] = size(temp);

            temp2 = img2_1(msearch1:msearch2,:);
            
            try 
                cc = normxcorr2(temp,temp2);
            catch
                cc = 0;
            end

            for a = (2*numbest+1):3*numbest
                [max_cc, imax] = max(cc(:));
                [ypeak, xpeak] = ind2sub(size(cc),imax(1));
                corr_offset = [ (ypeak-m2) (xpeak-n2) ];
                best(a,1) = cc(ypeak,xpeak);
                cc(ypeak,xpeak) = 0;
                best(a,2) = mblk1 - (msearch1 + corr_offset(1) - 1);
                best(a,3) = corr_offset(2) - 1;
                if (best(a,3) < 0)
                    best(a,1) = 0;
                end
            end
            clear temp cc img1sub

            %right
            temp = img2(mblk1:mblk2,(n0-ov+1):n0);
            [m2 n2] = size(temp);

            temp2 = img2_1(msearch1:msearch2,:);
            
            try 
                cc = normxcorr2(temp,temp2);
            catch
                cc = 0;
            end

            for a = (3*numbest+1):4*numbest
                [max_cc, imax] = max(cc(:));
                [ypeak, xpeak] = ind2sub(size(cc),imax(1));
                corr_offset = [ (ypeak-m2) (xpeak-n2) ];
                best(a,1) = cc(ypeak,xpeak);
                cc(ypeak,xpeak) = 0;
                best(a,2) = mblk1 - (msearch1 + corr_offset(1) - 1);
                best(a,3) = corr_offset(2) - (n0-ov+1);
                if (best(a,3) > 0)
                    best(a,1) = 0;
                end
            end
            clear temp temp2 cc img1sub
            
            best(logical(best(:,2)>(m0-1)),:) = [];
            best(logical(best(:,2)<-(m0-1)),:) = [];
            best(logical(best(:,3)>(n0-1)),:) = [];
            best(logical(best(:,3)<-(n0-1)),:) = [];
            best(logical(best(:,1)==0),:) = [];
            rel_reg(cnt,4) = 1;
            ratio = Inf*ones(1,size(best,1));
            
            temp2 = zeros(3*m0-2,3*n0-2);
            temp2(m0:(2*m0-1),n0:(2*n0-1)) = img2_1;
            msk2 = zeros(3*m0-2,3*n0-2);
            msk2(m0:(2*m0-1),n0:(2*n0-1)) = ones(m0,n0);

            for a = 1:size(best,1)
                temp = zeros(3*m0-2,3*n0-2);
                temp((m0+best(a,2)):(2*m0+best(a,2)-1),(n0+best(a,3)):(2*n0+best(a,3)-1)) = img2;
                msk = zeros(3*m0-2,3*n0-2);
                msk((m0+best(a,2)):(2*m0+best(a,2)-1),(n0+best(a,3)):(2*n0+best(a,3)-1)) = ones(m0,n0);

                mskx = msk.*msk2;
                sigma1 = std(temp(logical(mskx==1)));
                sigma2 = std(temp2(logical(mskx==1)));
                mu = mean(abs(temp(logical(mskx==1))-temp2(logical(mskx==1))));
                ratio(a) = mu^2/(sigma1*sigma2);
                
                if (sum(mskx(:)) < min_overlap)
                    ratio(a) = Inf;
                end
                clear temp msk
            end
            clear temp2 msk2

            if (~isempty(ratio))
                [mn ind] = min(ratio);
                rel_reg(cnt,1) = inc;
                rel_reg(cnt,2) = best(ind(1),2);
                rel_reg(cnt,3) = best(ind(1),3);

                test = zeros((3*m0-2),(3*n0-2),3);
                mx = max(img2_1(:));
                img2_1 = img2_1/mx;
                low_high = stretchlim(img2_1,[0.01 0.99]);
                img2_1 = 255*imadjust(img2_1,low_high,[0 1],1);
                test(m0:(2*m0-1),n0:(2*n0-1),1) = img2_1;
                mx = max(img2(:));
                img2_3 = img2/mx;
                low_high = stretchlim(img2_3,[0.01 0.99]);
                img2_3 = 255*imadjust(img2_3,low_high,[0 1],1);
                test((m0+best(ind(1),2)):(2*m0+best(ind(1),2)-1),(n0+best(ind(1),3)):(2*n0+best(ind(1),3)-1),3) = img2_3;
                clear img2_3
                if (usegui == 0)
                    imshow(uint8(test))
                elseif (usegui == 1)
                    handles = guidata(ArtRegister);
                    axes(handles.axes1)
                    imshow(uint8(test))
                    %set(handles.running,'String','Running ...');
                elseif (usegui == 2)
                    handles = guidata(HSIRegister);
                    axes(handles.axes1)
                    imshow(uint8(test))
                    set(handles.running,'String','Running ...');
                    drawnow;
                    pause(0.1)
                end
            else
                rel_reg(cnt,1) = inc;
                rel_reg(cnt,2) = 0;
                rel_reg(cnt,3) = 0;
                rel_reg(cnt,4) = 0;
            end
            pause(0.1)
        else
            rel_reg(cnt,1) = inc;
            rel_reg(cnt,2) = 0;
            rel_reg(cnt,3) = 0;
            rel_reg(cnt,4) = 1;
            first = 0;
        end
        img2_1 = img2;
        clear img2
        %rel_reg(cnt,:)
        pause(0.1)
    end
end
clear img2_1
rel_reg = rel_reg((rel_reg(:,1)~=0),:);

fn_full2 = [fn_path2 'offset_values.csv'];
csvwrite(fn_full2,rel_reg)

end

if (usegui == 0)
    reply = input('Does the alignment need to be corrected? Y/N [N]: ', 's');
    input('Press return to continue.');
else
    choice = questdlg('Does the alignment need to be corrected?','Fix Mosaic','Yes','No','No');
    switch choice
        case 'Yes'
            reply = 'Y';
        case 'No'
            reply = 'N';
    end
end
if isempty(reply)
    reply = 'N';
end
if (strcmpi(reply, 'Y') == 1)
    fix_mosaic
end
if (usegui == 2)
    handles = guidata(HSIRegister);
    set(handles.running,'String','');
    drawnow;
    pause(0.1)
end

function C = normxcorr2x(varargin)
%Mathworks function modified to force spatial correlation
%NORMXCORR2 Normalized two-dimensional cross-correlation.
%   C = NORMXCORR2(TEMPLATE,A) computes the normalized cross-correlation of
%   matrices TEMPLATE and A. The matrix A must be larger than the matrix
%   TEMPLATE for the normalization to be meaningful. The values of TEMPLATE
%   cannot all be the same. The resulting matrix C contains correlation
%   coefficients and its values may range from -1.0 to 1.0.
%
%   Class Support
%   -------------
%   The input matrices can be numeric. The output matrix C is double.
%
%   Example
%   -------
%   template = .2*ones(11); % make light gray plus on dark gray background
%   template(6,3:9) = .6;   
%   template(3:9,6) = .6;
%   BW = template > 0.5;      % make white plus on black background
%   figure, imshow(BW), figure, imshow(template)
%   % make new image that offsets the template
%   offsetTemplate = .2*ones(21); 
%   offset = [3 5];  % shift by 3 rows, 5 columns
%   offsetTemplate( (1:size(template,1))+offset(1),...
%                   (1:size(template,2))+offset(2) ) = template;
%   figure, imshow(offsetTemplate)
%   
%   % cross-correlate BW and offsetTemplate to recover offset  
%   cc = normxcorr2(BW,offsetTemplate); 
%   [max_cc, imax] = max(abs(cc(:)));
%   [ypeak, xpeak] = ind2sub(size(cc),imax(1));
%   corr_offset = [ (ypeak-size(template,1)) (xpeak-size(template,2)) ];
%   isequal(corr_offset,offset) % 1 means offset was recovered
%
%  See also CORRCOEF.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.12.4.10 $  $Date: 2006/12/10 20:03:20 $

%   Input-output specs
%   ------------------
%   T:    2-D, real, full matrix
%         logical, uint8, uint16, or double
%         no NaNs, no Infs
%         prod(size(T)) >= 2
%         std(T(:))~=0
%
%   A:    2-D, real, full matrix
%         logical, uint8, uint16, or double
%         no NaNs, no Infs
%         size(A,1) >= size(T,1)
%         size(A,2) >= size(T,2)
%
%   C:    double

[T, A] = ParseInputs(varargin{:});

%   We normalize the cross correlation to get correlation coefficients using the
%   definition of Haralick and Shapiro, Volume II (p. 317), generalized to
%   two-dimensions. 
%
%   Lewis explicitly defines the normalized cross-correlation in two-dimensions
%   in this paper (equation 2):
%   
%      "Fast Normalized Cross-Correlation", by J. P. Lewis, Industrial Light & Magic.
%      http://www.idiom.com/~zilla/Papers/nvisionInterface/nip.html
%
%   Our technical reference document on NORMXCORR2 shows how to get from
%   equation 2 of the Lewis paper to the code below.

xcorr_TA = xcorr2_fast(T,A);

[m n] = size(T);
mn = m*n;

local_sum_A = local_sum(A,m,n);
local_sum_A2 = local_sum(A.*A,m,n);

% Note: diff_local_sums should be nonnegative, but may have negative
% values due to round off errors. Below, we use max to ensure the
% radicand is nonnegative.
diff_local_sums = ( local_sum_A2 - (local_sum_A.^2)/mn );
denom_A = sqrt( max(diff_local_sums,0) ); 

denom_T = sqrt(mn-1)*std(T(:));
denom = denom_T*denom_A;
numerator = (xcorr_TA - local_sum_A*sum(T(:))/mn );

% We know denom_T~=0 from input parsing;
% so denom is only zero where denom_A is zero, and in 
% these locations, C is also zero.
C = zeros(size(numerator));
tol = 1000*eps( max(abs(denom(:))) );
i_nonzero = find(denom > tol);
C(i_nonzero) = numerator(i_nonzero) ./ denom(i_nonzero);

%-------------------------------
% Function  local_sum
%
function local_sum_A = local_sum(A,m,n)

% We thank Eli Horn for providing this code, used with his permission,
% to speed up the calculation of local sums. The algorithm depends on
% precomputing running sums as described in "Fast Normalized
% Cross-Correlation", by J. P. Lewis, Industrial Light & Magic.
% http://www.idiom.com/~zilla/Papers/nvisionInterface/nip.html

B = padarray(A,[m n]);
s = cumsum(B,1);
c = s(1+m:end-1,:)-s(1:end-m-1,:);
s = cumsum(c,2);
local_sum_A = s(:,1+n:end-1)-s(:,1:end-n-1);

%-------------------------------
% Function  xcorr2_fast
%
function cross_corr = xcorr2_fast(T,A)

T_size = size(T);
A_size = size(A);
outsize = A_size + T_size - 1;

% figure out when to use spatial domain vs. freq domain
conv_time = time_conv2(T_size,A_size); % 1 conv2
fft_time = 3*time_fft2(outsize); % 2 fft2 + 1 ifft2
%{
if (conv_time < fft_time)
    cross_corr = conv2(rot90(T,2),A);
else
    cross_corr = freqxcorr(T,A,outsize);
end
%}
sf = 0;
if (sf == 0)
    cross_corr = conv2(rot90(T,2),A);
end
if (sf == 1)
    cross_corr = freqxcorr(T,A,outsize);
end

%-------------------------------
% Function  freqxcorr
%
function xcorr_ab = freqxcorr(a,b,outsize)
  
% calculate correlation in frequency domain
Fa = fft2(rot90(a,2),outsize(1),outsize(2));
Fb = fft2(b,outsize(1),outsize(2));
xcorr_ab = real(ifft2(Fa .* Fb));


%-------------------------------
% Function  time_conv2
%
function time = time_conv2(obssize,refsize)

% time a spatial domain convolution for 10-by-10 x 20-by-20 matrices

% a = ones(10);
% b = ones(20);
% mintime = 0.1;

% t1 = cputime;
% t2 = t1;
% k = 0;
% while (t2-t1)<mintime
%     c = conv2(a,b);
%     k = k + 1;
%     t2 = cputime;
% end
% t_total = (t2-t1)/k;

% % convolution time = K*prod(size(a))*prod(size(b))
% % t_total = K*10*10*20*20 = 40000*K
% K = t_total/40000;

% K was empirically calculated by the commented-out code above.
K = 2.7e-8; 
            
% convolution time = K*prod(obssize)*prod(refsize)
time =  K*prod(obssize)*prod(refsize);


%-------------------------------
% Function  time_fft2
%
function time = time_fft2(outsize)

% time a frequency domain convolution by timing two one-dimensional ffts

R = outsize(1);
S = outsize(2);

% Tr = time_fft(R);
% K_fft = Tr/(R*log(R)); 

% K_fft was empirically calculated by the 2 commented-out lines above.
K_fft = 3.3e-7; 
Tr = K_fft*R*log(R);

if S==R
    Ts = Tr;
else
%    Ts = time_fft(S);  % uncomment to estimate explicitly
   Ts = K_fft*S*log(S); 
end

time = S*Tr + R*Ts;

% %-------------------------------
% % Function time_fft
% %
% function T = time_fft(M)

% % time a complex fft that is M elements long

% vec = complex(ones(M,1),ones(M,1));
% mintime = 0.1; 

% t1 = cputime;
% t2 = t1;
% k = 0;
% while (t2-t1) < mintime
%     dummy = fft(vec);
%     k = k + 1;
%     t2 = cputime;
% end
% T = (t2-t1)/k;


%-----------------------------------------------------------------------------
function [T, A] = ParseInputs(varargin)

iptchecknargin(2,2,nargin,mfilename)

T = varargin{1};
A = varargin{2};

iptcheckinput(T,{'logical','numeric'},{'real','nonsparse','2d','finite'},mfilename,'T',1)
iptcheckinput(A,{'logical','numeric'},{'real','nonsparse','2d','finite'},mfilename,'A',2)

checkSizesTandA(T,A)

% See geck 342320. If either A or T has a minimum value which is negative, we
% need to shift the array so all values are positive to ensure numerically
% robust results for the normalized cross-correlation.
A = shiftData(A);
T = shiftData(T);

checkIfFlat(T);

%-----------------------------------------------------------------------------
function B = shiftData(A)

B = double(A);

is_unsigned = isa(A,'uint8') || isa(A,'uint16') || isa(A,'uint32');
if ~is_unsigned
    
    min_B = min(B(:)); 
    
    if min_B < 0
        B = B - min_B;
    end
    
end

%-----------------------------------------------------------------------------
function checkSizesTandA(T,A)

if numel(T) < 2
    eid = sprintf('Images:%s:invalidTemplate',mfilename);
    msg = 'TEMPLATE must contain at least 2 elements.';
    error(eid,'%s',msg);
end

if size(A,1)<size(T,1) || size(A,2)<size(T,2) 
    eid = sprintf('Images:%s:invalidSizeForA',mfilename);
    msg = 'A must be the same size or larger than TEMPLATE.';
    error(eid,'%s',msg);
end

%-----------------------------------------------------------------------------
function checkIfFlat(T)

if std(T(:)) == 0
    eid = sprintf('Images:%s:sameElementsInTemplate',mfilename);
    msg = 'The values of TEMPLATE cannot all be the same.';
    error(eid,'%s',msg);
end
