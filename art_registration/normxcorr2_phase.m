function C = normxcorr2_phase(varargin)
%   Copyright 1993-2010 The MathWorks, Inc.
%   $Revision: 1.12.4.17 $  $Date: 2011/08/09 17:51:30 $
%      "Fast Normalized Cross-Correlation", by J. P. Lewis, Industrial Light & Magic.
%      http://www.idiom.com/~zilla/Papers/nvisionInterface/nip.html

[T, A] = ParseInputs(varargin{:});
T_size = size(T);
A_size = size(A);
outsize = A_size + T_size - 1;
[xcorr_TA T A] = freqxcorr(T,A,outsize);

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
tol = sqrt( eps( max(abs(denom(:)))) );
i_nonzero = find(denom > tol);
C(i_nonzero) = numerator(i_nonzero) ./ denom(i_nonzero);

% Another numerics backstop. If any of the coefficients are outside the
% range [-1 1], the numerics are unstable to small variance in A or T. In
% these cases, set C to zero to reflect undefined 0/0 condition.
C( ( abs(C) - 1 ) > sqrt(eps(1)) ) = 0;

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
% Function  freqxcorr
%
function [xcorr_ab Fa Fb] = freqxcorr(a,b,outsize)
  
% calculate correlation in frequency domain
Fa = fft2(rot90(a,2),outsize(1),outsize(2));
Fa = angle(Fa);
Fb = fft2(b,outsize(1),outsize(2));
Fb = angle(Fb);
CC = Fa + Fb;
CC = exp(1i*CC);
xcorr_ab = real(ifft2(CC));
Fa = fft2(rot90(a,2));
Fa = angle(Fa);
Fb = fft2(b);
Fb = angle(Fb);
Fa = real(ifft2(exp(1i*Fa)));
Fb = real(ifft2(exp(1i*Fb)));


%-----------------------------------------------------------------------------
function [T, A] = ParseInputs(varargin)

narginchk(2,2)

T = varargin{1};
A = varargin{2};

validateattributes(T,{'logical','numeric'},{'real','nonsparse','2d','finite'},mfilename,'T',1)
validateattributes(A,{'logical','numeric'},{'real','nonsparse','2d','finite'},mfilename,'A',2)

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
    error(message('images:normxcorr2:invalidTemplate'))
end

if size(A,1)<size(T,1) || size(A,2)<size(T,2) 
    error(message('images:normxcorr2:invalidSizeForA'))
end

%-----------------------------------------------------------------------------
function checkIfFlat(T)

if std(T(:)) == 0
    error(message('images:normxcorr2:sameElementsInTemplate'))
end
