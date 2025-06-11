% mi.m
% Author: Damon Conover
% Email: dconover@gwmail.gwu.edu
% Latest Revision: 18 July 2012
%
% Calculate mutual information between two images
%   img1 -> matrix 1 (0 -> 1)
%   img2 -> matrix 2 (0 -> 1)
%   bins -> number of bins
%   mi0 -> mutual information value
%
% Note: img1 and img2 must be coverted to double-precision values and
% normalized to be between 0 and 1

function [mi0] = mi(img1,img2,bins)

[m1 n1] = size(img1);

img1 = img1(:);
img2 = img2(:);
p1 = zeros(bins,1);
p2 = zeros(bins,1);
A1 = zeros(m1*n1,bins);
A2 = zeros(m1*n1,bins);

%for j = 0:(bins-1)
parfor j = 0:(bins-1)
    A1(:,j+1) = (img1 >= j/bins) & (img1 < (j+1)/bins);
    p1(j+1) = sum(A1(:,j+1));     %histogram 1
    A2(:,j+1) = (img2 >= j/bins) & (img2 < (j+1)/bins);
    p2(j+1) = sum(A2(:,j+1));     %histogram 2
end

pjoint = A1'*A2;
clear img1 img2 A1 A2

p1 = p1/sum(p1);   %probability 1
p1(p1==0) = [];
p2 = p2/sum(p2);   %probability 2
p2(p2==0) = [];
pjoint = pjoint/sum(pjoint(:)); %joint probability
pjoint = pjoint(:);
pjoint(pjoint==0) = [];

%calculate entropy of img1
H1 = -p1'*log2(p1);
clear p1

%calculate entropy of img2
H2 = -p2'*log2(p2);
clear p2

%calculate joint entropy of img1 and img2
H12 = -pjoint'*log2(pjoint);
clear pjoint

%calculate mutual information
mi0 = H1 + H2 - H12;
clear H1 H2 H12
