% fld.m
% Matlab script for ECE385

% Author: Damon Conover (G19699044)
% Email: dconover@gwmail.gwu.edu
% Latest Revision: 16 November 2008

% Use Fisher's Linear Discriminant to project data onto a line
%   xx1 -> multidimensional positive data
%   xx2 -> multidimensional negative data
%   x1 -> 1-D positive data
%   X2 -> 1-D negative data

function [x1,x2] = fld(xx1,xx2)

%xx1 = randn(16,3);
%xx2 = 1+ randn(16,3);

% sample means
mu1 = mean(xx1);
mu2 = mean(xx2);

% scatter matrices
[m0 n0] = size(xx1);
S1 = zeros(n0,n0);
for m = 1:m0
    S1 = S1 + (xx1(m,:)-mu1)'*(xx1(m,:)-mu1);
end
[m0 n0] = size(xx2);
S2 = zeros(n0,n0);
for m = 1:m0
    S2 = S2 + (xx2(m,:)-mu2)'*(xx2(m,:)-mu2);
end

% within-class scatter matrix
Sw = S1 + S2;

% between-class scatter matrix
Sb = (mu1-mu2)'*(mu1-mu2);

% Fisher's linear discriminant
w = pinv(Sw)*(mu1-mu2)';

% projection matrix
P = w*w'/(w'*w);

% project points onto w
[m0 n0] = size(xx1);
xx1p = [];
for m = 1:m0
    b = P*xx1(m,:)';
    xx1p = [xx1p b];
end
xx1p = xx1p';
[m0 n0] = size(xx2);
xx2p = [];
for m = 1:m0
    b = P*xx2(m,:)';
    xx2p = [xx2p b];
end
xx2p = xx2p';

% project points onto y-axis
[m0 n0] = size(xx1p);
wx = zeros(n0,1);
wx(1) = 1;
Px = wx*wx'/(wx'*wx);
xx1x = [];
for m = 1:m0
    b = Px*xx1p(m,:)';
    xx1x = [xx1x b];
end
xx1x = xx1x';
[m0 n0] = size(xx2p);
xx2x = [];
for m = 1:m0
    b = Px*xx2p(m,:)';
    xx2x = [xx2x b];
end
xx2x = xx2x';

% output
x1 = xx1x(:,1);
x2 = xx2x(:,1);

%{
figure(1)
plot3(xx1(:,1),xx1(:,2),xx1(:,3),'rx')
hold on
plot3(xx2(:,1),xx2(:,2),xx2(:,3),'bx')
hold on
plot3(xx1p(:,1),xx1p(:,2),xx1p(:,3),'r.')
hold on
plot3(xx2p(:,1),xx2p(:,2),xx2p(:,3),'b.')

figure(2)
plot(xx1x(:,1),xx1x(:,2),'r.')
hold on
plot(xx2x(:,1),xx2x(:,2),'b.')
%}
