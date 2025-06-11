%% Flat Field Automation
% flat_field.m
% Author: Damon Conover
% Email: dconover@gwmail.gwu.edu
% Latest Revision: 13 February 2012
%
% input:	white image (white.tif)
%			dark image (dark.tif)
%			last image number
%           data path or NULL
%           trial name of set (i.e., trial_)
% data:	set of image to be flat-fielded (trial_001.tif, trial_002.tif, ect.)
% output:	flat fielded images (one per input image, same name, but in a directory named out)
%
%   flat_field('white.tif','dark.tif',196,'c:\data\AIC\framelink_old\','framelink_')

function [] = flat_field(W,B,N,fn_path,trial)
if (nargin < 4 || nargin > 5)
    error('myApp:argChk', 'Wrong number of input arguments')
elseif (nargin == 4)
    trial = '';
end

if (strcmp(fn_path,'NULL') || strcmp(fn_path,'null'))
    fn_path = '';
end
fn_path2 = strcat(fn_path,'out/');

if (exist(fn_path2,'dir') == 7)
    rmdir(fn_path2,'s');
end
mkdir(fn_path2);

fn = strcat(fn_path,B);
black = double(imread(fn,'tif'));
fn = strcat(fn_path,W);
white = double(imread(fn,'tif'));

mx = 0;
mn = 2^16-1;
for i = 1:N
    if (N < 10)
        fn = [fn_path trial num2str(sprintf(['%01.0f'],i)) '.tif'];
    elseif (N < 100)
        fn = [fn_path trial num2str(sprintf(['%02.0f'],i)) '.tif'];
    elseif (N < 1000)
        fn = [fn_path trial num2str(sprintf(['%03.0f'],i)) '.tif'];
    elseif (N < 10000)
        fn = [fn_path trial num2str(sprintf(['%04.0f'],i)) '.tif'];
    end
    if (exist(fn,'file') == 2)
        img = double(imread(fn,'tif'));
        ff = (img - black)./(white - black);
        clear img
        mx0 = max(ff(:));
        if (mx0 > mx)
            mx = mx0;
        end
        mn0 = min(ff(:));
        if (mn0 < mn)
            mn = mn0;
        end
        clear ff
    end
end

for i = 1:N
    if (N < 10)
        fn = [fn_path trial num2str(sprintf(['%01.0f'],i)) '.tif'];
        fn2 = [fn_path2 trial num2str(sprintf(['%01.0f'],i)) '.tif'];
    elseif (N < 100)
        fn = [fn_path trial num2str(sprintf(['%02.0f'],i)) '.tif'];
        fn2 = [fn_path2 trial num2str(sprintf(['%02.0f'],i)) '.tif'];
    elseif (N < 1000)
        fn = [fn_path trial num2str(sprintf(['%03.0f'],i)) '.tif'];
        fn2 = [fn_path2 trial num2str(sprintf(['%03.0f'],i)) '.tif'];
    elseif (N < 10000)
        fn = [fn_path trial num2str(sprintf(['%04.0f'],i)) '.tif'];
        fn2 = [fn_path2 trial num2str(sprintf(['%04.0f'],i)) '.tif'];
    end
    if (exist(fn,'file') == 2)
        img = double(imread(fn,'tif'));
        ff = (img - black)./(white - black);
        clear img
        ff = (ff - mn)/(mx - mn);
        ff = (2^16-1)*ff;
        ff = uint16(ff);

        imwrite(ff,fn2,'tif','Compression','None');
        clear ff
    end
end
