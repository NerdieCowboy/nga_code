%function [] = rough_mosaicX(fn_path,trial,num1,num2,sc,pc,usegui)

warning('off','all')

trial = 'Mill';
fn_path = 'C:\damon\data\1942_9_62\';

N1 = 11;
N2 = 19;
sc = 0.5001;

mosaic_set = [];
m0 = [];
n0 = [];
for i = N1:N2
    fn = [fn_path trial '_' num2str(i) '.tif'];
    if (exist(fn,'file') == 2)
        mosaic_set = [mosaic_set i];
        img = imread(fn,'tif');
        img = imresize(img,sc,'bicubic');
        [m0(i-N1+1),n0(i-N1+1)] = size(img);
        fn_path2 = [fn_path '\sc_' trial];
        if (exist(fn_path2,'dir') ~= 7)
            mkdir(fn_path2);
        end
        fn = [fn_path 'sc_' trial '\sc_' trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
        imwrite(img,fn,'tif','Compression','None')
    end
end
sz = round(min(min(m0),min(n0))/5);
ov = 0;

corr_offset = zeros(length(mosaic_set),2);
first = 1;
cnt = 0;
for i = N1:N2
    fn_full2 = [fn_path 'sc_' trial '\sc_' trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
    if (exist(fn_full2,'file') == 2)
        cnt = cnt + 1;
        img = double(imread(fn_full2,'tif'));
        [m10 n10] = size(img);

        if (first == 0)
            corr_offset0 = [];
            mpnt = 1;
            mend = mpnt + sz - 1;
            npnt = 1;
            nend = npnt + sz - 1;
            inc = 1;
            statem = 0;
            staten = 0;
            forward = 1;
            while (statem < 2)
                while (staten < 2)
                    tmp = img(mpnt:mend,npnt:nend);
                    [m1,n1] = size(tmp);
                    cc = normxcorr2(tmp,img_prev);
                    [max_cc, imax] = max(cc(:));
                    [ypeak, xpeak] = ind2sub(size(cc),imax(1));
                    corr_offset0 = [corr_offset0;ypeak-m1-mpnt+1 xpeak-n1-npnt+1 max_cc];
                    clear tmp
                    
                    if (forward == 1)
                        npnt = nend - ov;
                        nend = npnt + sz - 1;
                        if (nend >= n10)
                            staten = staten + 1;
                            nend = n10;
                            npnt = n10 - sz + 1;
                        end
                    else
                        npnt = npnt - sz + ov;
                        nend = npnt + sz - 1;
                        if (npnt <= 1)
                            staten = staten + 1;
                            npnt = 1;
                            nend = sz;
                        end
                    end
                    inc = inc + 1;
                end
                if (forward == 1)
                    forward = 0;
                else
                    forward = 1;
                    npnt = 1;
                    nend = npnt + sz - 1;
                end
                staten = 0;
                mpnt = mend - ov;
                mend = mpnt + sz - 1;
                if (mend >= m10)
                    statem = statem + 1;
                    mend = m10;
                    mpnt = m10 - sz + 1;
                end
            end
            [mx ind] = max(corr_offset0(:,3));
            corr_offset(cnt,:) = corr_offset0(ind,1:2);
            %corr_offset(cnt,:) = round(median(corr_offset0(1:2)));
            %corr_offset(cnt,:) = mode(corr_offset0(1:2));
            
            [m1 n1] = size(img);
            [m2 n2] = size(img_prev);
            test = uint8(zeros(m2+2*(m1-1),n2+2*(n1-1),3));
            test(m1:(m1+m2-1),n1:(n1+n2-1),1) = img_prev;
            test((m1+corr_offset(cnt,1)):(m1+corr_offset(cnt,1)+m1-1),(n1+corr_offset(cnt,2)):(n1+corr_offset(cnt,2)+n1-1),3) = uint8(img);
            imshow(test)
            pause(1)
            clear test
        else
            first = 0;
        end
        img_prev = img;
        [m2 n2] = size(img_prev);
        clear img
        pause(0.1)
    end
end
clear img_prev

[mc nc] = size(corr_offset);
tmp = zeros(mc,4);
tmp(:,1) = mosaic_set';
tmp(:,2:3) = corr_offset;
tmp(:,4) = 1;
fn = [fn_path 'sc_' trial '\offset_values.csv'];
csvwrite(fn,tmp)
clear tmp

%{
    %{
    if (pc == 1)
        fn_full2 = [fn_path '\offset_values.csv'];
    elseif (pc == 0)
        fn_full2 = [fn_path '/offset_values.csv'];
    end
    csvwrite(fn_full2,rel_reg)
    %}







for i = num1:num2
    if (pc == 1)
        fn = [fn_path '\' trial '_' num2str(i) '.tif'];
    elseif (pc == 0)
        fn = [fn_path '/' trial '_' num2str(i) '.tif'];
    end
    if (exist(fn,'file') == 2)
        img = imread(fn,'tif');
        img = imresize(img,sc,'bicubic');
        if (pc == 1)
            fn = [fn_path '\sc_' trial '_' num2str(i) '.tif'];
        elseif (pc == 0)
            fn = [fn_path '/sc_' trial '_' num2str(i) '.tif'];
        end
        delete(fn)
    end
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
if (usegui == 3)
    handles = guidata(XrayRegister);
    set(handles.running,'String','');
    drawnow;
    pause(0.1)
end
%}