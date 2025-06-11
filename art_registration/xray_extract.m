%% Extract IR subimages
function [] = xray_extract(fn_path,trial0,j,img_ext,sc,sz,pc)

ov = round(sz/16);
if (~isempty(j))
    trial = [trial0 '_' num2str(j)];
else
    trial = trial0;
end
if (pc == 1)
    fn = [fn_path '\'  trial '.' img_ext];
    ir_fn = [fn_path '\' fn];
elseif (pc == 0)
    fn = [fn_path '/'  trial '.' img_ext];
    ir_fn = [fn_path '/' fn];
end
if (exist(ir_fn,'file') == 2)
info = imfinfo(ir_fn);
bd1 = info.BitsPerSample;
bd1 = bd1(1);

if (pc == 1)
   fn_path2 = [fn_path trial '\'];
elseif (pc == 0)
   fn_path2 = [fn_path trial '/'];
end
if (exist(fn_path2,'dir') ~= 7)
    mkdir(fn_path2);
end

img = imread(ir_fn);
img = imresize(img,sc,'bicubic');
[m1,n1,~] = size(img);

mpnt = 1;
mend = mpnt + sz - 1;
if (mend > m1)
    mend = m1;
    statem = 1;
else
    statem = 0;
end
npnt = 1;
nend = npnt + sz - 1;
if (nend > n1)
    nend = n1;
    staten = 1;
else
    staten = 0;
end
inc = 1;
forward = 1;
while (statem < 2)
    while (staten < 2)
        if (pc == 1)
            fn_full0 = [fn_path2 trial '_' num2str(sprintf('%03.0f',inc)) '.tif'];
        elseif (pc == 0)
            fn_full0 = [fn_path2 trial '_' num2str(sprintf('%03d',inc)) '.tif'];
        end

        tmp = img(mpnt:mend,npnt:nend,:);
        if (inc == 1)
            offset_values = [1 0 0 1];
        else
            offset_values = [offset_values;inc (mpnt-last(1)) (npnt-last(2)) 1];
        end
        last = [mpnt npnt];
        if (bd1 == 16)
            imwrite(uint16(tmp),fn_full0,'tif','Compression','None');
        elseif (bd1 == 8)
            imwrite(uint8(tmp),fn_full0,'tif','Compression','None');
        end
        clear tmp

        if (forward == 1)
            npnt = nend - ov;
            nend = npnt + sz - 1;
            if (nend >= n1)
                staten = staten + 1;
                nend = n1;
                npnt = n1 - sz + 1;
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
    if (mend >= m1)
        statem = statem + 1;
        mend = m1;
        mpnt = m1 - sz + 1;
    end
end
clear img

fn_full0 = [fn_path2 'offset_values.csv'];
csvwrite(fn_full0,offset_values);
clear offset_values
end
