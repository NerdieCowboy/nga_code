%% Extract IR subimages
function [] = ir_extractX(ir_fn,usegui,sc,sz,ov)

pc = ~ismac;
if (pc == 1)
    tmp = regexp(ir_fn,'\\');
elseif (pc == 0)
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

info = imfinfo(ir_fn);
bd1 = info.BitsPerSample;
bd1 = bd1(1);

if (usegui == 0)
    if (nargin < 5)
        reply = input('Do you wish to extract template images? Y/N [N]: ', 's');
        if isempty(reply)
            reply = 'N';
        end
    else
        reply = 'Y';
    end
elseif (usegui == 1)
    reply = 'Y';
end    
if (strcmpi(reply, 'Y') == 1)
    if (usegui == 0)
        if (nargin < 3)
            sc = input('What is the scale factor? [1]: ', 's');
            if isempty(sc)
                sc = '1';
            end
            sc = str2double(sc);
        end

        if (nargin < 4)
            sz = input('Set the block size in pixels [1024]: ', 's');
            if isempty(sz)
                sz = '1024';
            end
            sz = str2double(sz);
        end

        if (nargin < 5)
            ov = input('Set the overlap amount in pixels [128]: ', 's');
            if isempty(ov)
                ov = '128';
            end
            ov = str2double(ov);
        end
        disp('Template images are being extracted ...');
    elseif (usegui == 1)
        prompt = {'Set the block size in pixels: ','Set the overlap amount in pixels: '};
        dlg_title = 'Image segmentation settings';
        num_lines = 1;
        def = {'1024','128'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        sz = str2double(answer{1});
        ov = str2double(answer{2});
    end

    if (pc == 1)
        fn_path2 = [fn_path trial '\'];
    elseif (pc == 0)
        fn_path2 = [fn_path trial '/'];
    end
    if (exist(fn_path2,'dir') ~= 7)
        mkdir(fn_path2);
    end

    img = imread(ir_fn);
    [m1,n1,~] = size(img);

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
            if (pc == 1)
                fn_full0 = [fn_path2 trial '_' num2str(sprintf('%03.0f',inc)) '.tif'];
            elseif (pc == 0)
                fn_full0 = [fn_path2 trial '_' num2str(sprintf('%03d',inc)) '.tif'];
            end

            tmp = img(mpnt:mend,npnt:nend,:);
            tmp = imresize(tmp,sc);
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
    
    if (usegui == 0)
        %input('Delete all unwanted template images. Press return to continue.');
    elseif (usegui == 1)
        message = sprintf('Delete all unwanted template images.\n\nClick OK to continue');
        uiwait(msgbox(message));
        rough_mosaic(ir_fn,1)
    end
end
end
