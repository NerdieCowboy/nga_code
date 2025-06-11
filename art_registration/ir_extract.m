%% Extract IR subimages
function [] = ir_extract(ir_fn,tilt,usegui,sc)

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

if (usegui == 0)
    if (nargin < 4)
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
        if (nargin < 4)
            sc = input('What is the scale factor? [1]: ', 's');
            if isempty(sc)
                sc = '1';
            end
            sc = str2num(sc);
        end
        disp('Template images are being extracted ...');
    end
    
    if (pc == 1)
        fn_path2 = [fn_path trial '\'];
    elseif (pc == 0)
        fn_path2 = [fn_path trial '/'];
    end
    if (exist(fn_path2,'dir') ~= 7)
        mkdir(fn_path2);
    end

    %IR
    fn_full = [fn_path trial '.img'];
    if (exist(fn_full,'file') == 2)
        fid = fopen(fn_full,'r');
        header = fread(fid,2);  %header
        m0 = fread(fid,[1,1],'uint16',0,'ieee-le');
        n0 = fread(fid,[1,1],'uint16',0,'ieee-le');
        header = fread(fid,506);  %header
        %missing 512 Bytes
        clear header

        p1 = 0;
        while ~feof(fid)
            tmp = fread(fid,[m0,n0],'uint16',0,'ieee-le');
            if (~isempty(tmp))
                p1 = p1 + 1;
                if (pc == 1)
                    fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03.0f',p1)) '.tif'];
                elseif (pc == 0)
                    fn_full2 = [fn_path2 trial '_' num2str(sprintf('%03d',p1)) '.tif'];
                end
                if (tilt == 0)
                    tmp = tmp(end:-1:1,:);
                elseif (tilt == 1)
                    tmp = tmp(:,end:-1:1);
                elseif (tilt == 2)
                    tmp = tmp(end:-1:1,end:-1:1);
                elseif (tilt == 3)
                    tmp = tmp';
                end
                tmp = imresize(tmp,sc,'bicubic');
                imwrite(uint16(tmp),fn_full2,'tif','Compression','None');
                clear tmp
                pause(0.1)
            end
        end
        fclose(fid);
    end
    
    if (usegui == 0)
        %input('Delete all unwanted template images. Press return to continue.');
    elseif (usegui == 1)
        %message = sprintf('Delete all unwanted template images.\n\nClick OK to continue');
        %uiwait(msgbox(message));
        %rough_mosaic(ir_fn,1)
    end
end
end
