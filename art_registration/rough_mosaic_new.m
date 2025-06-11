trial = 'sc_ManSeated';
fn_path = '/Users/Shared/Syscat Dutch Regsitration/1937_1_77/sc_ManSeated/';
pc = 0;
sz = 1536;
ov = round(sz/16);

for i = 1:3
    last_pt = [0 0];
    offset_values = zeros(12,4);
    if (pc == 1)
        fn_path2 = [fn_path trial '_' num2str(i) '\'];
    elseif (pc == 0)
        fn_path2 = [fn_path trial '_' num2str(i) '/'];
    end
    if (exist(fn_path2,'dir') ~= 7)
        mkdir(fn_path2);
    end

    if (pc == 1)
        fn = [fn_path trial '_' num2str(i) '.tif'];
    elseif (pc == 0)
        fn = [fn_path trial '_' num2str(i) '.tif'];
    end
    if (exist(fn,'file') == 2)
        info = imfinfo(fn);
        bd1 = info.BitsPerSample;
        bd1 = bd1(1);

        img = imread(fn);
        [m1,n1] = size(img);

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
                    fn_full0 = [fn_path trial '_' num2str(i) '\' trial '_' num2str(i) '_' num2str(sprintf('%03.0f',inc)) '.tif'];
                elseif (pc == 0)
                    fn_full0 = [fn_path trial '_' num2str(i) '/' trial '_' num2str(i) '_' num2str(sprintf('%03.0f',inc)) '.tif'];
                end

                tmp = img(mpnt:mend,npnt:nend);
                offset_values(inc,:) = [inc mpnt-last_pt(1) npnt-last_pt(2) 1];
                last_pt = [mpnt npnt];
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
    end
    if (pc == 1)
        fn_csv = [fn_path trial '_' num2str(i) '\offset_values.csv'];
        csvwrite(fn_csv,offset_values)
    elseif (pc == 0)
        fn_csv = [fn_path trial '_' num2str(i) '/offset_values.csv'];
        csvwrite(fn_csv,offset_values)
    end
    
end
