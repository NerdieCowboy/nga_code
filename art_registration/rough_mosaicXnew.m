%function [] = rough_mosaicX(fn_path,trial,num1,num2,sc,pc,usegui)

fn_path = 'C:\data\1942_9_59';
trial = '1942_9_59';
num1 = 1;
num2 = 12;
sc = 0.30656;
pc = 1;
%usegui = 3;
usegui = 0;
sc1 = 1;

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
    warning('off','all')

    mosaic_set = [];
    m0 = [];
    n0 = [];
    for i = num1:num2
        if (pc == 1)
            fn = [fn_path '\' trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
        elseif (pc == 0)
            fn = [fn_path '/' trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
        end
        if (exist(fn,'file') == 2)
            mosaic_set = [mosaic_set i];
            img = imread(fn,'tif');
            img = imresize(img,sc,'bicubic');
            [m0(i),n0(i)] = size(img);
            if (pc == 1)
                fn = [fn_path '\sc_' trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
            elseif (pc == 0)
                fn = [fn_path '/sc_' trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
            end
            imwrite(img,fn,'tif','Compression','None')
        end
    end

    corr_offset = zeros(length(mosaic_set),4);
    first = 1;
    cnt = 0;
    for i = num1:num2
        if (pc == 1)
            fn_full2 = [fn_path '\sc_' trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
        elseif (pc == 0)
            fn_full2 = [fn_path '/sc_' trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
        end

        if (exist(fn_full2,'file') == 2)
            cnt = cnt + 1;
            temp = double(imread(fn_full2,'tif'));
            temp = imresize(temp,sc1,'bicubic');
            [mtemp,ntemp] = size(temp);
            mu_temp = mean(temp(:));
            sigma_temp = std(temp(:));
            temp = temp - mu_temp;

            if (first == 0)
                tmp1 = zeros(mref+2*(mtemp-1),nref+2*(ntemp-1));
                %tmp1 = sigma_temp*randn(mref+2*(mtemp-1),nref+2*(ntemp-1));
                tmp1(mtemp:((mtemp-1)+mref),ntemp:((ntemp-1)+nref)) = ref;
                
                tmp2 = zeros(mref+2*(mtemp-1),nref+2*(ntemp-1));
                %tmp2 = sigma_temp*randn(mref+2*(mtemp-1),nref+2*(ntemp-1));
                tmp2(1:mtemp,1:ntemp) = temp;
                
                %CC = conj(fft2(tmp1)).*fft2(tmp2);
                F1 = fft2(tmp1);
                F1 = angle(F1);
                F2 = conj(fft2(tmp2));
                F2 = angle(F2);
                CC = F1 + F2;
                clear F1 F2
                %CC = F1.*F2;
                CC = exp(1i*CC);
                %F1 = exp(1i*F1);
                %F2 = exp(1i*conj(F2));
                cc = real(ifft2(CC));
                %f1 = real(ifft2(F1));
                %f2 = real(ifft2(F2));
                clear CC
                
                [max_cc,imax] = max(cc,[],1);
                [~,xpeak] = max(max_cc);
                ypeak = imax(xpeak(1));
                corr_offset(cnt,1) = cnt;
                corr_offset(cnt,2:3) = [ypeak - mtemp, xpeak - ntemp];
                %ypeak - mtemp
                %xpeak - ntemp
                %{
                figure(100)
                cc(cc<=0) = 0;
                cc = log(cc+1);
                subplot(221),imagesc(tmp1)
                axis image
                subplot(222),imagesc(tmp2)
                axis image
                subplot(223),imagesc(cc)
                axis image
                %}
                test = uint8(zeros(mref+2*(mtemp-1),nref+2*(ntemp-1),3));
                test(mtemp:((mtemp-1)+mref),ntemp:((ntemp-1)+nref),1) = uint8(ref);
                test((mtemp+corr_offset(cnt,2)):(2*mtemp+corr_offset(cnt,2)-1),(ntemp+corr_offset(cnt,3)):(2*ntemp+corr_offset(cnt,3)-1),3) = uint8(temp);
                if (usegui == 0)
                    imshow(uint8(test))
                elseif (usegui == 3)
                    handles = guidata(XrayRegister);
                    axes(handles.axes1)
                    imshow(uint8(test))
                end
                clear test
            else
                first = 0;
                corr_offset(cnt,1) = cnt;
            end
            ref = temp;
            mu_ref = mu_temp;
            sigma_ref = sigma_temp;
            mref = mtemp;
            nref = ntemp;
            clear temp
            pause(0.1)
        end
    end
    clear img2_1

    if (pc == 1)
        fn_full2 = [fn_path '\offset_values.csv'];
    elseif (pc == 0)
        fn_full2 = [fn_path '/offset_values.csv'];
    end
    corr_offset(:,2:3) = round(corr_offset(:,2:3)/sc1);
    csvwrite(fn_full2,corr_offset)
end

for i = num1:num2
    if (pc == 1)
        fn = [fn_path '\' trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
    elseif (pc == 0)
        fn = [fn_path '/' trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
    end
    if (exist(fn,'file') == 2)
        img = imread(fn,'tif');
        img = imresize(img,sc,'bicubic');
        if (pc == 1)
            fn = [fn_path '\sc_' trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
        elseif (pc == 0)
            fn = [fn_path '/sc_' trial '_' num2str(sprintf('%03.0f',i)) '.tif'];
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

%end

