function [] = create_dark(fn_path,trial,N)
%fn_path = 'C:\data\blue_room_vnir\';
%trial = 'cube';
N1 = 1;
N2 = N;
fnh = [fn_path trial '_' num2str(N) '.hdr'];
fid = fopen(fnh);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
    [~,~,e] = regexp(line,'samples=','match','start','end');
    if (e>0)
        n0 = str2double(line((e+1):end));
    end
    [~,~,e] = regexp(line,'lines=','match','start','end');
    if (e>0)
        m0 = str2double(line((e+1):end));
    end
    [~,~,e] = regexp(line,'bands=','match','start','end');
    if (e>0)
        p0 = str2double(line((e+1):end));
    end
    [~,~,e] = regexp(line,'datatype=','match','start','end');
    if (e>0)
        datatype = str2double(line((e+1):end));
        if (datatype == 1)
            datatype = 'uint8';
            elseif (datatype == 4)
                datatype = 'single';
            elseif (datatype == 5)
                datatype = 'double';
            elseif (datatype == 12)
                datatype = 'uint16';
        end
    end
    [~,~,e] = regexp(line,'interleave=','match','start','end');
    if (e>0)
        interleave = line((e+1):end);
    end
    [~,~,e] = regexp(line,'byteorder=','match','start','end');
    if (e>0)
        byteorder = str2double(line((e+1):end));
        if (byteorder == 0)
            byteorder = 'ieee-le';
        elseif (byteorder == 1)
            byteorder = 'ieee-be';
        end
    end
    [~,~,e] = regexp(line,'headeroffset=','match','start','end');
    if (e>0)
        headeroffset = str2double(line((e+1):end));
    end
end
fclose(fid);

for j = N1:N2
    fn = [fn_path trial '_' num2str(j)];
    fn1 = [fn_path 'dark_' trial '_' num2str(j)];
    if (exist(fn,'file') == 2)
        cube = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Range',[1 4]});
        dark = mean(cube,3);
        clear cube
        
        Npoly = 1;
        [X,Y] = meshgrid(1:n0,1:m0);
        A = zeros(m0*n0,sum((1:(Npoly+1))));
        pnt = 1;
        for p = 0:Npoly
            for q = 0:(Npoly-p)
                A(:,pnt) = (X(:).^p).*(Y(:).^q);
                pnt = pnt + 1;
            end
        end
        pA = pinv(A);
        clear A
        coeff = pA*dark(:);
        clear pA
        dark0 = zeros(m0*n0,1);
        pnt = 1;
        for p = 0:Npoly
            for q = 0:(Npoly-p)
                dark0 = dark0 + coeff(pnt)*((X(:).^p).*(Y(:).^q));
                pnt = pnt + 1;
            end
        end
        %{
        plot3(X(:),Y(:),dark(:),'r.')
        hold on
        plot3(X(:),Y(:),dark0,'b.')
        hold off
        pause
        %}
        clear X Y dark
        
        dark0 = reshape(dark0,m0,n0);
        dark1 = repmat(dark0,[1 1 p0]);
        clear dark0
        if (strcmp(datatype,'uint8'))
            multibandwrite(uint8(dark1),fn1,interleave,[1,1,1],[m0,n0,p0],'offset',headeroffset,'machfmt',byteorder)
            elseif (strcmp(datatype,'single'))
                multibandwrite(single(dark1),fn1,interleave,[1,1,1],[m0,n0,p0],'offset',headeroffset,'machfmt',byteorder)
            elseif (strcmp(datatype,'double'))
                multibandwrite(double(dark1),fn1,interleave,[1,1,1],[m0,n0,p0],'offset',headeroffset,'machfmt',byteorder)
            elseif (strcmp(datatype,'uint16'))
                multibandwrite(uint16(dark1),fn1,interleave,[1,1,1],[m0,n0,p0],'offset',headeroffset,'machfmt',byteorder)
        end
        clear dark1
        copyfile([fn_path trial '_' num2str(j) '.hdr'],[fn_path 'dark_' trial '_' num2str(j) '.hdr'])
    end
end
    