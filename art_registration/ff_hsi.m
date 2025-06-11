function [] = ff_hsi(fn_path,trial,N,white1)
%fn_path = 'C:\data\blue_room_vnir\';
%trial = 'cube_1';
N1 = 1;
N2 = N;

% white
fnh = [fn_path white1 '.hdr'];
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
fn = [fn_path white1];
white = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);

% dark white
fnh = [fn_path 'dark_' white1 '.hdr'];
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
fn = [fn_path 'dark_' white1];
dark_white = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
den = (white - dark_white);
clear white dark_white

for j = N1:N2
    fnh = [fn_path trial '_' num2str(j) '.hdr'];
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

    fn = [fn_path trial '_' num2str(j)];
    fn1 = [fn_path trial '1_' num2str(j)];
    fn2 = [fn_path 'dark_' trial '_' num2str(j)];
    if (exist(fn,'file') == 2)
        dark = multibandread(fn2,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
        cube = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
        cube = (cube - dark);
        clear dark
        %cube = cube - min(cube(:));
        %den = den - max(den(:)) + 2^16 - 1;
        cube = cube./den;
        multibandwrite(single(cube),fn1,'bsq','machfmt','ieee-le')
        clear cube
        
        fn = [fn_path trial '1_' num2str(j) '.hdr'];
        fid = fopen(fn,'w');
        fprintf(fid,'ENVI\n');
        fprintf(fid,'description = {}\n');
        fprintf(fid,'samples = %u\n',n0);
        fprintf(fid,'lines   = %u\n',m0);
        fprintf(fid,'bands   = %u\n',p0);
        fprintf(fid,'header offset = 0\n');
        fprintf(fid,'file type = ENVI Standard\n');
        fprintf(fid,'data type = 4\n');
        fprintf(fid,'interleave = bsq\n');
        fprintf(fid,'sensor type = Unknown\n');
        fprintf(fid,'byte order = 0\n');
        fprintf(fid,'wavelength units = Unknown\n');
        fclose(fid);
    end
end
clear den

