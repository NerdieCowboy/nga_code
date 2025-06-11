clear all
close all
clc

% Required cube name format: <trial>_row<row#>
fn_path = 'D:\damon\20140610_conover_xNIR\';
trial = 'test_ff_DMC_One_125msec';
n1 = 1024;    %columns (pixels)
approx_ov = 64;    %desired overlap (pixels)
M = 2;  %number of rows

fn = [fn_path trial '_row' num2str(1)];
fnh = [fn '.hdr'];
lambdaon = 0;
lambda = [];
fid = fopen(fnh);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
    if (lambdaon == 0)
        [~,~,e] = regexp(line,'^samples=','match','start','end');
        if (e>0)
            n0 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^lines=','match','start','end');
        if (e>0)
            m0 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^bands=','match','start','end');
        if (e>0)
            p0 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^datatype=','match','start','end');
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
        [~,~,e] = regexp(line,'^interleave=','match','start','end');
        if (e>0)
            interleave = line((e+1):end);
        end
        [~,~,e] = regexp(line,'^byteorder=','match','start','end');
        if (e>0)
            byteorder = str2double(line((e+1):end));
            if (byteorder == 0)
                byteorder = 'ieee-le';
            elseif (byteorder == 1)
                byteorder = 'ieee-be';
            end
        end
        [~,~,e] = regexp(line,'^headeroffset=','match','start','end');
        if (e>0)
            headeroffset = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^wavelength={','match','start','end');
        if (e>0)
            lambdaon = 1;
            tmp1 = line((e+1):end);
            tmp2 = regexp(tmp1,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda = [lambda str2double(tmp2(i))];
                end
            end
        end
    end
    if (lambdaon == 1)
        [~,~,e] = regexp(line,'}$','match','start','end');
        if (e>0)
            lambdaon = 0;
            tmp1 = line(1:(e-1));
            tmp2 = regexp(tmp1,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda = [lambda str2double(tmp2(i))];
                end
            end
        else
            tmp2 = regexp(line,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda = [lambda str2double(tmp2(i))];
                end
            end
        end
    end
end
fclose(fid);

N = (n0 - approx_ov)/(n1 - approx_ov);
N = ceil(N);
ov = (n0 - N*n1)/(N + 1);
ov = ceil(ov);

fw = 1;
cnt = 1;
for m = 1:M
    fn = [fn_path trial '_row' num2str(m)];
    fnh = [fn '.hdr'];
    for n = 1:(n1-ov):n0
        if (fw == 1)
            n11 = n;
            n12 = n + n1 - 1;
            if (n12 > n0)
                n12 = n0;
                n11 = n12 - n1 + 1;
            end
        elseif (fw == 0)
            n12 = n0 - n + 1;
            n11 = n12 - n1 + 1;
            if (n11 < 1)
                n11 = 1;
                n12 = n11 + n1 - 1;
            end
        end
        [m cnt n11 n12]
        fn_out = [fn_path trial '_' num2str(cnt)];
        if (strcmp(datatype,'uint8'))
            cube = uint8(multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Column','Range',[n11 n12]}));
            multibandwrite(uint8(cube),fn_out,'bsq')
        elseif (strcmp(datatype,'single'))
            cube = single(multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Column','Range',[n11 n12]}));
            multibandwrite(single(cube),fn_out,'bsq')
        elseif (strcmp(datatype,'double'))
            cube = double(multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Column','Range',[n11 n12]}));
            multibandwrite(double(cube),fn_out,'bsq')
        elseif (strcmp(datatype,'uint16'))
            cube = uint16(multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Column','Range',[n11 n12]}));
            multibandwrite(uint16(cube),fn_out,'bsq')
        end
        clear cube

        fnh_out = [fn_out '.hdr'];
        fid1 = fopen(fnh_out,'w');
        fid = fopen(fnh,'r');
        while ~feof(fid)
            test = 0;
            line0 = fgetl(fid);
            line = line0;
            msk = isspace(line);
            line(msk==1) = '';
            [~,~,e] = regexp(line,'^samples=','match','start','end');
            if (e>0)
                test = 1;
                fprintf(fid1,'samples = %u\n',n1);
            end
            [~,~,e] = regexp(line,'^lines=','match','start','end');
            if (e>0)
                test = 1;
                fprintf(fid1,'lines = %u\n',m0);
            end
            [~,~,e] = regexp(line,'^bands=','match','start','end');
            if (e>0)
                test = 1;
                fprintf(fid1,'bands = %u\n',p0);
            end
            [~,~,e] = regexp(line,'^datatype=','match','start','end');
            if (e>0)
                test = 1;
                if (strcmp(datatype,'uint8'))
                    fprintf(fid1,'data type = 1\n');
                    elseif (strcmp(datatype,'single'))
                        fprintf(fid1,'data type = 4\n');
                    elseif (strcmp(datatype,'double'))
                        fprintf(fid1,'data type = 5\n');
                    elseif (strcmp(datatype,'uint16'))
                        fprintf(fid1,'data type = 12\n');
                end
            end
            [~,~,e] = regexp(line,'^interleave=','match','start','end');
            if (e>0)
                test = 1;
                fprintf(fid1,'interleave = bsq\n');
            end
            [~,~,e] = regexp(line,'^byteorder=','match','start','end');
            if (e>0)
                test = 1;
                fprintf(fid1,'byte order = 0\n');
            end
            [~,~,e] = regexp(line,'^headeroffset=','match','start','end');
            if (e>0)
                test = 1;
                fprintf(fid1,'header offset = 0\n');
            end
            if (test == 0)
                fprintf(fid1,[line0 '\n']);
            end
        end
        fclose(fid);
        fclose(fid1);

        cnt = cnt + 1;
    end
    if (fw == 1)
        fw = 0;
    elseif (fw == 0)
        fw = 1;
    end
end
