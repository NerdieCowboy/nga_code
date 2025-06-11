
fn_path = 'C:\data\conover_nir_vnir_xrf\';
trial = 'conover_xrf';
m1 = 34;
m2 = 312;
a0 = -0.0113;
a1 = 0.0137;
p1 = 1;
p2 = 1024;
lambda = a1*(p1:p2)+a0;

fnh = [fn_path trial '.hdr'];
fid = fopen(fnh);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
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
        datatype0 = str2double(line((e+1):end));
        if (datatype0 == 1)
            datatype = 'uint8';
        elseif (datatype0 == 4)
            datatype = 'single';
        elseif (datatype0 == 5)
            datatype = 'double';
        elseif (datatype0 == 12)
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
end
fclose(fid);

fn = [fn_path trial];
for m = m1:m2
    cube = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Row','Direct',m});
    sumspectra = sum(cube,2);
    sumspectra = sumspectra(:);
    fn_csv = [fn_path trial '_row' num2str(sprintf('%03.0f',m)) '.csv'];
    csvwrite(fn_csv,[lambda' sumspectra]);
end
