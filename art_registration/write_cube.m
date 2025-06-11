function [] = write_cube(cube,fn_in,fn_out,datatype)

[m0,n0,p0] = size(cube);

fnh_in = [fn_in '.hdr'];
fnh_out = [fn_out '.hdr'];
fid1 = fopen(fnh_out,'w');
fid = fopen(fnh_in,'r');
while ~feof(fid)
    test = 0;
    line0 = fgetl(fid);
    line = line0;
    msk = isspace(line);
    line(msk==1) = '';
    [~,~,e] = regexp(line,'^samples=','match','start','end');
    if (e>0)
        test = 1;
        fprintf(fid1,'samples = %u\n',n0);
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

if (strcmp(datatype,'uint8'))
    multibandwrite(uint8(cube),fn_out,'bsq','machfmt','ieee-le')
    elseif (strcmp(datatype,'single'))
        multibandwrite(single(cube),fn_out,'bsq','machfmt','ieee-le')
    elseif (strcmp(datatype,'double'))
        multibandwrite(double(cube),fn_out,'bsq','machfmt','ieee-le')
    elseif (strcmp(datatype,'uint16'))
        multibandwrite(uint16(cube),fn_out,'bsq','machfmt','ieee-le')
end
