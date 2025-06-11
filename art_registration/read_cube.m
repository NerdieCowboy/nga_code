function [cube,lambda,m0,n0,p0] = read_cube(fn,bsq,m10,m20,n10,n20,p10,p20)

lambdaon = 0;
lambda = [];
fnh = [fn '.hdr'];
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

if (strcmp(bsq,'bsq'))
    fn = [fn '.bsq'];
elseif (strcmp(bsq,'bil'))
    fn = [fn '.bil'];
elseif (strcmp(bsq,'bip'))
    fn = [fn '.bip'];
elseif (strcmp(bsq,'cube'))
    fn = [fn '.cube'];
end

if (isempty(p10))
    p10 = 1;
end
if (isempty(p20))
    p20 = p0;
end
if (isempty(m10))
    m10 = 1;
end
if (isempty(m20))
    m20 = m0;
end
if (isempty(n10))
    n10 = 1;
end
if (isempty(n20))
    n20 = n0;
end
if (strcmp(datatype,'uint8'))
    cube = uint8(multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Row','Range',[m10 m20]},{'Column','Range',[n10 n20]},{'Band','Range',[p10 p20]}));
elseif (strcmp(datatype,'single'))
    cube = single(multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Row','Range',[m10 m20]},{'Column','Range',[n10 n20]},{'Band','Range',[p10 p20]}));
elseif (strcmp(datatype,'double'))
    cube = double(multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Row','Range',[m10 m20]},{'Column','Range',[n10 n20]},{'Band','Range',[p10 p20]}));
elseif (strcmp(datatype,'uint16'))
    cube = uint16(multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Row','Range',[m10 m20]},{'Column','Range',[n10 n20]},{'Band','Range',[p10 p20]}));
end
