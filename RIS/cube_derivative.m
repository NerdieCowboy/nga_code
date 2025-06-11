function [] = cube_derivative(fn_path,trial,hsz,firstd)

%{
fn_path = 'C:\data\Pacino\RIS\';
trial = 'slope_slope_spatialavg_sub_sc_crop_cal_Pacino_ChristMajesty_vnir_xnir';
hsz = 3;
firstd = 0;
cube_derivative(fn_path,trial,hsz,firstd)
%}

% data cube
fn = [fn_path trial];
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
                [~,~,e1] = regexp(tmp2(i),'}$','match','start','end');
                if (e1{1} > 0)
                    tmp3 = cell2mat(tmp2(i));
                    tmp3 = tmp3(1:(e1{1}-1));
                    lambda = [lambda str2double(tmp3)];
                end
            end
            [~,~,e] = regexp(line,'}$','match','start','end');
            if (e>0)
                lambdaon = 0;
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
    clear msk tmp1 tmp2 line tmp3
end
fclose(fid);
lambda = lambda';

% data cube derivative
fn2 = [fn_path 'slope_' trial];
fn2h = [fn2 '.hdr'];
data3 = single(zeros(m0,n0,p0));
for p = (hsz+ceil(hsz/2)):(p0-hsz-floor(hsz/2))
    tmp = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Range',[p-hsz-floor(hsz/2) p+hsz+floor(hsz/2)]});
    if (firstd == 1)
        lt = mean(tmp(:,:,1:hsz),3);
        md = mean(tmp(:,:,(hsz+1):(2*hsz)),3);
        rt = mean(tmp(:,:,(2*hsz+1):(3*hsz)),3);
    else
        lt = tmp(:,:,ceil(hsz/2));
        md = tmp(:,:,(hsz + ceil(hsz/2)));
        rt = tmp(:,:,(2*hsz + ceil(hsz/2)));
    end
    data3(:,:,p) = single(rt - lt);
    clear lt md rt
    
    %tmpA = log(tmp.^(-1));
    %tmpA(isinf(tmpA)) = 1;
    %tmpA(isinf(-tmpA)) = -1;
    %mxtmpA = max(tmpA(:));
    %tmpA = tmpA/mxtmpA;
    %lt = mean(tmpA(:,:,1:hsz),3);
    %md = mean(tmpA(:,:,(hsz+1):(2*hsz)),3);
    %rt = mean(tmpA(:,:,(2*hsz+1):(3*hsz)),3);
    %dataA(:,:,p) = single(rt - lt);
    %clear tmp tmpA lt md rt
end
multibandwrite(single(data3(:,:,1)),fn2,'bsq',[1,1,1],[m0,n0,p0]);
%multibandwrite(single(dataA(:,:,1)),fn2A,'bsq',[1,1,1],[m0,n0,p0]);
for p = 1:p0
    multibandwrite(single(data3(:,:,p)),fn2,'bsq',[1,1,p],[m0,n0,p0]);
    %multibandwrite(single(dataA(:,:,p)),fn2A,'bsq',[1,1,p],[m0,n0,p0]);
end
%{
fid1 = fopen(fn2Ah,'w');
fprintf(fid1,'ENVI\n');
fprintf(fid1,'description = {}\n');
fprintf(fid1,'samples = %u\n',n0);
fprintf(fid1,'lines = %u\n',m0);
fprintf(fid1,'bands = %u\n',p0);
fprintf(fid1,'header offset = 0\n');
fprintf(fid1,'file type = ENVI Standard\n');
fprintf(fid1,'data type = 4\n');
fprintf(fid1,'interleave = bsq\n');
fprintf(fid1,'byte order = 0\n');
fprintf(fid1,'wavelength = {');
for i = 1:(p0-1)
    fprintf(fid1,'%f, ',lambda(i));
end
fprintf(fid1,'%f}\n',lambda(p0));
fclose(fid1);
%}
fnh_in = [fnh];
fnh_out = [fn2h];
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
        fprintf(fid1,'data type = 4\n');
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
