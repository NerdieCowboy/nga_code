close all
clc

fn_path = 'C:\data\Pacino\RIS\';
cube_fn = 'sub_sc_crop_cal_Pacino_ChristMajesty_vnir_xnir';
lib_test = 0;
lambdaon = 0;
lambda = [];
fnh = [cube_fn '.hdr'];
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
cube = multibandread([fn_path cube_fn],[m0,n0,p0],datatype,headeroffset,interleave,byteorder);

m = 105;
n = 153;
tmp = cube(m,n,:);
tmp = tmp(:);
figure,plot(lambda,tmp)
ind_test = ((features(:,1) == m) & (features(:,2) == n) & (features(:,3) == 13));
features(ind_test,:)

m = 114;
n = 16;
tmp = cube(m,n,:);
tmp = tmp(:);
figure,plot(lambda,tmp)
ind_test = ((features(:,1) == m) & (features(:,2) == n) & (features(:,3) == 13));
features(ind_test,:)

m = 112;
n = 156;
tmp = cube(m,n,:);
tmp = tmp(:);
figure,plot(lambda,tmp)
ind_test = ((features(:,1) == m) & (features(:,2) == n) & (features(:,3) == 13));
features(ind_test,:)
ind_test = ((features(:,1) == m) & (features(:,2) == n) & (features(:,3) == 23));
features(ind_test,:)
ind_test = ((features(:,1) == m) & (features(:,2) == n) & (features(:,3) == 16));
features(ind_test,:)

m = 48;
n = 16;
tmp = cube(m,n,:);
tmp = tmp(:);
figure,plot(lambda,tmp)
ind_test = ((features(:,1) == m) & (features(:,2) == n) & (features(:,3) == 12));
features(ind_test,:)