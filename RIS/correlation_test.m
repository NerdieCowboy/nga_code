clear all
close all
clc

cube_path = 'C:\data\Pacino\RIS\';
features_table_fn = 'C:\damon\dissertation\features_pigments_throughRenasaince.csv';
cube_fn = 'sub_sc_crop_cal_Pacino_ChristMajesty_vnir_xnir';
spectra_lib = 'resampled_EarlyRennaisance2';
spectra_names = {'lamp black', 'bone black', 'burnt umber', ...
    'red ochre', 'Yellow Ochre', 'gypsum', 'chalk', 'Madder lake', ...
    'Realgar', 'Malachite', 'Orpiment', 'Indigo', 'Azurite', ...
    'Red lead', 'vermilion', 'Green earth', 'verdigris', 'lead white', ...
    'ultramarine', 'Naples Yellow', 'Smalt', 'indian yellow', ...
    'Copper resinate', 'lead tin yellow', 'Van Dyke brown', ...
    'Carmine lake', 'cobalt blue', 'raw umber', 'Cd Red', ...
    'titanium white', 'Cd yellow', 'zinc white'};
derivative_fn = [cube_path 'slope_' cube_fn];
lib_derivative_fn = [cube_path 'slope_' spectra_lib];
wiggle = 5; %bands
blksz = [37 73 145];

% spectral library
fn = [cube_path spectra_lib];
lambdaon = 0;
lambda1 = [];
fnh = [fn '.hdr'];
fid = fopen(fnh);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
    if (lambdaon == 0)
        [~,~,e] = regexp(line,'^samples=','match','start','end');
        if (e>0)
            n01 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^lines=','match','start','end');
        if (e>0)
            m01 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^bands=','match','start','end');
        if (e>0)
            p01 = str2double(line((e+1):end));
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
                    lambda1 = [lambda1 str2double(tmp2(i))];
                end
                [~,~,e1] = regexp(tmp2(i),'}$','match','start','end');
                if (e1{1} > 0)
                    tmp3 = cell2mat(tmp2(i));
                    tmp3 = tmp3(1:(e1{1}-1));
                    lambda1 = [lambda1 str2double(tmp3)];
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
                    lambda1 = [lambda1 str2double(tmp2(i))];
                end
            end
        else
            tmp2 = regexp(line,'\,','split');
            for i = 1:length(tmp2)
                if (~isnan(str2double(tmp2(i))))
                    lambda1 = [lambda1 str2double(tmp2(i))];
                end
            end
        end
    end
    clear msk tmp1 tmp2 line tmp3
end
fclose(fid);
lambda1 = lambda1';
lib_der1 = multibandread(lib_derivative_fn,[m01,n01,p01],datatype,headeroffset,interleave,byteorder);

cc = zeros(m01,n01,m01,3);
for m1 = 1:m01
    s1 = lib_der1(m1,:);
    for m2 = m1:m01
        s2 = lib_der1(m2,:);
        for i = 1:length(blksz)
            for p = (ceil(blksz(i)/2)+wiggle):(n01-floor(blksz(i)/2)-wiggle)
                tmpcc = zeros(2*wiggle+1,1);
                inc = 1;
                for w = -wiggle:wiggle
                    tmp1 = s1((p-floor(blksz(i)/2)):(p+floor(blksz(i)/2)));
                    tmp1 = tmp1/std(tmp1);
                    tmp2 = s2((p-floor(blksz(i)/2)+w):(p+floor(blksz(i)/2)+w));
                    tmp2 = tmp2/std(tmp2);
                    tmpcc(inc) = sum(tmp1.*tmp2);
                    inc = inc + 1;
                end
                cc(m1,p,m2,i) = max(tmpcc);
            end
        end
    end
end
