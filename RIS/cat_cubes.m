clear all
close all
clc

fn_path = 'C:\data\Pacino\RIS\';
vnir_trial = 'cal_rot_crop_ff_vnir_pacino_cube_full';
nir_trial = 'sc_cal_xnir_pacino_cube_full';
trial = 'cal_Pacino_ChristMajesty_vnir_xnir';
N = 3;
p0_vnir_last = 255; %band
p0_nir_first = 62;

lambdaon = 0;
lambda = [];
fnh = [fn_path vnir_trial '.hdr'];
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
            p0_vnir = str2double(line((e+1):end));
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
lambda_vnir = lambda;

lambdaon = 0;
lambda = [];
fnh = [fn_path nir_trial '.hdr'];
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
            p0_nir = str2double(line((e+1):end));
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
lambda_nir = lambda;
clear lambda

fnh = [fn_path nir_trial '.hdr'];
fnh_out = [fn_path trial '1.hdr'];
fid1 = fopen(fnh_out,'w');
fid = fopen(fnh,'r');
bandnames_test = 0;
wavelength_test = 0;
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
        fprintf(fid1,'bands = %u\n',(p0_vnir_last + p0_nir - p0_nir_first + 1));
    end
    [~,~,e] = regexp(line,'^wavelength={','match','start','end');
    if (e>0)
        test = 1;
        wavelength_test = 1;
        fprintf(fid1,'wavelength = {%f, ',lambda_vnir(1));
        for i = 2:p0_vnir_last
            fprintf(fid1,'%f, ',lambda_vnir(i));
        end
        for i = p0_nir_first:(p0_nir-1)
            fprintf(fid1,'%f, ',lambda_nir(i));
        end
        fprintf(fid1,'%f}\n',lambda_nir(p0_nir));
    end
    [~,~,e] = regexp(line,'^bandnames={','match','start','end');
    if (e>0)
        test = 1;
        bandnames_test = 1;
    end
    [~,~,e] = regexp(line,'^interleave=','match','start','end');
    if (e>0)
        test = 1;
        fprintf(fid1,'interleave = bsq\n');
    end
    [~,~,e] = regexp(line,'^headeroffset=','match','start','end');
    if (e>0)
        test = 1;
        fprintf(fid1,'header offset = 0\n');
    end
    if ((test == 0) && (bandnames_test == 0) && (wavelength_test == 0))
        fprintf(fid1,[line0 '\n']);
    end
    if (wavelength_test == 1)
        [~,~,e] = regexp(line,'}$','match','start','end');
        if (e>0)
            wavelength_test = 0;
        end
    end
    if (bandnames_test == 1)
        [~,~,e] = regexp(line,'}$','match','start','end');
        if (e>0)
            bandnames_test = 0;
        end
    end
end
fclose(fid);
fclose(fid1);
copyfile(fnh_out,[fn_path trial '2.hdr']);

%{
fn = [fn_path vnir_trial];
std_vnir = single(zeros(p0_vnir,1));
tmp = single(zeros(m0,n0,N));
pnt = 1;
for p = 1:p0_vnir
    tmp(:,:,pnt) = single(multibandread(fn,[m0,n0,p0_vnir],'single',0,'bsq','ieee-le',{'Band','Direct',p}));
    tmp1 = tmp(:,:,pnt);
    tmp1 = tmp1(:);
    pnt = pnt + 1;
    if (pnt > N)
        pnt = 1;
    end
    std_vnir(p) = std(tmp1);
    clear tmp1
end
clear tmp

fn = [fn_path nir_trial];
std_nir = single(zeros(p0_nir,1));
tmp = single(zeros(m0,n0,N));
pnt = 1;
for p = 1:p0_nir
    tmp(:,:,pnt) = single(multibandread(fn,[m0,n0,p0_nir],'single',0,'bsq','ieee-le',{'Column','Range',[1 n0]},{'Row','Range',[1 m0]},{'Band','Direct',p}));
    tmp1 = tmp(:,:,pnt);
    tmp1 = tmp1(:);
    pnt = pnt + 1;
    if (pnt > N)
        pnt = 1;
    end
    std_nir(p) = std(tmp1);
    clear tmp1
end
clear tmp

plot(lambda_nir,std_nir,'b')
hold on
plot(lambda_vnir,std_vnir,'r')
hold off
%}

% gain for each pixel
fn_vnir = [fn_path vnir_trial];
fn_nir = [fn_path nir_trial];
A1 = single(ones(N,2));
cnt = 1;
for p = (p0_vnir_last - N):(p0_vnir_last + N)
    A1(cnt,1) = lambda_vnir(p);
    cnt = cnt + 1;
end
pA1 = pinv(A1);
A2 = single(ones(N,2));
cnt = 1;
for p = (p0_nir_first - N):(p0_nir_first + N)
    A2(cnt,1) = lambda_nir(p);
    cnt = cnt + 1;
end
pA2 = pinv(A2);
A3 = single(ones(N,1));
gn = single(zeros(m0,n0));
offs = single(zeros(m0,n0));
y1 = single(multibandread(fn_vnir,[m0,n0,p0_vnir],'single',0,'bsq','ieee-le',{'Band','Range',[(p0_vnir_last - N) (p0_vnir_last + N)]}));
y2 = single(multibandread(fn_nir,[m0,n0,p0_nir],'single',0,'bsq','ieee-le',{'Band','Range',[(p0_nir_first - N) (p0_nir_first + N)]}));
for m = 1:m0
    for n = 1:n0
        y11 = y1(m,n,:);
        y11 = y11(:);
        coeffs1 = pA1*y11;
        if (sum(y11)>0)
            y21 = y2(m,n,:);
            y21 = y21(:);
            coeffs2 = pA2*y21;
            y22 = A2*coeffs1;
            gn(m,n) = pinv(y21)*y22;
        end
    end
end
%{
% single gain
y2 = single(multibandread(fn_vnir,[m0,n0,p0_vnir],'single',0,'bsq','ieee-le',{'Band','Range',[p0_vnir_last (p0_vnir_last + N - 1)]}));
y1 = single(multibandread(fn_nir,[m0,n0,p0_nir],'single',0,'bsq','ieee-le',{'Band','Range',[(p0_nir_first - N + 1) p0_nir_first]}));
tmp = p0_vnir_last:(p0_vnir_last + N - 1);
tmp = reshape(tmp,[1 1 N]);
x2 = repmat(tmp,[m0 n0 1]);
clear tmp
tmp = (p0_nir_first - N + 1):p0_nir_first;
tmp = reshape(tmp,[1 1 N]);
x1 = repmat(tmp,[m0 n0 1]);
clear tmp

y1 = y1(:);
y2 = y2(:);
x1 = x1(:);
x2 = x2(:);
msk = (y1 > 0);
se = strel('square',15);
msk = imerode(msk,se);
y1(msk == 0) = [];
y2(msk == 0) = [];
x1(msk == 0) = [];
x2(msk == 0) = [];
A = single(ones(length(x1),2));
A(:,1) = x1;
coeffs = pinv(A)*y1;
A(:,1) = x2;
y21 = A*coeffs;
A = single(ones(length(x1),1));
A(:,1) = y2;
g = pinv(A)*y21;
%}
% adjust cubes
tmp = single(multibandread(fn_vnir,[m0,n0,p0_vnir],'single',0,'bsq','ieee-le',{'Column','Range',[1 n0]},{'Row','Range',[1 m0]},{'Band','Direct',1}));
%multibandwrite(g*tmp,[fn_path trial '1'],'bsq',[1 1 1],[m0 n0 (p0_vnir_last+p0_nir-(p0_nir_first-1))]);
multibandwrite(gn.*tmp,[fn_path trial '2'],'bsq',[1 1 1],[m0 n0 (p0_vnir_last+p0_nir-(p0_nir_first-1))]);
clear tmp
cnt = 1;
for p = 1:p0_vnir_last
    tmp = single(multibandread(fn_vnir,[m0,n0,p0_vnir],'single',0,'bsq','ieee-le',{'Column','Range',[1 n0]},{'Row','Range',[1 m0]},{'Band','Direct',p}));
    %multibandwrite(g*tmp,[fn_path trial '1'],'bsq',[1 1 cnt],[m0 n0 (p0_vnir_last+p0_nir-(p0_nir_first-1))]);
    multibandwrite(tmp,[fn_path trial '2'],'bsq',[1 1 cnt],[m0 n0 (p0_vnir_last+p0_nir-(p0_nir_first-1))]);
    cnt = cnt + 1;
    clear tmp
end
for p = p0_nir_first:p0_nir
    tmp = single(multibandread(fn_nir,[m0,n0,p0_nir],'single',0,'bsq','ieee-le',{'Column','Range',[1 n0]},{'Row','Range',[1 m0]},{'Band','Direct',p}));
    %multibandwrite(tmp,[fn_path trial '1'],'bsq',[1 1 cnt],[m0 n0 (p0_vnir_last+p0_nir-(p0_nir_first-1))]);
    multibandwrite(gn.*tmp,[fn_path trial '2'],'bsq',[1 1 cnt],[m0 n0 (p0_vnir_last+p0_nir-(p0_nir_first-1))]);
    cnt = cnt + 1;
    clear tmp
end
