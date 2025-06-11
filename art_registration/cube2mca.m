fn_path = 'C:\\data\\ApostlePaul';
fn = 'ap_top_right';

fnh = [fn_path '\\' fn '.hdr.'];
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
    clear line msk
end
fclose(fid);

fn1 = [fn_path '\\' fn];
if (exist(fn1,'dir') ~= 7)
    mkdir([fn1 '1']);
end
%cnt1 = 1;
ncores = feature('NumCores');
isOpen = matlabpool('size') > 0;
if (isOpen == 1)
    matlabpool close
end    
matlabpool(ncores)
for n = 1:n0
    parfor m = 1:m0
        cnt1 = n + (m-1)*n0;
        fn1 = [fn_path '\\' fn '1\\' fn '_' num2str(sprintf('%09d',cnt1))];
        tmp = single(multibandread([fn_path '\\' fn],[m0,n0,p0],'single',0,'bsq','ieee-le',{'Column','Direct',n},{'Row','Direct',m},{'Band','Range',[1 p0]}));
        tmp = tmp(:);
        tmp(isnan(tmp)) = 0;

        fid1 = fopen([fn1 '.mca'],'w');
        fprintf(fid1,['#F ' fn_path '\\' fn '.mca\n']);
        fprintf(fid1,'#D Mon Jan 07 00:00:01 2013\n');
        fprintf(fid1,'\n');
        fprintf(fid1,['#S 1 ' fn ' 1.1 Column 1\n']);
        fprintf(fid1,'#D Mon Jan 07 00:00:01 2013\n');
        fprintf(fid1,'#@MCA %%16C\n');
        fprintf(fid1,['#@CHANN ' num2str(p0) ' 0 ' num2str(p0-1) ' 1\n']);
        fprintf(fid1,'#@CALIB 1.345 0.0137 0\n');
        fprintf(fid1,'@A ');

        cnt = 1;
        for i = 1:p0
            fprintf(fid1,num2str(sprintf('%0.4f ',tmp(i))));
            cnt = cnt + 1;
            if (cnt > 16)
                cnt = 1;
                fprintf(fid1,' \\\n');
            end
        end
        fprintf(fid1,'\n');
        fclose(fid1);
        %cnt1 = cnt1 + 1;
    end
end
matlabpool close
