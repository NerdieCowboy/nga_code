features_table_fn = 'C:\damon\dissertation\features_pigments_throughRenasaince.csv';
pigment_family_col = 1;
xrf_cols = [14 113];

fid = fopen(features_table_fn);
features_table = textscan(fid,'%s','delimiter',',','EmptyValue',NaN);
fclose(fid);
tmp = reshape(features_table{1},[xrf_cols(2),33]);
tmp = tmp';
pigment_family = tmp(2:33,pigment_family_col);
clear tmp

sthresh = single(zeros(32,1));
score_test = score;
for i = 1:32
    pigment_family0 = pigment_family(i);
    family_msk = strcmp(pigment_family,pigment_family0);
    tmp3 = scoreL(:,i);
    tmp3(family_msk==1) = 0;
    [sthresh(i),~] = max(tmp3);
    if (sthresh(i) == 0)
        sthresh(i) = Inf;
    end
    tmp = score_test(:,:,i);
    tmp(tmp<=sthresh(i)) = 0;
    score_test(:,:,i) = tmp;
    clear tmp
end

%azurite
tmp = score_test(100,80,:);
tmp = tmp(:);
figure,plot(tmp,'b')
hold on
tmp2 = scoreL(13,:)';
plot(tmp2,'r')
hold on
tmp3 = scoreL(:,13);
plot(tmp3,'g')
hold off

%azurite
tmp = score_test(109,20,:);
tmp = tmp(:);
figure,plot(tmp,'b')
hold on
tmp2 = scoreL(13,:)';
plot(tmp2,'r')
hold on
tmp3 = scoreL(:,13);
plot(tmp3,'g')
hold off

%azurite
tmp = score_test(109,20,:);
tmp = tmp(:);
figure,plot(tmp,'b')
hold on
tmp2 = scoreL(14,:)';
plot(tmp2,'r')
hold on
tmp3 = scoreL(:,14);
plot(tmp3,'g')
hold off

% lead white
tmp = score_test(94,145,:);
tmp = tmp(:);
figure,plot(tmp,'b')
hold on
tmp2 = scoreL(18,:)';
plot(tmp2,'r')
hold on
tmp3 = scoreL(:,18);
plot(tmp3,'g')
hold off

% green
tmp = score_test(102,162,:);
tmp = tmp(:);
figure,plot(tmp,'b')
hold on
tmp2 = scoreL(23,:)';
plot(tmp2,'r')
hold on
tmp3 = scoreL(:,23);
plot(tmp3,'g')
hold off
