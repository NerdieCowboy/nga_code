%function [missing_penalty] = compute_missing_penalty(features_table_fn)

features_table_fn = 'C:\damon\dissertation\features_pigments_throughRenasaince.csv';

%features
pigment_family_col = 1;
pigment_col = 2;
color_col = 3;
reflectance_col = 4;
transition_edge_col = 5;
transition_edge_fwhm_col = 6;
transition_edge_asym_col = 7;
absorption_col = 8;
strong_absorption_center_col = 9;
strong_absorption_width_col = 10;
gradual_slope_col = 11;
interesting_col = 12;
organic_col = 13;
xrf_cols = [14 113];

% spectral features for each pigment
fid = fopen(features_table_fn);
features_table = textscan(fid,'%s','delimiter',',','EmptyValue',NaN);
fclose(fid);
tmp = reshape(features_table{1},[xrf_cols(2),33]);
tmp = tmp';
pigment_ris_table = tmp(2:33,reflectance_col:interesting_col);
clear tmp
[Lp,Lf] = size(pigment_ris_table);

% compute missing data penalty
missing_penalty = single(ones(Lp*Lf*100,5));
cnt = 1;
for i = 1:Lp
    for j = 1:Lf
        tmp = pigment_ris_table(i,j);
        tmp = tmp{1};
        if (~strcmp(tmp,'0'))
            tmp1 = regexp(tmp, '\;', 'split');
            Ltmp1 = length(tmp1);
            for k = 1:Ltmp1
                if (j == (reflectance_col-3))
                    pkdelta1 = mean(rpeak_list(:,5));
                    Lpk1 = length(rpeak_list(:,5));
                    missing_penalty(cnt,:) = [i j k pkdelta1 Lpk1];
                    cnt = cnt + 1;
                elseif (j == (transition_edge_col-3))
                    pkdelta1 = mean(edge_list(:,5));
                    Lpk1 = length(edge_list(:,5));
                    missing_penalty(cnt,:) = [i j k pkdelta1 Lpk1];
                    cnt = cnt + 1;
                    missing_penalty(cnt,:) = [i j+1 k pkdelta1 Lpk1];
                    cnt = cnt + 1;
                    missing_penalty(cnt,:) = [i j+2 k pkdelta1 Lpk1];
                    cnt = cnt + 1;
                elseif (j == (absorption_col-3))
                    pkdelta1 = mean(apeak_list(:,5));
                    Lpk1 = length(apeak_list(:,5));
                    missing_penalty(cnt,:) = [i j k pkdelta1 Lpk1];
                    cnt = cnt + 1;
                elseif (j == (strong_absorption_center_col-3))
                    tmp2 = pigment_ris_table(i,j+1);
                    tmp2 = tmp2{1};
                    tmp2 = regexp(tmp2, '\;', 'split');
                    val = single(str2double(tmp2));
                    missing_penalty(cnt,:) = [i j k val(k) 1];
                    cnt = cnt + 1;
                    missing_penalty(cnt,:) = [i j+1 k val(k) 1];
                    cnt = cnt + 1;
                elseif (j == (gradual_slope_col-3))
                    missing_penalty(cnt,:) = [i j k 1 1];
                    cnt = cnt + 1;
                elseif (j == (interesting_col-3))
                    missing_penalty(cnt,:) = [i j k 1 1];
                    cnt = cnt + 1;
                end
            end
        end
    end
end
missing_penalty(cnt:end,:) = [];
