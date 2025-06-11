sigma_score = single(zeros(Lp,1));
for i = 1:32
    tmp = score(:,:,i);
    tmp = tmp(:);
    sigma_score(i) = std(tmp);
    clear tmp
end

%azurite
tmp = score(100,80,:);
tmp = tmp(:);
conf = single(zeros(Lp,1));
for i = 1:Lp
    %tmp1 = scoreL(i,:)';
    tmp2 = scoreL(:,i);
    tmp3 = tmp.*tmp2;
    
    figure,stem(tmp,'bo')
    hold on
    stem(tmp2,'go')
    hold on
    stem(tmp3,'rx')
    hold on
    xlabel('pigment number')
    ylabel('score')
    axis([1 32 0 1])
    hold off
    
    [mx,ind] = max(tmp3);
    if (ind == i)
        tmp3(i) = 0;
        [mx0,~] = max(tmp3);
        conf(i) = abs(mx - mx0);
    else
        conf(i) = 0;
    end
end
figure,plot(conf)

%azurite mixture
tmp = score(109,20,:);
tmp = tmp(:);
conf = single(zeros(Lp,1));
for i = 1:Lp
    %tmp1 = scoreL(i,:)';
    tmp2 = scoreL(:,i);
    tmp3 = tmp.*tmp2;
    
    figure,stem(tmp,'bo')
    hold on
    stem(tmp2,'go')
    hold on
    stem(tmp3,'rx')
    hold on
    xlabel('pigment number')
    ylabel('score')
    axis([1 32 0 1])
    hold off
    
    [mx,ind] = max(tmp3);
    if (ind == i)
        tmp3(i) = 0;
        [mx0,~] = max(tmp3);
        conf(i) = abs(mx - mx0);
    else
        conf(i) = 0;
    end
end
figure,plot(conf)
