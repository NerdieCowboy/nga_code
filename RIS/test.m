%f = [0 1 0.5];
%f = [0.8 0.9 1];
%f = [0 0.1 0.2];
%f = rand(1,3);
mu = mean(f);
sigma = std(f);
f1 = (f-mu)/sigma;
[h,~] = hist(f1,1001);
nout = min(f1):0.001:max(f1);
cdf = zeros(1,length(h));
lim = [0 0];
lim1 = 0;
lim2 = 0;
for i = 1:length(h)
    cdf(i) = sum(h(1:i));
    if ((cdf(i) >= 0.01*sum(h)) && (lim1 == 0))
        lim(1) = nout(i);
        lim1 = 1;
    end
    if ((cdf(i) >= 0.99*sum(h)) && (lim2 == 0))
        lim(2) = nout(i);
        lim2 = 1;
    end
end
f1 = (f1 - lim(1))/(lim(2) - lim(1))/2;
f1(f1<0) = 0;
f1(f1>1) = 1;
f1
