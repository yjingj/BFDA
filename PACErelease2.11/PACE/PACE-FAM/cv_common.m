%  cross-validation to choose bandwidth choice for 1-d local polynomial
%  smoothing
function bw_opt = cv_common(kernel,nwe,npoly,nder,xin,yin)

nbw = 21;
bw = zeros(1,nbw-1);
rang = max(xin)-min(xin);
n = length(xin);
for i = 1:nbw-1
    bw(i) = 2.5*rang/n*(n/2.5)^((i-1)/(nbw-1));
end

pred = zeros(1,length(bw)); 
count = zeros(1,length(bw));
for k = 1:length(bw)
    for i = 1:length(xin)
        win = ones(1,length(xin));
        win(i) = 0;
        [invalid, temp] = lwls(bw(k),kernel,nwe,npoly,nder,xin,yin,win,xin(i));
        if invalid == 0
            pred(k) = pred(k)+(yin(i)-temp)^2;
            count(k) = count(k)+1;
        end
    end
end

ind = find(count/length(xin)>0.90);
cv = pred(ind)./count(ind);
bw = bw(ind);
bw_opt = bw(cv==min(cv));
