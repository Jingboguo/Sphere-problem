function [W3] = transition (k,t,W1)
% transition from n pair fail to n+1 pair fail with k out of n working
sw1 = size (W1);
i = sw1(1)-1;
n = sw1(2)/2;
W2 = zeros(sw1);
for l = 1 : 1 : i+1
    for m = 1 : 1 : n
        Wtemp = W1;
        if Wtemp(l, m) ~= 0
            Wtemp(l,m) = 0;
            Wtemp(l,m+n) = 0;
            if l ~= i + 1 && m == 1  % make sure the horizontal plane is consistant with VPs
                Wtemp((i+1), horizconsis(i,n,l)) = 0;
                Wtemp((i+1), horizconsis(i,n,l)+n) = 0;
            end
            if l == i + 1 && vpconsis (i, n, m)
                Wtemp(i/n*(m-1)+1,1) = 0;
                Wtemp(i/n*(m-1)+1,1+n) = 0;
            end
            flag1 = k_out_of_n(k, Wtemp);
            flag2 = sbalance(Wtemp,t);
            if flag1 && flag2
                W2 = [W2; Wtemp];
            end
        end
    end
end
W2 = W2(sw1(1)+1:end,:); % exclude the first zero rows
b = size(W2,1)/(i+1); % # of possible transitions
W3 = zeros(i+1,2*n,b); 
for j = 1 : 1 : b
    W3 (:,:,j)=W2(((j-1)*(i+1)+1):((j-1)*(i+1)+1+i),:);
end
