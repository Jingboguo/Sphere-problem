function [flag] = k_out_of_n (k,W1) % k out of n good system
sw1 = size(W1);
i = sw1(1)-1;
n = sw1(2)/2;
for l = 1 : 1 : i + 1
    if sum(W1(l,:)==0)/2 > (n-k) 
        flag = 0;
        break
    else
        flag = 1;
    end
end
        