function [flag] = k_out_of_n (k,W1)
sw1 = size(W1);
i = sw1(1)-1;
for l = 1 : 1 : i + 1
    if sum(W1(l,:)~=0)/2 > k 
        flag = 0;
        break
    end

end
        