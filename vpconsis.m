function [flag, l] = vpconsis (i, n, m) 
% This function calculate the (flag)th unit on horizontal plane should be 0
% when first pair on (l)th vertical plane is zero
% i: # of vp; n: # of pairs;  m: m th pair on horizontal plane
flag = 0;
l = 0;
ratio = n/i;
for j = 1 : 1 : i
    if m == (j-1)*ratio +1
        flag = 1;
        l = j;
        break
    end
end
end
