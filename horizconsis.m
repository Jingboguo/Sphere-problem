function [flag] = horizconsis (i, n, l) 
% This function calculate the (flag)th unit on horizontal plane should be 0
% when first pair on (l)th vertical plane is zero
% i: # of vp; n: # of pairs; l: vp m: units
ratio = n/i;
flag = (l-1)*ratio+1;
end
