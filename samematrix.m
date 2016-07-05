function [W2] = samematrix(W)
sw = size (W);
indexx = zeros(1);
for i = 1 : 1 : sw(3)
    for j = i+1 : 1 : sw(3)
        if W(:,:,i) == W(:,:,j)
            indexx = [indexx, j];
        end
    end
end
indexx = sort(indexx,'descend');
s = size(indexx);
for i = 1 : 1 : s
    W = W(:,:,[1:indexx(i)-1,indexx(i)+1:end]);
end
W2 = W
