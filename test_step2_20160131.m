%% STEP TWO : Check if all the vertical planes are balanced
h=0; % count the number of vertical planes
for i = 1 : 1 : m
    if VP (i, 2) == 0 % m+1 should be 2 if I solve the problem in first step
        break
    elseif VP (i, 2) ~= 0
        h = h +1;
    end
end
% vp is a matrix of grouped theta 
vp = ones (h , m);
for i = 1 : 1 : h
    for j = 1 : 1 : m
        vp (i, j) = U(VP(i, (j+1)) , 1);
    end   
end
% w is a matrix of grouped working status according to matrix vp
w = ones (h , m);
for i = 1 : 1 : h
    for j = 1 : 1 : m
        w (i, j) = W (VP (i, (j+1)));
    end
end
% Calculation of balance
M = zeros (h , m);
for i = 1 : 1 : h
    for j = 1 : 1 : m 
        M(i , j) = sin (vp (i, j)) .* (w(i, j));
    end
end
S = sum (M , 2);
for i = 1 : 1 : h
    if S(i) ~= 0
        disp('vertical plane', i, 'is not balanced')
    end
end