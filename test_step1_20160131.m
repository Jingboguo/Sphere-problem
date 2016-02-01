%% STEP ONE: Group all the points into different vertical planes
% Input: points coordinate U=[theta1,x1,y1;theta2,x2,y2;...thetan,xn,yn]; 
% working status W=[w1;w2;w3;...w4]
promt1 = 'What are the coordinates of the points';
U = input (promt1);
promt2 = 'What are the working situation of these points';
W = input (promt2);
% Assign coordinates into new vectors theta, X and Y
[m, n] = size (U);
theta = ones (m,1);
X = ones (m,1); 
Y = ones (m,1);
for i = 1 : 1 : m
    theta(i,1) = U (i, 1);
    X(i,1) = U (i, 2);
    Y(i,1) = U (i, 3);
end
% Group all the points into different veritcal planes
ratio = ones (m,1);
for i = 1 : 1 : m
    ratio(i,1) = X(i,1) ./ Y(i,1);
end

VP = zeros (m, (m+1));
VP (1, 1) = ratio (1,1); 
VP (1, 2) = 1;
for i = 2 : 1 : m
    for j = 1 : 1 : m
        if ratio (i,1) == VP (j, 1) && (VP (j, 2) ~= 0)
            for k = 2 : 1 : (m+1)
                if VP (j, k) == 0
                   VP (j, k) = i;
                   break %did not work, because the for loop still works
                end
            end
        elseif ratio (i,1) ~= VP (j, 1) && (VP (j, 1) == 0) && (VP (j, 2) == 0)
            VP (j, 1) = ratio (i,1);
            VP (j, 2) = i;
        end
    end
end