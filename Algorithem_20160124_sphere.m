%STEP ONE: Group all the points into different vertical planes
% Input: points coordinate U=[theta1,x1,y1;theta2,x2,y2;...thetan,xn,yn]; 
% working status W=[w1;w2;w3;...w4]
promt1 = 'What are the coordinates of the points';
U = input (promt1);
promt2 = 'What are the working situation of these points';
W = input (promt2);
% Assign coordinates into new vectors theta, X and Y
[m, n] = size (U);
theta = vector (m);
X = vector (m); 
Y = vector (m);
for i = 1 : 1 : m
    theta(i) = U (i, 1);
    X(i) = U (i, 2);
    Y(i) = U (i, 3);
end
% Group all the points into different veritcal planes
ratio = vector (m);
for i = 1 : 1 : m
    ratio(i) = X(i) ./ Y(i);
end

VP = matrix (m, (m+1));
VP (1, 1) = ratio (1); 
VP (1, 2) = 1;
for i = 2 : 1 : m
    for j = 1 : 1 : m
        if VP (j, 1) == ratio (i)
            for k = 2 : 1 : (m+1)
                if VP (j, k) == 0
                   VP (j, k) = i;
                   break
                end
            end
        elseif VP (j, 1) == 0
            VP (j, 1) = ratio (i);
            VP (j, 2) = i;
        end
    end
end

% STEP TWO : Check if all the vertical planes are balanced
% vp is a matrix of grouped theta 
h=0; % count the number of vertical planes
for i = 1 : 1 : m
    if VP (i, 2) == 0
        break
    elseif VP (i, 2) ~= 0
        h = h +1;
    end
end

vp = matrix (h , m);
for i = 1 : 1 : h
    for j = 2 : 1 : (m+1)
        vp (i, (j-1)) = U(1, VP(i, j));
    end   
end
% w is a matrix of grouped working status
w = matrix (h , m);
for i = 1 : 1 : h
    for j = 1 : 1 : m
        w (i, j) = W (VP (i, (j+1)));
    end
end
% Calculation of balance
M = matrix (h , m);
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

% STEP THREE: Calculate the moment difference for each pair of axis
% First, since units on the same vertical plane will be on the same axis,
% we extract one point from each vertical plane and save them in matrix P.
P = matrix (h, 2);
for i = 1 : 1 : h
    if VP(i,2) == 0
        break
    elseif VP(i,2) ~= 0
        P (i, 1) = X(VP(i, 2));
        P (i, 2) = Y(VP(i, 2));
    end
end
% Second, we calculate the angles for each axis
phi = matrix (h, m);
for i = 1 : 1 : h
    for j = 1 : 1 : m
        if P (i, 1) ~= U (j, 2) && (P (i, 2) ~= U (j, 3))
            phi(i, j) = arccos ((P(i, 1) .* U (j, 2) + P (i, 2) .* U (j, 3))./ (((P (i, 1) .^2 + P (i, 2) .^2)) .^ 0.5 .* (U (j, 2) .^2 .* U (j, 2) .^2).^ 0.5));
        end
    end
end
% make sure that the angles are within +- 90
for i = 1 : 1 : h
    for j = 1 : 1 : m
        if phi(i, j) > 45 || phi (i, j) < -45
            phi(i, j) = 0;
        end
    end
end
% Calculate the moment difference for each axis
MD1 = matrix (h, m);
for i = 1 : 1 : h
    for j = 1 : 1 : m
        MD1 (i, j) = sin (phi (i, j)) .* W (j);
    end
end
% Calculate the moment difference for perpendicular axis
% Solving equations for points' coordination in the perpendicular axis
P2 = matrix (h, 2);
syms x1 x2 y1 y2
for i = 1 : 1 : h
    x1 = P (i, 1);
    y1 = P (i, 2);
    P2 (i, :) = solve ([x1.*x2 + y1.*y2 == 0, x1^2 + y1^2 ==x2^2 + y2^2], [x2, y2]);
end
phi2 = matrix (h, m);
for i = 1 : 1 : h
    for j = 1 : 1 : m
        if P2 (i, 1) ~= U (j, 2) && (P2 (i, 2) ~= U (j, 3))
            phi2(i, j) = arccos ((P2(i, 1) .* U (j, 2) + P2 (i, 2) .* U (j, 3))./ (((P2 (i, 1) .^2 + P2 (i, 2) .^2)) .^ 0.5 .* (U (j, 2) .^2 .* U (j, 2) .^2).^ 0.5));
        end
    end
end
for i = 1 : 1 : h
    for j = 1 : 1 : m
        if phi2(i, j) > 45 || phi2 (i, j) < -45
            phi2(i, j) = 0;
        end
    end
end
MD2 = matrix (h, m);
for i = 1 : 1 : h
    for j = 1 : 1 : m
        MD2 (i, j) = sin (phi2 (i, j)) .* W (j);
    end
end
S1 = sum (MD1, 2);
S2 = sum (MD2, 2);
for i = 1 : 1 : h
    if S1 (i) == S2 (i)
        disp ('System is balanced according to ', i, 'th axis')
    else
        disp ('System is not balanced according to ', i, 'th axis')
    end
end


% STEP FOUR: Calulate the moment difference for all the other axis 
P3 = matrix (359 * h, 2);
syms x1 xnew y1 ynew
for i = 1 : 1 : h
    for angle = 1 : 1 : 359
        x1 = P (i, 1);
        y1 = P (i, 2);
        r = (x1^2 + x2^2)^0.5;
        P3 (i.*angle, :) = solve ([x1.*x2 + y1.*y2 == cos(angel) * r^2, xnew^2 + ynew^2 ==x1^2 + y2^2], [xnew, ynew]);
    end
end
        
            
    


            



        

        
        
        




    



                
            

    





    