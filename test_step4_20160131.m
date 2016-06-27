%% STEP FOUR: Calulate the moment difference for all the other axis 
P3 = ones (359 * h, 2);
syms xnew ynew
for i = 1 : 1 : h
    for j = 1 : 1 : 359
        x1 = P (i, 1);
        y1 = P (i, 2);
        [slov1 , slov2] = solve ([(x1 .* xnew + y1 .* ynew) ./ (x1 .^2 + y1 .^2) == cos(j .* pi ./180), xnew ^2 + ynew ^2 == x1 ^2 + y1 ^2], [xnew, ynew]);
        P3 (i.*j, 1) = slov1 (1);
        P3 (i.*j, 2) = slov2 (1);
    end
end
% Second, we calculate the angles for each axis
phi3 = zeros (359*h, m);
for i = 1 : 1 : 359*h
    for j = 1 : 1 : m
        phi3(i, j) = acos(((P3(i, 1) .* U (j, 2) + P3 (i, 2) .* U (j, 3))./ (((P3 (i, 1) .^2 + P3 (i, 2) .^2) .^ 0.5) .* (U (j, 2) .^2 + U (j, 3) .^2).^ 0.5))).*180./pi;
    end
end
% make sure that the angles are within +- 90
for i = 1 : 1 : 359*h
    for j = 1 : 1 : m
        if phi3(i, j) > 45 || phi3 (i, j) < -45
            phi3(i, j) = 0;
        end
    end
end
% Calculate the moment difference for each axis
MD3 = zeros (359*h, m);
for i = 1 : 1 : 359*h
    for j = 1 : 1 : m
        MD3 (i, j) = sin (phi3 (i, j)) .* W (j);
    end
end
% Calculate the moment difference for perpendicular axis
% Solving equations for points' coordination in the perpendicular axis
P4 = ones (359*h, 2);
syms x2 y2
for i = 1 : 1 : 359*h
    x1 = P3 (i, 1);
    y1 = P3 (i, 2);
    [slove1 , slove2] = solve ([x1.*x2 + y1.*y2 == 0, x1^2 + y1^2 ==x2^2 + y2^2], [x2, y2]);
    P4 (i , 1) = slove1 (1);
    P4 (i , 2) = slove2 (1);
end
phi4 = zeros (359*h, m);
for i = 1 : 1 : 359*h
    for j = 1 : 1 : m
        phi4(i, j) = acos(((P4(i, 1) .* U (j, 2) + P4 (i, 2) .* U (j, 3))./ (((P4 (i, 1) .^2 + P4 (i, 2) .^2) .^ 0.5) .* (U (j, 2) .^2 + U (j, 3) .^2).^ 0.5))).*180./pi;
    end
end
for i = 1 : 1 : 359*h
    for j = 1 : 1 : m
        if phi4(i, j) > 45 || phi4 (i, j) < -45
            phi4(i, j) = 0;
        end
    end
end
MD4 = zeros (359*h, m);
for i = 1 : 1 : 359*h
    for j = 1 : 1 : m
        MD4 (i, j) = sin (phi4 (i, j)) .* W (j);
    end
end
S3 = sum (MD3, 2);
S4 = sum (MD4, 2);
for i = 1 : 1 : 359*h
    if S3 (i) == 0 && S4 (i) == 0
        disp (['System is balanced according to ', num2str(i), 'th axis']);
    else
        disp (['System is not balanced according to ', num2str(i), 'th axis']);
    end
end
