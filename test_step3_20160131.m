%% STEP THREE: Calculate the moment difference for each pair of axis
% First, since units on the same vertical plane will be on the same axis,
% we extract one point from each vertical plane and save them in matrix P.
P = ones (h, 2);
for i = 1 : 1 : h
    if VP(i,2) == 0
        break
    elseif VP(i,2) ~= 0
        P (i, 1) = X(VP(i, 2));
        P (i, 2) = Y(VP(i, 2));
    end
end
% Second, we calculate the angles for each axis
phi = zeros (h, m);
for i = 1 : 1 : h
    for j = 1 : 1 : m
        %if P (i, 1) ~= U (j, 2) && (P (i, 2) ~= U (j, 3))
        phi(i, j) = acos(((P(i, 1) .* U (j, 2) + P (i, 2) .* U (j, 3))./ (((P (i, 1) .^2 + P (i, 2) .^2) .^ 0.5) .* (U (j, 2) .^2 + U (j, 3) .^2).^ 0.5))).*180./pi;
        %end
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
MD1 = zeros (h, m);
for i = 1 : 1 : h
    for j = 1 : 1 : m
        MD1 (i, j) = sin (phi (i, j)) .* W (j);
    end
end
% Calculate the moment difference for perpendicular axis
% Solving equations for points' coordination in the perpendicular axis
P2 = ones (h, 2);
syms x2 y2
for i = 1 : 1 : h
    x1 = P (i, 1);
    y1 = P (i, 2);
    [slo1 , slo2] = solve ([x1.*x2 + y1.*y2 == 0, x1^2 + y1^2 ==x2^2 + y2^2], [x2, y2]);
    P2 (i , 1) = slo1 (1);
    P2 (i , 2) = slo2 (1);
end
phi2 = zeros (h, m);
for i = 1 : 1 : h
    for j = 1 : 1 : m
        %if P2 (i, 1) ~= U (j, 2) && (P2 (i, 2) ~= U (j, 3))
        phi2(i, j) = acos(((P2(i, 1) .* U (j, 2) + P2 (i, 2) .* U (j, 3))./ (((P2 (i, 1) .^2 + P2 (i, 2) .^2) .^ 0.5) .* (U (j, 2) .^2 + U (j, 3) .^2).^ 0.5))).*180./pi;
        %end
    end
end
for i = 1 : 1 : h
    for j = 1 : 1 : m
        if phi2(i, j) > 45 || phi2 (i, j) < -45
            phi2(i, j) = 0;
        end
    end
end
MD2 = zeros (h, m);
for i = 1 : 1 : h
    for j = 1 : 1 : m
        MD2 (i, j) = sin (phi2 (i, j)) .* W (j);
    end
end
S1 = sum (MD1, 2);
S2 = sum (MD2, 2);
for i = 1 : 1 : h
    if S1 (i) == 0 && S2 (i) == 0
        disp (['System is balanced according to ', num2str(i), 'th axis']);
    else
        disp (['System is not balanced according to ', num2str(i), 'th axis']);
    end
end
