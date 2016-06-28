function [flag] = sbalance (n,i,Wtemp)
phi = 180/n; % angles on each planes
theta = 180/i ; % angles of projected vertical planes
% project all the units on the center line of each vertical plane. Z is a
% matrix that record the length of each unit to the center of vp.

Z = zeros (2*n , 1); 
for m = 1 : 1 : 2*n
    Z(m,1) = cos ((m-1)*phi*pi/180);
end
Z(abs(Z)<1e-5)=0;

% project every vertical plane on horizontal plane
M = zeros (i, 2);
for l = 1 : 1 : i
    M(l , 1) = cos ((360-(l-1)*theta)*pi/180);
    M(l , 2) = sin ((360-(l-1)*theta)*pi/180);
end
% coordinates of units on mth vertical plane
VP = zeros (2*n,2,i+1); %ith vertical plane, 2n units' coordinates
for m = 1 : 1 : i
    for l = 1 : 1 : 2*n %lth unit
        VP (l, 1, m) = M(m, 1) * Z(l);
        VP (l, 2, m) = M(m, 2) * Z(l);
    end
end
VP = round (VP, 5);
% horizontal plane
for l = 1 : 1 : 2*n
    VP(l,1,i+1) = cos ((360-(l-1)*phi)*pi/180);
    VP(l,2,i+1) = sin ((360-(l-1)*phi)*pi/180);
end
t = 1;
for count = 1: 1: t;
    Axis = zeros (2^count * i , 2); % x, y coordinates for ith axis and its perpendicular axis 
    Axisp = zeros (2^count * i , 2);
    Angle = zeros (2^count*i, 2*n, i+1);
    Anglep = zeros (2^count*i, 2*n, i+1);
    MD = zeros (i,2*n, 2^count*i); % ith vertical plane, 2^count*ith axis
    MDp = zeros (i,2*n, 2^count*i);
    h = zeros(i+1,2*n, 2^count*i); % sign matrix
    hp = zeros(i+1,2*n, 2^count*i);
    syms xnew ynew
    x1 = 1;
    y1 = 0;
    for l = 1 : 1 : 2^count*i %axis
        [slov1 , slov2] = solve([(x1 .* xnew + y1 .* ynew) ./ (x1 .^2 + y1 .^2) == cos((l-1)*180/i/2^count .* pi ./180), xnew ^2 + ynew ^2 == x1 ^2 + y1 ^2], [xnew, ynew]);
        Axis (l, 1) = slov1 (1); 
        Axis (l, 2) = slov2 (1);
        [slov1p , slov2p] = solve([(Axis (l, 1) .* xnew + Axis (l, 2) .* ynew) ./ ((Axis (l, 1)) .^2 + (Axis (l, 2)) .^2) == cos(90 .* pi ./180), xnew ^2 + ynew ^2 == (Axis (l, 1)) ^2 + (Axis (l, 2)) ^2], [xnew, ynew]);
        Axisp (l, 1) = slov1p(1);
        Axisp (l, 2) = slov2p(1);
    end
    Axis(abs(Axis)<1e-3)=0;
    Axisp(abs(Axisp)<1e-3)=0;
    %Check if the system is balanced according to these axis 
    % Angle1 is a matrix that stores all the angles. (i, 2*n, i) means the angle between 
    % (i)th axis in Axis and 2*n units on ith vetical plane 
    for j = 1: 1: 2^count*i
        for l = 1: 1: 2*n
            for m = 1: 1: i+1
                Angle(j, l, m) = real(acos(((Axis(j, 1) .* VP (l, 1, m) + Axis (j, 2) .* VP (l, 2, m))./ (((Axis (j, 1) .^2 + Axis (j, 2) .^2) .^ 0.5) .* (VP (l, 1, m) .^2 + VP (l, 2, m) .^2).^ 0.5))).*180./pi);
                Anglep(j, l, m) = real(acos(((Axisp(j, 1) .* VP (l, 1, m) + Axisp (j, 2) .* VP (l, 2, m))./ (((Axisp (j, 1) .^2 + Axisp (j, 2) .^2) .^ 0.5) .* (VP (l, 1, m) .^2 + VP (l, 2, m) .^2).^ 0.5))).*180./pi);
            end
        end
    end
    Angle(isnan(Angle))=0;
    Anglep(isnan(Anglep))=0;
    Angle(Angle> 90.001)= 0;
    Angle(Angle< -90.001)= 0;
    Anglep(Anglep> 90.001)= 0;
    Anglep(Anglep< -90.001)= 0;
    % Calculate the moment difference

    for j = 1 : 1 : 2^count*i % ith axis
        for l = 1 : 1 : 2*n 
            for m = 1: 1: i %vp
                h(m , l, j) = Axis(j, 1) * VP (l, 2, m)-Axis(j, 2) * VP (l, 1, m);
                hp(m, l, j) = Axisp(j, 1) * VP (l, 2, m)-Axisp(j, 2) * VP (l, 1, m);
                h(abs(h)<1e-3)=0;
                h(abs(h)<1e-3)=0;
                h (m , l, j) = sign(h (m , l, j));
                hp (m , l, j) = sign (hp (m , l, j));  
                MD (m, l, j) = h(m, l, j) * sin (Angle(j, l, m)*pi/180) .* Wtemp (m, l) *abs(Z(l, 1)) ;
                MDp (m, l, j) = hp(m, l, j) * sin (Anglep(j, l, m)*pi/180) .* Wtemp (m, l) *abs(Z(l, 1)) ;
            end
        end
    end   
    S = sum(MD, 2);
    Sp = sum(MDp, 2);
    S(abs(S)<1e-3)=0;
    Sp(abs(Sp)<1e-3)=0;
    br = 1;
    for j = 1 : 1 : 2^count*i
        if S (j) == 0 && Sp (j) == 0
            flag = 1;
            br = 2;
        else
            flag = 0;
        end
    end
    if br == 2,
        break
    end
end

end