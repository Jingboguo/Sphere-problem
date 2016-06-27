
promt1 = 'What are the values of k, n and i';
U = input (promt1);
promt2 = 'Please input the matrix of the configuration'; % should be redesigned
W = input (promt2);
k = U(1);
n = U(2);
i = U(3);
phi = 180/n; % angles on each planes
theta = 180/i ; % angles of projected vertical planes
% Check if the input is correct
% for j = 1 : 1 : i
%     for l= 1 : 1 : n
%         if 180 / i *(j-1) == 180 / n * (l-1)
%             flag = isequal (W(j, 1),W(j, (n+1)),W((i+1), l),W((i+1), (l+n)));
%             if  flag
%                 %disp('Inputs are correct');
%             else
%                 disp(['Inputs are not correct, check ', num2str(j), 'th vertical plane' ]);
%             end
%             break;
%         end
%     end
% end
%             

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
% % Creat axis from point (1,0), first check all the existed lines (projected vertical planes)
% % (2,2,i-1) means x and y coordinates of first and perpendicular axis of
% % (i-1)th pair of axis
% Axis1 = ones (2, 2, i); % x, y coordinates for ith axis and its perpendicular axis 
% syms xnew ynew
% x1 = 1;
% y1 = 0;
% for j = 1: 1: i
%     [slov1 , slov2] = solve([(x1 .* xnew + y1 .* ynew) ./ (x1 .^2 + y1 .^2) == cos((j-1)*180/i .* pi ./180), xnew ^2 + ynew ^2 == x1 ^2 + y1 ^2], [xnew, ynew]);
%     Axis1 (1, 1, j) = slov1 (1); 
%     Axis1 (1, 2, j) = slov2 (1);
%     [slov1p , slov2p] = solve([(Axis1 (1, 1, j) .* xnew + Axis1 (1, 2, j) .* ynew) ./ (Axis1 (1, 1, j) .^2 + Axis1 (1, 2, j) .^2) == cos(90 .* pi ./180), xnew ^2 + ynew ^2 == (Axis1 (1, 1, j)) ^2 + (Axis1 (1, 2, j)) ^2], [xnew, ynew]);
%     Axis1 (2, 1, j) = slov1p (1);
%     Axis1 (2, 2, j) = slov2p (1);
% end
% Axis1(abs(Axis1)<1e-3)=0;
% %Check if the system is balanced according to these axis 
% % Angle1 is a matrix that stores all the angles. (i, 2*n, i) means the angle between 
% % (i)th axis in Axis1 and 2*n units on ith vetical plane 
% Angle1 = ones (i, 2*n, i+1);
% Angle1p = ones (i, 2*n, i+1);
% for j = 1: 1: i
%     for l = 1: 1: 2*n
%         for m = 1: 1: i+1
%             Angle1(j, l, m) = real(acos(((Axis1(1, 1, j) .* VP (l, 1, m) + Axis1 (1, 2, j) .* VP (l, 2, m))./ ((((Axis1 (1, 1, j)) .^2 + (Axis1 (1, 2, j)) .^2) .^ 0.5) .* (VP (l, 1, m) .^2 + VP (l, 2, m) .^2).^ 0.5))).*180./pi);
%             Angle1p(j, l, m) = real(acos(((Axis1(2, 1, j) .* VP (l, 1, m) + Axis1 (2, 2, j) .* VP (l, 2, m))./ ((((Axis1 (2, 1, j)) .^2 + (Axis1 (2, 2, j)) .^2) .^ 0.5) .* (VP (l, 1, m) .^2 + VP (l, 2, m) .^2).^ 0.5))).*180./pi);
%         
%         end
%     end
% end
% Angle1(isnan(Angle1))=0;
% Angle1p(isnan(Angle1p))=0;
% Angle1(Angle1>90.001)= 0;
% Angle1(Angle1< -90.001)= 0;
% Angle1p(Angle1p > 90.001)= 0;
% Angle1p(Angle1p< -90.001)= 0;
% 
% % Calculate the moment difference %debug 06052016 abs(Z)should not be in
% % the calculation of horizontal plane
% MD1 = zeros (i+1,2*n, i); % (i+1)th vertical plane, ith axis
% MD1p = zeros (i+1,2*n, i);
% h = zeros (i+1, 2*n, i);
% hp = zeros (i+1, 2*n, i);
% for j = 1 : 1 : i % ith axis
%     for l = 1 : 1 : 2*n 
%         for m = 1: 1: i %vp 
%             h(m , l, j) = Axis1(1, 1, j) * VP (l, 2, m)-Axis1(1, 2, j) * VP (l, 1, m);
%             hp(m, l, j) = Axis1(2, 1, j) * VP (l, 2, m)-Axis1(2, 2, j) * VP (l, 1, m);
%             h(abs(h)<1e-3)=0;
%             hp(abs(hp)<1e-3)=0;
%             h (m , l, j) = sign(h (m , l, j));
%             hp (m , l, j) = sign (hp (m , l, j));
%             MD1 (m, l, j) = h(m , l, j)* sin (Angle1(j, l, m)*pi./180) .* W (m, l).* abs(Z(l,1)) ;
%             MD1p (m, l, j) = hp(m, l, j) * sin (Angle1p(j, l, m)*pi./180) .* W (m, l).* abs(Z(l,1));
%         end
%     end
% end 
% for j = 1 : 1 : i % ith axis
%     for l = 1 : 1 : 2*n 
%         h((i+1) , l, j) = Axis1(1, 1, j) * VP (l, 2, (i+1))-Axis1(1, 2, j) * VP (l, 1, (i+1));
%         hp((i+1), l, j) = Axis1(2, 1, j) * VP (l, 2, (i+1))-Axis1(2, 2, j) * VP (l, 1, (i+1));
%         h(abs(h)<1e-3)=0;
%         hp(abs(hp)<1e-3)=0;
%         h ((i+1), l, j) = sign(h ((i+1) , l, j));
%         hp ((i+1) , l, j) = sign (hp ((i+1) , l, j));
%         MD1 ((i+1), l, j) = h((i+1) , l, j)* sin (Angle1(j, l, (i+1))*pi./180) .* W ((i+1), l).* abs(Z(1,1)) ;
%         MD1p ((i+1), l, j) = hp((i+1), l, j) * sin (Angle1p(j, l, (i+1))*pi./180) .* W ((i+1), l).* abs(Z(1,1));
%     end
% end 
% MD1(abs(MD1)<1e-3)=0;
% MD1p(abs(MD1p)<1e-3)=0;
% S1 = sum (sum(MD1));
% S1p = sum (sum(MD1p));
% S1(abs(S1)<1e-3)=0;
% S1p(abs(S1p)<1e-3)=0;
% for j = 1 : 1 : i
%     if S1 (:,:,j) == 0 && S1p (:,:,j) == 0
%         disp (['System is balanced according to ', num2str(j), 'th axis']);
%     end
% end

% Now cut in half !!!question about solve equation (is the second equation correct???)
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
                MD (m, l, j) = h(m, l, j) * sin (Angle(j, l, m)*pi/180) .* W (m, l) *abs(Z(l, 1)) ;
                MDp (m, l, j) = hp(m, l, j) * sin (Anglep(j, l, m)*pi/180) .* W (m, l) *abs(Z(l, 1)) ;
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
            disp (['System is balanced according to ', num2str(j), 'th axis, when we cut ', num2str(count), ' times']);
            br = 2;
        end
    end
    if br == 2,
        break
    end
end
% Rebalance the system
S = round(S,3);
Sp = round(Sp,3);
MD = round(MD,3);
MDp = round(MDp,3);
W1 = W;
for count = 1: 1: t;
    for j = 1 : 1 : (2^count*i)/2 %axis
        for m = 1 : 1 : i %vp
            for l = 1 : 1 : 2*n %unit
                if l >6
                    id = l - 6;
                else
                    id = l + 6;
                end
                if (W(m,l) == 1) && (S(:,:,j) == MD(m, l, j)) && ((Sp(:,:,j) == MDp(m, l, j)) || (Sp(:,:,j) == MDp(m, id, j))) 
                    W2 = W;
                    W2(m, l) = 0;
                    W2(m, id) = 0;
                    W1 = [W1; W2];
                    disp (['Shut down ',num2str(l), 'th pair on ', num2str(m), 'th vertical plane make system balance according to ',num2str(j), 'th axis']);
                end
                
            end
            
        end 
    end
end

% rebalance the unbalanced system
% First, make a small vp matrix to denote the circle!
% vp(:,1) = VP(1,1,1:i);
% vp(:,2) = 1:1:i;
% vp = repmat(vp,3,1);
% vp = round(vp, 5);
% Axis = round(Axis, 5);
% Axis = repmat(Axis,2,1);
% vpp = ones(i/2,2,2^count*i);
% for j = 1 : 1 : 2^count*i %axis
%     for l = 1 : 1 : i
%         if isequal(Axis(j,1),vp(i+l,1))
%             for m =1 : 1 : i/2
%                 vpp(m, :, j) = [vp(i+l+m,2),vp(i+l-m,2)];
%             end
%         elseif isequal (Axis(j+1,1),vp(i+l,1))
%             for m =1 : 1 : i/2
%                 vpp(m, :, j) = [vp(i+l+m-1,2),vp(i+l-m,2)];
%             end
%             
%         end
%     end
% end
% this method has some problems... work on it later
% for j = 1 : 1 : 2^count*i
%     if S(:,:,j)~=0 || Sp(:,:,j)~=0
%         for m = 1 : 1 : i/2
%             if abs(sum(MD(vpp(m,1,j),:,j))) > abs(sum(MD(vpp(m,2,j),:,j)))
%                 for l = 1 : 1 : 2*n
%                     if MD(vpp(m,1,j),l,j)~=0 && MD(vpp(m,1,j),l,j) == round(sum(MD(vpp(m,1,j),:,j))+sum(MD(vpp(m,2,j),:,j)),3) && MDp(vpp(m,1,j),l,j)~=0 && MDp(vpp(m,1,j),l,j) == round(sum(MDp(vpp(m,1,j),:,j))+ sum(MDp(vpp(m,2,j),:,j)),3)
%                         disp (['Shut down ',num2str(l), 'th pair on ', num2str(vpp(m,1,j)), 'th vertical plane make system balance according to ',num2str(j), 'th axis']);
%                     end 
%                 end
%             elseif abs(sum(MD(vpp(m,1,j),:,j))) < abs(sum(MD(vpp(m,2,j),:,j)))
%                 for l = 1 : 1 : 2*n
%                     if MD(vpp(m,1,j),l,j)~=0  && MD(vpp(m,2,j),l,j) == round(sum(MD(vpp(m,1,j),:,j))+ sum(MD(vpp(m,2,j),:,j)),3) && MDp(vpp(m,1,j),l,j)~=0 && MDp(vpp(m,2,j),l,j) == round(sum(MDp(vpp(m,2,j),:,j))+ sum(MDp(vpp(m,2,j),:,j)),3)
%                         disp (['Shut down ',num2str(l), 'th pair on ', num2str(vpp(m,2,j)), 'th vertical plane make system balance according to ',num2str(j), 'th axis']);
%                     end 
%                 end
%             end
%         end
%     end
% end

% 




