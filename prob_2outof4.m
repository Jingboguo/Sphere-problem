% function that calculate the probability that i pairs failed in a working
% system
function y = prob_2outof4 (i, t)
    T01 = [1, 2, 1, 1];
    T12 = [2, 0, 0, 1, 0, 1, 0, 0; %1
           1, 1, 0, 0, 1, 0, 1, 0; %2
           0, 0, 0, 0, 0, 1, 2, 1; %3 shared
           0, 0, 0, 1, 2, 0, 0, 1];%4
    T23 = [0, 0, 1, 0, 0, 0, 0; %1
           0, 0, 0, 1, 2, 0, 0; %2
           1, 2, 0, 0, 0, 0, 0; %4
           0, 0, 2, 0, 0, 1, 0; %5 
           0, 0, 1, 1, 0, 0, 1; %7
           0, 0, 0, 0, 0, 1, 0; %10 shared
           0, 0, 0, 0, 0, 0, 1; %11 shared
           0, 0, 0, 0, 0, 1, 2]; %12 shared
    T34 = [2, 0, 0; %2
           1, 0, 1; %4
           0, 0, 0; %6
           0, 2, 0; %7
           0, 1, 1; %9
           0, 0, 0; %11 shared
           0, 0, 0]; %13 shared

    % Now define the reliability function of each unit (Suppose Normal distributed for each unit)
    mu = 50;
    sigma = 2;
    theta = 8;
    %F = @(t) expcdf(t, mu);
    %F = @(t) wblcdf(t, mu, sigma);
    F = @(t) normcdf(t, mu, theta);

    % Now define the relibility function of the system
    % number of transition path
    n = zeros (5, 1); 
    n(1) = 1; % 0 pairs failed
    n(2) = sum (T01, 2); % 1 pair failed
    n(3) = sum (T01*T12, 2); % 2 pairs failed 
    n(4) = sum (T01*T12*T23, 2); % 3 pairs failed 
    n(5) = sum (T01*T12*T23*T34, 2);  % 4 pairs failed  
    y = (n(i)/factorial(i-1))*(1-F(t))^(16-2*(i-1))*((1-(1-F(t))^2)^(i-1)); 
    
end

