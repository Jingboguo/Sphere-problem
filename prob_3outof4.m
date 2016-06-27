% function that calculate the probability that i pairs failed in a working
% system
function y = prob_3outof4 (i, t)
    T01 = [1, 2, 1, 1];
    T12 = [1, 0; %1
           0, 1; %2
           0, 0; %3 shared
           1, 2];%4
    % Now define the reliability function of each unit (Suppose Normal distributed for each unit)
    mu = 50;
    sigma = 2;
    theta = 8;
    %F = @(t) expcdf(t, mu);
    %F = @(t) wblcdf(t, mu, sigma);
    F = @(t) normcdf(t, mu, theta);

    % Now define the relibility function of the system
    % number of transition path
    n = zeros (3, 1); 
    n(1) = 1; % 0 pairs failed
    n(2) = sum (T01, 2); % 1 pair failed
    n(3) = sum (T01*T12, 2); % 2 pairs failed 
    y = (n(i)/factorial(i-1))*(1-F(t))^(16-2*(i-1))*((1-(1-F(t))^2)^(i-1)); 
    
end

