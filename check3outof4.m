    T01 = [1, 2, 1, 1];
    T12 = [1, 0; %1
           0, 1; %2
           0, 0; %3 shared
           1, 2];%4
    h = zeros (3, 1); 
    h(1) = 1; % 0 pairs failed
    h(2) = sum (T01, 2); % 1 pair failed
    h(3) = sum (T01*T12, 2); % 2 pairs failed 