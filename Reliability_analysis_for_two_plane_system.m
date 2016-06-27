% Now define the relibility function of the system
count = 80;
 % Probability that 3 out of 4 pairs working, system working
 P1 = zeros (3,count);
 for i = 1 : 1 : 3
     for t = 1 : 1 : count
         P1(i, t) = prob_3outof4 (i, t);
     end
 end
 R1 = sum (P1);
 mean1 = sum (R1);

 % Probability that 2 pairs working and system still maintain working
P2 = zeros (5, count);
for i = 1 : 1 : 5
    for t = 1 : 1 : count
        P2(i, t) = prob_2outof4 (i , t);
    end
end
R2 = sum (P2);
mean2 = sum (R2);

% Probability that 1 out of 4 pairs working, system working
 P3 = zeros (7,count);
 for i = 1 : 1 : 7
     for t = 1 : 1 : count
         P3(i, t) = prob_1outof4 (i, t);
     end
 end
R3 = sum (P3);
mean3 = sum (R3);

figure
X = 1:1:count;
plot (X, R1, ':g', X, R2,':b', X, R3, ':r' ) 
   
     