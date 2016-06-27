mu = 50;
sigma = 2;
theta = 8;
F = @(t) wblcdf(t, mu, sigma);
G = @(t) expcdf(t, mu);
M = @(t) normcdf(t, mu, theta);
count = 100;
P1 = zeros (1, count);
P2 = zeros (1, count);
P3 = zeros (1, count);
for t = 1 : 1 : count
    P1(1, t) = 1- F (t);
    P2(1, t) = 1- G (t);
    P3(1, t) = 1- M (t);
end
figure
X = 1:1:count;
plot (X, P1, ':r', X, P2,':b', X, P3, ':g') 

meanexp = sum (P2);
meanwbl = sum (P1);
meannorm = sum (P3);