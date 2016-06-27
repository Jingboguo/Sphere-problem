    mu = 900;
    sigma = 300; 
    M = @(t) normcdf(t, mu, sigma);
    N = ones (1000, 1);
    for t = 1: 1 :1000
        N (t) = 1 - M (t);
    end
    figure
    X = 1:1:1000;
    plot (X, N)