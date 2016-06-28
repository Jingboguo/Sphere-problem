function [W2] = transition (W1)
sw1 = size (W1);
i = sw1(1)-1;
n = sw1(2)/2;
W2 = zeros(sw1);
for l = 1 : 1 : i
    for m = 1 : 1 : n/2
        Wtemp = W1;
        if Wtemp(l, m) ~= 0
            Wtemp(l,m) = 0;
            Wtemp(l,m+n) = 0;
            if 180 / i *(l-1) == 180 / n * (m-1) % make sure the horizontal plane is consistant with VPs
                Wtemp((i+1), m) = 0;
                Wtemp((i+1), (m+n)) = 0;
            end
            W2 = [W2; Wtemp];
        end
    end
end
