% k out of n pairs with i vertical planes: G system
promt1 = 'What are the values of k, n and i';
U = input (promt1);
k = U(1);
n = U(2);
i = U(3);
W0 = ones(i+1,2*n);
W = zeros(i+1,2*n,n,i,(n-k)*i);
for l = 1 : 1 : (n-k)*i % how many pairs fail
    for j = 1 : 1 : i % note that we did not consider the horizontal plane yet
        for m = 1 : 1 : n
            if l == 1
                Wtemp = W0;
                Wtemp(j, m) = 0;
                Wtemp(j, m+n) = 0;
                if 180 / i *(j-1) == 180 / n * (m-1) % make sure the horizontal plane is consistant with VPs
                    Wtemp((i+1), m) = 0;
                    Wtemp((i+1), (m+n)) = 0;
                end
                % check if balanced 
                flag = sbalance(k,n,i,Wtemp);
                if flag == 1 % balanced
                    W(:,:,m,j,l) = Wtemp;
                end
            else if l ~= 1 && sum(all (W(:, :, m, j, l-1)~=0)) >0
                Wtemp = W(:, :, m, j, l-1);
                for x = 1 : 1 : i
                    for y = 1 : 1 : n
                        if Wtemp(x, y) ~= 0 && Wtemp(x, y+n) ~= 0
                            Wtemp(x, y) = 0;
                            Wtemp(x, y+n) = 0;
                            if 180 / i *(x-1) == 180 / n * (y-1) % make sure the horizontal plane is consistant with VPs
                                Wtemp((i+1), y) = 0;
                                Wtemp((i+1), (y+n)) = 0;
                            end
                            % check if balanced 
                            flag = sbalance(k,n,i,Wtemp);
                            if flag == 1 % balanced
                                W(:,:,m,j,l) = Wtemp;
                            end
                        end
                    end
                end
                end
            end
        end
    end
end
