function [flag] = srotate (W1,W2)
sw1 = size (W1);
sw2 = size (W2);
if sw1 ~= sw2
    flag = 0;
else
    Wtemp = W1;
    for i = 1 : 1 : sw1(1)-1
        Wtemp1 = Wtemp(1:end-1,:);
        Wtemp2 = Wtemp(end,:);
        Wtemp1 = circshift(Wtemp1,[-1,0]);
        Wtemp2 = circshift(Wtemp2,[0,-sw1(2)/2/(sw1(1)-1)]);
        Wtemp = [Wtemp1;Wtemp2];
        if all(W2(:) == Wtemp(:))
            flag = 1;
            break
        end       
    end
    flag = 0;
end
