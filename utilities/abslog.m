function [data] = abslog(data)
    

    z = find(data==0);
    neg = double(data>0);
    neg(neg==0)=-1;
    data = log(abs(data));
    data = data .* neg;
    data(z) = 0;