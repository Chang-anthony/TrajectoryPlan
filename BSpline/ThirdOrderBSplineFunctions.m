function N3 = ThirdOrderBSplineFunctions(k,t,x,SampleSize,periodic,N1,N2)

m1 = length(t);
m2 = length(x);
N3 = zeros(m1-3,m2);
d = 3;
    if (periodic == 1)
        k1 = 1;
        k2 = m1-3;
    else
        k1 = k-2;
        k2 = m1-k;
    end
    
i = k1;
j1 = (i-1)*SampleSize+1;
j2 = j1+3*SampleSize-1;

    if(periodic == 1)
        N3(i,j1:j2) = (x(j1:j2)-t(i))./(t(i+d-1)-t(i)).*N2(i,j1:j2)+(t(i+d)-x(j1:j2))./(t(i+d)-t(i+1)).*N2(i+1,j1:j2);
    else
        N3(i,j1:j2) = (t(i+d)-x(j1:j2))./(t(i+d)-t(i+1)).*N2(i+1,j1:j2);
    end

    for i = k1+1:k2-1
        j1 = (i-1)*SampleSize+1;
        j2 = j1+d*SampleSize-1;
        N3(i,j1:j2) = (x(j1:j2)-t(i))./(t(i+d-1)-t(i)).*N2(i,j1:j2)+(t(i+d)-x(j1:j2))./(t(i+d)-t(i+1)).*N2(i+1,j1:j2);
    end

i = k2;
j1 = (i-1)*SampleSize+1;
j2 = j1+d*SampleSize-1;
    if (periodic == 1)
        N3(i,j1:j2) = (x(j1:j2)-t(i))./(t(i+d-1)-t(i)).*N2(i,j1:j2)+(t(i+d)-x(j1:j2))./(t(i+d)-t(i+1)).*N2(i+1,j1:j2);
    else
        N3(i,j1:j2) = (x(j1:j2)-t(i))./(t(i+d-1)-t(i)).*N2(i,j1:j2);
    end
end