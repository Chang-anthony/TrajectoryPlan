function N = generalOrderBSplineFuntions(t,x,c,SampleSize,N1,k,periodic)
m1 = length(t);
m2 = length(x);
N = zeros(m1-c,m2);
N = N1;

    for j = 2:c
        d = j;
        if(periodic == 1)
            k1 = 1;
            k2 = m1-d;
        else
            k1 = k-(d-1);
            k2 = m1 -k;
        end

        i = k1;
        j1 = (i-1)*SampleSize+1;
        j2 = j1+d*SampleSize-1;

        if (periodic == 1)
            N(i,j1:j2) = (x(j1:j2)-t(i))./(t(i+d-1)-t(i)).*N(i,j1:j2)+(t(i+d)-x(j1:j2))./(t(i+d)-t(i+1)).*N(i+1,j1:j2);
        else
            N(i,j1:j2) = (t(i+d)-x(j1:j2))./(t(i+d)-t(i+1)).*N(i+1,j1:j2);
        end

        for i=k1+1:k2-1
            j1 = (i-1)*SampleSize+1;
            j2 = j1+d*SampleSize-1;
            N(i,j1:j2) = (x(j1:j2)-t(i))./(t(i+d-1)-t(i)).*N(i,j1:j2)+(t(i+d)-x(j1:j2))./(t(i+d)-t(i+1)).*N(i+1,j1:j2);
        end
        
        i = k2;
        j1 = (i-1)*SampleSize+1;
        j2 = j1+d*SampleSize-1;
        if (periodic == 1)
            N(i,j1:j2) = (x(j1:j2)-t(i))./(t(i+d-1)-t(i)).*N(i,j1:j2)+(t(i+d)-x(j1:j2))./(t(i+d)-t(i+1)).*N(i+1,j1:j2);
        else
            N(i,j1:j2) = (x(j1:j2)-t(i))./(t(i+d-1)-t(i)).*N(i,j1:j2);
        end

    end

end