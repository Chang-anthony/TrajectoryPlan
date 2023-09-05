function [t,Range] = UniformKnotVector(k,n,periodic,normalized)
    % k is the order of polynomial ,i.e. k=degree+1 , and
    % n+1 is the number of control points.
    t = zeros(1,n+k+1); % the length of knot vector is n+k+1 
    if (periodic==1)
        Range=n+k;
        t = 0:Range;
    else
        Range=n-k+2;
        for i=k+1:n+1
            t(i) = i-k;
        end
        t(n+2:n+k+1) = Range;
    end
    if (normalized == 1)
        t = t./Range; % the parameter range is 0<= t <= n+k
    end  
end