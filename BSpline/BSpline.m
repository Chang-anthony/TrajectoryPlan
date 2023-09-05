%B-Spline Knot Vector
clc;
clear;
n = 3;
k = 3;
c = 4; %c is the order of basis function we want to create
normalized=0;
periodic=1;
Bx = [0 3 6 9];
By = [0 9 3 6];

[t,Range] = UniformKnotVector(k,n,periodic,normalized);

SampleSzie = 500;
m = n+k+1;
% Sample size is number of division between the knot values
x = zeros(1,(n+k)*SampleSzie+1);
if (periodic == 1)
    x = 0:1/SampleSzie:Range;
else
    x((k-1)*SampleSzie+1:((n+1)*SampleSzie+1)) = 0:1/SampleSzie:Range;
    x(((n+1)*SampleSzie+2):end) = Range;
end
if (normalized == 1)
    x = x./Range;
end

N1 = FirstOrderBSplineFunctions(k,t,x,SampleSzie,periodic);

%plot(x,N1(1,:));
%for i=2:m-1
%   hold on
%    plot(x,N1(i,:));
%    hold off
%end

N2 = SecondOrderBSplineFunctions(k,t,x,SampleSzie,periodic,N1);

% plot(x,N2(1,:));
% for i=2:m-2
%    hold on
%     plot(x,N2(i,:));
%     hold off
% end

N3 = ThirdOrderBSplineFunctions(k,t,x,SampleSzie,periodic,N1,N2);

% plot(x,N3(1,:));
% for i=2:m-3
%    hold on
%     plot(x,N3(i,:));
%     hold off
% end


N = generalOrderBSplineFuntions(t,x,c,SampleSzie,N1,k,periodic);

% plot(x,N(1,:));
% for i=2:m-c
%    hold on
%     plot(x,N(i,:));
%     hold off
% end

Px = Bx*N(1:n+1,(k-1)*SampleSzie+1:(n+1)*SampleSzie);
Py = By*N(1:n+1,(k-1)*SampleSzie+1:(n+1)*SampleSzie);


plot(Bx,By);
hold on
plot(Px,Py);