clc;
clear;

Pmin = 0;
Pmax = 10; 
NInstance = 50; 
sigma = 2; %noise

%generator simple test data
X = [ones(NInstance,1) (Pmax-Pmin) * rand(NInstance,1)];
sort(X,1);
weight = rand(2,1);

Y = X * weight + rand(NInstance,1)*sigma;


%%fit use LSF method and get fit line predict yd
weight = inv(X'*X) * X' * Y;
yd = X * weight;

hold on
title('LSF Simple test');

%plot simple data
scatter(X(:,2),Y,"blue",'filled');
xlim([Pmin Pmax]);
ylim([min(Y)-sigma max(Y)+sigma]);

plotx = [min(X(:,2)) max(X(:,2))];
ploty = [min(yd) max(yd)];

%plot fit line
plot(plotx,ploty,'r');

hold off




