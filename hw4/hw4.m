clc;
clear;
addpath('./myclass');

line = SLine(314.15925,38,35,30);

hz = 1000; 
[acc,t] = line.Create_Acc_points(hz);
vel = line.Create_Vel_points(hz);
pos = line.Create_Pos_points(hz);

figure(1);
plot(t,acc,'R');
xlabel('T-axis');
ylabel('Acc-axis');
axis([0,10,-60,60]);
title('Acc-T  Sline');


figure(2);
plot(t,vel,'G');
xlabel('T-axis');
ylabel('Vel-axis');
axis([0,10,0,60]);
title('Vel-T  Sline');

figure(3);
plot(t,pos,'B');
xlabel('T-axis');
ylabel('Pos-axis');
%axis([0,10,0,400]);
title('Pos-T  Sline');



