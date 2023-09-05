clear;
clc;
addpath('./myclass');
%order 
%because knot is normalized vector
k = 3;
knot = [0 0 0 0.25 0.5 0.5 0.75 1 1 1];
weights = [1 0.5 0.5 1 0.5 0.5 1];
control_points = [[0,0];[0,50];[100,50];[100,0];[100,-50];[0,-50];[0,0]];

NURBS = NURBS(knot,control_points,weights,k);
T = 0.001;


line = SLine(314.15925,38,35,30);
% hz = 1000;
% vel = line.Create_Vel_points(hz);
% [vel,time] = line.Create_Vel_points_Time(T);
% figure(1);
% plot(t,vel);
% t = NURBS.Create_Interpolator_Vlist(T,vel,0);

t = NURBS.Get_Interpolator_param_list(line,T,0);

N1 = NURBS.Create_Basis_order_functions_paramlist(1,t);
figure(1);
plot(t,N1(1,:));
title('NURBS N1');

for i=2:length(NURBS.knot)
    hold on
    plot(t,N1(i,:));
    hold off
end

N2 = NURBS.Create_Basis_order_functions_paramlist(2,t);
figure(2);
plot(t,N2(1,:));
title('NURBS N2');
 
for i=2:length(NURBS.knot)-1
    hold on
    plot(t,N2(i,:));
    hold off
end

N3 = NURBS.Create_Basis_order_functions_paramlist(3,t);
figure(3);
plot(t,N3(1,:));
title('NURBS N3');
 
for i=2:length(NURBS.knot)-2
    hold on
    plot(t,N3(i,:));
    hold off
end

% 
% 
R1 = NURBS.Create_Rational_order_param_list(1,t);
figure(4);
plot(t,R1(1,:));
title('NURBS R1');

for i=2:length(NURBS.knot)
    hold on
    plot(t,R1(i,:));
    hold off
end

R2 = NURBS.Create_Rational_order_param_list(2,t);
figure(5);
plot(t,R2(1,:));
title('NURBS R2');
 
for i=2:length(NURBS.knot)-1
    hold on
    plot(t,R2(i,:));
    hold off
end

R3 = NURBS.Create_Rational_order_param_list(3,t);
figure(6);
plot(t,R3(1,:));
title('NURBS R3');
 
for i=2:length(NURBS.knot)-2
    hold on
    plot(t,R3(i,:));
    hold off
end

[Px,Py] = NURBS.Get_Control_points_value_param_list(t);

figure(7)
plot(t,Px);
title('NURBS X(u)');

figure(8)
plot(t,Py);
title('NURBS Y(u)');

figure(9)
plot(NURBS.P(:,1),NURBS.P(:,2));
hold on;
plot(Px,Py);
title('NURBS C(u)');





