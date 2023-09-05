clear;
clc;
addpath('./myclass');
%order 
%because knot is normalized vector
k = 3;
% circle
% knot = [0 0 0 0.25 0.5 0.5 0.75 1 1 1];
% weights = [1 0.5 0.5 1 0.5 0.5 1];
% control_points = [[0,0];[0,50];[100,50];[100,0];[100,-50];[0,-50];[0,0]];
% Heart1
% weights = [0.5 0.5 0.25 0.5 0.25 0.5 0.5];
% control_points = [[50,0];[25,50];[100,50];[100,0];[100,-50];[25,-50];[50,0]];


SampleSize = 500;

NURBS = NURBS(knot,control_points,weights,k);

t = NURBS.Create_Sample_param(SampleSize);

N1 = NURBS.Create_Basis_order_funtions(1,SampleSize+1);
% figure(1);
% plot(t,N1(1,:));
% title('NURBS N1');
% 
% for i=2:length(NURBS.knot)
%     hold on
%     plot(t,N1(i,:));
%     hold off
% end
% 
% N2 = NURBS.Create_Basis_order_funtions(2,SampleSize+1);
% figure(2);
% plot(t,N2(1,:));
% title('NURBS N2');
%  
% for i=2:length(NURBS.knot)-1
%     hold on
%     plot(t,N2(i,:));
%     hold off
% end
% 
% N3 = NURBS.Create_Basis_order_funtions(3,SampleSize+1);
% figure(3);
% plot(t,N3(1,:));
% title('NURBS N3');
%  
% for i=2:length(NURBS.knot)-2
%     hold on
%     plot(t,N3(i,:));
%     hold off
% end

% N = NURBS.Create_Basis_funtions(SampleSize+1);
% 
% figure(1);
% get_order1 =  N(1,:,:);
% N1 = permute(get_order,[2 3 1]);%交換維度 2維放到 3維 3維放到 2維
% 
% for j = 1:length(NURBS.knot)-k+1
%     hold on
%     plot(x,N1(j,:));
%     hold off
% end


R1 = NURBS.Create_Rational_order_functions(1,SampleSize+1);
% figure(4);
% plot(t,R1(1,:));
% title('NURBS R1');
% 
% for i=2:length(NURBS.knot)
%     hold on
%     plot(t,R1(i,:));
%     hold off
% end
% 
% R2 = NURBS.Create_Rational_order_functions(2,SampleSize+1);
% figure(5);
% plot(t,R2(1,:));
% title('NURBS R2');
%  
% for i=2:length(NURBS.knot)-1
%     hold on
%     plot(t,R2(i,:));
%     hold off
% end
% 
% R3 = NURBS.Create_Rational_order_functions(3,SampleSize+1);
% figure(6);
% plot(t,R3(1,:));
% title('NURBS R3');
%  
% for i=2:length(NURBS.knot)-2
%     hold on
%     plot(t,R3(i,:));
%     hold off
% end



[Px,Py] = NURBS.Get_Control_points_value(SampleSize);

% figure(7)
% plot(Px);
% title('NURBS X(u)');
% 
% figure(8)
% plot(Py);
% title('NURBS Y(u)');

figure(9)
plot(NURBS.P(:,1),NURBS.P(:,2));
hold on;
plot(Px,Py);
title('NURBS C(u)');

logical = NURBS.Can_Derivatives(0);



