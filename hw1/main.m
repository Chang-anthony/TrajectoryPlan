clc;
clear;
addpath('./myclass');

table = [[0,0,0,0,0];
        [5,5,0.7,6,0.8];
        [15,8,0,8,0];
        [30,9,0.25,3,-0.25];
        [40,10,0,2.5,0]];

hz = 1000; 
path = Path();
path.create_lines_table(table);

disp('R3 not normalization = ');
disp(path.lines(3).R);
disp('B3 not normalization = ');
disp(path.lines(3).B);
disp('B3_inv not normalization = ');
disp(inv(path.lines(3).B));
disp('A3 not normalization = ');
disp(path.lines(3).A);


[xt,yt] = path.get_time_point_value(35,4);
disp('xt yt on t = 35 =>');
%disp(xt);
%disp(yt);

%[xt,yt,t]= path.get_time_point_value_normalization(0.5,2);
%disp('on line2 and u = 0.5 => xt yt t = ');
disp(xt);
disp(yt);
%disp(t);

disp('R3 normalization = ');
disp(path.lines(3).R_n);
disp('B3  normalization = ');
disp(path.lines(3).B_n);
disp('B3_inv  normalization = ');
disp(inv(path.lines(3).B_n));
disp('A3  normalization')
disp(path.lines(3).A_n);

%[outx,outy,outt] = path.create_path_lines_point(hz);
[outx,outy,outt] = path.create_path_lines_point_normalization(hz);



figure(1);
plot(outt,outx,'R');
xlabel('T-axis');
ylabel('X-axis');
axis([0,50,0,15]);
title('X-T cubic line');


figure(2);
plot(outt,outy,'G');
xlabel('T-axis');
ylabel('Y-axis');
axis([0,50,0,15]);
title('Y-T cubic line');

figure(3);
plot(outx,outy,'B');
xlabel('X-axis');
ylabel('Y-axis');
axis([0,15,0,15]);
title('X-Y cubic line');

