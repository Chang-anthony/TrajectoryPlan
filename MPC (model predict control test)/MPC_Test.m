clear;
close all;
clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%

%%第一步，定義狀態空間矩陣

%%定義狀態矩陣 A , n x n 矩陣

A = [1 0.1; -1 2];

n= size (A,1);

%%定義輸入矩陣 B , n x p 矩陣

B = [ 0.2 1; 0.5 2];

p = size(B,2);

%%定義Q矩陣 , n x n 矩陣

Q=[100 0;0 1];

%%定義F矩陣 , n x n 矩陣

F=[100 0;0 1];

%%定義R矩陣, p x p 矩陣

R=[1 0 ;0 .1];

%%定義 step 數量 k

k_steps=100;

%%定義矩陣 X_k , n x k 矩陣

X_K = zeros(n,k_steps);

%%初始化狀態變量值 , n x 1 向量

X_K(:,1) =[20;-20];

%%定義輸入矩陣 U_K , p x k 矩陣
U_K=zeros(p,k_steps);

%%定義預測區間N

N=5;

%%Call MPC_Matrices 函數 求得 E,H 矩陣

[E,H]=MPC_Matrices(A,B,Q,R,F,N);


%%計算每步狀態變量的值

for k = 1:k_steps

%%求得U_K(:,k)

U_K(:,k) = Prediction(X_K(:,k),E,H,N,p);

%%計算K+1步時狀態變量的值

X_K(:,k+1)=(A*X_K(:,k)+B*U_K(:,k));

end



%%繪製狀態變量與輸入的變化

subplot(2, 1, 1);

hold;

for i =1 :size (X_K,1)

plot (X_K(i,:));

end

legend("x1","x2")

hold off;

subplot (2, 1, 2);

hold;

for i =1 : size (U_K,1)

plot (U_K(i,:));

end

legend("u1","u2")