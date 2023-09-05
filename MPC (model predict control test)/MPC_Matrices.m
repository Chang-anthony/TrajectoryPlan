function [E , H]=MPC_Matrices(A,B,Q,R,F,N)

n=size(A,1);% A 是 n x n 矩陣, 得到 n

p=size(B,2);% B 是 n x p 矩陣, 得到 p

%%%%%%%%%%%%

M=[eye(n);zeros(N*n,n)]; % 初始化 M 矩陣. M 矩陣是 (N+1)n x n的， 

% 它上面是 n x n 個 "I", 這一步先把下半部

% 分寫成 0 

C=zeros((N+1)*n,N*p); % 初始化 C 矩阵, 這一步令它有 (N+1)n x NP 個 0

%定義 M 和 C

tmp=eye(n); %定義一個 n x n 的 I 矩陣

%　更新M 和 C
for i=1:N % 循环，i 從 1到 N

rows =i*n+(1:n); %定義當前行數，從i x n開始，共n行 

 C(rows,:)=[tmp*B,C(rows-n, 1:end-p)]; %將 C 矩陣填滿

 tmp= A*tmp; %每一次將tmp左乘一次A 

 M(rows,:)=tmp; %將 C 矩陣寫滿

end



%定義 Q_bar 和 R_bar 

Q_bar = kron(eye(N),Q); %張量積

Q_bar = blkdiag(Q_bar,F);%Block diagonal matrix 塊對角矩陣

R_bar = kron(eye(N),R);%張量積


%計算 G , E , H

G=M'*Q_bar*M; % G: n x n

E=C'*Q_bar*M; % E: NP x n

H=C'*Q_bar*C+R_bar; % NP x NP 

end