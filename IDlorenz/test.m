clc;close all; clear;
load data_S5_0318.txt
data = load('data_S5_0318.txt');
friction = data(:,5)' ;
friction = lowp(friction,10,100,0.1,20,1000);
poscmd = data(:,3)';
velcmd = data(:,4)' ;
ts = 1e-3;
tt = length(friction)*ts;
t = linspace(ts,tt,tt/ts) ;

velcmd1 = velcmd;
friction = timeseries(friction,t);
poscmd = timeseries(poscmd,t);
velcmd = timeseries(velcmd,t);
Z0=0;Z1=0;tiltaZ0=0;tiltaZ1=0;
%k = 34.8 ;
k = 200;

%-------CSMC--------%
landa = 800; % no sig
k = 480;
% landa = 1200;
% k = 480;

%---------DFNTSMC-----%
% DFN_k1 = 0.8;   
% DFN_k2 = 1200;
% DFN_k3 = 190;
% DFN_k4 = 200;
% DFN_gamma = 1.8;
% DFN_landa = 250;
DFN_k1 = 0.8;     % no sig
DFN_k2 = 1200;
DFN_k3 = 190;
DFN_k4 = 200;
DFN_gamma = 1.8;
DFN_landa = 1000;
%%

lamda = 220;
KD = 2*lamda;
Kp = lamda^2;
M = 0.00018;



% tt = length(input)*ts ;
% t = linspace(ts,tt,tt/ts) ;
% sim('XXXX.slx');
sigma0=1.201; 
sigma1=0.1219;
sigma2=0.9758e-4;
st = 0.001
preZ0 = 0;
vs=2;
for i=1:length(velcmd1)
gv=(0.0293+(0.0305-0.0293)*exp(-power((velcmd1(i))/vs,2)));
dzdt=(velcmd1(i)-sigma0*abs(velcmd1(i))/gv*Z0);
Z0 = preZ0 + dzdt*st;
preZ0 = Z0;
F(i) = sigma0*Z0+sigma1*dzdt+sigma2*velcmd1(i);
end
% figure(1)
% plot(velcmd1, F)
% grid on

dfh = 0;fh =0;
sigmaF = 0;
%-----------
k1 = [5.741,2.1658,1.1768, 0.9294];
a = [0.299, 0.2813, 0.207, 0.211] ;
    for i=1:4
        if fh < a(i)*gv
            dfh = k1(i)*1;
            fh= fh + dfh*0.001; 
            sigmaF = sigmaF + fh;
        end
    end
function y = lowp(x,f1,f3,rp,rs,Fs)
   

    wp = 2*pi*f1/Fs;
    ws = 2*pi*f3/Fs;
    % ??切比雪夫濾波器
    [n,wn] = cheb1ord(wp/pi,ws/pi,rp,rs);
    [bz1,az1] = cheby1(n,rp,wp/pi);

    y = filter(bz1,az1,x);    % ?序列 x 濾波後得到的序列 y
end