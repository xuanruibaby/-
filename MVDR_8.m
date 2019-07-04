clear;
close all;
clc;

%产生零均值、方差为1的复高斯白噪声序列
N = 1000;
noise = (randn(1,N) + 1i*randn(1,N))/sqrt(2);
%产生带噪声的信号样本
signal1=exp(1i*0.5*pi*(0:N-1)+1i*2*pi*rand);%产生第一个信号
signal2=exp(-1i*0.3*pi*(0:N-1)+1i*2*pi*rand);%产生第二个信号
un=signal1+signal2+noise;%产生带噪声的信号

%计算自相关矩阵
M=8;  %自相关矩阵阶数
xs =zeros(M,N-M);
for k=1:N-M
    xs(:,k) = un(M+k-1:-1:k).'; %构造样本矩阵
end
R=xs*xs'/(N-M);   %按教材式（3.5.18）计算自相关矩阵

%计算MVDR谱
NF=2048;  %MVDR的扫描点数
step = 0.001;
w = -pi:step:pi;
Pmvdr = zeros(1,length(w));
for n=1:NF
    Aq=exp(-1i*2*pi*(n-1)/NF*(0:M-1)');
    Pmvdr(n)=1/(Aq'*inv(R)*Aq); %MVDR谱
end

Pmvdr = abs(Pmvdr);
Pmvdr_dB = 10*log10(Pmvdr/max(Pmvdr));
figure;
xx = w/(2*pi);
plot(xx,Pmvdr_dB);
xlabel('归一化频率/f');
ylabel('归一化功率谱/dB');






