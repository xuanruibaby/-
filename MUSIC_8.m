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

%自相关矩阵的特征值分解
[U,E]=svd(R);%矩阵U是特征矢量组成的矩阵，E是对角阵，
             %对角阵元素是由大到小排列的特征值
ev=diag(E);  %提取对角元素上的特征值

%利用AIC准则 和 MDL准则 估计信号源个数
AIC = zeros(1,M);
MDL = zeros(1,M);
for k = 1:M
    dec = prod(ev(k:M).^(1/(M-k+1)));%计算第一项中对数的自变量的分子
    nec = mean(ev(k:M));%计算第一项中对数的自变量的分母
    lnv = (dec/nec)^((M-k+1)*N);%计算第一项中对数的自变量
    AIC(k) = -2*log(lnv) +  2*(k-1)*(2*M-k+1);%计算AIC准则
    MDL(k) = -log(lnv) + (k-1)/2*(2*M-k+1)*log(N);%计算MDL准则
end
 [Amin1,K1] = min(AIC);%寻找使AIC准则最小的索引
 N1 = K1-1;%AIC准则信号源个数估计
 [Amin2,K2] = min(MDL);%寻找使MDL准则最小的索引
 N2 = K2-1;%MDL准则信号源个数估计

%计算MUSIC谱
En=U(:,N1+1:M);%噪声子空间的向量组成的矩阵
NF=2048;  %MUSIC的扫描点数
step = 0.001;
w = -pi:step:pi;
Pmusic = zeros(1,length(w));%为Pmusic预分配内存
for n=1:NF
    Aq=exp(-1i*2*pi*(n-1)/NF*(0:M-1)');
    Pmusic(n)=1/(Aq'*(En*En')*Aq); %Music谱
end

Pmusic = abs(Pmusic);
Pmusic_dB = 10*log10(Pmusic/max(Pmusic));
figure;
xx = w/(2*pi);
plot(xx,Pmusic_dB);
xlabel('归一化频率/f');
ylabel('归一化功率谱/dB');
