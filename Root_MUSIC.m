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
M = 16; %自相关矩阵阶数
K = 2; %信号源数目        
xs =zeros(M,N-M);
for k=1:N-M
    xs(:,k) = un(M+k-1:-1:k).'; %构造样本矩阵
end
R=xs*xs'/(N-M);   %按教材式（3.5.18）计算自相关矩阵

%自相关矩阵的特征值分解
[U,E]=svd(R);%矩阵U是特征矢量组成的矩阵，E是对角阵，
             %对角阵元素是由大到小排列的特征值
ev=diag(E);  %提取对角元素上的特征值

% Root-MUSIC 算法的实现
 G = U(:,3:M);         %噪声子空间的向量组成的矩阵
 Gr = G*G';
 co = zeros(2*M-1,1);  %初始化教材式（3.6.38)决定的2（M-1）次方程的系数
 for m=1:M
     co(m:m+M-1) = co(m:m+M-1)+Gr(M:-1:1,m);
 end                   %计算教材式（3.6.38）左边的多项式系数
 z = roots(co);        %多项式求根
 ph = angle(z)/(2*pi); %求所有跟的相位对应的归一化频率
 err = abs(abs(z)-1);  %求2（M-1）个根与单位圆之间的距离
 
 [t, index1] = sort(err); 
for i=1:K
    Rootmin(i) = z(index1(2*i-1));
    f(i) =  angle(Rootmin(i))/(2*pi); 
    %相邻每两个幅度对应同一个频率，间隔输出所有频率
end
display('Root-MUSIC算法中最接近单位圆的两个根:');
sort(Rootmin)
display('Root-MUSIC算法中最接近单位圆的两个根的相位对应的归一化频率:');
sort(f)

 
