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
M = 16;  %自相关矩阵阶数
K = 2;  %频率数目
xs =zeros(M,N-M);
for k=1:N-M
    xs(:,k) = un(M+k-1:-1:k).'; %构造样本矩阵
end
Rxx = xs(:,1:end-1)*xs(:,1:end-1)'/(N-M-1);%计算自相关矩阵Rxx
Rxy = xs(:,1:end-1)*xs(:,2:end)'/(N-M-1);  %计算互相关矩阵Rxy

%相关矩阵的特征值分解，寻找最小特征值
 [U,E] = svd(Rxx);%矩阵U是特征矢量组成的矩阵，E是对角阵，
                  %对角元素是由大到小排列的特征值
 ev = diag(E);    %提取对角元素上的特征值
 emin = ev(end);  %获取最小特征值
 
 %构造矩阵对{Cxx,Cxy}
 Z = [zeros(M-1,1),eye(M-1);0,zeros(1,M-1)];%构造矩阵Z
 Cxx = Rxx - emin*eye(M);%计算Cxx
 Cxy = Rxy - emin*Z;     %计算Cxy
 
 %矩阵队{Cxx,Cxy}的广义特征值分解
 [U,E] = eig(Cxx,Cxy);   %广义特征值分解
 z = diag(E);            %从对角矩阵E中提取特征值
 ph = angle(z)/(2*pi);   %求所有特征值的相位对应的归一化频率
 err = abs(abs(z)-1);    %与单位圆之间的距离
  [t, index2] = sort(err); 
for i=1:K
    Rootmin(i) = z(index2(2*i-1));
    f(i) =  angle(Rootmin(i))/(2*pi); 
    %相邻每两个幅度对应同一个频率，间隔输出所有频率
end
display('ESPRIT算法中最接近单位圆的两个根:');
sort(Rootmin)
display('ESPRIT算法中最接近单位圆的两个根的相位对应的归一化频率:');
sort(f)

