clear;
close all;
clc;

%�������ֵ������Ϊ1�ĸ���˹����������
N = 1000;
noise = (randn(1,N) + 1i*randn(1,N))/sqrt(2);
%�������������ź�����
signal1=exp(1i*0.5*pi*(0:N-1)+1i*2*pi*rand);%������һ���ź�
signal2=exp(-1i*0.3*pi*(0:N-1)+1i*2*pi*rand);%�����ڶ����ź�
un=signal1+signal2+noise;%�������������ź�

%��������ؾ���
M = 16;  %����ؾ������
K = 2;  %Ƶ����Ŀ
xs =zeros(M,N-M);
for k=1:N-M
    xs(:,k) = un(M+k-1:-1:k).'; %������������
end
Rxx = xs(:,1:end-1)*xs(:,1:end-1)'/(N-M-1);%��������ؾ���Rxx
Rxy = xs(:,1:end-1)*xs(:,2:end)'/(N-M-1);  %���㻥��ؾ���Rxy

%��ؾ��������ֵ�ֽ⣬Ѱ����С����ֵ
 [U,E] = svd(Rxx);%����U������ʸ����ɵľ���E�ǶԽ���
                  %�Խ�Ԫ�����ɴ�С���е�����ֵ
 ev = diag(E);    %��ȡ�Խ�Ԫ���ϵ�����ֵ
 emin = ev(end);  %��ȡ��С����ֵ
 
 %��������{Cxx,Cxy}
 Z = [zeros(M-1,1),eye(M-1);0,zeros(1,M-1)];%�������Z
 Cxx = Rxx - emin*eye(M);%����Cxx
 Cxy = Rxy - emin*Z;     %����Cxy
 
 %�����{Cxx,Cxy}�Ĺ�������ֵ�ֽ�
 [U,E] = eig(Cxx,Cxy);   %��������ֵ�ֽ�
 z = diag(E);            %�ӶԽǾ���E����ȡ����ֵ
 ph = angle(z)/(2*pi);   %����������ֵ����λ��Ӧ�Ĺ�һ��Ƶ��
 err = abs(abs(z)-1);    %�뵥λԲ֮��ľ���
  [t, index2] = sort(err); 
for i=1:K
    Rootmin(i) = z(index2(2*i-1));
    f(i) =  angle(Rootmin(i))/(2*pi); 
    %����ÿ�������ȶ�Ӧͬһ��Ƶ�ʣ�����������Ƶ��
end
display('ESPRIT�㷨����ӽ���λԲ��������:');
sort(Rootmin)
display('ESPRIT�㷨����ӽ���λԲ������������λ��Ӧ�Ĺ�һ��Ƶ��:');
sort(f)

