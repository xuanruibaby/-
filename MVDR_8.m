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
M=8;  %����ؾ������
xs =zeros(M,N-M);
for k=1:N-M
    xs(:,k) = un(M+k-1:-1:k).'; %������������
end
R=xs*xs'/(N-M);   %���̲�ʽ��3.5.18����������ؾ���

%����MVDR��
NF=2048;  %MVDR��ɨ�����
step = 0.001;
w = -pi:step:pi;
Pmvdr = zeros(1,length(w));
for n=1:NF
    Aq=exp(-1i*2*pi*(n-1)/NF*(0:M-1)');
    Pmvdr(n)=1/(Aq'*inv(R)*Aq); %MVDR��
end

Pmvdr = abs(Pmvdr);
Pmvdr_dB = 10*log10(Pmvdr/max(Pmvdr));
figure;
xx = w/(2*pi);
plot(xx,Pmvdr_dB);
xlabel('��һ��Ƶ��/f');
ylabel('��һ��������/dB');






