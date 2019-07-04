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

%����ؾ��������ֵ�ֽ�
[U,E]=svd(R);%����U������ʸ����ɵľ���E�ǶԽ���
             %�Խ���Ԫ�����ɴ�С���е�����ֵ
ev=diag(E);  %��ȡ�Խ�Ԫ���ϵ�����ֵ

%����AIC׼�� �� MDL׼�� �����ź�Դ����
AIC = zeros(1,M);
MDL = zeros(1,M);
for k = 1:M
    dec = prod(ev(k:M).^(1/(M-k+1)));%�����һ���ж������Ա����ķ���
    nec = mean(ev(k:M));%�����һ���ж������Ա����ķ�ĸ
    lnv = (dec/nec)^((M-k+1)*N);%�����һ���ж������Ա���
    AIC(k) = -2*log(lnv) +  2*(k-1)*(2*M-k+1);%����AIC׼��
    MDL(k) = -log(lnv) + (k-1)/2*(2*M-k+1)*log(N);%����MDL׼��
end
 [Amin1,K1] = min(AIC);%Ѱ��ʹAIC׼����С������
 N1 = K1-1;%AIC׼���ź�Դ��������
 [Amin2,K2] = min(MDL);%Ѱ��ʹMDL׼����С������
 N2 = K2-1;%MDL׼���ź�Դ��������

%����MUSIC��
En=U(:,N1+1:M);%�����ӿռ��������ɵľ���
NF=2048;  %MUSIC��ɨ�����
step = 0.001;
w = -pi:step:pi;
Pmusic = zeros(1,length(w));%ΪPmusicԤ�����ڴ�
for n=1:NF
    Aq=exp(-1i*2*pi*(n-1)/NF*(0:M-1)');
    Pmusic(n)=1/(Aq'*(En*En')*Aq); %Music��
end

Pmusic = abs(Pmusic);
Pmusic_dB = 10*log10(Pmusic/max(Pmusic));
figure;
xx = w/(2*pi);
plot(xx,Pmusic_dB);
xlabel('��һ��Ƶ��/f');
ylabel('��һ��������/dB');
