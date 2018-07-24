%% 1����ʼ��
clear;
n=4;
N=2^n;%ȷ���볤
Nr=N;
Np=N-Nr; 
K=8;%ȷ����Ϣλ�ĳ���
R=K/Nr;%ȷ������ 
k=20; %������֡��
maxii=1000000; %�����֡��
minii=1000; %��С����֡��
j=7;%ȷ����ѭ������
minebn0=0;
maxebn0=3;
be=zeros(1,j);
ber=zeros(1,j);
fe=zeros(1,j);
fer=zeros(1,j);

ebn0=linspace(minebn0,maxebn0,j);
SNR=ebn0+10*log10(2*R);

%% 2��������ѭ��
for i=1:j %��ͬ�������
    sigma = sqrt(1/10^(SNR(i)/10));
    ii=0;
    
%2.1 ȷ����Ϣλ����˹���ƣ�
[free_sets,frozen_sets,punctured_sets]=GAWRX_awgn_f(n,Nr,K,3);
    
while ((fe(i) <k) && (ii < maxii))||(ii<minii)    %֡ѭ������
    
%2.2 ������������
% u=randsrc(K,1,[0 1]);
u = ones(K,1);
uu=zeros(N,1);
uu(free_sets)=u;

%2.3 �������
x = encode(uu);
x = rvsl(x); %��֯
y0 = awgn(1-2*x, SNR(i));

%2.4 �������
%y0 = round(y0*10000)/10000;
p1=exp(-(y0-1).^2/(2*sigma^2));
pn1=exp(-(y0+1).^2/(2*sigma^2));
llr=log(p1./pn1);
llr(punctured_sets) = 16;
v =  SC_decoding_new_v7(llr,frozen_sets'-1,N,n,N-K);

%2.5 ����������/��֡��
diff=N-sum(v==uu);
be(i)=be(i)+diff;
ii=ii+1;
[ii be(i)]
if diff~=0
    fe(i)=fe(i)+1;
end
end
ber(i)=be(i)./(K*ii);
fer(i)=fe(i)/(ii);

semilogy(ebn0,ber,'b',ebn0,fer,'g');
grid on;
axis([0 maxebn0 0.0001 1]);
xlabel('Eb/N0(dB)');
ylabel('Bit/Frame Error Rate');
drawnow
end