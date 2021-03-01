clc;clear all; close all;
%%
data=load('gps_data.txt')
%x[n]=p[n]+w[n];
sigma_w=1;
%v[n+1]=u[n]+Tsa[n]
%p[n+1]=p[n]+Tsv[n]+Ts^2/2a[n]
Ts=1%s
%a[n+1]=a[n]+u[n]
%a[-1]=v[-1]=p[-1]=0 M[-1]=I

%%
%s[n]=As[n-1]+Bu[n]
%x[n]=Hs[n]+w[n]

%%
%sk[n|n-1]=Ask[n-1|n-1]
%M[n|n-1]=AM[n-1|n-1]A*+BQB*
%K[n]=M[n|n-1]H*(c+HM[n|n-1]H*)^{-1}
%sk[n|n]=sk[n|n-1]+K[n](x[n]-Hsk[n|n-1])
%M[n|n]=(I-K[n]H)M[n|n-1]

%% matrice
A=[[1, Ts,Ts^2/2];[0, 1, Ts];[0, 0, 1]];
B=[Ts^2/2,Ts,1];
H=[1,0,0];

%% vrednosti
sigma_u=0.00001;
Q=sigma_u;
C=sigma_w;
%% postavljanje pocetnih vrednosti
s=zeros(length(data),3,1);
M=zeros(length(data),3,3);
K=zeros(length(data),3,1);
s0=[0 0 0]';
M0=eye(3);


%% radjenje kalmana
preth_m=zeros(length(data),3,3);
for n=1:length(data)%+1
    if(n==1)
        s(n,:)=A*s0;
        M(n,:,:)=A*M0*A'+B*Q*B';
        pomocno=(C+H*reshape(M(n,:,:),3,[])*H');
        preth_m(n,:,:)=A*M0*A'+B*Q*B';
        K(n,:)=reshape(M(n,:,:),3,[])*H'*(pomocno)^(-1);
        s(n,:)=s(n,:)+K(n,:)*(data(n)-H*s(n,:)');
        M(n,:,:)=(eye(3)-K(n,:)'*H)*reshape(M(n,:,:),3,[]);
    else
        s(n,:)=A*s(n-1,:)';
        M(n,:,:)=A*reshape(M(n-1,:,:),3,[])*A'+B*Q*B';
        preth_m(n,:,:)=A*reshape(M(n-1,:,:),3,[])*A'+B*Q*B';
        K(n,:)=reshape(M(n,:,:),3,[])*H'*inv((C+H*M(n)*H'));
        s(n,:)=s(n,:)+K(n,:)*(data(n)-H*s(n,:)');
        M(n,:,:)=(eye(3)-K(n,:)'*H)*reshape(M(n,:,:),3,[]); 
    end
end
figure(1)
title('Estimacija s[n] za sigma_u')
subplot(3,1,1)
    hold all;
    plot(1:length(data),data);
    plot(s(:,1));
    title('polozaj za \sigma_u')
subplot(3,1,2)
    plot(s(:,2));
    title('brzina za \sigma_u')
subplot(3,1,3)
    plot(s(:,3));
    title('ubrzanje za \sigma_u')
  
saveas(figure(1),'sse_DZ2_2_1_sigma_u.png')

figure(4)
hold all;
plot(K(:,1));
plot(K(:,2));
plot(K(:,3));
title('Kalmanovo pojacanje')
legend('polozaj','brzina','ubrzanje');
saveas(figure(4),'sse_DZ2_2_1_kalman_1.png')


figure(7)
hold all;
plot(preth_m(:,1,1));
plot(preth_m(:,2,2));
plot(preth_m(:,3,3));
title('M[n|n-1] dijagnoala')
legend('(1,1)','(2,2)','(3,3)');
saveas(figure(7),'sse_DZ2_2_preth_poc1_1.png')



figure(12)
hold all;
plot(M(:,1,1));
plot(M(:,2,2));
plot(M(:,3,3));
title('M[n|n] dijagonala');
legend('(1,1)','(2,2)','(3,3)');
saveas(figure(12),'sse_DZ2_2_2_M_poc1_2.png')

%% sigma je 5* sigma_u
sigma_u_5=5*sigma_u;
Q=sigma_u_5;
preth_m=zeros(length(data),3,3);
for n=1:length(data)%+1
    if(n==1)
        s(n,:)=A*s0;
        M(n,:,:)=A*M0*A'+B*Q*B';
        pomocno=(C+H*reshape(M(n,:,:),3,[])*H');
        preth_m(n,:,:)=A*M0*A'+B*Q*B';
        K(n,:)=reshape(M(n,:,:),3,[])*H'*(pomocno)^(-1);
        s(n,:)=s(n,:)+K(n,:)*(data(n)-H*s(n,:)');
        M(n,:,:)=(eye(3)-K(n,:)'*H)*reshape(M(n,:,:),3,[]);
    else
        s(n,:)=A*s(n-1,:)';
        M(n,:,:)=A*reshape(M(n-1,:,:),3,[])*A'+B*Q*B';
        preth_m(n,:,:)=A*reshape(M(n-1,:,:),3,[])*A'+B*Q*B';
        K(n,:)=reshape(M(n,:,:),3,[])*H'*inv((C+H*M(n)*H'));
        s(n,:)=s(n,:)+K(n,:)*(data(n)-H*s(n,:)');
        M(n,:,:)=(eye(3)-K(n,:)'*H)*reshape(M(n,:,:),3,[]); 
    end
end
figure(2)
title('Estimacija s[n] za 5sigma_u')
subplot(3,1,1)
    hold all;
    plot(1:length(data),data);
    plot(s(:,1));
    title('polozaj za 5\sigma_u')
subplot(3,1,2)
    plot(s(:,2));
    title('brzina za 5\sigma_u')
subplot(3,1,3)
    plot(s(:,3));
    title('ubrzanje za 5\sigma_u')
  
saveas(figure(2),'sse_DZ2_2_1_sigma_u_5.png')

figure(5)
hold all;
plot(K(:,1));
plot(K(:,2));
plot(K(:,3));
title('Kalmanovo pojacanje')
legend('polozaj','brzina','ubrzanje');
saveas(figure(5),'sse_DZ2_2_1_kalman_2.png')

figure(8)
hold all;
plot(preth_m(:,1,1));
plot(preth_m(:,2,2));
plot(preth_m(:,3,3));
title('M[n|n-1] dijagnoala')
legend('(1,1)','(2,2)','(3,3)');
saveas(figure(8),'sse_DZ2_2_preth_poc1_2.png')

%% sigma je 1/5 sigma_u

sigma_kroz_5=(1/5)*sigma_u;
Q=sigma_kroz_5;
preth_m=zeros(length(data),3,3);
for n=1:length(data)%+1
    if(n==1)
        s(n,:)=A*s0;
        M(n,:,:)=A*M0*A'+B*Q*B';
        pomocno=(C+H*reshape(M(n,:,:),3,[])*H');
        preth_m(n,:,:)=A*M0*A'+B*Q*B';
        K(n,:)=reshape(M(n,:,:),3,[])*H'*(pomocno)^(-1);
        s(n,:)=s(n,:)+K(n,:)*(data(n)-H*s(n,:)');
        M(n,:,:)=(eye(3)-K(n,:)'*H)*reshape(M(n,:,:),3,[]);
    else
        s(n,:)=A*s(n-1,:)';
        M(n,:,:)=A*reshape(M(n-1,:,:),3,[])*A'+B*Q*B';
        preth_m(n,:,:)=A*reshape(M(n-1,:,:),3,[])*A'+B*Q*B';
        K(n,:)=reshape(M(n,:,:),3,[])*H'*inv((C+H*M(n)*H'));
        s(n,:)=s(n,:)+K(n,:)*(data(n)-H*s(n,:)');
        M(n,:,:)=(eye(3)-K(n,:)'*H)*reshape(M(n,:,:),3,[]); 
    end
end
figure(3)
title('Estimacija s[n] za 1/5sigma_u')
subplot(3,1,1)
    hold all;
    plot(1:length(data),data);
    plot(s(:,1));
    title('polozaj za 1/5\sigma_u')
subplot(3,1,2)
    plot(s(:,2));
    title('brzina za 1/5\sigma_u')
subplot(3,1,3)
    plot(s(:,3));
    title('ubrzanje za 1/5\sigma_u')
  
saveas(figure(3),'sse_DZ2_2_1_sigma_u_1_kroz_5.png')

figure(6)
hold all;
plot(K(:,1));
plot(K(:,2));
plot(K(:,3));
title('Kalmanovo pojacanje')
legend('polozaj','brzina','ubrzanje');
saveas(figure(6),'sse_DZ2_2_1_kalman_3.png')





figure(9)
hold all;
plot(preth_m(:,1,1));
plot(preth_m(:,2,2));
plot(preth_m(:,3,3));
title('M[n|n-1] dijagnoala')
legend('(1,1)','(2,2)','(3,3)');
saveas(figure(9),'sse_DZ2_2_preth_poc1_3.png')

%% promena pocetnih uslova
s=zeros(length(data),3,1);
M=zeros(length(data),3,3);
K=zeros(length(data),3,1);
s0=[5 5 5]';
M0=10.*eye(3);
Q=sigma_u;
preth_m=zeros(length(data),3,3);

for n=1:length(data)
    if(n==1)
        s(n,:)=A*s0;
        M(n,:,:)=A*M0*A'+B*Q*B';
        pomocno=(C+H*reshape(M(n,:,:),3,[])*H');
        preth_m(n,:,:)=A*M0*A'+B*Q*B';
        K(n,:)=reshape(M(n,:,:),3,[])*H'*(pomocno)^(-1);
        s(n,:)=s(n,:)+K(n,:)*(data(n)-H*s(n,:)');
        M(n,:,:)=(eye(3)-K(n,:)'*H)*reshape(M(n,:,:),3,[]);
    else
        s(n,:)=A*s(n-1,:)';
        M(n,:,:)=A*reshape(M(n-1,:,:),3,[])*A'+B*Q*B';
        preth_m(n,:,:)=A*reshape(M(n-1,:,:),3,[])*A'+B*Q*B';
        K(n,:)=reshape(M(n,:,:),3,[])*H'*inv((C+H*M(n)*H'));
        s(n,:)=s(n,:)+K(n,:)*(data(n)-H*s(n,:)');
        M(n,:,:)=(eye(3)-K(n,:)'*H)*reshape(M(n,:,:),3,[]); 
    end
end
close all;
figure(10)
hold all;
plot(K(:,1));
plot(K(:,2));
plot(K(:,3));
title('Kalmanovo pojacanje')
legend('polozaj','brzina','ubrzanje');
saveas(figure(10),'sse_DZ2_2_2_kalman_poc2.png')





figure(11)
hold all;
plot(preth_m(:,1,1));
plot(preth_m(:,2,2));
plot(preth_m(:,3,3));
title('M[n|n-1] dijagnoala')
legend('(1,1)','(2,2)','(3,3)');
saveas(figure(11),'sse_DZ2_2_2_preth_poc2_2.png')

figure(13)
hold all;
plot(M(:,1,1));
plot(M(:,2,2));
plot(M(:,3,3));
title('M[n|n] dijagonala');
legend('(1,1)','(2,2)','(3,3)');
saveas(figure(13),'sse_DZ2_2_2_M_poc2_2.png')











