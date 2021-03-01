clear all;close all;clc;
%konstante
%%

sigma_w=0.3;
fi=1;
N=100;
A=1;
f_0=0.2;



%% pravljenje merenja
x=zeros(1,N);
p_x_n_tetha=zeros(1,N);
K1=1/sqrt(2*pi*sigma_w^2);
for n=0:N-1
    x(n+1)=A*cos(2*pi*f_0*n+fi)+normrnd(0,sigma_w);
end
phi_arr=-5:0.1:5;
l_out=zeros(1,length(phi_arr));
l_prim_out=zeros(1,length(phi_arr));
l_sec_out=zeros(1,length(phi_arr));
%izracunavanje po fi izlaza
for i=1:length(l_out)
    l_out(i)=log_l(x,phi_arr(i));
    l_prim_out(i)=l_prim(x,phi_arr(i));
    l_sec_out(i)=l_sec(x,phi_arr(i));
end

%% plotovanje
figure(1)
plot(phi_arr,l_out);
xlabel('fi'),ylabel('l(fi)'),title('funkcija l(fi)');
grid on;
saveas(figure(1),'sse_ZAD_1_l_fi.png');

figure(2)
plot(phi_arr,l_prim_out);
xlabel('fi'),ylabel('l_{prim}(fi)'),title('prvi izvod funkcije l(fi)');
grid on;
saveas(figure(2),'sse_ZAD_1_l_fi_prim.png');

figure(3)
plot(phi_arr,l_sec_out);
xlabel('fi'),ylabel('l_{sec}(fi)'),title('drugi izvod funkcije l(fi)');
grid on;
saveas(figure(3),'sse_ZAD_1_l_fi_sec.png');


%% 2. vrednost procene fi_k sa kapom
pocetne=[1.4 0.7 1.6];
%formiranje signala
signali=zeros(5,N);
sigma_w_niz=[0.1 0.2 0.25 0.3 0.35];
n=1:N;
for a=1:5
    signali(a,:)=A*cos(2*pi*f_0*n+1)+randn(1,N)*sigma_w;
end
Nk=15;
kplot=1:Nk;
out_2=zeros(3,5,Nk);
for poc=1:3
    for sig=1:5 
        trenutni=signali(sig,:);
        out_2(poc,sig,:)=Newton_Raphson(trenutni,pocetne(poc),Nk);
    end
end
figure(4)
hold all;
%out_2(2,1,:)
plot(kplot,reshape(out_2(1,1,:),Nk,[]));
plot(kplot,reshape(out_2(1,2,:),Nk,[]));
plot(kplot,reshape(out_2(1,3,:),Nk,[]));
plot(kplot,reshape(out_2(1,4,:),Nk,[]));
plot(kplot,reshape(out_2(1,5,:),Nk,[]));
axis([1 15 0 2]);
saveas(figure(4),'sse_ZAD_1_2_razl1.png')


figure(5)
hold all;
%out_2(2,1,:)
plot(kplot,reshape(out_2(2,1,:),Nk,[]));
plot(kplot,reshape(out_2(2,2,:),Nk,[]));
plot(kplot,reshape(out_2(2,3,:),Nk,[]));
plot(kplot,reshape(out_2(2,4,:),Nk,[]));
plot(kplot,reshape(out_2(2,5,:),Nk,[]));
axis([1 15 0 2])
saveas(figure(5),'sse_ZAD_1_2_razl2.png')

figure(6)
hold all;
%out_2(2,1,:)
plot(kplot,reshape(out_2(3,1,:),Nk,[]));
plot(kplot,reshape(out_2(3,2,:),Nk,[]));
plot(kplot,reshape(out_2(3,3,:),Nk,[]));
plot(kplot,reshape(out_2(3,4,:),Nk,[]));
plot(kplot,reshape(out_2(3,5,:),Nk,[]));
axis([1 15 0 2])
saveas(figure(6),'sse_ZAD_1_2_razl3.png')

%%
%zadatak pod 3
Nr=1000;
fi_pocetno=0.5;
sigma_w_1=0.1
No=10;
%fisherova informacija
fisher=No*A^2/(2*sigma_w_1^2);
sigma_est=sqrt(1/fisher);
%formiranje signala
signali=zeros(Nr,No);
n=1:No;
close all;
for a=1:Nr
    signali(a,:)=A*cos(2*pi*f_0*n+1)+randn(1,No)*sigma_w_1;
end
%sa prethodnih grafika
Nk=5;
kplot=1:Nk;
out_3=zeros(Nr,Nk);

for sig=1:Nr 
    trenutni=signali(sig,:);
    out_3(sig,:)=Newton_Raphson(trenutni,fi_pocetno,Nk);
end
figure(6)
hold all;
%out_2(2,1,:)
for i=1:Nr
plot(1:Nk,reshape(out_3(i,:),Nk,[]));
end

space=0.8:0.01:1.2;

fgv=normpdf(space,1,sigma_est);
figure(7)
hold all;
histogram(out_3(:,Nk),15,'Normalization','pdf')
plot(space,fgv)
saveas(figure(7),'sse_DZ2_1_3_sigma_01.png');
%%

Nr=1000;
fi_pocetno=0.5;
sigma_w_2=1;
No=100;
%fisherova informacija
fisher=No*A^2/(2*sigma_w_2^2);
sigma_est=sqrt(1/fisher);
%formiranje signala
signali=zeros(Nr,No);
n=1:No;
close all;
for a=1:Nr
    signali(a,:)=A*cos(2*pi*f_0*n+1)+randn(1,No)*sigma_w_2;
end
%sa prethodnih grafika
Nk=5;
kplot=1:Nk;
out_3=zeros(Nr,Nk);

for sig=1:Nr 
    trenutni=signali(sig,:);
    out_3(sig,:)=Newton_Raphson(trenutni,fi_pocetno,Nk);
end
figure(6)
hold all;
%out_2(2,1,:)
for i=1:Nr
plot(1:Nk,reshape(out_3(i,:),Nk,[]));
end

space=0.5:0.01:1.6;

fgv=normpdf(space,1,sigma_est);
figure(8)
hold all;
histogram(out_3(:,Nk),15,'Normalization','pdf')
plot(space,fgv)
saveas(figure(8),'sse_DZ2_1_3_sigma_1.png');




%% definisanje log izraza l(fi)
function l_fi=log_l(x,fi_1)
    sigma_w=0.3;
    fi=1;
    A=1;
    f_0=0.2;
    sum_l_fi=0;    
    for n=1:length(x)
        sum_l_fi=sum_l_fi+(x(n)-A*cos(2*pi*f_0*(n)+fi_1))^2;
    end
    l_fi=-(length(x)/2)*log(2*pi*sigma_w^2)-(1/(2*sigma_w^2))*sum_l_fi; 
    
end

%% definisanje prvog izvda po fi l'(fi)
function l_prim_fi=l_prim(x,fi_1)
    sigma_w=0.3;
    A=1;
    f_0=0.2;
    sum_l_prvi_fi=0;
    for n=1:length(x)
        sum_l_prvi_fi=sum_l_prvi_fi+(x(n)*sin(2*pi*f_0*(n)+fi_1)-(A/2)*sin(4*pi*f_0*(n)+2*fi_1));
    end
    l_prim_fi=-(A)*sum_l_prvi_fi/(sigma_w^2);
end

%% definisanje drugog izvoda po fi l''(fi)
function l_sec_fi=l_sec(x,fi_1)
    sigma_w=0.3;
    A=1;
    f_0=0.2;
    sum_l_drugi_fi=0;
    for n=1:length(x)
        sum_l_drugi_fi=sum_l_drugi_fi+(1/2+(1/2)*x(n)*cos(4*pi*f_0*(n)+2*fi_1)-(cos(4*pi*f_0*(n)+2*fi_1)));
    end
    l_sec_fi=-(A)*sum_l_drugi_fi/(sigma_w^2);

    
end
%% NEwton-Raphson
function nr=Newton_Raphson(x,fi_poc,k)
    estimacija=zeros(1,k);
    estimacija(1)=fi_poc;
    for i=1:k-1
        step= l_prim(x,estimacija(i))./l_sec(x,estimacija(i));
        estimacija(i+1)=estimacija(i)-step;
    end
    estimacija
    nr=estimacija;
end



    
    
    
    