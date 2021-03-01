%Vrednosti parametara na osnovu formule i obrasca, Tijana Aleksic 2018/0455


P=2;Q=1;R=2;S=1;


%Zadatak 1  


%Tekst zadatka 
%Eksperiment izvlacenja papirica iz kutije. Na papiricima mogu da se nalaze
%brojevi 1 2 3 4 5 6. Verovatnoca da se dobije broj 2 je duplo veca nego
%verovatnoca da se dobije bilo koji od ostalih brojeva, a verovatnoca da je
%kutija prazna je 0.02



%izracunajmo prvo verovatnocu da je izvucen jedan papiric, ona je q
q=1-0.02;


%znamo da se dobije 2 je duplo veca sansa od ostalih papirica sto nam
%izjednacava vervoatnocu kao da dodamo jos jedan papiric sa dvojkom na sebi
%Od q uniformno su rasporedjeni 1 2 2 3 4 5 6  sto ih je 7 pa je vrv za
%jedan od papirica
%vrv1 je verovatnoca bilo kog osim 2 papirica a vrv2 je za dvojku
vrv1=q/7;
vrv2=2*vrv1;



%X nam predstavlja koje sve mogucnosti moze uzeti promenljiva X koja je
%izvlacenje papirica iz kutije pri cemu 0 znaci da nema papirica te nema
%papirica za izvuci

X=[0 1 2 3 4 5 6];


%ovo je kvadriranje matrice
X2=X.^2;

%ovo su verovatnoce koje odgovaraju datim vredostima u X, po tekstu zadatka
fgv=[0.02 vrv1 vrv2 vrv1 vrv1 vrv1 vrv1];

%radi lakseg predstavljanja funkcije raspodele 
F0=0;
F1=fgv(1);
F2=sum(fgv(1:2));
F3=sum(fgv(1:3));
F4=sum(fgv(1:4));
F5=sum(fgv(1:5));
F6=sum(fgv(1:6));
F7=sum(fgv(1:7));
%Funkcija raspodele promenljive X
F=[0 F1 F2 F3 F4 F5 F6 F7 1];



%po definiciji iz verovatnoce ovome je jednaka srednja vrednost
m=sum(X.*fgv);      

%izracunavamo prvo vaarijansu kao E(X^2)-(E(X))^2
varX=sum(X2.*fgv)-m^2;


%pravljenje slucajne promenljive sa datim verovatnocama
%funkcija rand daje random broj uniformno rasporedjen
%radnomi
n=1000; 
%kreiramo n razlicitih uniformnih brojeva
randomi = rand(1,n);

%na osnovu uniformne raspodele skaliramo na intervale 0:6
log1=(randomi>0.02 & randomi <=0.16);
log2=(randomi>0.16 & randomi <=0.44).*2;
log3=(randomi>0.44 & randomi <=0.58).*3;
log4=(randomi>0.58 & randomi <=0.72).*4;
log5=(randomi>0.72 & randomi <=0.86).*5;
log6=(randomi>0.86 & randomi <=1).*6;
%matrica sa brojevima od 0 do 6 koji su se dobili u  1000 puta
Rand0123456=log1+log2+log3+log4+log5+log6;

%ukupan broj dobijenih 0 1 2 3 4 5 6
sum0=sum(randomi<0.02);
sum1=sum(randomi>0.02 & randomi <=0.16);
sum2=sum(randomi>0.16 & randomi <=0.44);
sum3=sum(randomi>0.44 & randomi <=0.58);
sum4=sum(randomi>0.58 & randomi <=0.72);
sum5=sum(randomi>0.72 & randomi <=0.86);
sum6=sum(randomi>0.86 & randomi <=1);

a=Rand0123456;
%crtanje histograma
%figure,title('Histogram zadatak 1'),histogram(Rand0123456,'Normalization','probability');

%procena funkcija mase verovatnoce
Procenafgv=[sum0/n sum1/n sum2/n sum3/n sum4/n sum5/n sum6/n];
Fp0=0;
Fp1=Procenafgv(1);
Fp2=sum(Procenafgv(1:2));
Fp3=sum(Procenafgv(1:3));
Fp4=sum(Procenafgv(1:4));
Fp5=sum(Procenafgv(1:5));
Fp6=sum(Procenafgv(1:6));
Fp7=sum(Procenafgv(1:7));
ProcenaF=[0 Fp1 Fp2 Fp3 Fp4 Fp5 Fp6 Fp7 1];

%napomena za -1 i 7 : ubaceno je da bi se lepo video grafik "promene" u 0 i
%6 tj da se zna da je pre toga 0 a posle toga 1
xosa=[-1 0 1 2 3 4 5 6 7];
%F
figure(2);
hold on;
stairs(xosa,ProcenaF);stairs(xosa,F);
legend('eksperimentalno','analiticki');
title('Analiticka i eksperimentalna funkcija raspodele');
xlabel('vrednost papirica');
ylabel('F_x(vrednost papirica)');


analiticki=[0 0.02 0.14 0.28 0.14 0.14 0.14 0.14 0];
dobijeno=[0 sum0/n sum1/n sum2/n sum3/n sum4/n sum5/n sum6/n 0];
%fgv
figure;
hold on;
stem(xosa,analiticki);
stem(xosa,dobijeno);
legend('analiticki','eksperimentalno');
title('Funckija mase verovatnoce');
ylabel('p_x(vrednost)');
xlabel('vrednst papirica');

%matematicko ocekivanje eksperimenta
ProcenaM=sum(Rand0123456)/n;

%varijansa eksperimenta
ProcenaVarX=sum((Rand0123456-ProcenaM).^2)/(n-1);



%Zadatak 2 

% skiciranje grafika 
t=-0.8:0.001:1.4;
b=zeros(1,2201);
i=0;
for y=-0.8:0.001:1.4
    i=i+1;
    if(y<=-2/3)
        b(i)=0;
    end
    if(y>-2/3 && y<=0)
        b(i)=3/2*y+1;
    end
    if(y>0 && y<2/3)
        b(i)=0;
    end
    if(y>=2/3 && y<=4/3)
        b(i)=-3*y+4;
    end
    if(y>4/3)
        b(i)=0;
    end
end
%figure,plot(t,b)


%Tekst zadatka 
%Potrebno j generisati odbirke slucajne promenljice Y cija je fgv data u
%tabeli Q=1.
%a) izracunati vrednost realne konstante a.

% Iz racuna znamo da je a*1+2*a =1 odakle je 

a=2/3;



%c)Generisati N=10^5 odbiraka slucajne pormenljive Y i na osnovu njih
%proceniti pdgpvarajucu fgv koristeci histogram. Na istom grafiku prikazati
%i analiticku fgv iz tabele.
N2=10^5;
Rand2=rand(1,N2);
for x=1:N2
    if(Rand2(x)<=1/3 && Rand2(x)>=0)
        Rand2(x)=2*(sqrt(3*Rand2(x))-1)/3;
    elseif(Rand2(x)<=1 && Rand2(x)>1/3)
        Rand2(x)=(4-sqrt(6-6*Rand2(x)))/3;
    end
end

%figure,histogram(Rand2,50,'Normalization','pdf');
%hold all
%plot(t,b,'LineWidth',2)

%d)Analiticki odrediti matematicko ocekivanje i varijansu slucajne
%promenljive Y i uporediti ih sa eksperimentalnom procenjenom vrednoscu
%ovih parametara

%matematicko ocekivanje Y procena
ProcenaM2=sum(Rand2)/N2;

%varijansa Y
ProcenaVarY=sum((Rand2-ProcenaM2).^2)/(N2-1);


%Zadatak 3
N3=10^5;
sigma1=[3 2];
mu1=[0 0];

X1=randn(1,N3).*sigma1(1)+mu1(1);
X2=randn(1,N3).*sigma1(2)+mu1(2);
X=[X1;X2]';

%histogram ove dve promenljive
%figure,histogram2(X1,X2,20 ,'Normalization','pdf');
%xlabel ('x_1') ;  ylabel ('x_2') ;  zlabel ('f (x_1 , x_2)') ;
%title('Eksperimentalno dobijen histogram') ;

%analiticki  odredjena  fgv
xa1 = -3 * sigma1(1) : 0.3 : 3 * sigma1(1) ;
xa2 = -3 * sigma1(2) : 0.3 : 3 * sigma1(2) ;
[Xa1,Xa2] = meshgrid (xa1 , xa2 ) ;
Xanaliticki = [Xa1;Xa2]';
F = 1/(2*pi*sigma1(1)*sigma1(2) )*exp (-( Xa1.^2)*0.5/( sigma1(1) ^2)-(Xa2.^2)*0.5/( sigma1(2) ^2) ) ;

%figure(2),surf(xa1 , xa2 ,F) ;
%xlabel ('x_1') ;  ylabel ('x_2') ;  zlabel ('f (x_1 , x_2)') ;
%title('Analiticki  odredjena  fgv') ;


A=[0.611 -0.4;0 0.5];
b=[-2 1]';
Y=A*X'+b;
m1exp=sum(Y(1,:))/N3;
m2exp=sum(Y(2,:))/N3;


varY1 = sum((Y(1,:)-m1exp).^2)/(N3-1);
varY2 = sum((Y(2,:)-m2exp).^2)/(N3-1);
covY1Y2 = sum((Y(1,:)-m1exp).*(Y(2,:)-m2exp))/(N3-1);

A1=[0.4 0.8;0 0.5];
A2=[2/3 0 ;0 0.5];
Y1=A1*X'+b;
Y2=A2*X'+b;

%figure(3) ;
%plot (Y(1,:),Y(2,:) ,'x') ;
%xlabel ('Y_1') ;ylabel ('Y_2') ;
%title('Odbirci  vektora Y kada  je  \rho =-0 .4') ;

%figure(4) ;
%plot (Y1(1,:),Y1(2,:) ,'x') ;
%xlabel ('Y2_1') ;ylabel ('Y2_2') ;
%title('Odbirci  vektora Y_2 kada  je  \rho =0.8') ;

%figure(5) ;
%plot (Y2(1,:),Y2(2,:) ,'x') ;
%xlabel ('Y1_1') ;ylabel ('Y1_2') ;
%title('Odbirci  vektora Y_1 kada  je  \rho =0') ;


