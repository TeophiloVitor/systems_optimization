clear
// Develop
clear
//
//x1=-32.768:0.5:32.768;
x1=-30:0.75:30;
//    x =-10:0.25:10;
//    y =-10:0.25:10;
//    x=x1;
//    y=x1;
//
x2=x1;
//
x=x1;
y=x2;
p1=x1;
p2=x2;
//
[X1 X2]=meshgrid(x1,x2);
a=20;
b=0.1;
c=0.5*%pi;
//
F10=-a*exp(-b*sqrt((X1.^2+X2.^2)/2))-exp((cos(c*X1)+cos(c*X2))/2)+a+exp(1);
//
//    J=0.5-(sin^2(sqrt(x^2+y^2))-0.5)/(1+0.001*(x^2+y^2))^2;
[n nx]=size(x);
[n ny]=size(y);
for i=1:nx
for j=1:ny
    z(i,j)=0.5+((sin(sqrt(x(i)^2+y(j)^2)))^2-0.5)/(1+0.001*(x(i)^2+y(j)^2))^2;
end
end
//plot3d(x,y,z)
//surf(X1,X2,F10)
Fobj=F10+z+a*cos(X1/30);
//
title('Ackleys path f10 + Schaeffer','fontsiz',4)
surf(X1,X2,Fobj)
//
//
npart=8;
itermax=20;
C1max=2.0;
C1min=0.5;
C1=2.0;
C2=2.0;
W=0.5;
Wmax=0.9;
Wmin=0.4;
acermax=1.0;
acermin=0.2;
//
x1max=max(p1);
x1min=min(p1);
x2max=x1max;
x2min=x1min;
//
//
// Nuvem inicial
        Y0 = grand(2,npart, "unf", x1min, x1max);
        Ybest=Y0;
        Ymean=[mean(Y0(1,:));mean(Y0(2,:))];
//
//
iter=0;
x1=Y0(1,:);
x2=Y0(2,:);
z12part=-a*exp(-b*sqrt((x1.^2+x2.^2)/2))-exp((cos(c*x1)+cos(c*x2))/2)+a+exp(1);
//
z12=0.5*ones(x1)-((sin(sqrt(x1.^2+x2.^2))).^2-0.5*ones(x1))./(ones(x1)+0.1*(x1.^2+x2.^2)).^2;
//plot3d(x,y,z)
//surf(X1,X2,F10)
z12=z12part+z12+a*cos(x1/30);
//
///////////////////////////////// Elitismo: garantia da média //////////
[WORST Posworst]=max(z12);
Y0(:,Posworst)=Ymean;
//
//x1=Y0(1,:);
//x2=Y0(2,:);
//z12=-a*exp(-b*sqrt((x1.^2+x2.^2)/2))-exp((cos(c*x1)+cos(c*x2))/2)+a+exp(1);
/////////////////////////////////////////////////////////////////////////
//
z12hist=z12;
zbest=z12;
[BEST Posbest]=min(z12);
histbest=BEST;
Gbest=Ybest(:,Posbest);
x1hist=Y0(1,:);
x2hist=Y0(2,:);
vel0=zeros(2,npart);
velk=vel0;
deltav=vel0;
//
////////////////////////// while /////
while iter<=itermax
//
pause
clf()
//
subplot(211)
title('Ackleys path f10 + Schaeffer','fontsiz',4)
plot3d(p1,p2,F10)
//
t1=Y0(1,:);
t2=Y0(2,:);
z1part=-a*exp(-b*sqrt((t1.^2+t2.^2)/2))-exp((cos(c*t1)+cos(c*t2))/2)+a+exp(1);
//
z2part=0.5*ones(t1)-((sin(sqrt(t1.^2+t2.^2))).^2-0.5*ones(t1))./(ones(t1)+0.1*(t1.^2+t2.^2)).^2;
zpart=z1part+z2part;
//
for i=1:npart
    P(i,:)=t1;
    Q(i,:)=t2;
    zpq(i,:)=zpart;
end
//
mesh(P,Q,zpq,'edgeco','b','marker','o','markersiz',9,'markeredg','red','linestyle','none');
//
//
    r1=rand();
    r2=rand();
    W=Wmax-iter*(Wmax-Wmin)/itermax;
    acer=acermax-iter*(acermax-acermin)/itermax;
    C1=C1max-iter*(C1max-C1min)/itermax;
for m=1:npart
    r1=rand();
    r2=rand();
    dist2best(m)=norm(Gbest-Y0(:,m));
    if m <> Posbest then
//
    velk(:,m)=W*(vel0(:,m)-0.5*acer*deltav(:,m))+r1*C1*(Ybest(:,m)-Y0(:,m))+r2*C2*(Gbest-Y0(:,m));
    //+r3*C3*velbest;
//
    else
//    velk(:,m)=(itermax-iter)*Wmin*vel0(m)/itermax;
    velk(:,m)=W*(vel0(:,m)+0.5**0.8*acer*deltav(:,m))+r1*C1*(Ybest(:,m)-Y0(:,m));
//    velk(:,Posbest)=[0 0]';
    end
end
//////////////////////////Atualização////////////
Y1test=Y0+velk;
//for k=1:npart
//    if norm(Gbest-Y1test(:,k))<dist2best(k) then
//        Y1(:,k)=Y1test(:,k);
//    else
//        Y1(:,k)=Y0(:,k);
//    end
//end
//Y1(:,Posbest)=Y0(:,Posbest)+velk(:,Posbest);
//////////////////////////////////////////////////////////
Y1=Y1test;
//////////////////////////////////////////////////////////
/////////////////////// Testa Limites ///////
for j=1:npart
        if Y1(1,j)>x1max then
            Y1(1,j)=x1max
            velk(1,j)=0;
        end
        if Y1(1,j)<x1min then
            Y1(1,j)=x1min
            velk(1,j)=0;
        end
        if Y1(2,j)>x2max then
            Y1(2,j)=x2max
            velk(2,j)=0;
        end
        if Y1(2,j)<x2min then
            Y1(2,j)=x2min
            velk(2,j)=0;
        end
end
/////////////////////////////////////////////////
Ymean=[mean(Y1(1,:));mean(Y1(2,:))];
t1=Y1(1,:);
t2=Y1(2,:);
//x1hist=[x1hist;t1];
//x2hist=[x2hist;t2];
z12=-a*exp(-b*sqrt((t1.^2+t2.^2)/2))-exp((cos(c*t1)+cos(c*t2))/2)+a+exp(1);
//z12hist=[z12hist;z12];
//
///////////////////////////////////////////////////////
[WORST Posworst]=max(z12);
//Y1(:,Posworst)=Ymean;
t1=Y1(1,:);
t2=Y1(2,:);
x1hist=[x1hist;t1];
x2hist=[x2hist;t2];
z12=-a*exp(-b*sqrt((t1.^2+t2.^2)/2))-exp((cos(c*t1)+cos(c*t2))/2)+a+exp(1);
//
z2part=0.5*ones(t1)-((sin(sqrt(t1.^2+t2.^2))).^2-0.5*ones(t1))./(ones(t1)+0.1*(t1.^2+t2.^2)).^2;
z12=z12+z2part;
//
z12hist=[z12hist;z12];
////////////////////////////////////////////////////////
for m=1:npart
    if z12(m)<zbest(m) then
        Ybest(:,m)=Y1(:,m);
        zbest(m)=z12(m);
    end
end
//
[BEST Posbest]=min(zbest);
if histbest > BEST then
    histbest=BEST;
    Gbest=Ybest(:,Posbest);
    //velk(:,Posbest)=[0 0]';
    //velbest=velk(:,Posbest);
end
//Gbest=Ybest(:,Posbest);
Y0=Y1;
deltav=velk-vel0;
vel0=velk;
//
iter=iter+1;
//
//
end
//
subplot(212)
title('Trajetória da melhor partícula','fontsiz',4)
//plot3d(t,t,z)
xx=x1hist(:,Posbest);
yy=x2hist(:,Posbest);
zz=z12hist(:,Posbest);
[nlin uuu]=size(xx);
for i=1:nlin
for j=1:nlin
 Pxx(i,j)=xx(i);
 Qyy(i,j)=yy(i);
 Zpq(i,j)=zz(i);
end
end
comet3d(Pxx,Qyy,Zpq,0.01)
