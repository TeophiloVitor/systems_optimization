clear
// Ackley's path f10
clear
//
fa=sqrt(5)/2-0.5;
//
xq=[0 -1 -2 0 3 -2 -9 3];
yq=[0 1 0 -2 1 6 -2 -15];
xqmax=max(xq);
xqmin=min(xq);
yqmax=max(yq);
yqmin=min(yq);
//
//x1=-32.768:0.5:32.768;
x1=-50:1:50;
    x =-10:0.25:10;
    y =-10:0.25:10;
//    x=x1;
//    y=x1;
//
x2=x1;
//
x=x1;
y=x2;
p1=x1;
p2=p1;
//
[X1 X2]=meshgrid(x1,x2);
a=500;
b=0.1;
c=0.5*%pi;
//
F10=-a*exp(-b*sqrt((X1.^2+X2.^2)/2))-exp((cos(c*X1)+cos(c*X2))/2)+a+exp(1);
//
Fobj=F10;
//
title('Ackleys path f10','fontsiz',4)
surf(X1,X2,Fobj)
//
pause
//
npart=8;
itermax=10;
C1max=2.0;
C1min=0.5;
C1=2.0;
C2=2.0;
//W=0.5;
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
//
///////////////////////////////// Elitismo: garantia da média //////////
//[WORST Posworst]=max(z12);
//Y0(:,Posworst)=Ymean;
//
x1=Y0(1,:);
x2=Y0(2,:);
z12=-a*exp(-b*sqrt((x1.^2+x2.^2)/2))-exp((cos(c*x1)+cos(c*x2))/2)+a+exp(1);
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
sleep(1000)
clf()
//
subplot(211)
title('Ackleys path f10','fontsiz',4)
plot3d(p1,p2,F10)
//
t1=Y0(1,:);
t2=Y0(2,:);
zpart=-a*exp(-b*sqrt((t1.^2+t2.^2)/2))-exp((cos(c*t1)+cos(c*t2))/2)+a+exp(1);
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
    acer=0.0;
    C1p=C1max-iter*(C1max-C1min)/itermax;
    C1m=C1max-C1p;
    C1m=0;
for m=1:npart
    for j=1:2
        r1=rand();
        r2=rand();
    //
//        if r1 > 0.5 then
//            r1=fa;
//        else
//            r1=1-fa;
//        end
//    //
//        if r2 > 0.5 then
//            r2=fa;
//        else
//            r2=1-fa;
//        end
    //
//    dist2best(m)=norm(Gbest-Y0(:,m));
//
    if m <> Posbest then
//
    velk(:,m)=W*(vel0(:,m)-0.5*acer*deltav(:,m))+r1*C1*(Ybest(:,m)-Y0(:,m))+r2*C2*(Gbest-Y0(:,m));
    //+r3*C3*velbest;
    velk(:,m)=W*(vel0(:,m)-0.5*acer*deltav(:,m))+r1*(C1p*(Ybest(:,m)-Y0(:,m))+C1m*(Ymean-Y0(:,m)))+r2*C2*(Gbest-Y0(:,m));
//
    else
//    velk(:,m)=(itermax-iter)*Wmin*vel0(m)/itermax;
    velk(:,m)=W*(vel0(:,m)+0.5**0.8*acer*deltav(:,m))+r1*C1*(Ybest(:,m)-Y0(:,m));
//    velk(:,Posbest)=[0 0]';
    end
end
//
    if norm(velk(:,m))>0.8*norm(Y0(:,m)) then
        velk(:,m)=(0.8*velk(:,m)/norm(velk(:,m)))*norm(Y0(:,m));
    end
//
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
//
///////////////////////////////////////////////////////
t1=Y1(1,:);
t2=Y1(2,:);
x1hist=[x1hist;t1];
x2hist=[x2hist;t2];
z12=-a*exp(-b*sqrt((t1.^2+t2.^2)/2))-exp((cos(c*t1)+cos(c*t2))/2)+a+exp(1);
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
end
//
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
