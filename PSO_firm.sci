clear
t=linspace(-2*%pi,2*%pi,2*30);
function z=my_surface(x,y),z=x*sin(x)^2*cos(y)+0.7*y^2,endfunction
//
//
z=feval(t,t,my_surface);
title('z=x*sin(x)^2*cos(y)+0.7*y^2','fontsiz',4)
plot3d(t,t,z,'edgeco','none')
//
//
npart=8;
itermax=10;
C1max=2.0;
C1min=1.5;
C1=2.0;
C2=2.0;
W=0.5;
Wmax=0.9;
Wmin=0.6;
acermax=1.0;
acermin=0.2;
//
x1max=2*%pi;
x1min=-2*%pi;
x2max=x1max;
x2min=x1min;
//
//
// Nuvem inicial
        Y0 = grand(2,npart, "unf", x1min, x1max);
        Ybest=Y0;
        Ymean=[mean(Y0(1,:));mean(Y0(2,:))];
//
x1hist=Y0(1,:);
x2hist=Y0(2,:);
//
iter=0;
x1=Y0(1,:);
x2=Y0(2,:);
z12=(x1.*sin(x1).^2).*cos(x2)+0.7*x2.^2;
//
///////////////////////////////// Elitismo: garantia da média //////////
[WORST Posworst]=max(z12);
Y0(:,Posworst)=Ymean;
//
x1=Y0(1,:);
x2=Y0(2,:);
z12=(x1.*sin(x1).^2).*cos(x2)+0.7*x2.^2;
/////////////////////////////////////////////////////////////////////////
//
z12hist=z12;
zbest=z12;
[BEST Posbest]=min(z12);
histbest=BEST;
Gbest=Ybest(:,Posbest);
vel0=zeros(2,npart);
velk=vel0;
deltav=vel0;
//
////////////////////////// while /////
//
while iter<=itermax
//
pause
clf()
//
z=feval(t,t,my_surface);
subplot(211)
title('z=x*sin(x)^2*cos(y)+0.7*y^2','fontsiz',4)
plot3d(t,t,z)
//
t1=Y0(1,:);
t2=Y0(2,:);
zpart=(t1.*sin(t1).^2).*cos(t2)+0.7*t2.^2;
for i=1:npart
    P(i,:)=t1;
    Q(i,:)=t2;
    zpq(i,:)=zpart;
end
//
mesh(P,Q,zpq,'edgeco','b','marker','o','markersiz',9,'markeredg','red');
//
//
    r1=rand();
    r2=rand();
    W=Wmax-iter*(Wmax-Wmin)/itermax;
    C1=C1max-iter*(C1max-C1min)/itermax;
    acer=acermax-iter*(acermax-acermin)/itermax;
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
    velk(:,m)=W*(vel0(:,m)+0.5*0.8*acer*deltav(:,m))+r1*C1*(Ybest(:,m)-Y0(:,m));
//    velk(:,Posbest)=[0 0]';
    end
end
Y1=Y0+velk;
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
//////////////////////////////////////////////
Ymean=[mean(Y1(1,:));mean(Y1(2,:))];
x1=Y1(1,:);
x2=Y1(2,:);
//x1hist=[x1hist;x1];
//x2hist=[x2hist;x2];
z12=(x1.*sin(x1).^2).*cos(x2)+0.7*x2.^2;
//z12hist=[z12hist;z12];
/////////////////////////////////////////////////
[WORST Posworst]=max(z12);
Y1(:,Posworst)=Ymean;
x1=Y1(1,:);
x2=Y1(2,:);
x1hist=[x1hist;x1];
x2hist=[x2hist;x2];
z12=(x1.*sin(x1).^2).*cos(x2)+0.7*x2.^2;
z12hist=[z12hist;z12];
////////////////////////////////////////////////////////
//
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
//clf()
//
//z=feval(t,t,my_surface);
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
