clear all
//plot Lista PNLR - problema 2
x=3:-0.1:1;
p1=-1:0.2:2;
p2=-1:0.2:2;
[P1 P2] = meshgrid(p1,p2);
//objfunc=-2*P1-P2;
objfunc = P1.^2+P2.^2;
minfunc=min(objfunc);
maxfunc=max(objfunc);
//
y1=(x-1).*sqrt(x-1);
//y1=(x-1).*sqrt(x-1);
y2=-(x-1).*sqrt(x-1);
//y2=-(x-1).*sqrt(x-1);
//y=2*x;
//
[m n]=size(x);
[m np]=size(p1);
Z=zeros(n,n);
Y1=zeros(n,n);
Y2=zeros(n,n);
X=zeros(n,n);
for i=1:n
    Z(i,:)=x;//*(maxfunc-minfunc);
    if i<=n then
        X(i,:)=x;
        Y1(i,:)=y1;
        Y2(i,:)=y2;
     end
end

W=Z';

surf(X,Y1,W,'facecolor','blue')
surf(X,Y2,W,'facecolor','blue')
surf(P1,P2,objfunc,'facecolor','yellow')
set(gca(),"grid",[1 1 1])
