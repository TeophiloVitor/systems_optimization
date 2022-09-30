clear
x=-%pi/2:0.05*%pi:%pi/2;
y=x;
x0=-1;
y0=1;
tol=1e-2;
[m n]=size(x);
//
for j=1:n
   for i=1:n
   z(i,j)=x(i)^2-sin(y(j))+1;
   end
   end
//
plot3d(x,y,z)
viter=[x0;y0];
for k=1:4
//
gradf=[2*x0;-cos(y0)];
hessf=[2 0;0 sin(y0)];
alfa(k)=gradf'*gradf/(gradf'*hessf*gradf);
vn=[x0;y0]-alfa(k)*gradf;
//vn=[x0;y0]-inv(hessf)*gradf;
viter=[viter vn];
if norm(vn-[x0;y0])<tol then
    disp('Convergiu com '+string(k)+' iterações')
end
//
x0=vn(1);
y0=vn(2);
//
end


