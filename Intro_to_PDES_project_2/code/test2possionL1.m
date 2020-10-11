M=513;
N=513;
L=2;

xx = linspace(0,1,M);
yy = linspace(0,L,N);
f=zeros(N,M);
u=f;
for n=1:N
    for m=1:M
        x=xx(1,m);
        y=yy(1,n);
     f(n,m)=((y-2)*y*x*(x+3)+(x-1)*x*(y*y-6*y+6))*exp(x-y);
     u(n,m)=(x-1)*x*y*(y-2)*exp(x-y);
    end
end

%
u1=FullMGL1(f,L);

z=u1-u;

 figure(1)
mesh(xx,yy,u)
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right
title('GS solution with multigrid of (9) with N,M=513')

 figure(2)
mesh(xx,yy,u1)
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right
title('exact solution of (9) with N,M=65')


 figure(3)
contour(xx,yy,z)
colorbar
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right
title('error of GS solution with multigrid of (9) with N,M=513')

figure(4)
mesh(xx,yy,z)
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right
title('error of GS solution with multigrid of (9) with N,M=65')