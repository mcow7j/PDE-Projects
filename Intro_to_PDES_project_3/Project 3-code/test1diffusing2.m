N=129;
M=129;
L=2;
Niter=20000;
k1=0.0005;
d=1;

t=0;
uit=uexact(N,M,L,t);
u=uit;

xx=linspace(0,1,M);
yy=linspace(0,L,N);

figure(1)
mesh(xx,yy,uit)
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right
title('Initial conditions for (6), with L=2,M=129,N=129,k=0.0005')



B1x=zeros(N,M);
B1y=zeros(N,M);

unum=cell(Niter,1);
uexac=cell(Niter,1);

for i=1:Niter
    Q=Qmakerdtest1(uit,N,M,L,(t+0.5*k1));
    A=Acalculator1dhomo(uit,L); 
    f=(-uit-0.5*k1*d*A-k1*Q);
    u = FullMGL2timedbound2(uit,f,L,k1,d,B1x,B1y);
    uit=u;
    t=t+k1;
    u1=uexact(N,M,L,t);
 if mod(i,100)==0
    unum{i}=u;
    uexac{i}=u1;
 else
    
    
    
%     pause(0.01)
%         figure(2)
%     mesh(xx,yy,u)
%    xlabel('x','FontSize',18)
%     ylabel('y','FontSize',18)
%     zlabel('u','FontSize',18)
%       grid off
%     lighting phong
%     camlight headlight
%    camlight right
%  %title('numerical method with N,M=65,t=0.05,k=0.0005')
% %  
%  figure(3)
%  mesh(xx,yy,u1)
%  xlabel('x','FontSize',18)
%  ylabel('y','FontSize',18)
%  zlabel('u','FontSize',18)
%  grid off
% lighting phong
%  camlight headlight
%   camlight right
% 
%   figure(4)
%   contour(xx,yy,abs(u-u1),40)
%   colorbar
%   xlabel('x','FontSize',18)
%   ylabel('y','FontSize',18)
%   zlabel('u','FontSize',18)
%   grid off
%   lighting phong
%   camlight headlight
%   camlight right
% title('error of numerical method with N,M=65,t=0.05,k=0.0005')
end
end% t
% 
%  figure(2)
%   mesh(xx,yy,u)
%   xlabel('x','FontSize',18)
%   ylabel('y','FontSize',18)
%   zlabel('u','FontSize',18)
%   grid off
%   lighting phong
%   camlight headlight
%   camlight right
%   title('numerical method with N,M=257,t=1.00,k=0.0005')
%   
%   figure(4)
%  contour(xx,yy,abs((u-u1)),20)
%  colorbar
%  xlabel('x','FontSize',18)
%  ylabel('y','FontSize',18)
%  zlabel('u','FontSize',18)
%  grid off
%  lighting phong
%  camlight headlight
% camlight right
%  title('error of numerical method with N,M=257,t=1.00,k=0.0005')
  