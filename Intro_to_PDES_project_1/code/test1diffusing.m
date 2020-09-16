N=257;
M=257;
L=2;
Niter=2000;
k1=0.0005;


t=0;
uit=uexact4(N,M,L,t);
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




for i=1:Niter
    Q=Qmaker4(uit,N,M,L,(t+0.5*k1));
    A=Acalculator1(uit,L,k1);
    f=-A-k1*Q;
    u = FullMGL2time(uit,f,L,k1);
    uit=u;
    t=t+k1;
    u1=uexact4(N,M,L,t);
%  pause(0.01)       
%        figure(2)
%    mesh(xx,yy,u)
%   xlabel('x','FontSize',18)
%    ylabel('y','FontSize',18)
%    zlabel('u','FontSize',18)
%      grid off
%    lighting phong
%    camlight headlight
%   camlight right
 %title('numerical method with N,M=65,t=0.05,k=0.0005')
%  
% figure(3)
% mesh(xx,yy,u1)
% xlabel('x','FontSize',18)
% ylabel('y','FontSize',18)
% zlabel('u','FontSize',18)
% grid off
% lighting phong
% camlight headlight
%  camlight right

%  figure(4)
%  contour(xx,yy,abs(u-u1),20)
%  colorbar
%  xlabel('x','FontSize',18)
%  ylabel('y','FontSize',18)
%  zlabel('u','FontSize',18)
%  grid off
%  lighting phong
%  camlight headlight
%  camlight right
% title('error of numerical method with N,M=65,t=0.05,k=0.0005')
end
% t
% 
 figure(2)
  mesh(xx,yy,u)
  xlabel('x','FontSize',18)
  ylabel('y','FontSize',18)
  zlabel('u','FontSize',18)
  grid off
  lighting phong
  camlight headlight
  camlight right
  title('numerical method with N,M=257,t=1.00,k=0.0005')
  
  figure(4)
 contour(xx,yy,abs((u-u1)),20)
 colorbar
 xlabel('x','FontSize',18)
 ylabel('y','FontSize',18)
 zlabel('u','FontSize',18)
 grid off
 lighting phong
 camlight headlight
camlight right
 title('error of numerical method with N,M=257,t=1.00,k=0.0005')
  