N=129;
M=129;
L=2;
Niter=20000;
k1=0.0005;
d=1;

t=0;
uit=uexact3(N,M,L,t);
u=uit;
u3=uit;
uit3=u;
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
title('Initial conditions for (7), with L=2,M=129,N=129,k=0.0005')


unumm=cell(Niter,1);
unumFM=cell(Niter,1);
uexac=cell(Niter,1);

Bx=zeros(2,M);
By=zeros(2,N);
By(2,1:N)=ones(1,N)*(5*exp(-t));

for i=1:Niter
    Q=Qmakerdtest3(uit,N,M,L,(t+0.5*k1));
     A=Acalculator1d7(uit,L,Bx,By); 
    By(2,1:M)=ones(1,M)*(5*exp(-(t+k1)));
    f=(-uit-0.5*k1*d*A-k1*Q);
    u=FullMGL2timedbound7(uit,f,L,k1,d,Bx,By);
     u3=MultigridVL1timedbound7(uit3,f,L,k1,d,Bx,By);
%     u = GSL1timedbound7(uit,f,300,L,k1,d,Bx,By);
    uit=u;
    uit3=u3;
    t=t+k1;
    u1=uexact3(N,M,L,t);
 
   if mod(i,10)==0
    unumm{i}=u;
    unumFM{i}=u3;
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
 %title('numerical method with N,M=65,t=0.05,k=0.0005')
%  
%  figure(3)
%  mesh(xx,yy,u1)
%  xlabel('x','FontSize',18)
%  ylabel('y','FontSize',18)
%  zlabel('u','FontSize',18)
%  grid off
% lighting phong
%  camlight headlight
%   camlight right

%   figure(4)
%   contour(xx,yy,abs(u1-u),20)
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
end
% t
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
  