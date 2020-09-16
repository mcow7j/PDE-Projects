

a=0.07;
b=0.95;
gamma=160;

N=129;
M=129;
L=2;
Niter=2000;
k1=0.0005;
du=1;
dv=10;
t=0;

uit=(a+b)*ones(N,M);
vit=(b\((a+b)^2))*ones(N,M);
B1x=zeros(N,M);
B1y=zeros(N,M);
xx=linspace(0,1,M);
yy=linspace(0,L,N);
[xmat,ymat]=meshgrid(xx,yy);
uit=uit+0.05*cos(20*xmat);
vit=vit+0.05*cos(ymat*20);







figure(1)
mesh(xx,yy,uit)
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right




  for i=1:2000
           QoldU=QoldmakerUderiv(uit,vit,gamma,N,M);
      QU=QmakerUlin(uit,vit,QoldU,gamma,a,N,M);
      AU=Acalculator1dhomo(uit,L);
    fu=(-uit-0.5*k1*du*AU-k1*QU);
    u =FullMGL2timedbound2lin(uit,QoldU,fu,L,k1,du);
     
     
     
     QoldV=QoldmakerVderiv(u,gamma,N,M);
      QV=QmakerVlin(u,vit,QoldV,gamma,b,N,M);
      AV=Acalculator1dhomo(vit,L);
      fv=(-vit-0.5*k1*dv*AV-k1*QV);
    v =FullMGL2timedbound2lin(vit,QoldV,fv,L,k1,dv);
%      
%      
     vit=v;
    uit=u;
      t=t+k1;
%    
    
  
%     figure(3)
%   contour(xx,yy,u,1)
%   colorbar
%   xlabel('x','FontSize',18)
%   ylabel('y','FontSize',18)
%   zlabel('u','FontSize',18)
%   grid off
%   lighting phong
%   camlight headlight
%   camlight right
%     
%   figure(4)
%   contour(xx,yy,v,2)
%   colorbar
%   xlabel('x','FontSize',18)
%   ylabel('y','FontSize',18)
%   zlabel('u','FontSize',18)
%   grid off
%   lighting phong
%   camlight headlight
%   camlight right
%     
% %     
%      pause(0.0001)
%   figure(2)
%      mesh(xx,yy,u)
%     xlabel('x','FontSize',18)
%      ylabel('y','FontSize',18)
%     zlabel('u','FontSize',18)
%        grid off
%      lighting phong
%      camlight headlight
%    camlight right   
% %     
%  
%  figure(3)
%     mesh(xx,yy,v)
%    xlabel('x','FontSize',18)
%     ylabel('y','FontSize',18)
%     zlabel('u','FontSize',18)
%       grid off
%     lighting phong
%     camlight headlight
%    camlight right   
%    
%    
%    
%  figure(4)
%     mesh(xx,yy,v-u)
%    xlabel('x','FontSize',18)
%     ylabel('y','FontSize',18)
%     zlabel('u','FontSize',18)
%       grid off
%     lighting phong
%     camlight headlight
%    camlight right   
% %    
%    
  end

 for i=1:Niter
%     
%     QoldU=QoldmakerUderiv(uit,vit,gamma,N,M);
%     QU=QmakerUlin(uit,vit,QoldU,gamma,a,N,M);
%     AU=Acalculator1dhomo(uit,L);
%     fu=(-uit-0.5*k1*du*AU-k1*QU);
%     u =FullMGL2timedbound2lin(uit,QoldU,fu,L,k1,du);
%     
%     
%    
%     QoldV=QoldmakerVderiv(u,gamma,N,M);
%     QV=QmakerVlin(u,vit,QoldV,gamma,b,N,M);
%     AV=Acalculator1dhomo(vit,L);
%     fv=(-vit-0.5*k1*dv*AV-k1*QV);
%     v =FullMGL2timedbound2lin(vit,QoldV,fv,L,k1,dv);
%     
%     
%     
%     vit=v;
%     uit=u;
%     t=t+k1;
%   
%     
%   
% %     figure(3)
% %   contour(xx,yy,u,1)
% %   colorbar
% %   xlabel('x','FontSize',18)
% %   ylabel('y','FontSize',18)
% %   zlabel('u','FontSize',18)
% %   grid off
% %   lighting phong
% %   camlight headlight
% %   camlight right
% %     
% %   figure(4)
% %   contour(xx,yy,v,2)
% %   colorbar
% %   xlabel('x','FontSize',18)
% %   ylabel('y','FontSize',18)
% %   zlabel('u','FontSize',18)
% %   grid off
% %   lighting phong
% %   camlight headlight
% %   camlight right
% %     
% %     
     pause(0.001)
  figure(2)
     mesh(xx,yy,u)
    xlabel('x','FontSize',18)
     ylabel('y','FontSize',18)
     zlabel('u','FontSize',18)
       grid off
     lighting phong
     camlight headlight
    camlight right   
%     
%  
% %  figure(3)
% %     mesh(xx,yy,v)
% %    xlabel('x','FontSize',18)
% %     ylabel('y','FontSize',18)
% %     zlabel('u','FontSize',18)
% %       grid off
% %     lighting phong
% %     camlight headlight
% %    camlight right   
% %    
% %    
% %    
% %  figure(4)
% %     mesh(xx,yy,v-u)
% %    xlabel('x','FontSize',18)
% %     ylabel('y','FontSize',18)
% %     zlabel('u','FontSize',18)
% %       grid off
% %     lighting phong
% %     camlight headlight
% %    camlight right   
%    
%    
%    
% 
%     t
 end


udone=u;
vdone=v;

