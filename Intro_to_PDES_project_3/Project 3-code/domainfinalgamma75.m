a=0.07;
b=0.95;
gamma=75;
k1=0.0005;
du=1;
dv=10;
t=0;
babes=20;

N1=65;
M1=65;
N2=65;
M2=65;
L1=4;
L2=0.5;


hy1=L1/(N1-1);
hy2=L2/(N2-1);
hx2=1/(M2-1);


Bxu1=zeros(2,M1);
Byu1=zeros(2,N1);
Bxv1=zeros(2,M1);
Byv1=zeros(2,N1);
Bxu2=zeros(2,M2);
Byu2=zeros(2,N2);
Bxv2=zeros(2,M2);
Byv2=zeros(2,N2);

Bny=L2*(N2-1)/L1+1;
Bnx=11;

z=L1/L2;
w=log2(z);






xxin1=linspace(0,1,M1);
xxin2=xxin1+(1-(Bnx-1)*hx2);
yyin1=linspace(0,L1,N1);
yyin2=linspace(0,L2,N2);
xx2=linspace(1,2,M1-Bnx+1);


[xmat1,ymat1]=meshgrid(xxin1,yyin1);
[xmat2,ymat2]=meshgrid(xxin2,yyin2);



uit1=(a+b)*ones(N1,M1) ;
vit1=(b\((a+b)^2))*ones(N1,M1);
uit2=(a+b)*ones(N2,M2) ;
vit2=(b\((a+b)^2))*ones(N2,M2);
uit1=uit1+0.05*cos(ymat1*2*pi)*cos(xmat1*2*pi);
vit1=vit1+0.05*cos(ymat1*2*pi)*cos(xmat1*2*pi);
uit2=uit2+0.05*cos(ymat2*2*pi)*cos(xmat2*2*pi);
vit2=vit2+0.05*cos(ymat2*2*pi)*cos(xmat2*2*pi);


% figure(1)
% mesh(xxin1,yyin1,uit1)
% xlabel('x','FontSize',18)
% ylabel('y','FontSize',18)
% zlabel('u','FontSize',18)
% grid off
% lighting phong
% camlight headlight
% camlight right
% title('Initial conditions of u for N,M=129,L=2')
% 
% figure(2)
% mesh(xxin2,yyin2,uit2)
% xlabel('x','FontSize',18)
% ylabel('y','FontSize',18)
% zlabel('u','FontSize',18)
% grid off
% lighting phong
% camlight headlight
% camlight right
% title('Initial conditions of u for N,M=129,L=2')
% 
% figure(3)
% grid off
% mesh(xxin1,yyin1,uit1)
% hold on
% mesh(xxin2,yyin2,uit2)
% hold off
% xlabel('x','FontSize',18)
% ylabel('y','FontSize',18)
% zlabel('u','FontSize',18)
% lighting phong
% camlight headlight
% camlight right
% title('Initial conditions of u for N,M=129,L=2')

%xxin=linspace(0,2-(Bnx-1)*hx2,N1+N2-Bnx);

ufin1 = cell(babes,1);
ufin1{1} = uit1;
vfin1 = cell(babes,1);
vfin1{1} = vit1;
ufin2 = cell(babes,1);
ufin2{1} = uit2;
vfin2 = cell(babes,1);
vfin2{1} = vit2;







 i=1;
 

for  loop1=1:1000
    
   %solve first rectangle 
    [u1,v1]=turdomlin(a,b,gamma,N1,M1,L1,k1,uit1,vit1,Bxu1,Byu1,Bxv1,Byv1);
   
    %find boundaries data from first solution 
    gu=(u1(Bny+1,M1-Bnx+1:M1)-u1(Bny-1,M1-Bnx+1:M1))/(2*hy1);
    bodu=(u1(1:Bny,M1-Bnx)'-u1(1:Bny,M1-Bnx+2)')*(M1-1)*0.5;
    gv=(v1(Bny+1,M1-Bnx+1:M1)-v1(Bny-1,M1-Bnx+1:M1))/(2*hy1);
    bodv=(v1(1:Bny,M1-Bnx)'-v1(1:Bny,M1-Bnx+2)')*(M1-1)*0.5;
    
    %interpolate boundaries so data matches dimensions
    for gk=1:w
        bodu=interpolatebound(bodu);
        bodv=interpolatebound(bodv);
    end
    
   %alter boundaries
    Byu2(1,1:N2)=bodu(1,1:N2);
    Byv2(1,1:N2)=bodv(1,1:N2);
    Bxu2(2,1:Bnx)=gu;
    Bxv2(2,1:Bnx)=gv;
    
%find 2nd rectangles solution 
    [u2,v2]=turdomlin(a,b,gamma,N2,M2,L2,k1,uit2,vit2,Bxu2,Byu2,Bxv2,Byv2);

    %adjust solution 1 boundaries
    
    Byu1(2,1:Bny)=(u2(1:z:N2,Bnx+1)'-u2(1:z:N2,Bnx-1)')*(M1-1)*0.5;
    Byv1(2,1:Bny)=(v2(1:z:N2,Bnx+1)'-v2(1:z:N2,Bnx-1)')*(M1-1)*0.5;    
end 

%update initial conditions after 1st time step

 vit1=v1;
 uit1=u1;
 vit2=v2;
 uit2=u2;

 Niter=1;
 
 for timtimes=1:babes
 for iter1=1:2000
for  loop2=1:10
    
    
    [u1,v1]=turdomlin(a,b,gamma,N1,M1,L1,k1,uit1,vit1,Bxu1,Byu1,Bxv1,Byv1);
    
    gu=(u1(Bny+1,M1-Bnx+1:M1)-u1(Bny-1,M1-Bnx+1:M1))/(2*hy1);
    bodu=(u1(1:Bny,M1-Bnx)'-u1(1:Bny,M1-Bnx+2)')*(M1-1)*0.5;
    gv=(v1(Bny+1,M1-Bnx+1:M1)-v1(Bny-1,M1-Bnx+1:M1))/(2*hy1);
    bodv=(v1(1:Bny,M1-Bnx)'-v1(1:Bny,M1-Bnx+2)')*(M1-1)*0.5;
    
    
    for gk=1:w
        bodu=interpolatebound(bodu);
        bodv=interpolatebound(bodv);
    end
    
    
    
    Byu2(1,1:N2)=bodu(1,1:N2);
    Byv2(1,1:N2)=bodv(1,1:N2);
    
    Bxu2(2,1:Bnx)=gu;
    Bxv2(2,1:Bnx)=gv;
    
    
    
    [u2,v2]=turdomlin(a,b,gamma,N2,M2,L2,k1,uit2,vit2,Bxu2,Byu2,Bxv2,Byv2);

    
    
    Byu1(2,1:Bny)=(u2(1:z:N2,Bnx+1)'-u2(1:z:N2,Bnx-1)')*(M1-1)*0.5;
    Byv1(2,1:Bny)=(v2(1:z:N2,Bnx+1)'-v2(1:z:N2,Bnx-1)')*(M1-1)*0.5;
    
    
    
    
    
end 

%update solutions and time
vit1=v1;
 uit1=u1;
 vit2=v2;
 uit2=u2;
 


 end
 
 Niter=Niter+1;
ufin1{Niter} = u1;
% vfin1{Niter} = v1;
ufin2{Niter} = u2;
% vfin2{Niter} = v2;



 end
 

 
 
 
%  u3=u2(1:N2,Bnx:M2);
 
 
%  
%  pause(0.01)
% figure(7)
% grid off
% mesh(xxin1,yyin1,u1)
% hold on
% mesh(xx2,yyin2,u3)
% hold off
% xlabel('x','FontSize',18)
% ylabel('y','FontSize',18)
% zlabel('u','FontSize',18)
% lighting phong
% camlight headlight
% camlight right
% title('Initial conditions of u for N,M=129,L=2')
 






