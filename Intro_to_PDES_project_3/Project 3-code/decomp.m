a=0.07;
b=0.95;
gamma=160;
Niter=100;
k1=0.005;
du=1;
dv=10;
t=0;


N1=65;
M1=65;
N=65;
M=65;
L1=16;
uit1=(a+b)*ones(N1,M1) ;
vit1=(b\((a+b)^2))*ones(N1,M1);

N2=65;
M2=65;
L2=1;
hy1=L1/(N1-1);
hy2=L2/(N2-1);
uit2=(a+b)*ones(N2,M2) ;
vit2=(b\((a+b)^2))*ones(N2,M2);
xx=linspace(0,1,M1);
yy1=linspace(0,L1,N1);
yy2=linspace(0,L2,N2);


% xx=linspace(0,1,M);
% yy=linspace(0,L,N);
% [xmat,ymat]=meshgrid(xx,yy);
% uit=uit+0.1*cos(xmat*20);
% vit=vit+0.1*cos(ymat*20);


Bxu1=zeros(2,M1);
Byu1=zeros(2,N1);
Bxv1=zeros(2,M1);
Byv1=zeros(2,N1);
Bxu2=zeros(2,M2);
Byu2=zeros(2,N2);
Bxv2=zeros(2,M2);
Byv2=zeros(2,N2);

Bny=L2*(N2-1)/L1+1;
Bnx=15;

z=L1/L2;
w=log2(z);

Num=(N1-1)*z+1;

uit3=zeros(Num,M1+M2-Bnx);
vit3=zeros(Num,M1+M2-Bnx);
xx3=1:Num;
yy3=1:M1+M2-Bnx;

ufin1 = cell(Niter,1);
ufin1{1} = uit1;
vfin1 = cell(Niter,1);
vfin1{1} = vit1;
ufin2 = cell(Niter,1);
ufin2{1} = uit2;
vfin2 = cell(Niter,1);
vfin2{1} = vit2;

ufin3 = cell(Niter,1);
ufin3{1} = uit2;
vfin3 = cell(Niter,1);
vfin3{1} = vit2;


 i=1;
 
    for  il=1:100
    [u1,v1]=turdomlin(a,b,gamma,N1,M1,L1,k1,uit1,vit1,Bxu1,Byu1,Bxv1,Byv1);
    vit1=v1;
    uit1=u1;
    
    
    
   
    Bxu2(2,1:Bnx)=1:Bnx;
     Bxv2(2,1:Bnx)=1:Bnx;
    
    gu=(-u1(Bny+1,M1-Bnx:M1)+u1(Bny-1,M1-Bnx:M1))/2*hy1;
    bodu=(u1(1:Bny,Bnx+1)'-u1(1:Bny,Bnx-1)')*(M-1)*0.5;
     gv=(-v1(Bny+1,M1-Bnx:M1)+v1(Bny-1,M1-Bnx:M1))/2*hy1;
    bodv=(v1(1:Bny,Bnx+1)'-v1(1:Bny,Bnx-1)')*(M-1)*0.5;
    
    solu=u1;
    solv=v1;
    
    
    for gk=1:w
        bodu=interpolatebound(bodu);
        bodv=interpolatebound(bodv);
         solu=interpolateuni(solu);
        solv=interpolateuni(solv);
    end
    
     uit3(1:Num,1:M1)=solu; 
    vit3(1:Num,1:M1)=solv; 
    
    
    Byu2(1,1:N2)=bodu(1,1:N2);
    Byv2(1,1:N2)=bodv(1,1:N2);
    
    [u2,v2]=turdomlin(a,b,gamma,N2,M2,L2,k1,uit2,vit2,Bxu2,Byu2,Bxv2,Byv2);
    vit2=v2;
    uit2=u2;
    
     uit3(1:N2,M1+1:M2+M1-Bnx)=u2(1:N2,Bnx+1:M2);
     vit3(1:N2,M1+1:M2+M1-Bnx)=v2(1:N2,Bnx+1:M2);
    Byu1(2,1:Bny)=(u2(1:z:N2,Bnx-1)'-u2(1:z:N2,Bnx+1)')*(M-1)*0.5;
    Byv1(2,1:Bny)=(v2(1:z:N2,Bnx-1)'-v2(1:z:N2,Bnx+1)')*(M-1)*0.5;
    end 
   
%     ufin1{i+1} = u1;
%     vfin1{i+1} = v1;
%     ufin2{i+1} = u2;
%     vfin2{i+1} = u2;


for i=1:Niter
    
    for ikj=1:40
    for  il=1:8
    [u1,v1]=turdomlin(a,b,gamma,N1,M1,L1,k1,uit1,vit1,Bxu1,Byu1,Bxv1,Byv1);
    vit1=v1;
    uit1=u1;
    
    
    
   
    Bxu2(2,1:Bnx)=1:Bnx;
     Bxv2(2,1:Bnx)=1:Bnx;
    
    gu=(-u1(Bny+1,M1-Bnx:M1)+u1(Bny-1,M1-Bnx:M1))/2*hy1;
    bodu=(u1(1:Bny,Bnx+1)'-u1(1:Bny,Bnx-1)')*(M-1)*0.5;
     gv=(-v1(Bny+1,M1-Bnx:M1)+v1(Bny-1,M1-Bnx:M1))/2*hy1;
    bodv=(v1(1:Bny,Bnx+1)'-v1(1:Bny,Bnx-1)')*(M-1)*0.5;
    
    solu=u1;
    solv=v1;
    
    
    for gk=1:w
        bodu=interpolatebound(bodu);
        bodv=interpolatebound(bodv);
         solu=interpolateuni(solu);
        solv=interpolateuni(solv);
    end
    
     uit3(1:Num,1:M1)=solu; 
    vit3(1:Num,1:M1)=solv; 
    
    
    Byu2(1,1:N2)=bodu(1,1:N2);
    Byv2(1,1:N2)=bodv(1,1:N2);
    
    [u2,v2]=turdomlin(a,b,gamma,N2,M2,L2,k1,uit2,vit2,Bxu2,Byu2,Bxv2,Byv2);
    vit2=v2;
    uit2=u2;
    
    
    Byu1(2,1:Bny)=(u2(1:z:N2,Bnx-1)'-u2(1:z:N2,Bnx+1)')*(M-1)*0.5;
    Byv1(2,1:Bny)=(v2(1:z:N2,Bnx-1)'-v2(1:z:N2,Bnx+1)')*(M-1)*0.5;
    end
    end 
    
    solu=u1;
    solv=v1;
    
    for gk2=1:w
         solu=interpolateuni(solu);
        solv=interpolateuni(solv);
    end
    
    uit3(1:Num,1:M1)=solu; 
    vit3(1:Num,1:M1)=solv; 
    
    uit3(1:N2,M1+1:M2+M1-Bnx)=u2(1:N2,Bnx+1:M2);
    vit3(1:N2,M1+1:M2+M1-Bnx)=v2(1:N2,Bnx+1:M2);
    
    ufin3{i+1} = uit3;
    vfin3{i+1} = vit3;
    
    
%     ufin1{i+1} = u1;
%     vfin1{i+1} = v1;
%     ufin2{i+1} = u2;
%     vfin2{i+1} = v2;
end
