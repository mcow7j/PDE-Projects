gam=10;
stab=zeros(1,gam);
ufin=cell(gam,1);
vfin=cell(gam,1);


 for lol=1:gam
a=0.07;
b=0.95;
gamma=20+10*lol;

N=129;
M=129;
L=2;
k1=0.001;
du=1;
dv=10;
t=0;

u=(a+b)*ones(N,M) ;
v=(b\((a+b)^2))*ones(N,M);
B1x=zeros(N,M);
B1y=zeros(N,M);
xx=linspace(0,1,M);
yy=linspace(0,L,N);
[xmat,ymat]=meshgrid(xx,yy);
u=u+0.05*cos(ymat*20);
v=v+0.05*cos(ymat*20);
uit3=8*u;
vit3=8*v;

 while( max(max(abs(u-uit3)))+max(max(abs(v-vit3)))>0.00005 && t<10)
    uit3=u;
    vit3=v;
    for i=1:50
    
    vit=v;
    uit=u;   
    
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
    t=t+k1;
    end
           
 end  
 stab(1,lol)=t;  
 ufin{lol}=v;
vfin{lol}=u;
end
 