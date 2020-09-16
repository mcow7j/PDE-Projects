A=1;
u=linspace(-0.5,2.5,10000);
y=A.*u.*(1-u).*(2-u);
w=u*0;
yy=linspace(-2,2,1000);
xx=yy*0;
plot(u,y,'r',u,w,'bl',xx,yy,'bl')
xlabel('u')
ylabel('u_t')
title('phase plane analysis for A=1, ignoring diffusive terms')