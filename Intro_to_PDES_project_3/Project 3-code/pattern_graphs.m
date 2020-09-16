amen=amen+1;
u=ufin{amen};
v=vfin{amen};
stab2(1,amen)
gamma=20+10*amen
figure(1)
mesh(xx,yy,u)
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right
title('Initial conditions 4 of u for N,M=129,L=2')


figure(2)
mesh(xx,yy,v)
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right
title('Initial conditions 4 of v for N,M=129,L=2')