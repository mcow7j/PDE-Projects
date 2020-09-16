amen=amen+1;
u=ufin{amen};
v=vfin{amen};
t=stab2(1,amen);
gamma=20+10*amen;

figure(1)
mesh(xx,yy,u)
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right
title(['Full Multigrid solution of u for N,M=129,L=2,t=' num2str(t) ' \gamma=' num2str(gamma)])


figure(2)
mesh(xx,yy,v)
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right
title(['Full Multigrid solution of v for N,M=129,L=2,t=' num2str(t) ' \gamma=' num2str(gamma)])