%inputs needed to create raw data, k being desired time step, N+1 being
%number of points used for r, iter being the number of time steps desired.
%eta relevant to intial equation


N=100;
iter=100000;
eta=0.04;
k=0.0002/(3*eta-0.0002);
[a,b]=rawdatafindiff(N,k,iter,eta);

%i is the time step at which we want the graph to be produced





[r,theta] = meshgrid(0:2/N:2,-pi:pi/20:pi);
x = r.*cos(theta);
y = r.*sin(theta);


i=3/k-mod(3/k,1)+1;
z=a(i,1:N+1).*cos(theta)-b(i,1:N+1).*sin(theta);
figure(1)
contourf(x,y,z,20,'r');
colorbar
title('graph showing contours of \phi at t=3, for \eta=0.04')
axis image



i=6/k-mod(6/k,1)+1;
z=a(i,1:N+1).*cos(theta)-b(i,1:N+1).*sin(theta);
figure(2)
contourf(x,y,z,20,'r');
colorbar
title('graph showing contours of \phi at t=6, for \eta=0.04')
axis image

i=9/k-mod(9/k,1)+1;
z=a(i,1:N+1).*cos(theta)-b(i,1:N+1).*sin(theta);
figure(3)
contourf(x,y,z,20,'r');
colorbar
title('graph showing contours of \phi at t=9, for \eta=0.04')
axis image

i=iter;
z=a(i,1:N+1).*cos(theta)-b(i,1:N+1).*sin(theta);
figure(4)
contourf(x,y,z,20,'r');
colorbar
title('graph showing contours of \phi for t very large, for \eta=0.04')
axis image

t=0:k:k*(iter);
t=t';
figure(5)
plot(t,a(1:iter+1,26),t,b(1:iter+1,26),t,a(1:iter+1,76),t,b(1:iter+1,76))
title('Value of a(t) and b(t) for particular values of r with \eta=0.04')
xlabel('time(seconds))')
ylabel('value a(t)and b(t)')
legend('value of a at r=0.5','value of b at r=0.5','value of a at r=1.5','value of b at r=1.5')





