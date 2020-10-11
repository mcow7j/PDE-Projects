
%inputs for polar graph :

%time step at which to start shoopwing polar graph
%can skip extra iterations when showing graphs
%pause time of graphs

iter_1=1;
iterstep=50;
t1=0.001;



%inputs needed to create raw data, k being desired time step, N+1 being
%number of points used for r, iter being the nuber of time steps desired.
%eta relevant to intial equation
N=100;
iter=20000;
eta=0.01;
k=0.0001/(3*eta-0.0001);
%inputs for one value a and b against time
%r_i1 and r_i2 are postiotns of r_i used

%r_i1=N/4;
%r_i2=(3*N/2);


%following code produces raw data

%sets up basica variables
r=0:(2/N):2;
a=zeros(iter+1,N+1);
b=zeros(iter+1,N+1);
a(1,1:N+1)=r;
a(2:iter+1,N+1)=2*ones(iter,1);


%s inictaes which entry of r for which r>1, producing q at which it happens
s=mod(N,2);
if s==0 
    q=(N/2)+1;
else
    q=(N+1)/2;
end

%set up coefficents matrix
c=zeros(3,N+1);
for i=2:N
c(1:3,i)=[k*eta*((N/2)^2-(N/(4*r(1,i))));(1-k*eta*(0.5*(N^2)+(1/(r(1,i)^2)))); k*eta*((N/2)^2+(N/(4*r(1,i))))]   ;
end

%j represents time step 
%i represents postion with respect to r e.g i=10 represent r_10

for j=1:iter
    for i=2:q
       a(j+1,i)=(a(j,i-1:i+1)*c(1:3,i))+k*b(j,i);
       b(j+1,i)=(b(j,i-1:i+1)*c(1:3,i))-k*a(j,i);
    end
    %represents how w changes so w=0 for r>1
    for i=q+1:N
        a(j+1,i)=(a(j,i-1:i+1)*c(1:3,i));
        b(j+1,i)=(b(j,i-1:i+1)*c(1:3,i));
    end
end

%produces a and b, matrices with row j representing value of a(or b) at time j*k
%and column i presenting its value at r_i. 


%a and b graph plot

%t=0:k:k*(iter);
%t=t';
%figure(1)
%plot(t,a(1:iter+1,r_i1),t,b(1:iter+1,r_i1),t,a(1:iter+1,r_i2),t,b(1:iter+1,r_i2))
%title('Value of a(t) and b(t) for particular values of r')
%xlabel('time(seconds))')
%ylabel('value a(t)and b(t)')
%legend('value of a at r=0.5','value of b at r=0.5','value of a at r=1.5','value of b at r=1.5')

%the following code plots the polar graph of raw data
for i=iter_1:iterstep:iter+1
[r,theta] = meshgrid(0:2/N:2,-pi:pi/20:pi);

z=a(i,1:N+1).*cos(theta)-b(i,1:N+1).*sin(theta);

x = r.*cos(theta);

y = r.*sin(theta);



pause(t1) 
figure(3)

contourf(x,y,z,20,'r');
colorbar
% Hide the POLAR function data and annotations

%set(h,'Visible','on')

% Turn off axes and set square aspect ratio

axis off
%colorbar
axis image
end