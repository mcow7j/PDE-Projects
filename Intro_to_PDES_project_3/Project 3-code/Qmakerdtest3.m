function Q=Qmakerdtest3(uu,N,M,L,t)

xx=linspace(0,1,M);
yy=linspace(0,L,N);
Q=zeros(N,M);

for m=1:M
for n=1:N
x=xx(1,m);
y=yy(1,n);
u=uu(n,m);
Q(n,m)=(3*y^2-y^3-x^3-x^2-6*y-6*x+4)*exp(-t);
end
end

end