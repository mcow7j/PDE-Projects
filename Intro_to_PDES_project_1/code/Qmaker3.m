function Q=Qmaker3(uu,N,M,L,t,A)

xx=linspace(0,1,M);
yy=linspace(0,L,N);
Q=zeros(N,M);

for m=1:M
for n=1:N
x=xx(1,m);
y=yy(1,n);
u=uu(n,m);
Q(n,m)=A*u*(1-u)*(2-u);
end
end

end