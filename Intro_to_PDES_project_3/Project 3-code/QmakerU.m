function Q=QmakerU(uu,vv,gamma,a,N,M)

Q=zeros(N,M);

for m=1:M
for n=1:N
u=uu(n,m);
v=vv(n,m);
Q(n,m)=gamma*(a-u+(u^2)*v);
end
end

end