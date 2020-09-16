function Q=QoldmakerUderiv(uu,vv,gamma,N,M)

Q=zeros(N,M);

for m=1:M
for n=1:N
u=uu(n,m);
v=vv(n,m);
Q(n,m)=gamma*(u*v-0.5);
end
end

end