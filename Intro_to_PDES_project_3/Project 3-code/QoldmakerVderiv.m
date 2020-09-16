function Q=QoldmakerVderiv(uu,gamma,N,M)

Q=zeros(N,M);

for m=1:M
for n=1:N
u=uu(n,m);
Q(n,m)=(-0.5)*gamma*(u*u);
end
end

end