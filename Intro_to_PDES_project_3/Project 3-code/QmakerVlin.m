function Q=QmakerVlin(uu,vv,Qold,gamma,b,N,M)

Q=zeros(N,M);

for m=1:M
for n=1:N
u=uu(n,m);
v=vv(n,m);
Q(n,m)=gamma*(b-(u^2)*v);
end
end
Q=Q-Qold.*vv;
end