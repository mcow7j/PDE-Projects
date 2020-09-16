function u=uexact(N,M,L,t);

xx=linspace(0,1,M);
yy=linspace(0,L,N);
u=zeros(N,M);



for n=1:N
    for m=1:M
        x=xx(1,m);
        y=yy(1,n);
        u(n,m)=(x-1)*x*y*(y-L)*exp(x-y)*exp(-t);
    end
end

end

