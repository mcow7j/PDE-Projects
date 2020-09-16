function u=uexact3(N,M,L,t);

xx=linspace(0,1,M);
yy=linspace(0,L,N);
u=zeros(N,M);



for n=1:N
    for m=1:M
        x=xx(1,m);
        y=yy(1,n);
        u(n,m)=(y^3-3*y^2+x^3+x^2)*exp(-t);
    end
end

end

