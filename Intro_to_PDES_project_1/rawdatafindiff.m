function [a,b]=rawdatafindiff(N,k,iter,eta)
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

for j=1:iter
    for i=2:q
       a(j+1,i)=(a(j,i-1:i+1)*c(1:3,i))+k*b(j,i);
       b(j+1,i)=(b(j,i-1:i+1)*c(1:3,i))-k*a(j,i);
    end
    for i=q+1:N
        a(j+1,i)=(a(j,i-1:i+1)*c(1:3,i));
        b(j+1,i)=(b(j,i-1:i+1)*c(1:3,i));
    end
end

end