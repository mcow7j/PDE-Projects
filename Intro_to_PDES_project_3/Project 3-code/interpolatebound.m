function fine=interpolatebound(cor)

[m,n]=size(cor);
r=n-1;

fine=zeros(m,2*r+1);

fine(1:m,1:2:2*r+1)=cor(1:m,1:n);
fine(1:m,2:2:2*r)=0.5*cor(1:m,2:n)+0.5*cor(1:m,1:r);




end