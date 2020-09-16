function u = GSL1timedbound8(uinit,Qold,f,Niter,L,k1,d,Bx,By);


%------------------------------------------%
% Takes Niter Gauss-Seidel iterations for  %
% the discrete Poisson equation in form    %
% Au=b, on unit square starting from uinit.%
%------------------------------------------%


%L rectangle in y,m direction
%k is timestep
%f is o.5kQ-A

  [N,M]=size(f);
  %dx = 1/(m-1);
  %dy = 1/(n-1);
  M1=(M-1);
  N1=(N-1)/L;
  M2 = (M-1)*(M-1);
  N2 = ((N-1)*(N-1))/(L^2);
  NMtotal = (1+k1*d*(M2 + N2)); 

  % initialization
  unew = uinit;
  
  % iteration
  for k=1:Niter
      
      
      
      for n=1
          for m=2:M-1
              unew(n,m) = (0.5*k1*d*(Bx(1,m)*N1*2+N2*2*unew(n+1,m)+M2*(unew(n,m+1)+unew(n,m-1)))-f(n,m))/(NMtotal-k1*Qold(n,m));
          end
      end
      for n=N
          for m=2:M-1
              unew(n,m) = (0.5*k1*d*(Bx(2,m)*N1*2+N2*2*unew(n-1,m)+M2*(unew(n,m+1)+unew(n,m-1)))-f(n,m))/(NMtotal-k1*Qold(n,m));
          end
      end
      
      for m=1
          for n=2:N-1
              unew(n,m) = (0.5*k1*d*(N2*(unew(n+1,m)+unew(n-1,m))+By(1,n)*M1*2+M2*2*unew(n,m+1))-f(n,m))/(NMtotal-k1*Qold(n,m));
          end
      end
      for m=M
          for n=2:N-1
              unew(n,m) = (0.5*k1*d*(N2*(unew(n+1,m)+unew(n-1,m))+By(2,n)*M1*2+M2*2*unew(n,m-1))-f(n,m))/(NMtotal-k1*Qold(n,m));
          end
      end
      
       unew(1,1)=(k1*d*(Bx(1,1)*N1+N2*unew(2,1)+By(1,1)*M1+M2*unew(1,2))-f(1,1))/(NMtotal-k1*Qold(1,1));
      unew(N,M)=(k1*d*(Bx(2,M)*N1+N2*unew(N-1,M)+By(2,N)*M1+M2*unew(N,M-1))-f(N,M))/(NMtotal-k1*Qold(N,M));
      unew(1,M)=(k1*d*(Bx(1,M)*N1+N2*unew(2,M)+By(2,1)*M1+M2*unew(1,M-1))-f(1,M))/(NMtotal-k1*Qold(1,M));
      unew(N,1)=(k1*d*(Bx(2,1)*N1+N2*unew(N-1,1)+By(1,N)*M1+M2*unew(N,2))-f(N,1))/(NMtotal-k1*Qold(N,1));
      
      
      
    for n=2:N-1
      for m=2:M-1
          % For Gauss-Seidel, overwrite u and use calculated new values
	unew(n,m) = (0.5*k1*d*(N2*(unew(n+1,m)+unew(n-1,m)) + M2*(unew(n,m+1)+unew(n,m-1))) - f(n,m))/(NMtotal-k1*Qold(n,m));
      end
    end
  end

  u = unew;