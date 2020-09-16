function A=Acalculator1dhomo(uit,L);


[N,M]=size(uit);
  %dx = 1/(m-1);
  %dy = 1/(n-1);
%   M1=(M-1);
%   N1=(N-1)/(L);
  M2 = (M-1)*(M-1);
  N2 = ((N-1)*(N-1))/(L^2);
  NMtotal = (M2 + N2)*2; 

  % initialization
  A = zeros(N,M);
  
      A(1,1)=N2*2*uit(2,1)+M2*2*uit(1,2)-uit(1,1)*NMtotal;
      A(N,M)=N2*2*uit(N-1,M)+M2*2*uit(N,M-1)-uit(N,M)*NMtotal;
      A(1,M)=N2*2*uit(2,M)+M2*2*uit(1,M-1)-uit(1,M)*NMtotal;
      A(N,1)=N2*2*uit(N-1,1)+M2*2*uit(N,2)-uit(N,1)*NMtotal;
      
      
      for n=1
          for m=2:M-1
              A(n,m) = N2*2*uit(n+1,m)+M2*(uit(n,m+1)+uit(n,m-1))- uit(n,m)*NMtotal;
          end
      end
      for n=N
          for m=2:M-1
              A(n,m) = N2*2*uit(n-1,m)+M2*(uit(n,m+1)+uit(n,m-1))-uit(n,m)*NMtotal;
          end
      end
      
      for m=1
          for n=2:N-1
              A(n,m) = N2*(uit(n+1,m)+uit(n-1,m))+M2*2*uit(n,m+1)-uit(n,m)*NMtotal;
          end
      end
      for m=M
          for n=2:N-1
              A(n,m) = N2*(uit(n+1,m)+uit(n-1,m))+M2*2*uit(n,m-1)-uit(n,m)*NMtotal;
          end
      end
  
  % iteration
    for n=2:N-1
      for m=2:M-1
          % calulates old values of u
	A(n,m) = N2*(uit(n+1,m)+uit(n-1,m)) + M2*(uit(n,m+1)+uit(n,m-1))-NMtotal*uit(n,m);
      end
    end
    end
