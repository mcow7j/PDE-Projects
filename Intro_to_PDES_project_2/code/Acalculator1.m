function A=Acalculator1(uit,L,k);


[N,M]=size(uit);
  %dx = 1/(m-1);
  %dy = 1/(n-1);
  M2 = (M-1)*(M-1);
  N2 = ((N-1)*(N-1))/(L^2);
  NMtotal = (1-k*(M2 + N2)); 

  % initialization
  A = zeros(N,M);
  
  % iteration
    for n=2:N-1
      for m=2:M-1
          % calulates old values of u
	A(n,m) = 0.5*k*(N2*(uit(n+1,m)+uit(n-1,m)) + M2*(uit(n,m+1)+uit(n,m-1)))+NMtotal*uit(n,m);
      end
    end
    end
