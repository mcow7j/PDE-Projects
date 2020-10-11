function u = GSL1time(uinit,f,Niter,L,k1);

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
  M2 = (M-1)*(M-1);
  N2 = ((N-1)*(N-1))/(L^2);
  NMtotal = (1+k1*(M2 + N2)); 

  % initialization
  unew = uinit;
  
  % iteration
  for k=1:Niter
    for n=2:N-1
      for m=2:M-1
          % For Gauss-Seidel, overwrite u and use calculated new values
	unew(n,m) = (0.5*k1*(N2*(unew(n+1,m)+unew(n-1,m)) + M2*(unew(n,m+1)+unew(n,m-1))) - f(n,m))/NMtotal;
      end
    end
  end

  u = unew;