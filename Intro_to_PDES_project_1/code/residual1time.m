   function res = residual1time(u,f,L,k1) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %
% Residuals: What is f-Au over the grid (1/m, 1/n)? %
%                                                   %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [N,M]=size(f);
  %dx = 1/(m-1);
  %dy = 1/(n-1);
  M2 = (M-1)*(M-1);
  N2 = ((N-1)*(N-1))/(L^2);
  NMtotal = (1+k1*(M2 + N2));
  res = zeros(size(u)); 
  
  for n=2:N-1
   for m=2:M-1
   res(n,m) = f(n,m)+u(n,m)*NMtotal -(u(n,m+1)+u(n,m-1))*M2*k1*0.5 -(u(n+1,m)+u(n-1,m))*N2*k1*0.5; 
   end
   end