   function res = residual1timedboundhomolin(u,Qold,f,L,k1,d) 

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
  NMtotal = (1+k1*d*(M2 + N2));
  res = zeros(size(u)); 
  

   res(1,1)=u(1,1)*(NMtotal-Qold(1,1)*k1)+f(1,1)-k1*d*(N2*u(2,1)+M2*u(1,2));
   res(N,M)=u(N,M)*(NMtotal-Qold(N,M)*k1)+f(N,M)-k1*d*(N2*u(N-1,M)+M2*u(N,M-1));
   res(1,M)=u(1,M)*(NMtotal-Qold(1,M)*k1)+f(1,M)-k1*d*(N2*u(2,M)+M2*u(1,M-1));
   res(N,1)=u(N,1)*(NMtotal-Qold(N,1)*k1)+f(N,1)-k1*d*(N2*u(N-1,1)+M2*u(N,2));
      
      
      for n=1
          for m=2:M-1
              res(n,m) = u(n,m)*(NMtotal-Qold(n,m)*k1)-0.5*k1*d*(N2*2*u(n+1,m)+M2*(u(n,m+1)+u(n,m-1)))+f(n,m);
          end
      end
      for n=N
          for m=2:M-1
              res(n,m) = u(n,m)*(NMtotal-Qold(n,m)*k1)-0.5*k1*d*(N2*2*u(n-1,m)+M2*(u(n,m+1)+u(n,m-1)))+f(n,m);
          end
      end
      
      for m=1
          for n=2:N-1
              res(n,m) = u(n,m)*(NMtotal-Qold(n,m)*k1)-0.5*k1*d*(N2*(u(n+1,m)+u(n-1,m))+M2*2*u(n,m+1))+f(n,m);
          end
      end
      for m=M
          for n=2:N-1
              res(n,m) = u(n,m)*(NMtotal-Qold(n,m)*k1)-0.5*k1*d*(N2*(u(n+1,m)+u(n-1,m))+M2*2*u(n,m-1))+f(n,m);
          end
      end
      
      
      
  
  for n=2:N-1
   for m=2:M-1
   res(n,m) = f(n,m)+u(n,m)*(NMtotal-Qold(n,m)*k1) -(u(n,m+1)+u(n,m-1))*M2*k1*0.5*d -(u(n+1,m)+u(n-1,m))*N2*k1*0.5*d; 
   end
  end
   end