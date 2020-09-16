function Uout = MultigridVL1timedboundhomolin(Uin,Qold,f,L,k1,d); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multigrid V-cycle                                                %
%                                                                  %
% Calls itself recursively                                         %
% Homogeneous Dirichlet boundary conditions                        %
% square grid. Dimension N=min(m,n) has to be of the form 2^k + 1  %
%                                                                  %
% Uin:  initial estimate. Right-hand side (n x m)-matrix           %
%                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = size(f,1); 
    m = size(f,2);  

% if we are at the coarsest level take 10 GS iterations - should be enough
    if ((n==3)|(m==3)) 
      Uout = GSL1timedboundhomolin(Uin,Qold,f,4,L,k1,d); 
    else
% otherwise begin the cycle from fine to coarsest
%
%   Start by smoothing input with 10 GS iterations - could be too many
      Usmooth = GSL1timedboundhomolin(Uin,Qold,f,4,L,k1,d); 
%
%   compute the residuals 
      res  = residual1timedboundhomolin(Usmooth,Qold,f,L,k1,d); 

%     and restrict the residual to a coarser grid, half the size
      reshalf = restrict2(res); 

%     Now  solve the error equation A(error)=residulal on the next grid
%      Do this by calling this same routine recursively!
      err = MultigridVL1timedboundhomolin(zeros(size(reshalf)),Qold,reshalf,L,k1,d); 

%     Now interpolate the course error onto finer grid and add to smoothed
      Usmooth = Usmooth + interpolate2(err); 

%     Finally, smooth out any new high-frequency error (post-smoothing) 
      Uout = GSL1timedboundhomolin(Usmooth,Qold,f,4,L,k1,d);

% This completes a Multigrid V-cycle. If we want, we can call it again
    end

