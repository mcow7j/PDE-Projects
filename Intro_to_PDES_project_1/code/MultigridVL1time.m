function Uout = MultigridVL1time(Uin,f,L,k1); 

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
      Uout = GSL1time(Uin,f,10,L,k1); 
    else
% otherwise begin the cycle from fine to coarsest
%
%   Start by smoothing input with 10 GS iterations - could be too many
      Usmooth = GSL1time(Uin,f,10,L,k1); 
%
%   compute the residuals 
      res  = residual1time(Usmooth,f,L,k1); 

%     and restrict the residual to a coarser grid, half the size
      reshalf = restrict1(res); 

%     Now  solve the error equation A(error)=residulal on the next grid
%      Do this by calling this same routine recursively!
      err = MultigridVL1time(zeros(size(reshalf)),reshalf,L,k1); 

%     Now interpolate the course error onto finer grid and add to smoothed
      Usmooth = Usmooth + interpolate1(err); 

%     Finally, smooth out any new high-frequency error (post-smoothing) 
      Uout = GSL1time(Usmooth,f,10,L,k1);

% This completes a Multigrid V-cycle. If we want, we can call it again
    end

