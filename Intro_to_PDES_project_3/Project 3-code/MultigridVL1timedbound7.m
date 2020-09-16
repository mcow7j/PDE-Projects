function Uout = MultigridVL1timedbound7(Uin,f,L,k1,d,Bx,By); 

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
      Uout = GSL1timedbound7(Uin,f,10,L,k1,d,Bx,By); 
    else
% otherwise begin the cycle from fine to coarsest
%
%   Start by smoothing input with 10 GS iterations - could be too many
      Usmooth = GSL1timedbound7(Uin,f,10,L,k1,d,Bx,By); 
%
%   compute the residuals 
      res  = residual1timedbound7(Usmooth,f,L,k1,d,Bx,By); 

%     and restrict the residual to a coarser grid, half the size
      reshalf = restrict2(res); 

%     Now  solve the error equation A(error)=residulal on the next grid
%      Do this by calling this same routine recursively!
      err = MultigridVL1timedboundhomo(zeros(size(reshalf)),reshalf,L,k1,d); 

%     Now interpolate the course error onto finer grid and add to smoothed
      Usmooth = Usmooth + interpolate2(err); 

%     Finally, smooth out any new high-frequency error (post-smoothing) 
      Uout = GSL1timedbound7(Usmooth,f,10,L,k1,d,Bx,By);

% This completes a Multigrid V-cycle. If we want, we can call it again
    end

