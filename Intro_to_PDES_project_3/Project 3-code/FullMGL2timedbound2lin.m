     function u = FullMGL2timedbound2lin(uin,Qold,f,L,k1,d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% full geometric multigrid                                        %
%         Solves delsqr(u)=f on a grid given by size of f         %                                                        %
% (homogeneous Dirichlet boundary conditions)                     %
% (dimension N has to be of the form 2^klevel + 1)                %
%                                                                 %
% f:    right-hand side (n x m)-matrix                            %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n=size(f,1);
    m=size(f,2);
 
%   determine maximum number of grid levels
    klevel=round(log(min(n,m)-1)/log(2));

%    As we will solve on all grids, we need to calculate the
%    restrictions of the right-hand side to all grid levels
%    ff(1) is finest grid,clea ff(klevel) coarsest
    ff = cell(klevel,1);
    ff{1} = f;
    gg = cell(klevel,1);
    gg{1} = uin;
    for k=2:klevel
      ff{k} = restrict2(ff{k-1});
      gg{k} = restrict2(gg{k-1});
    end

    
%   solve the equation at the coarsest level
    f0 = ff{klevel};
%    u0 = initu(zeros(size(f0,1),size(f0,2)));
 %   u0=zeros(size(f0,1),size(f0,2));
 u0=gg{klevel};
    u0 = GSL1timedboundhomolin(u0,Qold,f0,10,L,k1,d);

%   loop over all higher levels (coarser grids) using V-cycles
    for k=klevel-1:-1:1

%     interpolate solution to next-higher grid level
%      u1 = initu(prolong(u0));
      u1=interpolate2(u0);
      f1 = ff{k};

%     perform one multigrid V-cycle
      u0=MultigridVL1timedboundhomolin(u1,Qold,f1,L,k1,d);

%  figure(k)
%  [nn,mm] = size(u0);
%  xx = linspace(0,1,mm);
%  yy = linspace(0,1,nn);
% mesh(xx,yy,u0)
% xlabel('x','FontSize',18)
% ylabel('y','FontSize',18)
% zlabel('u','FontSize',18)
% grid off
% lighting phong
% camlight headlight
% camlight right
% pause(1)

    end

%   output solution
    u=u0;

