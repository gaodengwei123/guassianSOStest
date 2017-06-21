% /*! @funnel_composition.m
% *************************************************************************
% <PRE>
% file.name       : funnel_composition.m
% related files   :
% function&ablity : 
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       : 
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/03/2     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output: rho is the factor for shrink
% *************************************************************************
function rho = funnel_composition(sys,V,OL_time,OL_state, P2, rho2)
x = msspoly('x',4);

%% compete the funnel V1
taus = linspace(0,OL_time(end),20);
sys.breaks = taus;
V.getbreak = taus;
statefun = spline(OL_time,OL_state');
sys.FunTraj = polyniminalTrajectory(statefun);
V.x0 = sys.FunTraj;
P1 = ppval(V.Spp,taus(end)); % the first P
% move the ROA to orign
V1 = x'*P1*x;
V2 = x'*P2*x;
% options.x0 = [0;0;0;0];
% xfun1 = getLevelSet(x,V1,options);
% xfun2 = getLevelSet(x,V2,options);
% plot3(xfun1(1,:),xfun1(2,:),repmat(.5,1,size(xfun1,2)),'b','LineWidth',1);
% plot3(xfun2(1,:),xfun2(2,:),repmat(.5,1,size(xfun2,2)),'r','LineWidth',1);

%% V1 is in V2 coordinate projection ==============================
% no_plot_dims=1:length(options.x0);  no_plot_dims(options.plotdims)=[];
x0 = [0;0;0;0];
cyclicIdx = [3 4];      % cyclic  coordinate
nonCyclicIdx = [1 2];   % noncyclic  coordinate
V1 = subs(V1,x(cyclicIdx),x0(cyclicIdx));
V2 = subs(V2,x(cyclicIdx),x0(cyclicIdx));
x = x(nonCyclicIdx);
P1 = double(0.5*diff(diff(V1,x)',x));
P2 = double(0.5*diff(diff(V2,x)',x));
V1 = x'*P1*x;
V2 = x'*P2*x;
%% ================================================================
tic
prog = spotsosprog;
prog = prog.withIndeterminate(x);
Lmonom = monomials(x,0:2);
% rho1 is free constraint which we want to know
[prog,rho] = prog.newFree(1);
[prog,L] = prog.newSOSPoly(Lmonom);
prog = prog.withSOS(L*(rho2 - V2)-(rho-V1));
option = spot_sdp_default_options();
sol = prog.minimize(-rho,@spot_mosek,option);
rho = double(sol.eval(rho));
toc



end







