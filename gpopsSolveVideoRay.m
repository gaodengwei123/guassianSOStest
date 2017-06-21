% /*! @gpopsSolveVideoRay.m
% *************************************************************************
% <PRE>
% file.name       : gpopsSolveVideoRay.m
% related files   :
% function&ablity : gpops to Solve car trajectory
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/8/29     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:
%
% output:
% *************************************************************************
% state
function [output, gpopsHistory] = gpopsSolveVideoRay(guessTime, guessstate, guessControl, x0, xT)
% -----------------------
% Optimization Problem
% -----------------------
% The problem solved here is given as follows :
% Minimize
% t_f , J
% subject to the dynamic constraints
% dx/dt = u*cos(psi) - v*sin(psi)
% dy/dt = u*sin(psi) + v*cos(psi)
% dpsi/dt = r
% du/dt = -r*v+Xu*u/m+f1
% dv/dt = r*u+Yv*v/m+f2
% dr/dt = Nr/Izz*r+f3
% with the boundary conditions

%===================================================%
%------------- Data Required by Problem ------------%
%===================================================%
global CONSTANS
auxdata.m = CONSTANS.m; 
auxdata.Xu = CONSTANS.Xu;
auxdata.Yv = CONSTANS.Yv; 
auxdata.Nr = CONSTANS.Nr; 
auxdata.Izz = CONSTANS.Izz;
auxdata.rG = CONSTANS.rG;

t0 = 0;
xmin = min(guessstate(:,1));
xmax = max(guessstate(:,1));
ymin = min(guessstate(:,2));
ymax = max(guessstate(:,2));
phimin = -pi/2;
phimax = pi/2;
dxmin = min(guessstate(:,4));
dxmax = max(guessstate(:,4));
dymin = min(guessstate(:,5));
dymax = max(guessstate(:,5));
dphimin = -CONSTANS.momentSat;
dphimax = CONSTANS.momentSat;

param_min = [];
param_max = [];
path_min = [];
for i = 1:size(CONSTANS.obstacle,1)
    path_min = [path_min;(CONSTANS.obstacleRadius_multiplier*CONSTANS.obstacleRadius)^2];
end
path_max = xmax^2*ones(1,size(CONSTANS.obstacle,1))';
event_min = [];
event_max = [];
duration_min = [];
duration_max = [];

% limits.meshPoints = [-1 +1];
% limits.nodesPerInterval = 10;
limits.time.min = [0 3];
limits.time.max = [0 3*max(guessTime)];
limits.control.min(1,:) = -CONSTANS.controlSat;
limits.control.max(1,:) = CONSTANS.controlSat;
% if the Y axis of videoRay have thrust or not
if CONSTANS.underactuation == 1
    limits.control.min(2,:) = -0.0;
    limits.control.max(2,:) =  0.0;
else
    limits.control.min(2,:) = -CONSTANS.controlSat;
    limits.control.max(2,:) =  CONSTANS.controlSat;
end
limits.control.min(3,:) = -CONSTANS.momentSat;
limits.control.max(3,:) =  CONSTANS.momentSat;

limits.state.min(1 ,:) = [ x0(1) xmin xT(1)];
limits.state.max(1 ,:) = [ x0(1) xmax xT(1)];
limits.state.min(2 ,:) = [ x0(2) ymin xT(2)];
limits.state.max(2 ,:) = [ x0(2) ymax xT(2)];
limits.state.min(3 ,:) = [ x0(3) -pi/2 xT(3)];
limits.state.max(3 ,:) = [ x0(3) pi/2 xT(3)];
limits.state.min(4 ,:) = [ x0(4) -CONSTANS.speed xT(4)];
limits.state.max(4 ,:) = [ x0(4) CONSTANS.speed xT(4)];
limits.state.min(5 ,:) = [ x0(5) -CONSTANS.speed xT(5)];
limits.state.max(5 ,:) = [ x0(5) CONSTANS.speed xT(5)];
limits.state.min(6 ,:) = [ x0(6) -CONSTANS.momentSat xT(6)];
limits.state.max(6 ,:) = [ x0(6) CONSTANS.momentSat xT(6)];

limits.parameter.min = param_min ;
limits.parameter.max = param_max ;
limits.path.min = path_min ;
limits.path.max = path_max ;
limits.event.min = event_min ;
limits.event.max = event_max ;
limits.duration.min = duration_min ;
limits.duration.max = duration_max ;

%===================================================%
%--------------- Set Up Initial Guess --------------%
%===================================================%
guess.time = guessTime;
guess.state(: ,1) = guessstate(:,1);
guess.state(: ,2) = guessstate(:,2);
guess.state(: ,3) = guessstate(:,3);
guess.state(: ,4) = guessstate(:,4);
guess.state(: ,5) = guessstate(:,5);
guess.state(: ,6) = guessstate(:,6);
guess.control(: ,1) = guessControl(:,1);
guess.control(: ,2) = guessControl(:,2);
guess.control(: ,3) = guessControl(:,3);
guess.parameter = [];

% ===================================================%
% --------------- Set Up for bounds -----------------%
% ===================================================%
% bounds.finalstate.lower = 
% ===================================================%
% --------------- Set Up for Solver -----------------%
% ===================================================%
setup.name = 'VideoRay-Problem';
setup.funcs.cost = 'VideoRay_controlCost';
setup.funcs.dae = 'VideoRay_controlDae';
setup.limits = limits ;
setup.guess = guess ;
setup.auxdata = auxdata;
% setup.bounds = bounds;
setup.derivatives = 'finite-difference';
setup.direction = 'increasing';
setup.autoscale = 'on';
setup.mesh.tolerance = 1e-3;
setup.mesh.maxiterations = 10;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;

% ===================================================%
% ------ Solve Problem and Extract Solution ---------%
% ===================================================%
[output, gpopsHistory] = gpops(setup ) ;