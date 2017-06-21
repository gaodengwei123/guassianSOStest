% /*! @gpopsSolveCar.m
% *************************************************************************
% <PRE>
% file.name       : gpopsSolveCar.m
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
function [output, gpopsHistory] = gpopsSolveCar(guessTime,guessstate,guessControl,x0,xT)
% -----------------------
% Optimization Problem
% -----------------------
% The problem solved here is given as follows :
% Minimize
% t_f , J
% subject to the dynamic constraints
% dx / dt = v * cos ( psi )
% dy / dt = v * sin ( psi )
% dpsi / dt = ddpsi
% ddpsi / dt = u
% with the boundary conditions
% --------------------------------------------------
global CONSTANS
CONSTANS.PRMflag = 0;
t0 = 0;
xmin = min(guessstate(:,1));
xmax = max(guessstate(:,1));
ymin = min(guessstate(:,2));
ymax = max(guessstate(:,2));
phimin = -pi;
phimax = pi;
dphimin = -pi;
dphimax = pi;

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
limits.time.min = [0 t0];
limits.time.max = [0 1.5*max(guessTime)];
limits.control.min = -CONSTANS.controlSat;
limits.control.max = CONSTANS.controlSat;

limits.state.min = [ x0(1) xmin-10 xT(1)
                    x0(2) ymin-10 xT(2)
                    x0(3) phimin xT(3)
                    x0(4) dphimin xT(4)];
limits.state.max = [ x0(1) xmax+10 xT(1)
                    x0(2) ymax+10 xT(2)
                    x0(3) phimax xT(3)
                    x0(4) dphimax xT(4)];

limits.parameter.min = param_min ;
limits.parameter.max = param_max ;
limits.path.min = path_min ;
limits.path.max = path_max ;
limits.event.min = event_min ;
limits.event.max = event_max ;
limits.duration.min = duration_min ;
limits.duration.max = duration_max ;

guess.time = guessTime;
guess.state(: ,1) = guessstate(:,1);
guess.state(: ,2) = guessstate(:,2);
guess.state(: ,3) = guessstate(:,3);
guess.state(: ,4) = guessstate(:,4);
guess.control(: ,1) = guessControl;
guess.parameter = [];

setup.name = 'CarControl-Problem';
setup.funcs.cost = 'Car_controlCost';
setup.funcs.dae = 'Car_controlDae';
setup.limits = limits ;
setup.guess = guess ;
setup.derivatives = 'finite-difference';
setup.direction = 'increasing';
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-3;
setup.mesh.maxiterations = 10;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;
[output, gpopsHistory] = gpops(setup) ;