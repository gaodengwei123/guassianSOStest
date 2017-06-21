function [output, gpopsHistory] = PRM_gpopsSolveCar(sys,x0,xT)
global CONSTANS
CONSTANS = sys.INPUTS;
CONSTANS.state_goal = xT;
CONSTANS.PRMflag = 1;
t0 = 0;
xmin = min(x0(1),xT(1));
xmax = max(x0(1),xT(1));
ymin = min(x0(2),xT(2));
ymax = max(x0(2),xT(2));
phimin = -pi;
phimax = pi;
dphimin = -pi;
dphimax = pi;

param_min = [];
param_max = [];
path_min = [];
path_max = [];

event_min = [];
event_max = [];
duration_min = [];
duration_max = [];

% limits.meshPoints = [-1 +1];
% limits.nodesPerInterval = 10;
limits.time.min = [0 t0];
limits.time.max = [0 10];
limits.control.min = -pi;
limits.control.max = pi;

limits.state.min (1 ,:) = [ x0(1) xmin-10 xT(1)];
limits.state.max (1 ,:) = [ x0(1) xmax+10 xT(1)];
limits.state.min (2 ,:) = [ x0(2) ymin-10 xT(2)];
limits.state.max (2 ,:) = [ x0(2) ymax+10 xT(2)];
limits.state.min (3 ,:) = [ x0(3) phimin xT(3)];
limits.state.max (3 ,:) = [ x0(3) phimax xT(3)];
limits.state.min (4 ,:) = [ x0(4) dphimin xT(4)];
limits.state.max (4 ,:) = [ x0(4) dphimax xT(4)];

limits.parameter.min = param_min ;
limits.parameter.max = param_max ;
limits.path.min = path_min ;
limits.path.max = path_max ;
limits.event .min = event_min ;
limits.event.max = event_max ;
limits.duration.min = duration_min ;
limits.duration.max = duration_max ;

guess.time = [0; 3];
guess.state(: ,1) = [x0(1); xT(1)];
guess.state(: ,2) = [x0(2); xT(2)];
guess.state(: ,3) = [x0(3); xT(3)];
guess.state(: ,4) = [x0(4); xT(4)];
guess.control(: ,1) = [0;0];
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
[output, gpopsHistory] = gpops(setup);


end




