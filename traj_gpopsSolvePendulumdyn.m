function traj = traj_gpopsSolvePendulumdyn(sys,Time,Point)
global CONSTANS
CONSTANS.Q = sys.Q;
CONSTANS.R = sys.R;
x0 = Point(:,1);
xT = Point(:,2);
CONSTANS.state_goal = xT;
controlBound = 3;
xmin = -pi;
xmax = pi;
ymin = -2*pi;
ymax = 2*pi;

limits.time.min = [0 Time(1)];
limits.time.max = [0 Time(2)];
limits.control.min = -controlBound;
limits.control.max = controlBound;

limits.state.min (1 ,:) = [ 0 xmin pi];
limits.state.max (1 ,:) = [ 0 xmax pi];
limits.state.min (2 ,:) = [ 0 ymin 0];
limits.state.max (2 ,:) = [ 0 ymax 0];


limits.parameter.min = [] ;
limits.parameter.max = [] ;
limits.path.min = [] ;
limits.path.max = [] ;
limits.event.min = [] ;
limits.event.max = [] ;
limits.duration.min = [] ;
limits.duration.max = [] ;
guess.time = [0; Time(2)];
guess.state(: ,1) = [x0(1); xT(1)];
guess.state(: ,2) = [x0(2); xT(2)];
guess.control(: ,1) = [0;0];
    
    
setup.name = 'pendulumdyn_name';
setup.funcs.cost = 'pend_controlCost';
setup.funcs.dae = 'pend_controlDae';
setup.limits = limits ;
setup.guess = guess ;
setup.derivatives = 'finite-difference';
setup.checkDerivatives = 0;
setup.direction = 'increasing';
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-2;
setup.mesh.maxiterations = 4;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;
tic
[output, ~] = gpops(setup);
toc

traj.state = output.solutionPlot.state;
traj.control = output.solutionPlot.control;
traj.time = output.solutionPlot.time;
traj.cost = output.cost;
end




