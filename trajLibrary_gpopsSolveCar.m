function traj = trajLibrary_gpopsSolveCar(sys,Time,Point,waypoints)
global CONSTANS
CONSTANS = sys.INPUTS;
x0 = Point(1,:);
xT = Point(2,:);
CONSTANS.state_goal = xT;
CONSTANS.PRMflag = 1;
controlBound = 1000;

xmin = min(x0(1),xT(1));
xmax = max(x0(1),xT(1));
ymin = min(x0(2),xT(2));
ymax = max(x0(2),xT(2));
phimin = -2*pi;
phimax = 2*pi;
dphimin = -controlBound;
dphimax = controlBound;

limits.time.min = [0 Time(1)];
limits.time.max = [0 Time(2)];
limits.control.min = -controlBound;
limits.control.max = controlBound;

limits.state.min (1 ,:) = [ x0(1) xmin-5 xT(1)];
limits.state.max (1 ,:) = [ x0(1) xmax+5 xT(1)];
limits.state.min (2 ,:) = [ x0(2) ymin-5 xT(2)];
limits.state.max (2 ,:) = [ x0(2) ymax+5 xT(2)];
limits.state.min (3 ,:) = [ x0(3) phimin xT(3)];
limits.state.max (3 ,:) = [ x0(3) phimax xT(3)];
limits.state.min (4 ,:) = [ x0(4) dphimin xT(4)];
limits.state.max (4 ,:) = [ x0(4) dphimax xT(4)];

limits.parameter.min = [] ;
limits.parameter.max = [] ;
limits.path.min = [] ;
limits.path.max = [] ;
limits.event.min = [] ;
limits.event.max = [] ;
limits.duration.min = [] ;
limits.duration.max = [] ;
if ~isempty(waypoints)
    guess.time = [0; waypoints(:,1);Time(2)];
    guess.state(: ,1) = [x0(1); waypoints(:,2); xT(1)];
    guess.state(: ,2) = [x0(2); waypoints(:,3); xT(2)];
    guess.state(: ,3) = [x0(3); waypoints(:,4); xT(3)];
    guess.state(: ,4) = [x0(4); waypoints(:,5); xT(4)];
    guess.control(: ,1) = [0;waypoints(:,6);0];
else
    guess.time = [0; Time(2)];
    guess.state(: ,1) = [x0(1); xT(1)];
    guess.state(: ,2) = [x0(2); xT(2)];
    guess.state(: ,3) = [x0(3); xT(3)];
    guess.state(: ,4) = [x0(4); xT(4)];
    guess.control(: ,1) = [0;0];
end
guess.parameter = [];

setup.name = 'CarControl-Problem';
setup.funcs.cost = 'Car_controlCost';
setup.funcs.dae = 'Car_controlDae';
setup.limits = limits ;
setup.guess = guess ;
setup.derivatives = 'finite-difference';
setup.checkDerivatives = 0;
setup.direction = 'increasing';
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-6;
setup.mesh.maxiterations = 10;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;
[output, ~] = gpops(setup);

traj.state = output.solutionPlot.state;
traj.control = output.solutionPlot.control;
traj.time = output.solutionPlot.time;
traj.cost = output.cost;
end




