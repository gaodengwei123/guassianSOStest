clear all
close all
clc
sys = CarDynamic();

figure(1)
hold on
point = [0 0 0 0 
         1 2 2 0];
Time = [0 1];
waypoints  = addwaypoints(sys, point,5);
output = trajLibrary_gpopsSolveCar(sys, Time, point, waypoints);
plot(output.solutionPlot.state(:,1),output.solutionPlot.state(:,2),'r')

% check open loop
time = output.solutionPlot.time;
control = output.solutionPlot.control;
Spp = spline(time,control);
u = @(t) ppval(Spp,t);
odefun = @(t,x)sys.dynamics(t,x,u(t));
sol = ode45(odefun,[time(1),time(end)],[0;0;0;0]);
ts = linspace(time(1),time(end),100);
state = deval(sol,ts);
plot(state(1,:),state(2,:),'k')
% figure
% plot(control)