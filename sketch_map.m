% sketch_map
clc
clear all
close all
figure
hold on
x=[1 2 2 1];
y=[2 2 4 4];
fill(x,y,'k')
sys = CarDynamic();
INPUTS = sys.INPUTS;
InitialState = [0 0 0 0]'; % real state: add Gaussian noise at initial state
invtime  = [0 2];
point = [0 0 pi/3 0
         4 4 pi/2 0];
traj = trajLibrary_gpopsSolveCar(sys, invtime, point, []);
ts = traj.time;
Time = [ts(1) ts(end)];
control = spline(ts,traj.control);
u = @(t)ppval(control,t)+rand;
Stochsys = StochasticSystem(sys,100*INPUTS.Qk,INPUTS.Rk);
odefun = @(t,x)Stochsys.dynamics(t,x,u(t),(2*rand(4,1)-1));
sol = ode45(odefun,Time,point(1,:));
% draw real trajectory
taus = linspace(Time(1),Time(2),100);
state = deval(sol,taus);
plot(state(1,:),state(2,:),'b')





point = [0 0 pi/3 0
         1.5 2 pi/6 0];
traj = trajLibrary_gpopsSolveCar(sys, invtime, point, []);
ts = traj.time;
Time = [ts(1) ts(end)];
control = spline(ts,traj.control);
u = @(t)ppval(control,t)+rand;
Stochsys = StochasticSystem(sys,100*INPUTS.Qk,INPUTS.Rk);
odefun = @(t,x)Stochsys.dynamics(t,x,u(t),(2*rand(4,1)-1));
sol = ode45(odefun,Time,point(1,:));
% draw real trajectory
taus = linspace(Time(1),Time(2),100);
state = deval(sol,taus);
plot(state(1,:),state(2,:),'b')



point = [0 0 pi/3 0
         3 0 -pi/3 0];
traj = trajLibrary_gpopsSolveCar(sys, invtime, point, []);
ts = traj.time;
Time = [ts(1) ts(end)];
control = spline(ts,traj.control);
u = @(t)ppval(control,t)+rand;
Stochsys = StochasticSystem(sys,100*INPUTS.Qk,INPUTS.Rk);
odefun = @(t,x)Stochsys.dynamics(t,x,u(t),(2*rand(4,1)-1));
sol = ode45(odefun,Time,point(1,:));
% draw real trajectory
taus = linspace(Time(1),Time(2),100);
state = deval(sol,taus);
plot(state(1,:),state(2,:),'b')



