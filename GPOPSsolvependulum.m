clear all
close all
clc
dbstop if error
%%
Pendsys = Pendulumdyn();
A = @(x)([0  1; -9.8*cos(x(1)) -0.1]);
B = [0;1];
Q = diag([1,1]);R=1;
range = [-pi -2.5*pi pi 2.5*pi];  % [x1l x2l x1u x2u]


sys.A = A;
sys.B = @(x)[0;1];
sys.f = Pendsys;
sys.Q = Q;
sys.R = R;
sys.uppercontrol = 3;
Num = 5000;

xT = [pi;0];
x0 = [0;0];
Time = [0;30];
Point = [x0 xT];
traj = traj_gpopsSolvePendulumdyn(sys,Time,Point);


