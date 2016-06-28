%-------------------------------------%
% BEGIN: function local_controlDae.m %
%-------------------------------------%
function dae = local_controlDae(sol)
global INPUTS
t = sol.time ;
x = sol.state(:,1);
y = sol.state(:,2);
psi = sol.state(:,3);
dpsi = sol.state(:,4);
u = sol.control;

xdot = INPUTS.speed*cos(psi);
ydot = INPUTS.speed*sin(psi);
psidot = dpsi;
psi2dot = u;
path = [];
for i = 1: size(INPUTS.obstacle,1)
    currentPath = (x- INPUTS.obstacle(i,1)).^2+(y- INPUTS.obstacle(i,2)).^2;
    path = [path currentPath];
end
dae = [xdot ydot psidot psi2dot path];


