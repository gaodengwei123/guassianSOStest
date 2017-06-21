%-------------------------------------%
% BEGIN: function local_controlDae.m %
%-------------------------------------%
function dae = Car_controlDae(sol)
global CONSTANS
t = sol.time ;
x = sol.state(:,1);
y = sol.state(:,2);
psi = sol.state(:,3);
dpsi = sol.state(:,4);
u = sol.control;

xdot = CONSTANS.speed*cos(psi);
ydot = CONSTANS.speed*sin(psi);
psidot = dpsi;
psi2dot = u;
path = [];
if CONSTANS.PRMflag == 0
    for i = 1: size(CONSTANS.obstacle,1)
        currentPath = (x- CONSTANS.obstacle(i,1)).^2+(y- CONSTANS.obstacle(i,2)).^2;
        path = [path currentPath];
    end
end
dae = [xdot ydot psidot psi2dot path];


