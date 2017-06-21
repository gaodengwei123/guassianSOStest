%-------------------------------------%
% BEGIN: function local_controlDae.m %
%-------------------------------------%
function dae = pend_controlDae(sol)
% t = sol.time ;
x = sol.state(:,1);
y = sol.state(:,2);

u = sol.control;

xdot = y;
ydot = u-0.1*y-9.8*sin(x);

path = [];
dae = [xdot ydot path];


