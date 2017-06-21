function [Mayer, Lagrange]= pend_controlCost(sol)
global CONSTANS
% tf = sol.terminal.time;
% t0 = sol.initial.time;
% t = sol.time;
xf = sol.initial.state;
x = sol.state;%-CONSTANS.state_goal;
u = sol.control;

Q = CONSTANS.Q;
R = CONSTANS.R;
Mayer = dot(xf,Q*xf);
Lagrange = dot(x,x*Q',2)+dot(u,u*R',2);

end
