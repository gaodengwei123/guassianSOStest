function [Mayer, Lagrange]= Car_controlCost(sol)
global CONSTANS
tf = sol.terminal.time;
t0 = sol.initial.time;
t = sol.time;
xf = sol.initial.state;
% x = sol.state(: ,1)-CONSTANS.state_goal(1);
% y = sol.state(: ,2)-CONSTANS.state_goal(2);
% psi = sol.state(: ,3)-CONSTANS.state_goal(3);
% dpsi = sol.state(: ,4)-CONSTANS.state_goal(4);
x = sol.state;%-CONSTANS.state_goal;

u = sol.control;

Q = CONSTANS.Q(1);
R = CONSTANS.R(1);
Mayer = dot(xf,Q*xf);%t'*t;%tf *10^1;
Lagrange = dot(x,x*Q',2)+dot(u,u*R',2);
end
