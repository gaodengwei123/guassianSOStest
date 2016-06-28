function [Mayer , Lagrange]= local_controlCost(sol)
global INPUTS
tf = sol.terminal.time ;
t = sol.time ;
x = sol.state(: ,1)-INPUTS.state_goal(1);
y = sol.state(: ,2)-INPUTS.state_goal(2);
psi = sol.state(: ,3)-INPUTS.state_goal(3);
dpsi = sol.state(: ,4)-INPUTS.state_goal(4);
u = sol.control;
Mayer = tf *10^1;
Q = INPUTS.Q(1);
R = INPUTS.R(1);
Lagrange = Q(1,1)*x.^2+Q(2,2)*y.^2+Q(3,3)*psi.^2+Q(4,4)*dpsi.^2+R(1,1)*u.^2;
end
