function [Mayer, Lagrange]= Car_controlCost_Phases(sol)
global CONSTANS
tf = sol.terminal.time;
t = sol.time ;

x = sol.state(: ,1)-CONSTANS.state_goal_phase(sol.phase,1);
y = sol.state(: ,2)-CONSTANS.state_goal_phase(sol.phase,2);
psi = sol.state(: ,3)-CONSTANS.state_goal_phase(sol.phase,3);
dpsi = sol.state(: ,4)-CONSTANS.state_goal_phase(sol.phase,4);
u = sol.control;
Mayer = tf *10^1;
Q = CONSTANS.Q(1);
R = CONSTANS.R(1);
Lagrange = Q(1,1)*x.^2+Q(2,2)*y.^2+Q(3,3)*psi.^2+Q(4,4)*dpsi.^2+R(1,1)*u.^2;
end
