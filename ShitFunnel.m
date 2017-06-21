% /*! @ShitFunnel.m
% *************************************************************************
% <PRE>
% file.name       : ShitFunnel.m
% related files   :
% function&ablity :
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/03/6     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output:
% *************************************************************************
function Ftuple = ShitFunnel(funnel, traj, node)
% t0 = node.t0;
x_now = node.pos;
% shift funnel along cycile coordinate
deltaX = [x_now - traj.state(1,1:2), 0, 0];
OL_time = traj.time;
OL_state = traj.state;
shift_state = OL_state+repmat(deltaX,size(OL_state,1),1);


% shift trajectory
statefun = spline(OL_time,shift_state');
Traj = polyniminalTrajectory(statefun);
% shift funnel
funnel.V = funnel.V.inFrame(Traj);  % 0 ~ t
N = length(funnel.V.getbreak);
% time = linspace(Traj.tspan(1),Traj.tspan(end),N);
% funnel.V = funnel.V.updatebreaks(time);
% shift feedback gain
K = funnel.gain;
A = funnel.A;
B = funnel.B;
u0 = funnel.u0;
Ftuple = FunnelTuple();
Ftuple = Ftuple.setFunnelTuple(funnel.V, Traj, A, B, K, u0, traj.cost);
% plot(Ftuple);
end









