% /*! @gao_EKF.m
% *************************************************************************
% <PRE>
% file.name       : gao_EKF.m
% related files   :
% function&ablity :
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/1/24    1.00                     build

% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input :

% output:
% *************************************************************************
function [dert_x,PP] = gao_EKF(sys,t,P0,xdot,Zk1,Obs)
Ac =  sys.A(t) + sys.B(t)*sys.K(t);
B = sys.B(t);
K = sys.K(t);
INPUTS = sys.INPUTS;
Q = INPUTS.Qk;
R = INPUTS.Rk;
m = length(INPUTS.noise_v);
n = length(xdot);
Qt = Q;
temp
Hk1 = [eye(m),zeros(m,sys.getNumInput)];    % linearlize the measurement function
Xk = zeros(n,1);

[dert_x,P1] = kfilter(Ac, Xk, Qt, Hk1, Zk1, R, P0, temp, n, Obs);

% error dyanmic
Bk = B*temp;
Ak = eye(n)+Ac*temp+Bk*K;
beta1 = sqrt(trace(Q))/(sqrt(trace(-Bk*K*P0))+sqrt(trace(Q)));
Qmao = -Bk*K*P0/(1-beta1)+Q/beta1;
beta = sqrt(trace(Qmao))/(sqrt(trace(Ak*P0*Ak'))+sqrt(trace(Qmao)));
PP = double(Ak*P0/(1-beta)*Ak'+Qmao/beta);






