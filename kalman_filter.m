% /*! @kalman_filter.m
% *************************************************************************
% <PRE>
% file.name       : kalman_filter.m
% related files   :     
% function&ablity : 
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         : 
% --------------------------------------------------------------------------------
% record of modify : 
% date          version     name         content 
% 2015/10/29    1.00                     build
% 2016/3/2        2.00                     rebuild 3D
% </PRE>
% ********************************************************************************
% 
% * right(c)
% 
% *************************************************************************
% input : x10:prior state estimate;real_Z: the real observed value
%        Q: process noise covariance; R:measurement noise covariance;
%        P0:prior variance;U:force and moment;temp:simulantion step.i:step number 

% output:dert_x: EKF update; x1: the state of next step;P1:posteriori estimate covariance
% *************************************************************************
function [dert_x,x1,P1] = kalman_filter(x10,C,D,inv_MS,J_n,ROV,real_Z,Q,R,P0,temp,i)

% observation
Obs =0; % measurement mark
Hk1 = [eye(6) zeros(6)]; % measurement matrix
if rem(i-1,1/temp) == 0  % measurement update(every sec)
    Zk1 = real_Z - [x10(1) x10(2) x10(3) x10(4) x10(5) x10(6)]';
    Obs = 1;
else
    Zk1 = [0;0;0;0;0;0];
end

% error state
Xk = zeros(12,1); % reset the state error
G = eye(12);
Qt=G*Q*G';

%% linearliation
W = ROV.m*9.8;
B = W;
phi = x10(4);
sita = x10(5);
g15 = (W-B)*(1-1/6*sita^2+1/120*sita^4);
g24 = -(W-B)*cos(sita)*(1-1/6*phi^2+1/120*phi^4);
g35 = -(W-B)*(-sita/2+1/24*sita^3)*cos(phi);
g44 = (ROV.zG*W-ROV.zB*B)*cos(sita)*(1-1/6*phi^2+1/120*phi^4);
g45 = -(ROV.yG*W-ROV.yB*B)*(-sita/2+1/24*sita^3)*cos(phi);
g54 = (ROV.xG*W-ROV.xB*B)*cos(sita)*(-phi/2+1/24*phi^3);
g55 = (ROV.zG*W-ROV.zB*B)*(1-1/6*sita^2+1/120*sita^4);
g64 = -(ROV.xG*W-ROV.xB*B)*cos(sita)*(1-1/6*phi^2+1/120*phi^4);
g65 = -(ROV.yG*W-ROV.yB*B)*(1-1/6*sita^2+1/120*sita^4);
Gn = [zeros(6,3) [0 g15 0;g24 0 0;0 g35 0;g44 g45 0;g54 g55 0;g64 g65 0]];

Ac = [zeros(6)     J_n
     -inv_MS*Gn   -inv_MS*(C+D)];

[dert_x,P1] = kfilter(Ac, Xk, Qt, Hk1, Zk1, R, P0, temp, 12, Obs);
x1 = x10+dert_x;






