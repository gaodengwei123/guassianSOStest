% /*! @Funnel_filter.m
% *************************************************************************
% <PRE>
% file.name       : Funnel_filter.m
% related files   :
% function&ablity :
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/3/3      3.00                     kalman filter
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input :

% output:
% *************************************************************************
function [SigmaF,LambdaF] = Funnel_filter(sys,Time, Sigma0, Lambda0)
xtraj = sys.FunTraj;
INPUTS = sys.INPUTS;

n = sys.getNumStates;
Sigma{1} = Sigma0;           % convariance of estimate state x_hat
Lambda{1} = Lambda0;          % convariance of x_tilde = x - x_bar
Hk = eye(size(INPUTS.Rk,1),n);  % measurement Jacobian matrix

SampleNum = length(Time);
for i = 1:SampleNum-1
    A = sys.A(Time(i));
    BK = sys.B(Time(i))*sys.K(Time(i));
    % measurement condition: near obstacle <15

    Obs = NearSensor(sys.INPUTS.field,xtraj.eval(Time(i)));
    Tkf = Time(i+1)-Time(i);
    [Sigma{i+1},Lambda{i+1}] = kfilter2(A, BK, INPUTS.Qk, Hk, INPUTS.Rk, Sigma{i}, Lambda{i}, Tkf, Obs);
%     [Sigma{i+1},Lambda{i+1}] = SMFfilter(A, BK, INPUTS.Qk, Hk, INPUTS.Rk, Sigma{i}, Lambda{i}, Tkf, Obs);

end
SigmaF = Sigma{i+1};
LambdaF = Lambda{i+1};
end






