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
% 2017/2/23     3.00                     kalman filter
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input :

% output:
% *************************************************************************
function [Sigma,Lambda,Obs] = Gao_filter(sys,Time)
xtraj = sys.FunTraj;
INPUTS = sys.INPUTS;
% INPUTS.Rk = [INPUTS.Rk zeros(3,1);zeros(1,3) 0.00001];
n = sys.getNumStates;
Sigma{1} = INPUTS.Qk;           % convariance of estimate state x_hat
Lambda{1} = INPUTS.Qk;          % convariance of x_tilde = x - x_bar
Hk = eye(size(INPUTS.Rk,1),n);  % measurement Jacobian matrix
checkstate = xtraj.eval(Time);
SampleNum = length(Time);
for i = 1:SampleNum-1
    A = sys.A(Time(i));
    BK = sys.B(Time(i))*sys.K(Time(i));
    % measurement condition: near obstacle <15
    for j=1:1%size(INPUTS.obstacle,1)
        if norm(checkstate(1:2,i)'-INPUTS.obstacle(j,:))<15
            Obs(i) = 1;
        else
            Obs(i) = 0;
        end
    end
    Tkf = Time(i+1)-Time(i);
    [Sigma{i+1},Lambda{i+1}] = kfilter2(A, BK, INPUTS.Qk, Hk, INPUTS.Rk, Sigma{i}, Lambda{i}, Tkf, Obs(i));
%     [Sigma{i+1},Lambda{i+1}] = SMFfilter(A, BK, INPUTS.Qk, Hk, INPUTS.Rk, Sigma{i}, Lambda{i}, Tkf, Obs(i));
    if any(eig(Lambda{i})<0)
        mm = 0;
    end
    if any(eig(Sigma{i})<0)
        nn = 0;
    end
end
end






