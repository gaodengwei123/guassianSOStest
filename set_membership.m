% /*! @set_membership.m
% *************************************************************************
% <PRE>
% file.name       : set_membership.m
% related files   :
% function&ablity : set_membership filter
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/5/22     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input : 
%
% output:
% *************************************************************************
function error_boundary = set_membership(f,x0,Q,R,Z,P0)
x10=f(x0);

% linearlized noise Qbar
Xk = [x0-sqrt(diag(P0)),x0+sqrt(diag(P0))];
syms v1 v2 v3 v4 v5 v6
S = f(0,[v1, v2, v3, v4, v5, v6],zeros(3,1));
% hessian matrix
Hf1 = hessian(S(1),[v1, v2, v3, v4, v5, v6]);
Hf2 = hessian(S(2),[v1, v2, v3, v4, v5, v6]);
Hf3 = hessian(S(3),[v1, v2, v3, v4, v5, v6]);
Hf4 = hessian(S(4),[v1, v2, v3, v4, v5, v6]);
Hf5 = hessian(S(5),[v1, v2, v3, v4, v5, v6]);
Hf6 = hessian(S(6),[v1, v2, v3, v4, v5, v6]);
% jacobi matrix
Phi = jacobi(S,[v1, v2, v3, v4, v5, v6]);
%% ¸³Öµ

%%
XR2 = 0.5*diag(Xk'*[Hf1,Hf2,Hf3,Hf4,Hf5,Hf6]*Xk);
Qbar = 2*XR2^2;
betaQ = sqrt(trace(Q))/(sqrt(trace(Qbar))+sqrt(trace(Q)));
Qhat = Qbar/(1-betaQ)+Q/betaQ;
beta = sqrt(trace(Qhat))/(sqrt(trace(Phi*P0*Phi'))+sqrt(trace(Qhat)));
P10 = Phi*P0/(1-beta)*Phi'+Qbar/beta;

%% measurement update
h10 = f(x10);%Ô¤²â
H = 0;%Ô¤²â

pm = max(H*P10*H');
rm = max(R);
rou = sqrt(rm)/(sqrt(pm)+sqrt(rm));
invW = inv(H*P10/(1-rou)*H+R/rou);
K = P10/(1-rou)*H'*invW;
x1 = x10+K*(Z-h10);
P1bar = P10/(1-rou)-P10/(1-rou)*H'*invW*H*P10/(1-rou);
dert = 1-(Z-h10)'*invW*(Z-h10);
P1 = dert*P1bar;
miu = 0;% MIT




