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
% input : state10: feedbackcontrol next step;state0:optimal
% estimate;nominal_state: open loop nonminal state;P0:state covariance
% U: feedback control
% output:
% *************************************************************************
function [dert_X,P1,PP] = set_membership(sys,A,B,K,state10,state0,P0,ts,t_step,n,U)
INPUTS = sys.INPUTS;
Q = INPUTS.Qk;
R = INPUTS.Rk;
%% ESMF
% interval of time step
dert = (ts(t_step+1)-ts(t_step));
% linearlized noise Qbar
Xk = [state0-sqrt(diag(P0)),state0+sqrt(diag(P0))];
syms v1 v2 v3 v4
% discretization system
S = [v1 v2 v3 v4]'+sys.f(0,[v1, v2, v3, v4],U)*dert;
% Hessian matrix for each state
Hf1 = subs(hessian(S(1),[v1, v2, v3, v4]),[v1, v2, v3, v4],state0');
Hf2 = subs(hessian(S(2),[v1, v2, v3, v4]),[v1, v2, v3, v4],state0');
Hf3 = subs(hessian(S(3),[v1, v2, v3, v4]),[v1, v2, v3, v4],state0');
Hf4 = subs(hessian(S(4),[v1, v2, v3, v4]),[v1, v2, v3, v4],state0');
% state transform matrix
x = msspoly('x',n);
Phi = eye(n)+double(subs(diff(sys.f(0,x,U),x),x,state0))*dert;

% Linearization & recursion Linearization & recursion
XR2 = 0.5*repmat(diag((Xk(:,1)+(Xk(:,2)-Xk(:,1)).*rand(n,1))'-state0'),1,n)*[Hf1;Hf2;Hf3;Hf4]*(Xk(:,1)+(Xk(:,2)-Xk(:,1)).*rand(n,1)-state0);% random 
Qbar = 2*diag(XR2).^2;
% set membership filter
betaQ = sqrt(trace(Q))/(sqrt(trace(Qbar))+sqrt(trace(Q)));
% Qhat = double(Qbar/(1-betaQ)+Q/betaQ)
Qhat = Q;
betak = sqrt(trace(Qhat))/(sqrt(trace(Phi*P0*Phi'))+sqrt(trace(Qhat)));
P10 = Phi*P0/(1-betak)*Phi'+Qhat/betak;

%error dynamic
Bk = B*dert;
Ak = eye(n)+A*dert+Bk*K;
%% measurement update
if rem(t_step,5)==0
    H = [eye(3),zeros(3,1)];% linearlize the measurement function
    Z = state10(1:3)+INPUTS.noise_v.*(2*rand(3,1)-1);
    
    [~,pmS,~] = svd(H*P10*H',0);
    [~,rmS,~] = svd(R,0);
    pm = pmS(1,1);
    rm = rmS(1,1);
    
    rou = sqrt(rm)/(sqrt(pm)+sqrt(rm));
    inv_W = inv(H*P10/(1-rou)*H'+R/rou);
    Lk1 = P10/(1-rou)*H'*inv_W;
    dert_X = Lk1*(Z-H*state10);
    P1bar = P10/(1-rou)-P10/(1-rou)*H'*inv_W*H*P10/(1-rou);
    dertk1 = 1-(Z-H*state10)'*inv_W*(Z-H*state10);
    P1 = double(dertk1*P1bar);
else
    P1=double(P10);
    dert_X = zeros(n,1);
end
    % error dyanmic
    beta1 = sqrt(trace(Q))/(sqrt(trace(-Bk*K*P0))+sqrt(trace(Q)));
    Qmao = -Bk*K*P0/(1-beta1)+Q/beta1;
%     dert_X = Lk1*(H*sqrt(diag(P0)).*(2*rand(3,1)-1)+INPUTS.noise_v.*(2*rand(3,1)-1));
    beta = sqrt(trace(Qmao))/(sqrt(trace(Ak*P0*Ak'))+sqrt(trace(Qmao)));
    PP = double(Ak*P0/(1-beta)*Ak'+Qmao/beta);

end





