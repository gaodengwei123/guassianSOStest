% /*! @RRTPropagate.m
% *************************************************************************
% <PRE>
% file.name       : RRTPropagate.m
% related files   :
% function&ablity :
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/03/07    1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output:
% *************************************************************************
function x_son = RRTPropagate(x_parent, x_son, sys)
% add parents of the new state
x_son.path = [x_parent.path; x_son.Tag];
ts = x_son.Ftuple.breaks;

sys.K = x_son.Ftuple.K;              % time varying feedback gain for local equilibrium
sys.FunTraj = x_son.Ftuple.traj;    
sys.A =  x_son.Ftuple.A;
sys.B =  x_son.Ftuple.B;

Sigma0 = x_parent.Sigma;
Lambda0 = x_parent.Lambda;
% expand the EKF
[Sigma,Lambda] = Funnel_filter(sys,ts,Sigma0,Lambda0);
% update the convariance of new node
x_son.Sigma = Sigma;
x_son.Lambda = Lambda;
x0 = x_son.Ftuple.traj.eval(ts);
x_son.t0 = x_parent.t0+ts(end);
x_son.cost = x_parent.cost+x_son.Ftuple.cost;
%% ================================================================
pdim =  [1 0 0 0;0 1 0 0]';
ELambda = projection(ellipsoid(x0(:,end),Lambda),pdim);
ESigma = projection(ellipsoid(x0(:,end),Sigma),pdim);
% Covariance = Lambda;
% ELambda = error_ellipse(x0(end,1:2)',Covariance(1:2,1:2));
% Covariance = Sigma;
% ESigma = error_ellipse(x0(end,1:2)',Covariance(1:2,1:2));
% drawnow
% plot(ELambda,'r');
% plot(ESigma,'y');

% plot SOS funnel
% x_son.FunnelHandle = plot(x_son.Ftuple);
x_son.FunnelHandle =[];
V = x_son.Ftuple.V;
boundPointMat = PdiffVE(V,ESigma,ts(end),[1;2]);
% plot3(boundPointMat(1,:),boundPointMat(2,:),repmat(.2,1,size(boundPointMat,2)),'g');
if isempty(boundPointMat)
    x_son = [];
    return
end
x_son.Ftuple.Pdiff = boundPointMat;
% check bound is inside ellipsoid or not    
a = isinternal(ELambda,boundPointMat,'i'); 

if ~all(a==0)
    x_son = [];
    return
end

end

