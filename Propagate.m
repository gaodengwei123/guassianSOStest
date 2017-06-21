% /*! @Propagate.m
% *************************************************************************
% <PRE>
% file.name       : Propagate.m
% related files   :
% function&ablity :
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/02/28    1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output:
% *************************************************************************
function New = Propagate(POP, node, GaoNodes, sys)
%% choose the traj=================================================

n = sys.getNumStates;
Vec = vec2mat([GaoNodes.trajLibrary(:).mark],3);
path = POP.path;
if length(path)==1
    Trajvec = [0,path(end),node];
else
    Trajvec = [path(end-1:end),node];
end
index = find(ismember(Vec,Trajvec,'rows'),1);
if isempty(index)  % can not be linked
    New = [];
    return
end
%% ================================================================

x0 = GaoNodes.trajLibrary(index).state;
u0 = GaoNodes.trajLibrary(index).control;
ts = GaoNodes.trajLibrary(index).time;
sys.PlotObj.solutionPlot.state = x0;

sys = sys.updateNominal(ts,x0,u0);
sys = sys.timecalculation(20);
sys = sys.tv_poly_linearize;
%% ===================time-varying Riccati equation=================
Qf = diag([1 1 10 10]); % this is must smaller than ROA
[tv,Ss] = tv_lqr_riccati(sys.Maxinterval,sys.A,sys.B,sys.INPUTS.Q,sys.INPUTS.R,Qf);
Spp = spline(tv,Ss);
S = @(t) ppval(Spp,t);
K = @(t) -inv(sys.INPUTS.R(t))*sys.B(t)'*S(t);
sys.K = K;          % time varying feedback gain for local equilibrium

Sigma0 = POP.Sigma;
Lambda0 = POP.Lambda;
% expand the EKF
[Sigma,Lambda] = Funnel_filter(sys,sys.breaks,Sigma0,Lambda0);
New.Sigma = Sigma;
New.Lambda = Lambda;
New.node = node;
New.funnel = GaoNodes.FunnleLibrary(index);
New.traj = GaoNodes.trajLibrary(index);
New.cost = GaoNodes.trajLibrary(index).cost;
%% ================================================================

Lambdachol = chol(Lambda+eps*eye(n));    % passive
Lambda = Lambdachol'*Lambdachol;
Sigmachol = chol(Sigma+eps*eye(n));
Sigma = Sigmachol'*Sigmachol;
pdim =  [1 0 0 0;0 1 0 0]';
ELambda = projection(ellipsoid(x0(end,:)',Lambda),pdim);
ESigma = projection(ellipsoid(x0(end,:)',Sigma),pdim);
% Covariance = Lambda;
% ELambda = error_ellipse(x0(end,1:2)',Covariance(1:2,1:2));
% Covariance = Sigma;
% ESigma = error_ellipse(x0(end,1:2)',Covariance(1:2,1:2));

drawnow
% plot(ELambda,'r');
% plot(ESigma,'y');

% plot SOS funnel
V = GaoNodes.FunnleLibrary(index).V;
V.x0 = sys.FunTraj;
% options.plotdims=[1 2];
% options.inclusion='slice';
% options.x0 = sys.FunTraj;
% VFrame = V.inFrame(sys.FunTraj);
% plot_myFunnel(sys,VFrame,options);

boundPointMat = PdiffVE(V,ESigma,sys.breaks(end),[1;2]);
% plot3(boundPointMat(1,:),boundPointMat(2,:),repmat(.2,1,size(boundPointMat,2)),'g');
if isempty(boundPointMat)
    New = [];
end
%%
rho = GaoNodes.FunnleLibrary(index).rho;
fact = linspace(2,1,20); % factor to enlage funnel for compution error
rho = (fact.*rho')';
V = V.updateV(foh(V.getbreak,1./rho'));
% ESigma = projection(ellipsoid(x0(end,:)',diag([10 10 10 10])*Sigma),pdim);% situation 2
ESigma = projection(ellipsoid(x0(end,:)',diag([7 7 7 7])*Sigma),pdim);% situation 3
% plot(ELambda,'y')
boundPointMat2 = PdiffVE(V,ESigma,sys.breaks(end),[1;2]);
New.Pdiff.P1 = boundPointMat;
New.Pdiff.P2 = boundPointMat2;




%%
a = isinternal(ELambda,boundPointMat,'i'); % check bound is inside ellipsoid or not

if ~all(a==0)
    New = [];
end

end

