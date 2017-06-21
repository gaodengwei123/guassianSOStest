function [Pini,rhoini] = plotFunnelLibary(sys,funnel,OL_time,OL_state,ShrinkRho,j)
V = funnel.V;
rho = funnel.rho;
rho = rho*ShrinkRho;
fact = linspace(2,1,20); % factor to enlage funnel for compution error
rho = (fact.*rho')';
taus = linspace(0,OL_time(end),20);
options.plotdims=[1 2];
options.inclusion='slice';
sys.breaks = taus;
V.getbreak = taus;
statefun = spline(OL_time,OL_state');
sys.FunTraj = polyniminalTrajectory(statefun);
V.x0 = sys.FunTraj;
Pini = ppval(V.Spp,0); % the first P
rhoini = rho(1);     % the first rho
V = V.updateV(foh(taus,1./rho'));

options.x0 = sys.FunTraj;
VFrame = V.inFrame(sys.FunTraj);
sys.PlotObj.solutionPlot.state = OL_state;
plot_myFunnel(sys,VFrame,options);



%%
fact = linspace(0.1,1,20); % factor to enlage funnel for compution error
rho = (fact.*rho')';
V = V.updateV(foh(taus,1./rho'));
options.x0 = V.x0.eval(0);
xfun0 = getLevelSet(V.x,V.getPoly(0),options);
h0 = plot(xfun0(1,:),xfun0(2,:),'y','linewidth',2);
% 
% fact = linspace(2,1,20); % factor to enlage funnel for compution error
% rho = (fact.*rho')';
% V = V.updateV(foh(taus,1./rho'));
% options.x0 = V.x0.eval(0);
% xfun0 = getLevelSet(V.x,V.getPoly(0),options);
% h1 = plot(xfun0(1,:),xfun0(2,:),'g','linewidth',2);




%%
if j==1
    
    %     Sigma = diag([0.3 0.3]);% situation 1
    %     Sigma = diag([0.6 0.6]);% situation 2
    Sigma = diag([0.9 0.9]);% situation 3
    ESigma = ellipsoid([8.5;8.5],Sigma);
    boundPointMat = PdiffVE(VFrame,ESigma,sys.breaks(1),[1;2]);
    plot(boundPointMat(1,:),boundPointMat(2,:),'g');
end
end






