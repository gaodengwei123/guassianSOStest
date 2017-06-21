function FunnleLibrary = FunnleLibrary_primitives(sys,trajLibrary)
% use spot package calculate the funnel in msspoly variable
t = msspoly('t');
u = msspoly('u',sys.getNumInput);
x = msspoly('x',sys.getNumStates);
sys.p_t = t;
sys.p_u = u;
sys.p_x = x;
sys = sys.updateNominal(trajLibrary.time,trajLibrary.state,trajLibrary.control);
sys = sys.timecalculation(20);
sys = sys.tv_poly_linearize; % linearize the sys
%% ===================time-varying Riccati equation=================
Qf = diag([1 1 10 10]);
[tv,Ss] = tv_lqr_riccati(sys.Maxinterval,sys.A,sys.B,sys.INPUTS.Q,sys.INPUTS.R,Qf);
Spp = spline(tv,Ss);
S = @(t) ppval(Spp,t);
K = @(t) -inv(sys.INPUTS.R(t))*sys.B(t)'*S(t);
sys.K = K;          % time varying feedback gain for local equilibrium
% build V function
V = V_function(sys.p_t,sys.p_x,Spp,sys.breaks);

%% Options for funnel computation
options = struct();
options.n = sys.getNumStates;               % demination of state
options.m = sys.getNumInput;                % demination of control
options.rho0_tau = 10;                      % Determine initial guess for rho
options.max_iterations = 3;                 % Maximum number of iterations to run for
options.stability = false;                  % true implies that we want exponential stability
options.converged_tol = 1e-3;               % Tolerance for checking convergence
options.lyap_parameterization = 'rho';      % Lyapunov upper symbol
options.rho0 = 1;                           % Initial "guessed" rho
options.controller_deg = 1;                 % Degree of polynomial controller to search for
options.clean_tol = 1e-3;                   % tolerance for cleaning small terms
options.backoff_percent = 5;                % 1 percent backing off
options.degL1 = 4;                          % Chosen to do degree matching deg(Vdot{1},x);
options.degLu = options.controller_deg - 1;
options.degLu1 = 2;
options.degLu2 = 2;
options.degLup = 2;
options.degLum = 2;
options.saturations = 0;
options.plot_rho = 1;
%% slove the funnel of system
% sampling-based method to solve Funnel
tic
[V,rho] = myReachabilityFunnel(sys,V,sys.breaks,options);
toc
FunnleLibrary.V = V;
FunnleLibrary.rho = rho;
V = V.updateV(foh(sys.breaks,1./rho'));
% plot funnel
options.plotdims=[1 2];
options.inclusion='slice';
options.x0 = sys.FunTraj;
VFrame = V.inFrame(sys.FunTraj);
sys.PlotObj.solutionPlot.state = trajLibrary.state;
plot_myFunnel(sys,VFrame,options);
% shrink the funnel according to Obstacles
% Vtrim = ShrinkFunnelwithObstacles(sys,VFrame);


end









