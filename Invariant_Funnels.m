% /*! @Invariant_Funnels.m
% *************************************************************************
% <PRE>
% file.name       : Invariant_Funnels.m
% related files   :
% function&ablity : Calculate Invariant_Funnels(reference:Invariant Funnels around Trajectories using Sum-of-Squares Programming)
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/4/26     1.00        dengwei      build
% 2016/11/11    2.00        dengwei      add sample-based ROA
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:xx:nominal trajectory state,u0:nominal open loop
% control;Maxinterval:max time
%
% output:
% *************************************************************************
function [Vtrim, sys] = Invariant_Funnels(sys)

INPUTS = sys.INPUTS;

if 0
    load('my_testdata.mat')
else
    % use spot package calculate the funnel in msspoly variable
    t = msspoly('t');
    u = msspoly('u',sys.getNumInput);
    x = msspoly('x',sys.getNumStates);
    sys.p_t = t;
    sys.p_u = u;
    sys.p_x = x;
    %% ============linearlize the system near every points ===========
    sys = sys.timecalculation(50);% calcluate the time points
    taus = sys.breaks;
    sys = sys.tv_poly_linearize; % linearize the sys
    A = sys.A;
    B = sys.B;
    Qf = diag([1 1 10 10]); % this is must smaller than ROA
    % compute the region of attraction for target aero
    % V = regionOfAttraction_levelSetMethodYalmip(Qf)
    
    %% ===================time-varying Riccati equation=================
    [tv,Ss] = tv_lqr_riccati(sys.Maxinterval,A,B,INPUTS.Q,INPUTS.R,Qf);
    Spp = spline(tv,Ss);
    S = @(t) ppval(Spp,t);
    K = @(t) -inv(INPUTS.R(t))*B(t)'*S(t);
    sys.K = K;          % time varying feedback gain for local equilibrium
    %         Ac = @(t) A(t) + B(t)*K(t);
    %         Q0  = @(t) (INPUTS.Q(t) + S(t)*B(t)*inv(INPUTS.R(t))*B(t)'*S(t));
    %% ===============lypunove function===================
    % xT is the center of the goal region.
    % S0 defines the goal region (x-xT)'S0(x-xT) < = 1.
    % We use an estimate of the  basin of attraction for infinite horizon LQR, but choose your own.
    %         [taus1,Pstv] = tv_lyapunov(sys.Maxinterval,@(t)Ac(t),Q0,Qf);
    %         taus2 = flipud(taus1);
    %         Pstv0 = flip(Pstv,3);
    %         Pmes = spline(taus2,Pstv0); % x-x0 is expand to feedback
    
end
%% =======create funnel Lyapunov function ============
V = V_function(sys.p_t,sys.p_x,Spp,taus);

%% Options for funnel computation
options = struct();
options.n = sys.getNumStates;               % demination of state
options.m = sys.getNumInput;                % demination of control
options.rho0_tau = 10;                      % Determine initial guess for rho
options.max_iterations = 6;                 % Maximum number of iterations to run for
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
% regionOfAttraction
% V = regionOfAttraction_levelSetMethodYalmip(V0,f,A,ts,options);
if sys.SetVariable.ROAflag == 1
%     sampling-based method to solve Funnel%
    [V,rho] = SampleReachabilityFunnel(sys,V,taus,options);
else
%     sampling-based method to solve Funnel
    [V,rho] = myReachabilityFunnel(sys,V,taus,options);
end
disp('The funnel have done');
%%
% load('datafunnel0122.mat','rho')
save datafunnel20170210
V = V.updateV(foh(taus,1./rho'));
%%
figure
options.plotdims=[1 2];
options.inclusion='slice';
options.x0 = sys.FunTraj;
VFrame = V.inFrame(sys.FunTraj);
plot_myFunnel(sys,VFrame,options);
Vtrim = ShrinkFunnelwithObstacles(sys,VFrame);

%% Gaussian programming
% Vtrim = GMPLdatabase(sys.FunTraj,Vtrim);
%% =============end==============





