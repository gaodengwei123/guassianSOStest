clear all
close all
clc
dbstop if error


t = msspoly('t');
u = msspoly('u',1);
x = msspoly('x',4);

sys = Acrobot();
Q = diag([10,10,1,1]); R = 1;
[A,B] = linearize(sys,0,sys.xG,sys.uG);
[K,S] = lqr(full(A),full(B),Q,R);


sys = feedbackSys(sys,K,B);

V0 = V_function(t,x,S,[]);
options.method={'levelset'};
V0.x0 = sys.sys.xG;
V = regionOfAttraction(sys,V0,options);
V = extractV_Function(V,x,V0.x0);

options = struct();
options.controller_deg = 1; % Degree of polynomial controller to search for
options.max_iterations = 10; % Maximum number of iterations (3 steps per iteration)
options.converged_tol = 1e-3; % Tolerance for checking convergence
options.rho0 = 0.001; % Initial guess for rho
options.clean_tol = 1e-6; % tolerance for cleaning small terms
options.backoff_percent = 5; % 1 percent backing off
options.degL1 = options.controller_deg + 1; % Chosen to do degree matching
options.degLu = options.controller_deg - 1; 
options.degLu1 = 2; 
options.degLu2 = 2;
options.degLup = 2;
options.degLum = 2; 

options.controller_deg = 1;
% sys = setInputLimits(sys,-20,20);
options.max_iterations = 10;
sys = Acrobot();
[c,Vmax] = maxROAFeedbackcontrol(sys,K,V,options);

%% plot
Vmax = inFrame(Vmax,sys.xG);
V = inFrame(V,sys.xG);
% x1 x2
figure
hold on
axHandle = gca;
xlabel(axHandle,'x1');
ylabel(axHandle,'x2');
options.color = 'b';
options.x0 = sys.xG;
plotROA(Vmax,options)
options.color = 'r';
plotROA(V,options)

% x1 x3
figure
hold on
axHandle = gca;
xlabel(axHandle,'x1');
ylabel(axHandle,'x3');
options.color = 'b';
options.x0 = sys.xG;
options.plotdims = [1 3];
options.inclusion = 'projection';
plotROA(Vmax,options)

options.color = 'r';
plotROA(V,options)

