%% ==============ground Vehicle model=============================
clear all
clc
global INPUTS

eps = 1e-10;
timenow = 0;
% simulantion total step
end_step = 300;
% simulantion step
temp = 1;
n = 4;   % Number of States
m = 1; % Number of Inputs
%% ======================dynamic model============================
% pvar x1 x2 x3 x4 u t;
% x = [x1,x2,x3,x4]';
% Xu = 20.55;
% Yv = 5.0326;
% Nr = 0.05;
% Izz = 0.02532;
% B = @(t)[zeros(3);eye(3)];
%% canstant set
INPUTS.speed = 10;
INPUTS.obstacle = [10 30
                   30 40
                   25 20];
INPUTS.obstacleRadius = 5;

% show the equation in local_controlDae
f = @(t,x,u)[INPUTS.speed*pcos(x(3))
             INPUTS.speed*psin(x(3))
             x(4)
             u];
% inital state node
U_int(:,1) = 0;   % control law

% start point
OL_state(:,1) = [0;0;0;0];
% terminal point
INPUTS.state_goal = [50;40;1;0];
INPUTS.limit_x1_min = 0;
INPUTS.limit_x2_min = 0;
INPUTS.limit_x3_min = -pi;
INPUTS.limit_x4_min = -pi;
INPUTS.umin = -pi;
INPUTS.limit_x1_max = 50;
INPUTS.limit_x2_max = 50;
INPUTS.limit_x3_max = pi;
INPUTS.limit_x4_max = pi;
INPUTS.umax = pi;
% cost function Q and R also be used in local_controlCost
INPUTS.Q = @(t)eye(n);
INPUTS.R = @(t)10*eye(m);

%% =================open loop nominal trajectory===================
[OL_state,OL_control,OL_time] = optimal_local_programming(OL_state(:,1),INPUTS.state_goal);
statefun1 = spline(OL_time,OL_state');
statefun = @(t) ppval(statefun1,t);
controlfun1 = spline(OL_time,OL_control');
controlfun = @(t) ppval(controlfun1,t);
Maxinterval = [0 max(OL_time)]; % cost time
[rho0, Ppp, upp, K0] = Invariant_Funnels(f,statefun,controlfun,Maxinterval);


%% ===================end =======================================

%% ==========================guss learning processing=======================
GMPLprogama
%% ===================assignment claculate==================================
for i=1:1000
    % extend kalman filter
    %     [dert_x,x1,P1] = kalman_filter(x10,C,D,inv_MS,J_n,ROV,real_Z,Q,R,P0,temp,i);
    %     set_membership
    
    
end


% plot figure
figure
th = linspace(-pi,pi,100);
ell = S0(1:2,1:2)^(-1/2)*[cos(th) ; sin(th)];
fill(ell(1,:),ell(2,:),0.8*ones(1,3))

figure
plot3(statedint(1,:),statedint(2,:),statedint(3,:),'r')
[C,h] = pcontour(V,0.093);  % plot ROA
hold on
plot3(x0(1,:),x0(2,:),x0(3,:))
grid on
figure
plot(U_body(1,:))
hold on
plot(U_body(2,:),'y')
plot(U_body(3,:),'r')
plot(U_body(4,:),'k')
save g2
save par_state2
save control_gain2
