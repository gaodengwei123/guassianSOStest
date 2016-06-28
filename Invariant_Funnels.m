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
% 2016/4/26     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:
%
% output:
% *************************************************************************
function [rho0, Ppp, upp, K0] = Invariant_Funnels(f0,xx,u0,Maxinterval)
global INPUTS
Q = INPUTS.Q(1);
R = INPUTS.R(1);
n = length(xx(0)); m = length(u0(0));

%% Guessian guess programming
% xT = [INPUTS.limit_x1_min+(INPUTS.limit_x1_max-INPUTS.limit_x1_min)*rand(1)
%       INPUTS.limit_x2_min+(INPUTS.limit_x2_max-INPUTS.limit_x2_min)*rand(1)
%       INPUTS.limit_x3_min+(INPUTS.limit_x3_max-INPUTS.limit_x3_min)*rand(1)
%       INPUTS.limit_x4_min+(INPUTS.limit_x4_max-INPUTS.limit_x4_min)*rand(1)];
% uT = xT(4);

% Maxinterval = [0 max(taus)]; % cost time
[A,B] = tv_poly_linearize(f0,xx,u0);
Qf = 10*Q;
[ts,Ss] = tv_lqr_riccati(Maxinterval,A,B,INPUTS.Q,INPUTS.R,Qf);
Spp = spline(ts,Ss);
S = @(t) ppval(Spp,t);
K = @(t) inv(INPUTS.R(t))*B(t)'*S(t);
Ac = @(t) A(t) - B(t)*K(t);
Q0  = @(t) (INPUTS.Q(t) + S(t)*B(t)*inv(INPUTS.R(t))*B(t)'*S(t));
% xT is the center of the goal region.
% S0 defines the goal region (x-xT)'S0(x-xT) < = 1.
% We use an estimate of the  basin of attraction for infinite horizon LQR, but choose your own.
[taus1,Ps] = tv_lyapunov(Maxinterval,@(t) Ac(t),Q0,Qf);
taus = flipud(taus1);
N = length(taus);
Ppp = interp1(taus,reshape(permute(Ps,[3 1 2]),N,n*n),'linear', 'pp');
upp = spline(taus,u0(taus'));
% 
% integral the state with feedback
state(:,1) = xx(0);
stateOL(:,1) = xx(0);
for i = 1:length(taus)-1
    [~,state1] = ode45(@(t,x)f0(t,x,u0(t)-K(t)*(state(:,i)-xx(taus(i)))),taus(i:i+1),state(:,i));
    [~,stateOL1] = ode45(@(t,x)f0(t,x,u0(t)),taus(i:i+1),stateOL(:,i));
    state(:,i+1) = state1(end,:)';
    stateOL(:,i+1) = stateOL1(end,:)';
    Ufeedback1(i,:) = u0(taus(i))-K(taus(i))*(state(:,i)-xx(taus(i)));
end
Ufeedback1(i+1,:) = u0(taus(i+1))-K(taus(i+1))*(state(:,i+1)-xx(taus(i+1)));
Ufeedback0 = spline(taus,Ufeedback1');
Ufeedback = @(t) ppval(Ufeedback0,t);
figure(1)
plot(state(1,:),state(2,:),'*b')
hold on 
plot(stateOL(1,:),stateOL(2,:),'k')
% feedback control and real state
xpp = spline(taus,state);
upp = spline(taus,Ufeedback(taus));

% rou = GMPLdatabase();
% Pp = @(t) ppval(Ppp,t);
% Up = @(t) ppval(upp,t);

% GMPLprogama
% S0 = 1.01*S0/rho0;

%% ================ellipsoid calculate============================
Max_t = length(taus); % time step
ts = flipud(taus);
Ps = roundn(Ps,-6); % precision adjust
figure(1)
hold on
pB1 = [1 0 0 0
    0 1 0 0
    0 0 1 0]';
plot_elliposoid = 1;
if plot_elliposoid ==1
    for i=1:Max_t
        EE(i) = ellipsoid(xx(taus(i)),Ps(:,:,Max_t-i+1));
        % projection the ellipsoid to different basis
        EE1(i) = projection(EE(i),pB1);
        plot(EE1(i))
        drawnow;
        frame(i)=getframe(gcf);
    end
    writegif('test1.gif',frame,0.1);
end

c = 3;
rhot = exp(c*(taus-max(taus))/(max(taus)-min(taus)));
rhopp = interp1(taus,rhot,'linear','pp');
% lagrange_multipliers



% % not bad direction(dimensions are less than 3)
% RA1 = minkdiff(EE2,EE1);
% RA2 = minkdiff(FE2,FE1);
% RA = [RA1;RA2];
% 
% for i=1:size(RA,2)
%     EEE = minkdiff_ia(E2,E1,RA(:,i));
%     EEF = minkdiff_ea(E2,E1,RA(:,i));
% end


