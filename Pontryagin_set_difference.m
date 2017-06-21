% /*! @Pontryagin_set_difference.m
% *************************************************************************
% <PRE>
% file.name       : Pontryagin_set_difference.m
% related files   :
% function&ablity : Calculate the Pontryagin set difference from sensor noise and Funnel arround nominal trajectory
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/7/27     1.00                     build
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
function current_ell = Pontryagin_set_difference(sys,V)
% eps = 1e-2;
% just for videoRay model
xtraj = sys.FunTraj;
utraj = sys.FunContr;
taus = sys.breaks;

INPUTS = sys.INPUTS;
% integral the state with feedback
n = sys.getNumStates;
m = sys.getNumInput;
N = sys.timeNum;

% dynamic integral
figure
hold on
plotDynNum = 1;    % Monte Carlo simulation
SampleNum = 100;
ts = linspace(taus(1),taus(end),SampleNum);  % set sample time
for i = 1:plotDynNum
    InitialState = xtraj.eval(0)+INPUTS.noise_w.*(2*rand(n,1)-1); % real state: add Gaussian noise at initial state
    [~, ~, ~, ~, feedbacksys] = TaylorExpansion(sys,3);           % feedback sys
    Stochsys = StochasticSystem(feedbacksys,INPUTS.Qk,INPUTS.Rk);
    odefun = @(t,x)Stochsys.dynamics(t,x,zeros(sys.getNumInput,1), INPUTS.noise_w.*(2*rand(n,1)-1));%, INPUTS.noise_v.*(2*randn(Numoutput,1)-1)
    sol = ode45(odefun,[taus(1) taus(end)],InitialState);
    % draw real trajectory
    state = deval(sol,ts);
    plot(state(1,:),state(2,:),'k')
    drawnow
end
% draw nominal trajectory
nominal_traj = xtraj.eval(ts);
plot(nominal_traj(1,:),nominal_traj(2,:))

%% kalman filter
kalmanflag = 1;
if kalmanflag == 1
  [Sigma,Lambda,Obs] = Gao_filter(sys,ts);
end

%% plot elliposid
for j = 1:size(INPUTS.obstacle,1)
    OBSEEE(j) = ellipsoid(INPUTS.obstacle(j,:)',diag(15^2*ones(2,1)));
    plot(OBSEEE(j),'k')
end
for i = 1:SampleNum
    pdim =  [1 0 0 0;0 1 0 0]';
    EEEE(i) = projection(ellipsoid(nominal_traj(:,i),Lambda{i}),pdim);
    EEE(i) = projection(ellipsoid(nominal_traj(:,i),Sigma{i}),pdim);
%     Covariance = Lambda{i};
%     EEE(i) = error_ellipse(nominal_traj(1:2,i),Covariance(1:2,1:2));
%     Covariance = Sigma{i};
%     EEEE(i) = error_ellipse(nominal_traj(1:2,i),Covariance(1:2,1:2));

    drawnow
    plot(EEE(i),'r');
    plot(EEEE(i),'b');

end

% plot SOS funnel
options.plotdims=[1 2];
options.inclusion='slice';
options.x0 = sys.FunTraj;
VFrame = V.inFrame(sys.FunTraj);
plot_myFunnel(sys,VFrame,options);

for i=1:N  
    boundPointMat = PdiffVE(V,EEE(i),sys.breaks(i),[1;2]);
    Savebound{i} = boundPointMat;
    plot3(boundPointMat(1,:),boundPointMat(2,:),repmat(.2,1,size(boundPointMat,2)),'g');
end

%% ===============shrink the funnel based on noise=================
Vtrim =  ShrinkFunnelwithnoise(V,Savebound);
options.color = [0 1 0];
plot_myFunnel(sys,Vtrim,options);

% rho = ones(length(taus),1);
% for i = 1:length(taus)
%     if i<10
%         rho(i) = 0.8*rho(i);
%     end
% end
%
% V = V.updateV(foh(taus,1./rho'));
% figure
% options.plotdims=[1 2];
% options.inclusion='slice';
% options.x0 = sys.FunTraj;
% VFrame = V.inFrame(sys.FunTraj);
% plot_myFunnel(sys,VFrame,options);
%% ===============


%% error ellipse computation
for i = 1:length(ts)
    RealxEll{i} = ellipsoid(state,Pk);
    NomxEll{i} = ellipsoid(nominal_traj,Perror);
end

% Bk = B*temp;
% Ak = eye(n)+Ac*temp+Bk*K;
% beta1 = sqrt(trace(Q))/(sqrt(trace(-Bk*K*P0))+sqrt(trace(Q)));
% Qmao = -Bk*K*P0/(1-beta1)+Q/beta1;
% beta = sqrt(trace(Qmao))/(sqrt(trace(Ak*P0*Ak'))+sqrt(trace(Qmao)));
% PP = double(Ak*P0/(1-beta)*Ak'+Qmao/beta);

%% ================ellipsoid calculate============================
% projection the ellipsoid to the first two dimension
if sys.mark ==1
    pB1 = [1 0 0 0
        0 1 0 0]';
    Ps = V.S;
    figure(10)
    show(sys.PlotObj.map)
    xlabel('x(m)');
    ylabel('y(m)');
    title('')
    %     ell_measure1 = ellipsoid(INPUTS.obstacle(1,1:2)',1.5*diag([9*INPUTS.obstacleRadius^2,9*INPUTS.obstacleRadius^2]));
    %     ell_measure2 = ellipsoid(INPUTS.obstacle(2,1:2)',1.5*diag([9*INPUTS.obstacleRadius^2,9*INPUTS.obstacleRadius^2]));
    %     ell_measure3 = ellipsoid(INPUTS.obstacle(3,1:2)',1.5*diag([9*INPUTS.obstacleRadius^2,9*INPUTS.obstacleRadius^2]));
    %     ell_measure4 = ellipsoid(INPUTS.state_goal(1:2),1.5*diag([9*INPUTS.obstacleRadius^2,9*INPUTS.obstacleRadius^2]));
    hold on
    %     plot(ell_measure1,'b.')
    %     plot(ell_measure2,'b.')
    %     plot(ell_measure3,'b.')
    %     plot(ell_measure4,'b.')
    
elseif sys.mark ==2
    pB1 = [1 0 0 0 0 0
        0 1 0 0 0 0]';
    Ps = V.S;
    figure(10)
    show(sys.PlotObj.prmSimple.Map)
    xlabel('x(m)');
    ylabel('y(m)');
    title('')
    hold on
    VideoRay_measure = ellipsoid([5;0],diag([25,25]));
    plot(VideoRay_measure,'b.')
end

% V = V.shrinkV(repmat(1.025,V.timeNum,1)); % test the car funnels

% plot funnel
options.plotdims=[1 2];
options.inclusion='slice';
plot_myFunnel(V,options);

if INPUTS.plot_elliposoid ==2  % use elliposoid toolbox to calculate
    muli = 1;
    for i=1:N
        i
        % ellipsoid around nominal trajectory(ROA)
        EEE(i) = ellipsoid(xtraj(taus(i)),diag(diag(Ps{i})));
        %         collision_found = collision_check_point(sys,EEE(i));   % check the safty of path planning
        %         if collision_found
        %             current_ell = EEE(i);
        %             return
        %         end
        EEE1(i) = projection(EEE(i),pB1);
        EE(i) = ellipsoid(state(:,i),diag(diag(sqrt(Perror{i}))));
        % projection the ellipsoid to different basis
        EE1(i) = move2origin(projection(EE(i),pB1));
        % Pontryagin set difference
        if EEE1(i)>EE1(i)
            [~, boundPointMat] = minkdiff(EEE1(i),EE1(i));
        else
            warning('out of ROA')
            current_ell = EEE(i);
            return
        end
        %         plot(boundPointMat(1,:),boundPointMat(2,:))
        patch([boundPointMat(1,:)],[boundPointMat(2,:)], [1 0 0])
        
        drawnow;
        
        % [Poly_ROA, Rho_ROA] = ellipsoid2Polyhedron(nom_x(taus(i)),V.x,V.V{i});
        % %% this is calculate by MPT Polyhedron
        % options.x0 = nom_x(taus(i));
        % Poly_ROA = getLevelSet(V.x,V.V{i},options);
        % P = Polyhedron('V',Poly_ROA');
        % % P.plot('color','b')
        % Per = Perror{i};
        % Per(Per<0)=0;
        % Verror = V.x'*(sqrt(Per)^-1)*V.x;
        % options.x0 = zeros(n,1);
        % addorigin = nom_x(taus(i));
        % Poly_error1 = getLevelSet(V.x,Verror,options);
        % Poly_error = Poly_error1+repmat(addorigin(1:2),1,size(Poly_error1,2));
        % S = Polyhedron('V',Poly_error');
        % % S.plot('color','k')
        % Q = P-S;
        % Q.plot('color','g');
    end
end
hold off
current_ell = ellipsoid(xtraj(taus(end)),diag(diag(Ps{end})));
% plot figure
sys.PlotObj.state = state;
sys.PlotObj.stateOL = nominal_traj;
sys.PlotObj.Max_t = Max_t;
sys.PlotObj.Esetdiff = Eea;
sys.PlotObj.Eerror = EE1;
