% /*! @optimal_local_programming.m
% *************************************************************************
% <PRE>
% file.name       : optimal_local_programming.m
% related files   :
% function&ablity :gpops-II to solve NLP problem
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       : 
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/05/7     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:f0:the state function;xT&uT:the state and control law,Q&R:paramaeters of cost

% output:

% [ts,Ss] = tv_lqr_riccati(tspan,A,B,Q,R,S0);
% *************************************************************************
function [openloopstate,control,time] = optimal_local_programming(x0,xT)
global INPUTS
%% =======================gpops-II============================begin
% range of time & states & control output (guess)
t0 = 0; % tf = norm(xT-x0)/10;

% -----------------------
% Optimization Problem
% -----------------------
% The problem solved here is given as follows :
% Minimize
% t_f , J
% subject to the dynamic constraints
% dx / dt = v * cos ( psi )
% dy / dt = v * sin ( psi )
% dpsi / dt = ddpsi
% ddpsi / dt = u
% with the boundary conditions
% --------------------------------------------------
%% =============== initial trajectory guess =======================
%build map and guess initial points
[guessTime,guessstate,guessControl] = Optguesspoint(x0(1:2),xT(1:2));
%%
tic
% state
xmin = min(guessstate(:,1));
xmax = max(guessstate(:,1));
ymin = min(guessstate(:,2));
ymax = max(guessstate(:,2));
phimin = -pi;
phimax = pi;
dphimin = INPUTS.umin;
dphimax = INPUTS.umax;

param_min = [];
param_max = [];
path_min = [INPUTS.obstacleRadius^2;INPUTS.obstacleRadius^2;INPUTS.obstacleRadius^2];
path_max = xmax^2*ones(1,size(INPUTS.obstacle,1))';
event_min = [];
event_max = [];
duration_min = [];
duration_max = [];

% limits.meshPoints = [-1 +1];
% limits.nodesPerInterval = 10;
limits.time.min = [0 t0];
limits.time.max = [0 1.5*max(guessTime)];
limits.control.min = INPUTS.umin;
limits.control.max = INPUTS.umax;

limits.state.min (1 ,:) = [ x0(1) xmin-10 xT(1)];
limits.state.max (1 ,:) = [ x0(1) xmax+10 xT(1)];
limits.state.min (2 ,:) = [ x0(2) ymin-10 xT(2)];
limits.state.max (2 ,:) = [ x0(2) ymax+10 xT(2)];
limits.state.min (3 ,:) = [ x0(3) phimin xT(3)];
limits.state.max (3 ,:) = [ x0(3) phimax xT(3)];
limits.state.min (4 ,:) = [ x0(4) dphimin xT(4)];
limits.state.max (4 ,:) = [ x0(4) dphimax xT(4)];

limits.parameter.min = param_min ;
limits.parameter.max = param_max ;
limits.path.min = path_min ;
limits.path.max = path_max ;
limits.event .min = event_min ;
limits.event.max = event_max ;
limits.duration.min = duration_min ;
limits.duration.max = duration_max ;

guess.time = guessTime;
guess.state(: ,1) = guessstate(:,1);
guess.state(: ,2) = guessstate(:,2);
guess.state(: ,3) = guessstate(:,3);
guess.state(: ,4) = guessstate(:,4);
guess.control(: ,1) = guessControl;
guess.parameter = [];

setup.name = 'localControl-Problem';
setup.funcs.cost = 'local_controlCost';
setup.funcs.dae = 'local_controlDae';
setup.limits = limits ;
setup.guess = guess ;
setup.derivatives = 'finite-difference';
setup.direction = 'increasing';
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-3;
setup.mesh.maxiterations = 10;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;
[output, gpopsHistory] = gpops(setup ) ;

toc
% assignment the value 
solutionPlot = output.solutionPlot ;
openloopstate = solutionPlot.state;
control = solutionPlot.control(: ,1);
time = solutionPlot.time;
if gpopsHistory(1,end).output.SNOPT_info == 1
    disp('Optimality conditions satisfied .')
    disp('Acceptable Course ')
    disp(['Total Cost of Course:',sprintf('%3.4f', gpopsHistory(1 ,1).output.cost )])
else if gpopsHistory(1,end). output.SNOPT_info == 3
        disp('Optimization sucessful , but best accuracy not acheived .')
        disp('Course is usable , but optimizer had difficulty.')
        disp(['Total Cost of Course : ',sprintf('%3.4 f',gpopsHistory(1 ,1). output.cost )]);
    else if gpopsHistory(1,end). output.SNOPT_info == 41
            disp(' Optimization experienced numerical difficulties')
            
            disp('Current point cannot be improved , though solution has been provided .')
            disp('Use of solution possible , though more analysis is required .')
            disp([ 'Total Cost of Course : ',sprintf('%3.4f', gpopsHistory(1 ,1).output.cost )]);
        else
            disp(' Unsucessfull Optimization .')
            disp('Ignore Output ')
        end
    end
end

% plot figure
hold on
plot(solutionPlot.state(: ,1) , solutionPlot.state(: ,2),'r');
axis equal
title('Optimized local control')
xlabel('x- position')
ylabel('y- position')

figure(2)
plot(solutionPlot.time , solutionPlot.control(: ,1))
title('Heading Derivative in Degrees ')



