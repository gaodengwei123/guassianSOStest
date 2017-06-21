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
function [sys,openloopstate,control,time] = optimal_local_programming(sys,x0,current_ell)
global CONSTANS
xT = sys.INPUTS.state_goal;
CONSTANS = sys.INPUTS;
CONSTANS.PRMflag = 1;
% ellipsoid major axis
[currentq,currentQ] = double(current_ell);
majorAxis  = sqrt(min(diag(currentQ)));
% adjust the multiplier to make the path safer
if CONSTANS.obstacleRadius_multiplier<majorAxis
    CONSTANS.obstacleRadius_multiplier = majorAxis;
end
%% =======================gpops-II============================begin
% range of time & states & control output (guess)
 % tf = norm(xT-x0)/10;


%% =============== initial trajectory guess =======================
%build map and guess initial points
[sys,guessTime,guessstate,guessControl] = Optguesspoint(sys,x0(1:2),xT(1:2));
%%
tic
switch (sys.mark)
    case 1  % car
        [output, gpopsHistory] = gpopsSolveCar(guessTime,guessstate,guessControl,x0,xT);
    case 2  % videoRay
        [output, gpopsHistory] = gpopsSolveVideoRay(guessTime,guessstate,guessControl,x0,xT);
end
toc
% assignment the value 
solutionPlot = output.solutionPlot ;
openloopstate = solutionPlot.state;
control = solutionPlot.control;
time = solutionPlot.time;


%%
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
sys.PlotObj.solutionPlot = solutionPlot;




