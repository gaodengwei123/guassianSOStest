% /*! @courseOptguesspoint.m
% *************************************************************************
% <PRE>
% file.name       : courseOptguesspoint.m
% related files   :
% function&ablity : build 2D map and use PRM to find a reasonable path
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/06/16    1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output:
% *************************************************************************
function [sys,nuTime,guesspoint,nuControl] = Optguesspoint(sys,startLocation,endLocation)
if 1
    INPUTS = sys.INPUTS;
    obstacle = INPUTS.obstacle;
    obstacleRadius = INPUTS.obstacleRadius;
    %% build map
    x0 = startLocation(1);
    y0 = startLocation(2);
    xf = endLocation(1);
    yf = endLocation(2);
    % plotMin = min([x0 ,xf ,y0 ,yf , obstacle(:,1)', obstacle(:,2)']);
    plotMin = 0;
    plotMax = max([x0 ,xf ,y0 ,yf , obstacle(:,1)', obstacle(:,2)']);
    resolution = 0.1;% interval between map points
    map_x_length = plotMax-plotMin;
    map_y_length = plotMax-plotMin;
    map_x_grid = map_x_length/resolution;
    map_y_grid = map_y_length/resolution;
    simpleMap = zeros(map_x_grid,map_y_grid);
    for current_x_grid = 1:map_x_grid
        for current_y_grid = 1:map_y_grid
            node = [current_x_grid current_y_grid]*resolution;
            if any(sqrt(sum(((repmat(node,sys.getNumObstacles,1) - obstacle(:,1:2)).^2)'))<obstacleRadius)
                simpleMap(map_y_grid-current_y_grid,current_x_grid) = 1;% I don't known why it must  over turn
            end
        end
    end
%     sys = sys.BuildMap(simpleMap,1/resolution); % build 3D map
    %% PRM algrithm
%     rngState = rng;
    map = robotics.BinaryOccupancyGrid(simpleMap,1/resolution);

    prmSimple = robotics.PRM(map,100); % 100 PRM nodes in map
    prmSimple.ConnectionDistance=10;
    path = findpath(prmSimple,startLocation',endLocation');

    sys.PlotObj.prmSimple = prmSimple;
    sys.PlotObj.path = path;
    sys.PlotObj.simpleMap = simpleMap;
    sys.PlotObj.map = map;
    options.plotPRM = 1;
%     options.plotMap = 0;
    PlotFigure(sys,options);
    
    %% judge the trajectory is good or not
    reply = input('Do you want continue with current guess and radar field ? Y/N [Y]: ', 's');
    if strcmp(reply ,'Y') || strcmp(reply ,'y')
%         close all;
    end
    if ~ strcmp(reply ,'Y') && ~ strcmp(reply ,'y')
        error('Run Cancelled ')
    end
    %% calculate the state for open loop with GPOPS
    [nuTime, guesspoint, nuControl ] = OptGuessEnhancer(sys, INPUTS, path);
    
    % save guessdata
else
    load('guessdata')
end


