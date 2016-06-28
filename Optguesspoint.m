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
function [nuTime,guesspoint,nuControl] = Optguesspoint(startLocation,endLocation)
global INPUTS
obstacle = INPUTS.obstacle;
obstacleRadius = INPUTS.obstacleRadius;
%% build map
x0 = startLocation(1);
y0 = startLocation(2);
xf = endLocation(1);
yf = endLocation(2);
plotMin = min([x0 ,xf ,y0 ,yf , obstacle(:,1)', obstacle(:,2)']);
plotMax = max([x0 ,xf ,y0 ,yf , obstacle(:,1)', obstacle(:,2)']);
resolution = 0.1;% interval between map points
map_x_length = plotMax-plotMin;
map_y_length = plotMax-plotMin;
map_x_grid = map_x_length/resolution;
map_y_grid = map_y_length/resolution;
simpleMap = zeros(map_x_grid,map_y_grid);
% three obstacles
for current_x_grid = 1:map_x_grid
    for current_y_grid = 1:map_y_grid
        node = [current_x_grid current_y_grid]*resolution;
        if norm(node-obstacle(1,1:2))< obstacleRadius||norm(node-obstacle(2,1:2))< obstacleRadius||norm(node-obstacle(3,1:2)) < obstacleRadius
            simpleMap(map_y_grid-current_y_grid,current_x_grid) = 1;% I do known why it must  over turn
        end
    end
end

%% PRM algrithm
map = robotics.BinaryOccupancyGrid(simpleMap,1/resolution);
prmSimple = robotics.PRM(map,100); % 100 PRM nodes in map
figure(1)
show(prmSimple)
path = findpath(prmSimple,startLocation',endLocation'); 
show(prmSimple)

% map_root = [0 0];
% planner = prm;
% planner = planner.path_planning(simpleMap, map_root, 1/resolution, startLocation', endLocation');
% %% Plot
% % Save and Load Map Image
% figure(99); clf(99); set(gcf,'Visible', 'off'); hold on;
% contour(flipud(planner.map)); %%
% box off; axis off;
% set(gca,'position',[0 0 1 1],'units','normalized');
% print('map','-dpng');
% p = figure(1);
% clf(1);
% img = imread('map.png');
% hold on;
% box on; axis on;
% set(gcf, 'position', [50 50 length(map(1,:))/length(map(:,1))*500 500]);
% imagesc([planner.map_root(1) planner.map_resolution*planner.map_x_grid+planner.map_root(1)],...
%         [planner.map_root(2) planner.map_resolution*planner.map_y_grid+planner.map_root(2)], img);
% % Plot Edges
% gplot(planner.adjacencyMatrix,planner.nodeList,'g'); %drawnow; 
% % Plot Nodes
% scatter(planner.nodeList(:,1), planner.nodeList(:,2), 'k.');
% %% Plot Path
% plot(start(1),start(2),'k.','MarkerSize',30);
% plot(goal(1),goal(2),'r.','MarkerSize',30);
% if ~isempty(planner.path)
%     plot(planner.path(:,1), planner.path(:,2),'k','LineWidth',1.5);
% else
%     disp('No Path Found!');
% end


%% judge the trajectory is good or not 
reply = input('Do you want continue with current guess and radar field ? Y/N [Y]: ', 's');
if ~ strcmp(reply ,'Y') && ~ strcmp(reply ,'y')
    error('Run Cancelled ')
end
%% calculate the state for open loop
[nuTime, guesspoint, nuControl ] = OptGuessEnhancer(INPUTS, path );


