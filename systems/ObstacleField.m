classdef ObstacleField
    
    properties
        obstacles;               % cell array of obstacle instances
        number_of_obstacles = 0  % number of obstacles currently in the field
        Map
        resolution
        range
        sensor
    end
    
    methods
        function obj = GenerateObstacles(obj,vertex,OBS_Position,NUM_OBSTACLES)
            %% --------------generate obs randomly----------------------------
            if isempty(OBS_Position)
                x_low = min(vertex(:,1));
                x_high = max(vertex(:,1));
                z_low = min(vertex(:,2));
                z_high = max(vertex(:,2));
                OBS_Position = [rand(1,NUM_OBSTACLES)*(x_high-x_low)+x_low; rand(1,NUM_OBSTACLES)*(z_high-z_low)+z_low]';
            end
            if size(OBS_Position,1)~=NUM_OBSTACLES
                error('number of OBS must same with pos')
            end
            %%----------------------------------------------------------------------
            OBS_MIN_radius = 3;
            OBS_MAX_radius = 3;
            OBS_MIN_EDGE = 4;
            OBS_MAX_EDGE = 5;
            %%
            % figure
            % hold on
            for i = 1:NUM_OBSTACLES
                obs_num_edges = floor((OBS_MAX_EDGE-OBS_MIN_EDGE)*rand+OBS_MIN_EDGE);
                radius = (OBS_MAX_radius-OBS_MIN_radius)*rand+OBS_MIN_radius;
                
                obj= AddObstacle(obj,OBS_Position(i,:),radius,obs_num_edges);
            end
        end
        
        % add obstacles
        function obj = AddObstacle(obj, OBS_Position, radius, obs_num_edges)
            % add polygon obstacles
            obs_angles = 2*pi*rand(1,obs_num_edges);
            obs_coords = [OBS_Position(1)+radius*cos(obs_angles)
                OBS_Position(2)+radius*sin(obs_angles)];
            obj.number_of_obstacles = obj.number_of_obstacles+1;
            obj.obstacles{obj.number_of_obstacles} =obs_coords(:,convhull(obs_coords(1,:)',obs_coords(2,:)'));
        end
        
        function obj = buildQRMap(obj, range, resolution)
            %             resolution: interval between map points
            map_x_length = range(3);
            map_y_length = range(4);
            map_x_grid = map_x_length/resolution;
            map_y_grid = map_y_length/resolution;
            simpleMap = zeros(map_x_grid,map_y_grid);
            X = 1:map_x_grid;
            Y = 1:map_y_grid;
            [xx,yy]=meshgrid(Y,X); % grid in simplemap
            x = xx*resolution; y = yy*resolution;% position in simplemap
            
            % if point is inside polygon map will set to 1
            for i = 1:obj.number_of_obstacles
                simpleMap = InsideOBS(obj, simpleMap, x(:),y(:),obj.obstacles{i});
            end
            obj.Map = flipud(simpleMap);
            obj.range = range;
            obj.resolution = resolution;
        end
        
        function obj = addsensors(obj,pos,radius)
            % add obsever in mapfield
            % pos: is the postion of sensor
            % radius£º is the radius of sensor
            sensorNum = size(pos,1);
            for i = 1:sensorNum
                obj.sensor{i}.pos = pos(i,:);
                obj.sensor{i}.radius = radius;
            end
        end
        
        function axHandle = show(obj, varargin)
            axHandle = gca;
            hold on
            xlimits = [obj.range(1) obj.range(3)];
            ylimits = [obj.range(2) obj.range(4)];

            xlabel(axHandle,'X[meters]');
            ylabel(axHandle,'Y[meters]');
            axHandle.XLim = xlimits;
            axHandle.YLim = ylimits;
            axHandle.Visible = 'on';
            for i = 1:obj.number_of_obstacles
                obstacle = obj.obstacles{i};
                xv = obstacle(1,:)';
                yv = obstacle(2,:)';
                fill(xv,yv,'k')
            end
            % plot sensors if exist
            if ~isempty(obj.sensor)
                for i = 1:length(obj.sensor)
                    plot(obj.sensor{i}.pos(1),obj.sensor{i}.pos(2),'b.','markersize',20);
                    circle(obj,i);
                end
            end
            
        end
        
        function flag = NearSensor(obj,pos)
            % judge the pos is near sensors
            flag = 0;
            for i = 1:length(obj.sensor)
                if norm(pos(1:2)' - obj.sensor{i}.pos)<obj.sensor{i}.radius
                    flag = 1;
                end
            end
        end
        
        function obj = Map2Obstacle(obj)
            obj.number_of_obstacles =3;
            obj.obstacles{1} = [14.3692	14.9533
15.5374	15.8879
14.0187	19.6262
14.486	19.9766
16.0047	16.8224
17.8738	18.2243
19.1589	16.7056
15.6542	13.4346
14.3692	14.9533]';
            obj.obstacles{2} = [16.9393	23.5981
18.2243	24.6495
17.072	24.7664
17.072	26.8692
18.8084	26.8692
18.9252	25.1168
25.5841	28.8551
25.3505	30.6075
20.5607	30.6075
20.5607	31.1916
26.285	31.3084
26.285	27.8037
27.6869	24.6495
19.1589	19.5093
16.9393	23.5981]';
            obj.obstacles{3} = [32.5935	30.4907
35.9813	31.7757
37.0327	28.3879
33.7617	26.8692
32.5935	30.4907]';
        end
        
        function axHandle = showMap(obj, varargin)
            if nargin>1
                Resolution = varargin{1};
            else
                Resolution = 1/obj.resolution;
            end
            axHandle = gca;
            hold on
            if isempty(axHandle)
                axHandle = newplot;
            end
            cmap = colormap('gray');
            % Flip color map to make unoccupied (free) cells white
            cmap = flip(cmap);
            
            % Display the grid map
            imghandle = imshow(obj.Map, 'Parent', axHandle,'InitialMagnification', 'fit');
            colormap(axHandle, cmap)
            title(axHandle, 'ObstacleField');
            xlimits = [obj.range(1) obj.range(3)];
            ylimits = [obj.range(2) obj.range(4)];
            correction = 1/(2*Resolution);
            
            % Set XData and YData
            if (abs(xlimits(1)-xlimits(2)+2*correction) < eps)
                % Special case when there is only one cell
                imghandle.XData = [xlimits(1), xlimits(2)];
            else
                imghandle.XData = [xlimits(1)+correction, xlimits(2)-correction];
            end
            
            if (abs(ylimits(1)-ylimits(2)+2*correction) < eps)
                imghandle.YData = [ylimits(2), ylimits(1)];
            else
                imghandle.YData = [ylimits(2)-correction, ylimits(1)+correction];
            end
            
            xlabel(axHandle,'X[meters]');
            ylabel(axHandle,'Y[meters]');
            
            % Set the axes
            set(axHandle, 'YDir','normal');
            grid(axHandle, 'off');
            % Set XLim and YLim
            axHandle.XLim = xlimits;
            axHandle.YLim = ylimits;
            %
            %                         Make axes visible
            axHandle.Visible = 'on';
            
            % plot sensors if exist
            if ~isempty(obj.sensor)
                for i = 1:length(obj.sensor)
                    plot(obj.sensor{i}.pos(1),obj.sensor{i}.pos(2),'b.','markersize',20);
%                     circle(obj,i);
                end
            end
            
        end
        
    end
    
    methods (Access = private)
        function simpleMap = InsideOBS(obj, simpleMap,xq,yq,obstacle)
            %             jugdge if the nodes are in side the obstacle polygon
            xv = obstacle(1,:)';
            yv = obstacle(2,:)';
            [in, on] = inpolygon(xq,yq,xv,yv);
            simpleMap(in) = 1;
            simpleMap(on) = 1;
        end
        
        function h = circle(obj,i)
            x0 = obj.sensor{i}.pos(1);
            y0 = obj.sensor{i}.pos(2);
            r = obj.sensor{i}.radius;
            theta=0:pi/50:2*pi;
            x=x0+r*cos(theta);
            y=y0+r*sin(theta);
            h = plot(x,y,'b--');
        end
    end
    
    
    
    
    
    
end
