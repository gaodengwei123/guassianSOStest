classdef CarDynamic<GaoSystem
    methods
        function obj = CarDynamic()
            mark = 1;                           % 1 is the Car 2-D model
            obj = obj@GaoSystem(mark,4,1);
            obj = SetParameter(obj);   % set parameters for the vechicle
            obj.getNumObstacles = size(obj.INPUTS.obstacle,1);   % Number of Obstacle
            obj = setInputLimits(obj,-inf,inf);
        end
        
        function obj = SetParameter(obj)
            INPUTS.state_goal = [50;40;1;0];                    % terminal point
            INPUTS.state_initial = zeros(obj.getNumStates,1);   %initilal points
            INPUTS.speed = 10;
            INPUTS.controlSat = pi;     % control saturation
            
            INPUTS.obstacle = [10 30
                30 40
                25 20];
            INPUTS.obstacleRadius = 4;
            INPUTS.obstacleRadius_multiplier = 2.1;
            
            % INPUTS.state_goal = [5;9;1;0];
            
            % cost function Q and R also be used in local_controlCost
            INPUTS.Q = @(t)diag([1 1 10 10]);
            INPUTS.R = @(t)10*eye(obj.getNumInput);
            INPUTS.filter_flag = 2;         %1:set membership filter,2:Kalman filter
            INPUTS.plot_elliposoid = 2;     %1:plot elliposoid in the path;0:don't plot elliposoid
            % noise set
            INPUTS.noise_w = 0.1*[1,1,0.1,0.01]';
            INPUTS.noise_v = 0.1*[1,1,0.1]';
            INPUTS.filter_P0 = diag(INPUTS.noise_w.^2);
            INPUTS.Qk = diag(INPUTS.noise_w.^2);
            INPUTS.Rk = diag(INPUTS.noise_v.^2);
            obj.INPUTS = INPUTS;
        end
        
        function [xdot,df,d2f,d3f ]= dynamics(obj,t,x,u)
            % x = [ x; y; theta; thetadot ]
            % u = [ thetaddot ]
            %             Car_controlDae must chnage when change the dynamic
            INPUTS = obj.INPUTS;
            theta = x(3);
            xdot = [ INPUTS.speed * cos(theta);
                INPUTS.speed * sin(theta);
                x(4);
                u(1) ];
            
            if (nargout>1)
                [df,d2f,d3f]= dynamicsGradients(obj,t,x,u,nargout-1);
            end
        end
        % systme output
        function y = output(obj,t,x,u)
            y = x;
        end
        
        function [z,H,Obs] = measurement(obj,t,x,u)
            if ~isa(obj.INPUTS.field,'ObstacleField')
                error('add field before get measurement')
            else
                Obs = NearSensor(obj.INPUTS.field,x(1:2));
            end
            
            if Obs
                z = x([1;2;3]);    % x y theta
                H = [eye(3) zeros(3,1)];
            else
                z=[];
                H=[];
            end
        end
        
        function A = dynamicA(obj, t, x, u)
            [~,df,~,~] = dynamics(obj,t,x,u);
            A = full(df(:,2:5));
        end
        function B = dynamicB(obj, t, x, u)
            [~,df,~,~] = dynamics(obj,t,x,u);
            B = full(df(:,6));
        end
    end
    
    methods (Access=private)
        function [df,d2f,d3f]= dynamicsGradients(obj, t, a3, u, order)
            % Symbol table:
            a3_3=a3(3);
            
            df = sparse(4,6);
            df(1,4) = -obj.INPUTS.speed*sin(a3_3);
            df(2,4) = obj.INPUTS.speed*cos(a3_3);
            df(3,5) = 1;
            df(4,6) = 1;
            
            % d2f
            if (order>=2)
                d2f = sparse(4,36);
                d2f(1,22) = -obj.INPUTS.speed*cos(a3_3);
                d2f(2,22) = -obj.INPUTS.speed*sin(a3_3);
            else
                d2f=[];
            end
            
            % d3f
            if (order>=3)
                d3f = sparse(4,216);
                d3f(1,130) = obj.INPUTS.speed*sin(a3_3);
                d3f(2,130) = -obj.INPUTS.speed*cos(a3_3);
            else
                d3f=[];
            end
            if (order>=4)
                error('add Gradients here')
            end
        end
    end
end