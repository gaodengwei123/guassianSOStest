classdef FunnelTuple
    
    properties
        V           % is a Vfunction class with lyapunv function
%         rho         % the radius of funnel in different time
        traj        % the spline open-loop trajectory
        K           % feedback gain
        u0          % open-loop control
        cost        % cost of this trajectory
        breaks      % time
        A           % coefficient of sys after linearization
        B           % coefficient of sys after linearization
        Pdiff       % % Pontryagin difference
    end
    
    methods
        function obj = setFunnelTuple(obj,V,traj,A,B,K,u0,cost)
            obj.V = V;          % V function
%             obj.rho = rho;      % 1*N vec
            obj.traj = traj;    % polyniminalTrajectory
            obj.K = K;          % function handle
            obj.u0 = u0;        % polyniminalTrajectory
            obj.cost = cost;    % scalar
            obj.breaks = obj.V.getbreak;
            obj.A = A;
            obj.B = B;
        end
        
        function h = plot(obj)
            sys.PlotObj.solutionPlot.state = obj.traj.eval(obj.V.getbreak)';
            options.x0 = obj.traj;
%             ShrinkRho = funnel_composition(sys,V,OL_time,OL_state, P2, rho2);
            h = plot_myFunnel(sys,obj.V,options);
        end
        
        function flag = isbigger(F1,F2)
            if isempty(F1.Pdiff)||isempty(F2.Pdiff)
                error('do not have Pontryagin difference')
            end
            xq = F1.Pdiff(1,:);
            yq = F1.Pdiff(2,:);
            xv = F2.Pdiff(1,:);
            yv = F2.Pdiff(2,:);
            [in, ~] = inpolygon(xq,yq,xv,yv);
            flag1 = sum(in);
            flag = ~flag1;
        end
        
    end
    methods(Access = private)

    end
end