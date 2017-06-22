classdef GaoNodes
    
    properties
        node            % position of vertex
        Tag             % tag of vertex
        link            % have edge with this vertex
        linkTag         % linked vertex Tag
        linkNum         % linked vertex number
        linkweights     % linked vertex cost
        trajLibrary     % trajectories in each phrase with open-loop control
        FunnleLibrary   % funnel library with each trajLibrary
        Graph           % my graph
        gains           % feedback gain
    end
    
    methods
        function obj = getEdges(obj,path)
            hold on
            dist = path.ConnectionDistance;   % connected distance
            nodes = path.Nodes;               % position of nodes
            idx = rangesearch(nodes,nodes,dist,'Distance','euclidean');
            for i = 1:path.NumNodes
                nextnode = idx{i};
                nextnode(1) = [];             % remove the first node which is itself
                Next = nodes(nextnode,:);   % position of linked nodes
                NextNum = length(nextnode);
                weights = zeros(1,NextNum);
                for j = 1:NextNum
                    weights(j) = norm(nodes(i,:) - Next(j,:));
                    line([nodes(i,1);Next(j,1)],[nodes(i,2);Next(j,2)],'LineStyle',':','linewidth',1)
                end
                Nodes(i) = buildNodes(obj,nodes(i,:),i,Next,nextnode,length(nextnode),weights);
                plot(Nodes(i).node(1),Nodes(i).node(2),'k.','markersize',5)
%                 text(Nodes(i).node(1),Nodes(i).node(2),num2str(i))
            end
            obj.Graph = Nodes;
        end
        
        function obj = Calculate_trajLibrary(obj,sys,GaoPath, start, target)
%             obj.trajLibrary = trajectory_primitives(obj, GaoPath, sys, start, target);% consider about xT
            obj.trajLibrary = trajectory_primitives2(obj, GaoPath, sys, start, target);% consider about xFrom
        end
        
        function obj = RRT_trajLibrary(obj,sys,delta)
            obj.trajLibrary = RRT_trajectory_primitives(sys, delta);% consider about xFrom
        end
        
        
        function obj = Calculate_FunnelLibrary(obj,sys)
            if isempty(obj.trajLibrary)
                error('must compete trajLibrary first')
            end
            for i = 1:1%length(obj.trajLibrary)
                Funnel(i) = FunnleLibrary_primitives(sys,obj.trajLibrary(i));
            end
            obj.FunnleLibrary = Funnel;
        end
        
        function obj = RRT_FunnelLibrary(obj,sys)
            if isempty(obj.trajLibrary)
                error('must compete trajLibrary first')
            end
            for i = 1:length(obj.trajLibrary)
                Funnel(i) = FunnleLibrary_primitives(sys,obj.trajLibrary(i));
            end
            obj.FunnleLibrary = Funnel;
        end
        
        function obj = FeedbackGain(obj,sys)
%             figure
%             axHandle = gca;
%             hold on
            for i = 1:length(obj.FunnleLibrary)
                x0 = obj.trajLibrary(i).state;
                u0 = obj.trajLibrary(i).control;
                ts = obj.trajLibrary(i).time;
                sys = sys.updateNominal(ts,x0,u0);
                sys.Maxinterval = [ts(1) ts(end)];
                N = length(obj.FunnleLibrary(1).rho);
                sys = sys.timecalculation(N);
                sys = sys.tv_poly_linearize;
%                 %% ===================time-varying Riccati equation=================
%                 Qf = diag([1 1 10 10]); % this is must smaller than ROA
%                 [tv,Ss] = tv_lqr_riccati(sys.Maxinterval,sys.A,sys.B,sys.INPUTS.Q,sys.INPUTS.R,Qf);
%                 Spp = spline(tv,Ss);
%                 S = @(t) ppval(Spp,t);
%                 gain{i} = @(t) -inv(sys.INPUTS.R(t))*sys.B(t)'*S(t);
               obj.FunnleLibrary(i).gain = @(t) -inv(sys.INPUTS.R(t))*[0 0 0 1]*ppval(obj.FunnleLibrary(i).V.Spp,t);
               obj.FunnleLibrary(i).A = sys.A;
               obj.FunnleLibrary(i).B = sys.B;
               % update LF
               V = obj.FunnleLibrary(i).V ;
               rho = obj.FunnleLibrary(i).rho;
               fact = linspace(2,1,20); % factor to enlage funnel for compution error
               rho = (fact.*rho')';
               obj.FunnleLibrary(i).V = V.updateV(foh(V.getbreak,1./rho'));
               obj.FunnleLibrary(i).rho = [];
               
               % spline trajectory
               statefun = spline(ts,x0');
               Traj = polyniminalTrajectory(statefun);
               % spline control
               controlfun = spline(ts,u0');
               u0 = polyniminalTrajectory(controlfun);
               obj.FunnleLibrary(i).Traj = Traj;
               obj.FunnleLibrary(i).u0 = u0;
               
               %% plot
               
%                sys.PlotObj.solutionPlot.state = obj.FunnleLibrary(i).Traj.eval(obj.FunnleLibrary(i).V.getbreak)';
%                options.x0 = obj.FunnleLibrary(i).Traj;
%                %             ShrinkRho = funnel_composition(sys,V,OL_time,OL_state, P2, rho2);
%                if i==1
%                V = obj.FunnleLibrary(i).V.inFrame(Traj);
%                h = plot_myFunnel(sys,V,options);
%                end
%                h0=plot(sys.PlotObj.solutionPlot.state(: ,1), sys.PlotObj.solutionPlot.state(: ,2),'k');
%                title(axHandle, 'Trajectories and Funnel Libraries');
%                xlimits = [-2 6];
%                ylimits = [-6 6];
%                xlabel(axHandle,'X[meters]');
%                ylabel(axHandle,'Y[meters]');
               
            end

        end
        
    end
    methods(Access = private)
        function Nodes = buildNodes(obj,nodes,tag,Next,nextnode,linknum,weights)
            Nodes.node =  nodes;             % position of vertex
            Nodes.Tag = tag;                   % tag of vertex
            Nodes.link = Next;               % have edge with this vertex
            Nodes.linkTag = nextnode';        % linked vertex Tag
            Nodes.linkNum = linknum;         % linked vertex number
            Nodes.linkweights = weights';     % linked vertex cost
        end
    end

end