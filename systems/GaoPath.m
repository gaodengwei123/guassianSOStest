classdef GaoPath
    
    properties
        field               % Gao map
        NumNodes            % nodes number for random roadmap
        ConnectionDistance  % Connection Distance
        path                % path
        cost                % cost of path
        Nodes               % vertex in field
    end
    
    methods
        function obj = setfield(obj,field, numnodes, ConnectionDistance)
            %¡¡numnodes£ºthe number of vertex
            % ConnectionDistance: connect distance
            if ~isa(field,'ObstacleField')
                error('field must be ObstacleField')
            end
            obj.field = field;
            obj.ConnectionDistance = ConnectionDistance;  % Connection Distance
            obj.NumNodes = numnodes;
            obj = generateNodes(obj, field, numnodes);
        end
        
        function obj = findpath(obj, GaoNodes, sys, start, target, method)
            % Find path
            switch method
                case 'Astar'
                    Findpath = BRM_gao_node(GaoNodes,start,target);
                case 'BFS'
                    Findpath = BFS_gao_node(GaoNodes,start,target);
                case 'DFS'
                    Findpath = DFS_gao_node(GaoNodes,start,target);
                case 'FunnelAstar'
                    Findpath = BFS_funnel_node(GaoNodes, sys, start,target);
                case 'FunnelRRT'
                    Findpath = RRT_funnel_node(GaoNodes, sys, start,target);
                otherwise
                    error('choose a method')
            end
            obj.path = Findpath.path;
            obj.cost = Findpath.cost;
        end
    end
    
    methods(Access = private)
        function obj = generateNodes(obj, field, numnodes)
            % PRM nodes are ranom node
            % obj.Nodes = generateprmNode(numnodes, field, obj.ConnectionDistance);
            obj.Nodes = NodefromBRM();
        end
    end
end