% /*! @BRM_gao_node.m
% *************************************************************************
% <PRE>
% file.name       : BRM_gao_node.m
% related files   :
% function&ablity : 
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       : 
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/02/28    1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:Nodes:GaoNodes;  start: start node;   target: goal node

% output: 
% BRM_gao_node is the A* version bfs
% *************************************************************************
function Path = BRM_gao_node(GaoNodes,start,target)
if ~isa(GaoNodes,'GaoNodes')
    error('Nodes must be GaoNodes')
end
Nodes = GaoNodes.Graph;
close_set = start; %11 15
visit(1).node = start;         % add first node to queue
visit(1).path = start;         % add path to queue
visit(1).weights = 0;
visit(1).TargetCost = norm(Nodes(start).node - Nodes(target).node);
visit(1).StartCost = 0;
visit(1).HeurCost = visit(1).TargetCost;
while ~isempty(visit)                  % judge queue is empty or not
    % -----------POP lawest cost node in Queue---------------
    Heurdist = [visit(:).HeurCost];
    [~,i] = min(Heurdist);
    POP.node = visit(i).node;
    POP.path = visit(i).path;
    POP.c = visit(i).StartCost;
    visit(i) = [];
    % --------------------------------------------
    % if arrive the target we need to record path
    if POP.node == target                      
        Path.path = POP.path;
        Path.cost = POP.c;
        break    
    end
    % --------------------------------------------------
    % -----------find----------------
    for g = 1:Nodes(POP.node).linkNum
        Newnode = Nodes(POP.node).linkTag(g);
        if ~any(close_set == Newnode)         % new element do not reach current path
            close_set = [close_set, Newnode];
            % calculate the Heuristic distance for each link
            TargetCost = norm(Nodes(Newnode).node - Nodes(target).node);
            StartCost = POP.c + Nodes(POP.node).linkweights(g);
            HeurCost = StartCost + TargetCost;%
            
            % update the node and path and cost in visit
            visit(end+1).node = Newnode;
            visit(end).path = [POP.path, Newnode]; % add new node in path
            visit(end).HeurCost = HeurCost;
            visit(end).TargetCost = TargetCost;
            visit(end).StartCost = StartCost;
        end
    end

    
end
end






