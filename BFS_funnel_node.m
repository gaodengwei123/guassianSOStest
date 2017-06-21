% /*! @BFS_funnel_node.m
% *************************************************************************
% <PRE>
% file.name       : BFS_funnel_node.m
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
% BFS_funnel_node is the A* version bfs with considering the funnel and
% P-subs
% *************************************************************************
function Path = BFS_funnel_node(GaoNodes, sys, start,target)
if ~isa(GaoNodes,'GaoNodes')
    error('Nodes must be GaoNodes')
end
Nodes = GaoNodes.Graph;
visit(1).node = start;         % add first node to queue
visit(1).path = start;         % add path to queue
visit(1).weights = 0;
visit(1).TargetCost = norm(Nodes(start).node - Nodes(target).node);
visit(1).StartCost = 0;
visit(1).HeurCost = visit(1).TargetCost;
visit(1).funnel = [];
visit(1).traj = [];
visit(1).cost = [];
visit(1).Sigma = 10*sys.INPUTS.Qk; 
visit(1).Lambda = 10*sys.INPUTS.Qk; 
visit(1).Pdiff = [];
Path = [];
while ~isempty(visit)                  % judge queue is empty or not
    % -----------POP lawest cost node in Queue---------------
    Heurdist = [visit(:).HeurCost];
    [~,i] = min(Heurdist);
    POP.node = visit(i).node;
    POP.path = visit(i).path;
    POP.c = visit(i).StartCost;
    POP.funnel = visit(i).funnel;
    POP.traj = visit(i).traj;
    POP.cost = visit(i).cost;
    POP.Sigma = visit(i).Sigma;
    POP.Lambda = visit(i).Lambda;
    POP.Pdiff = visit(i).Pdiff;
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
        if ~any(POP.path == Newnode)         % new element do not exist current path
            % calculate the Heuristic distance for each link
            TargetCost = norm(Nodes(Newnode).node - Nodes(target).node);
            StartCost = POP.c + Nodes(POP.node).linkweights(g);
            HeurCost = StartCost + TargetCost;%
            
            % update the node and path and cost in visit
            New = Propagate(POP, Newnode, GaoNodes, sys);
            if ~isempty(New)
                visit(end+1).node = New.node;
                visit(end).path = [POP.path, New.node]; % add new node in path
                visit(end).HeurCost = HeurCost;
                visit(end).TargetCost = TargetCost;
                visit(end).StartCost = StartCost;
                visit(end).funnel = [POP.funnel,New.funnel];
                visit(end).traj = [POP.traj,New.traj];
                visit(end).cost = New.cost;
                visit(end).Sigma = New.Sigma;
                visit(end).Lambda = New.Lambda;
                visit(end).Pdiff = [POP.Pdiff;New.Pdiff];
            end
        end
    end
    
    
end
if isempty(Path)
    error('do not find a feasible path')
end

for j = 1:length(POP.Pdiff)-1
    Pb1 = POP.Pdiff(j).P1;
%     Pb2 = POP.Pdiff(j).P2;
    plot(Pb1(1,:),Pb1(2,:),'g','linewidth',2); % entrance
%     plot(Pb2(1,:),Pb2(2,:),'y','linewidth',2); % out
end


end






