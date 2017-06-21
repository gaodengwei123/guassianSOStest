% /*! @BFS_gao_node.m
% *************************************************************************
% <PRE>
% file.name       : BFS_gao_node.m
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
% *************************************************************************
function Path = BFS_gao_node(GaoNodes,start,target)
if ~isa(GaoNodes,'GaoNodes')
    error('Nodes must be GaoNodes')
end
Nodes = GaoNodes.Graph;
top = 1;                         % top of stack
k = 1;
Queue(top).node = start;         % add first node to queue
Queue(top).path = start;         % add path to queue
Queue(top).weights = 0;
while top ~= 0                   % judge queue is empty or not
    % -----------POP node in Queue---------------
    POP.node = Queue(1).node;
    POP.path = Queue(1).path;
    POP.c = Queue(1).weights;
    Queue(1) = [];
    top=top-1;
    % --------------------------------------------
    % if arrive the target we need to record path
    if POP.node == target                      
        Path{k}.path = POP.path;
        Path{k}.cost = POP.c;
        k = k+1;
        continue    % find other possible path
    end
    % --------------------------------------------
    % -----------add new elment in queue----------------
    for g = 1:Nodes(POP.node).linkNum
        if ~any(POP.path == Nodes(POP.node).linkTag(g)) % new element do not reach by path
            Newnode = Nodes(POP.node).linkTag(g);
            % update the node and path and cost in end queue
            Queue(end+1).node = Newnode;
            top = top+1;
            Queue(end).path = [POP.path, Newnode]; % add new node in path
            Queue(end).weights = POP.c+Nodes(POP.node).linkweights(g); % update cost
        end
    end
end
end






