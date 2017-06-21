% /*! @DFS_gao_node.m
% *************************************************************************
% <PRE>
% file.name       : DFS_gao_node.m
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
function Path = DFS_gao_node(GaoNodes,start,target)
if ~isa(GaoNodes,'GaoNodes')
    error('Nodes must be GaoNodes')
end
Nodes = GaoNodes.Graph;
top = 1;                         % top of stack
k = 1;
stack(top).node = start;         % add first node to stack
stack(top).path = start;
stack(top).weights = 0;
while top ~= 0                  % judge stack is empty or not
    if stack(top).node == target
        Path{k}.path = stack(top).path;
        Path{k}.cost = stack(top).weights;
        k = k+1;
        stack(top) = [];
        top=top-1;
        if top == 0 
            break
        end
        continue;
    end
    i = stack(top).node;
    toppath = stack(top).path;
    topweights = stack(top).weights;
    stack(top) = [];
    top=top-1;

    for g = 1:Nodes(i).linkNum
        if ~any(toppath==Nodes(i).linkTag(g))
            top=top+1;
            stack(top).node = Nodes(i).linkTag(g);
            if ~isfield(stack(top),'path')
                stack(top).path = toppath;
            end
            if ~isfield(stack(top),'weights')
                stack(top).weights = topweights;
            end
            stack(top).path = [stack(top).path, stack(top).node];
            stack(top).weights = stack(top).weights+Nodes(i).linkweights(g);
        end
    end
end
end






