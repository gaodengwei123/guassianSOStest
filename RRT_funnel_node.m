% /*! @RRT_funnel_node.m
% *************************************************************************
% <PRE>
% file.name       : RRT_funnel_node.m
% related files   :
% function&ablity :
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/03/6     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:Nodes:GaoNodes; start: is the position of start node; target: is
% the position of target node

% output:
% *************************************************************************
function FindPath = RRT_funnel_node(GaoNodes, sys, start,target)
if ~isa(GaoNodes,'GaoNodes')
    error('Nodes must be GaoNodes')
end
M = 10000;
i = 1;
Vertex = struct([]);
Vertex(1).pos = start;                % the position of vertex
Vertex(1).Tag = 1;
% Vertex(1).Graph = struct([]);       % graph from start to current vertex
% Vertex(1).Path = [];
Vertex(1).Ftuple = FunnelTuple();
Vertex(1).path = [];
Vertex(1).Sigma = sys.INPUTS.Qk;
Vertex(1).Lambda = sys.INPUTS.Qk;
Vertex(1).t0 = 0;
Vertex(1).cost = 0;
Vertex(1).FunnelHandle = [];
FindPath = [];
randHandle = [];
% add the Queue
top = 1;
Queue = struct([]);
% this region can not reach
% AvoidRegin = [0 0;35.5 0;0 30;25 50;0 50;0 0];
AvoidRegin = [0 0;40 0;0 30;0 0];
AvoidRegin1 = [0 40;15 50;0 50;0 40];
hold on
% random sample
Randomp = haltonset(2);
Randomp = scramble(Randomp,'RR2');
X0 = net(Randomp,M+1);

rj = 2;
while i<M                 % judge queue is empty or not
%     delete(randHandle)
    % random sample
    
    x_rand = [80*X0(rj,1) 50*X0(rj,2)];
    rj = rj+1;
    [insideAVoid, ~] = inpolygon(x_rand(1),x_rand(2),AvoidRegin(:,1),AvoidRegin(:,2));
    [insideAVoid2, ~] = inpolygon(x_rand(1),x_rand(2),AvoidRegin1(:,1),AvoidRegin1(:,2));
    if insideAVoid==1
        x_rand = [80 50*rand(1)];
    end
    if insideAVoid2==1
        x_rand = target;
    end
%     randHandle = plot(x_rand(1),x_rand(2),'k*');
    Vpos = vec2mat([Vertex(:).pos],2); % pos of all vertex
    [idx_near,~]=knnsearch(Vpos, x_rand); % this will return the first one in The Tree(may be more near state)
    x_near = Vertex(idx_near);  % find the nearest node in the tree
    
    % if any predecessor of target in tree
    x_grow = repmat(target,9,1)+[-5*ones(9,1) linspace(-4,4,9)'];
    % find the path
    IndexTarget = ismember(Vpos,x_grow,'rows');
    if any(IndexTarget)
        if i>300  % after sufficient optimal
            % find the best predecessor
            bestFather = BestPredecessor(Vertex, IndexTarget);
%             FindPath = [Bestparent.parent;target];
            FindPath = plotpath(Vertex, bestFather, GaoNodes, target);
            % here need a shift funnel----
%             Pathplot = plot(FindPath(:,1),FindPath(:,2),'g','LineWidth',3);
            % -----------------------------
            continue
        end
    end
    % growing a vertex in random
    x_new.pos = x_near.pos+[5,randi([-4,4])];
    x_new.Tag = i+1;
    
    % find the dyanmic trajectory in trajectory library
    MarkIndex = x_new.pos(2) - x_near.pos(2);
    Mark = [GaoNodes.trajLibrary(:).mark];
    index = find(Mark==MarkIndex);
    Ftuple = ShitFunnel(GaoNodes.FunnleLibrary(index),GaoNodes.trajLibrary(index), x_near);
    x_new.Ftuple = Ftuple;
    [x_new,flag]= CheckFunnelcollision(x_new, sys.INPUTS.field);
    % if traj is collision with obs continue
    if flag
        continue
    end
    
    NewVertex = RRTPropagate(x_near, x_new, sys);
    % we can guarantee the states are inside funnel or not
    if isempty(NewVertex)
        continue
    end
    
    % push postion in queue
    Queue = [Queue, x_near];
    
    % find an optimal trajectory from parents£¨update Tree£©
    while ~isempty(Queue)
        % -----------POP node in Queue---------------
        POP = Queue(1);
        Queue(1) = [];
        top=top-1;
        % find all the nearest predecessors of newnode this is a Near
        % choice like RRT*
        x_Parents = repmat(NewVertex.pos,9,1)+[-5*ones(9,1) linspace(-4,4,9)'];
        IndexinTree = ismember(Vpos,x_Parents,'rows');
        fatherInedx = find(IndexinTree==1);  % the index of newnode's father in the tree
        
        % maybe more than 9 fathers since the difference path
        for j = 1:length(fatherInedx)
            x_father = Vertex(fatherInedx(j));
            
            % skip the same parent
            Path1 = x_father.path;
            Path2 = POP.path;
            if length(Path1)==length(Path2)||~all(x_father.parent == POP.parent)
                % find index in funnel library
                LI =  NewVertex.pos(2)-x_father.pos(2) + 5;
                % find a shorter path
                if x_father.cost + GaoNodes.trajLibrary(LI).cost < NewVertex.cost
                    % find a better link
                    x_link.pos = NewVertex.pos;
                    x_link.Tag = NewVertex.Tag;
                    Ftuple = ShitFunnel(GaoNodes.FunnleLibrary(LI),GaoNodes.trajLibrary(LI), x_father);
                    x_link.Ftuple = Ftuple;
                    [x_link,flag]= CheckFunnelcollision(x_link, sys.INPUTS.field);
                    % if traj is collision with obs continue
                    if flag
                        continue
                    end
                    
                    x_link = RRTPropagate(x_father, x_link, sys);
                    % we canguarantee the states are inside funnel or not
                    if ~isempty(x_link)%&&isbigger(x_link.Ftuple,NewVertex.Ftuple)  % if this PD funnel is also bigger than exist
                        % replace old node
                        delete(NewVertex.FunnelHandle)
                        NewVertex = x_link;
                        % push into queue
                        Queue(end+1) = x_father;
                    end
                end
            end
        end
    end
    
    i=i+1;
    % growing the tree
    Vertex(i).pos = NewVertex.pos;
    Vertex(i).Tag = NewVertex.Tag;
    Vertex(i).Ftuple = NewVertex.Ftuple;
    Vertex(i).path = NewVertex.path;
    Vertex(i).Sigma = NewVertex.Sigma;
    Vertex(i).Lambda = NewVertex.Lambda;
    Vertex(i).cost = NewVertex.cost;
    drawnow
%     plot(NewVertex.pos(1),NewVertex.pos(2),'.r')
%     plot(NewVertex.Ftuple)
end

if isempty(FindPath)
    error('do not find a feasible path')
end
end










