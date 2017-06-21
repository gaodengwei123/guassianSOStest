function trajlibrary = trajectory_primitives(GaoGraph, GaoPath, sys, first_vertex, end_vertex)

% PRM in map and trajectory from one vertex to another
Nodes = GaoGraph.Graph;     % graph
vertex = GaoPath.Nodes;     % vertex
trajlibrary = [];
hold on
%% the first vertex================================================
x0 = Nodes(first_vertex).node;
for i = 1:Nodes(first_vertex).linkNum
    xTTag = Nodes(first_vertex).linkTag(i);
    xT = Nodes(first_vertex).link(i,:);
    for j = 1:Nodes(xTTag).linkNum
        xTo = Nodes(xTTag).link(j,:);
        xToTag = Nodes(xTTag).linkTag(j);
        [Time, point,Flag_feasible] = trajNode(sys,x0,xT,[],xTo);
        if Flag_feasible==1  % if this tarjectory is feasible
            traj = trajLibrary_gpopsSolveCar(sys, Time, point, []);
            traj.mark = [first_vertex,xTTag,xToTag];
            trajlibrary = [trajlibrary, traj];
            plot(traj.state(:,1),traj.state(:,2))
            drawnow
        end
    end
end
%% =============================================================
for i = 1:length(Nodes)
    if i == first_vertex
        continue
    end
    x0 = Nodes(i).node;
    x0tag = i;
    % x0--->xT
    for k = 1:Nodes(i).linkNum
        xT = Nodes(i).link(k,:);
        xTTag = Nodes(i).linkTag(k);
        if xTTag == end_vertex
            xTo = [];
            [Time, point,Flag_feasible] = trajNode(sys,x0,xT,1,xTo);
            if Flag_feasible==1  % if this tarjectory is feasible
                traj = trajLibrary_gpopsSolveCar(sys, Time, point, []);
                traj.mark = [x0tag,xTTag,0];
                trajlibrary = [trajlibrary, traj];
                plot(traj.state(:,1),traj.state(:,2))
                drawnow
            else
                continue
            end
        else
            % if not the endvertex
            for t = 1:Nodes(xTTag).linkNum
                xToTag = Nodes(xTTag).linkTag(t);
                if xToTag == x0tag % cannot be a circle
                    continue
                else
                    xTo = Nodes(xTTag).link(t,:);
                end
                [Time, point,Flag_feasible] = trajNode(sys,x0,xT,1,xTo);
%                 if norm(xTo-Nodes(end_vertex).node)>norm(xT-Nodes(end_vertex).node) % xto must near to the end
%                     Flag_feasible = 0;
%                 end
                if Flag_feasible==1  % if this tarjectory is feasible
                    traj = trajLibrary_gpopsSolveCar(sys, Time, point, []);
                    traj.mark = [x0tag,xTTag,xToTag];
                    trajlibrary = [trajlibrary, traj];
                    plot(traj.state(:,1),traj.state(:,2))
                    drawnow
                else
                    continue
                end
            end
        end
    end
end
end







