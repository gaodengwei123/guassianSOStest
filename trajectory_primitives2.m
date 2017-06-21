function trajlibrary = trajectory_primitives2(GaoGraph, GaoPath, sys, first_vertex, end_vertex)

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
    
    [Time, point,Flag_feasible] = trajNode2(sys,x0,xT,[],[]);
    if Flag_feasible==1  % if this tarjectory is feasible
        traj = trajLibrary_gpopsSolveCar(sys, Time, point, []);
        traj.mark = [0,first_vertex,xTTag];
        trajlibrary = [trajlibrary, traj];
        plot(traj.state(:,1),traj.state(:,2))
        drawnow
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
        for t = 1:Nodes(i).linkNum
            xfromTag = Nodes(i).linkTag(t);
            if xfromTag~=xTTag    % donot the target
                xform = Nodes(i).link(t,:);
                [Time, point,Flag_feasible] = trajNode2(sys,x0,xT,xform,[]);
                if Flag_feasible==1  % if this tarjectory is feasible
                    traj = trajLibrary_gpopsSolveCar(sys, Time, point, []);
                    traj.mark = [xfromTag,x0tag,xTTag];
                    trajlibrary = [trajlibrary, traj];
                    plot(traj.state(:,1),traj.state(:,2))
                    drawnow
                else
                    continue
                end
            else
                continue
            end
        end
    end
end







