% /*! @RRT_trajectory_primitives.m
% *************************************************************************
% <PRE>
% file.name       : RRT_trajectory_primitives.m
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
% input:Nodes:GaoNodes;

% output:
% *************************************************************************
function trajlibrary = RRT_trajectory_primitives(sys, delta)
trajlibrary = [];
figure
hold on
Time = [0 1.5*delta/sys.INPUTS.speed];
N = 9;
list = [-4 -3 -2 -1 0 1 2 3 4];
for i = 1:N
    a = list(i);
    point = [0      0       0       0
             delta  a       0       0];

    traj = trajLibrary_gpopsSolveCar(sys, Time, point, []);
    traj.mark = a;
    trajlibrary = [trajlibrary, traj];
    plot(traj.state(:,1),traj.state(:,2))
    drawnow
end
end







