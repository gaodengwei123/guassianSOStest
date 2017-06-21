% /*! @generateprmNode.m
% *************************************************************************
% <PRE>
% file.name       : generateprmNode.m
% related files   :
% function&ablity : generate prm Nodes
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/2/25     1.00        dengwei      build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:
%
% output:
% *************************************************************************
function Nodes = generateprmNode(numnodes, field, Distance)

obj.nodeList = zeros(numnodes-2, 2);
obj.map_x_length = obj.map_x_axis(2) - obj.map_x_axis(1);
obj.map_y_length = obj.map_y_axis(2) - obj.map_y_axis(1);
while (obj.numNodesinmap <= (obj.numNodes-2))
    newNode = rand(1,2);
    newNode(1) = newNode(1) * obj.map_x_length + obj.map_root(1);
    newNode(2) = newNode(2) * obj.map_y_length + obj.map_root(2);
    if ~onObstacle(obj,newNode) % not on obstacle add new node to the list
        obj.nodeList(obj.numNodesinmap, :) = newNode(1:2);
        obj.numNodesinmap = obj.numNodesinmap + 1;
    end
end
obj.nodeList = [obj.initialPose; obj.nodeList];
obj.nodeList = [obj.nodeList; obj.goalPose];




end























