% /*! @collision_check_point.m
% *************************************************************************
% <PRE>
% file.name       : collision_check_point.m
% related files   : ellipose toolbox 3.00
% function&ablity : check the ellipse convergence is collision or not
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/07/24     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output:

% *************************************************************************

function collision_found = collision_check_point(sys,current_ell0)
INPUTS = sys.INPUTS;
if sys.mark == 1
    pB1 = [1 0 0 0
        0 1 0 0]';
elseif sys.mark == 2
    pB1 = [1 0 0 0 0 0
        0 1 0 0 0 0]';
end
current_ell = projection(current_ell0,pB1);

% [q,~] = double(current_ell);
% px=q(1);
% py=q(2);
% returning 1 indicates collision risk
% returning 0 indicates no collision risk

collision_found = 0;

num_obstacles = size(INPUTS.obstacle,1);

for i = 1:num_obstacles
    obstacles =  ellipsoid(INPUTS.obstacle(i,:)',diag([INPUTS.obstacleRadius^2,INPUTS.obstacleRadius^2]));
    collision_found = intersect(current_ell,obstacles);
    %     if distance(current_ell, obstacles)<=1.5;
    
end
