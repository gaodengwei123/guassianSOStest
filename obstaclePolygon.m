% /*! @obstaclePolygon.m
% *************************************************************************
% <PRE>
% file.name       : obstaclePolygon.m
% related files   :
% function&ablity : make the obstacle to polygon 
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/8/1      1.00                     build
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
function obs_coords = obstaclePolygon(obs_pos,obs_lengths)
    obs_angles = 0:pi/6:2*pi;  % every 30deg will have a value to make obs
    obs_coords=[cos(obs_angles).*obs_lengths + obs_pos(1); sin(obs_angles).*obs_lengths + obs_pos(2)];
    obs_coords=obs_coords(:,convhull(obs_coords(1,:)',obs_coords(2,:)'));
end