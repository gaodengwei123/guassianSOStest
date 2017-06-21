% /*! @addwaypoints.m
% *************************************************************************
% <PRE>
% file.name       : addwaypoints.m
% related files   :
% function&ablity : to guess the state in each vertexes
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/02/24     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output:

% *************************************************************************
function waypoint = addwaypoints(sys, vertex, N)
% N is the number of points between the way
speed = sys.INPUTS.speed;

nuGuess = vertex(1,:);
nuTime = 0;
dist = norm(vertex(2,:)-vertex(1,:));
dertD = dist/(N+1);
for i = 1:N+1
    deltaTime = dertD/speed;
    nuTime = [nuTime;nuTime(end)+deltaTime];
    h = atan2(vertex(2,2) -vertex(1 ,2) ,vertex(2 ,1) -vertex(1 ,1));
    nuGuess = [nuGuess; nuGuess(end,1)+cos(h)*dertD, nuGuess(end,2)+sin(h)*dertD, h, 0];
end
nuControl = zeros(length(nuTime),1);
waypoint = [nuTime, nuGuess, nuControl];
% reserve the points between initilal and target
waypoint(1,:) = [];
waypoint(end,:) = [];


end








