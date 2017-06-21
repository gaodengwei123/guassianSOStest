% calculate the trajectory in library
%
% ---xfrom------                 ---xto------
% ----xfrom------x0------------xT---xto------
% -----xfrom-----                 ---xto------
%
function [Time, point, Flag] = trajNode(sys,x0,xT,xFrom,xTo)
Flag = 1;
speed = sys.INPUTS.speed;
dist = norm(x0(1:2)-xT(1:2));
deltaTime = dist/speed;
Time = 2*[0;deltaTime];

if isempty(xFrom)
    h0 = sys.INPUTS.state_initial(3);      % initial 
    dtheta0 = sys.INPUTS.state_initial(4);
else
    h0 = atan2((xT(2) - x0(2)),(xT(1) - x0(1))); % x0 theta [-pi pi]
    dtheta0 = 0;
end

if isempty(xTo)
    hT = atan2((xT(2) - x0(2)),(xT(1) - x0(1)));
    dthetaT = 0;
else
    hT = atan2((xTo(2) - xT(2)),(xTo(1) - xT(1)));     % xT theta [-pi pi]
    dthetaT = 0;
    line1 = xT - x0;
    line2 = xTo - xT;
    angle = acos(dot(line1,line2)/norm(line1)/norm(line2));
    if angle > pi/2
        Flag = 0;
    end
end
% modify the angle which can not be over pi
if abs(hT-h0)>pi
    h0 = sign(hT)*2*pi+h0;
end

point = [x0(1), x0(2), h0, dtheta0
         xT(1), xT(2), hT, dthetaT];
end


