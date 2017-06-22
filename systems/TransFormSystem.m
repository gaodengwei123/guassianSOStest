% /*! @TransFormSystem.m
% *************************************************************************
% <PRE>
% file.name       : TransFormSystem.m
% related files   :
% function&ablity : TransForm System to another coordinate 
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/5/4      1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output:
% *************************************************************************
function g = TransFormSystem(xdot,V)
% from xdot = f(t,x,u) to zdot = g(t,z,u) via z=c(t,x) and x=d(t,z)
 % zdot = dcdt(t,d(t,z)) + dcdx(t,d(t,z)*f(t,d(t,z),u)

c = V.x - V.x0;
d = V.x + V.x0;
x = V.x;

g = diff(c,x)*xdot;
g = subss(g,x,d);

end

