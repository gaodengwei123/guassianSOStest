% /*! @rechability_set.m
% *************************************************************************
% <PRE>
% file.name       : rechability_set.m
% related files   :
% function&ablity :
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/7/7      1.00                     build
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
function [setout,Xout ]= rechability_set(A,B,U,T,x0,n,xest,x10,P0)
global INPUTS
X0 = ellipsoid(x0,P0);
X10 = ellipsoid(x10,P0);
% PBB = [1 0 0 0; 0 1 0 0]';
% PBB1 = [0 0 1 0;0 0 0 1]';
% Xps0 = projection(X0,PBB);
% Xps1 = projection(X10,PBB);
% Xps01 = projection(X0,PBB1);
% Xps11 = projection(X10,PBB1);
% RA0 = minkdiff(Xps0, Xps1);
% RA1 = minkdiff(Xps01, Xps11);
L0 = x0;
% isbaddirection(Xps0,Xps1,RA0)
G = eye(n);
lsys = linsys(A,B,U,G,INPUTS.noise_w);

rs = reach(lsys,X0,L0,T);
BB = [1 0 0 0; 0 1 0 0]';
ps = projection(rs,BB);
endrs = cut(rs,T(2));
setout.center = get_center(rs);
% setout.shape = 
Xout = ellipsoid(xest+setout.center,setout.shape);

end
