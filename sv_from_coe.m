% /*! @sv_from_coe.m
% *************************************************************************
% <PRE>
% file.name       : sv_from_coe.m
% related files   :
% function&ablity : build system
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :six elements transform to position and veclcity
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/10/5     2.00                     add miu
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
function [r,v]=sv_from_coe(h,e,RA,incl,w,TA,miu)


rp=(h^2/miu)*(1/(1+e*cos(TA)))*(cos(TA)*[1;0;0]+sin(TA)*[0;1;0]);
vp=(miu/h)*(-sin(TA)*[1;0;0]+(e+cos(TA))*[0;1;0]);
R3_W=[cos(RA) sin(RA) 0;-sin(RA) cos(RA) 0;0 0 1];
R1_i=[1 0 0;0 cos(incl) sin(incl);0 -sin(incl) cos(incl)];
R3_w=[cos(w) sin(w) 0;-sin(w) cos(w) 0;0 0 1];
Q_pX=R3_W'*R1_i'*R3_w';
r=Q_pX*rp;
v=Q_pX*vp;

end



