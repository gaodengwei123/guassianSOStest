% /*! @TaylorApproxSys.m
% *************************************************************************
% <PRE>
% file.name       : TaylorApproxSys.m
% related files   :
% function&ablity : build system extend by taylor approximate
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/5/3      1.00        gaodengwei   build
% TaylorApprox system as ploynominal 
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
function xdot = TaylorApproxSys(sys,x0,u0,order)
obj.x0 = x0;
obj.u0 = u0;
obj.order = order;
obj.p_x = msspoly('x',sys.getNumStates);
obj.p_u = msspoly('u',sys.getNumInput);
obj.sys = sys;
p_xu = [obj.p_x; obj.p_u];
xdot = build_poly(@sys.dynamics,0,obj.x0,obj.u0,obj.order,p_xu);

end



function p=build_poly(fun,t,x0,u0,order,p_xu)
nX=length(x0);
nU=length(u0);
xu0=[x0;u0];
xu=TaylorVar.init(xu0,order);
x0=xu(1:nX);
u0=xu(nX+(1:nU));
p=getmsspoly(fun(t,x0,u0),p_xu-xu0);
end




