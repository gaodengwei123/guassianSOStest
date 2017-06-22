% /*! @Invariant_Funnels.m
% *************************************************************************
% <PRE>
% file.name       : Invariant_Funnels.m
% related files   :
% function&ablity : Calculate Invariant_Funnels(reference:Invariant Funnels around Trajectories using Sum-of-Squares Programming)
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% Linearize polynomial dynamics about a nominal trajectory.
%
% [A,B] = tv_poly_linearize(xdot,x0,u0)
%
% Input:
%   xdot  -- xdot(t,x,u) takes (1-by-1)x(n-by-1)x(m-by-1) msspolys
%                returns n-by-1 msspoly.
%   x0    -- x0(t) function n-by-1 double to linearize around.
%   u0    -- u0(t) function m-by-1 double to linearize around.
% Output:
%   A     -- A(t) n-by-n Jacobian of xdot w.r.t. x at x0(t).
%   B     -- B(t) n-by-m Jacobian of xdot w.r.t. u at u0(t).
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/7/24     1.00        dengwei      change the msspoly2sym.m from drake, since msspoly is from spot-master package
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
function sys = tv_poly_linearize(sys,varargin)
x0 = varargin{1};
t = sys.p_t;
x = sys.p_x;

if size(x0(0),2) ~= 1, error('x0 must be a column.'); end

if varargin < 3 % no control input
    if varargout > 1
        error(['u0 must be provided to output B ' ...
            'matrix.']); 
    end
    xdot = sys.f(t,x);
%     A = @(t)double(subs(dfdx,[tt;xx],[t;x0(t)]));
    A = @(t)double(subs(diff(xdot,x),x,x0(t)));
else
    u0 = varargin{2};
    if size(u0(0),2) ~= 1, error('u0 must be a column.'); end
    u = sys.p_u;

    xdot = sys.f(t,x,u);
    A = @(t)double(subs(diff(xdot,x),x,x0(t)));
    B = @(t)double(subs(diff(xdot,u),u,u0(t)));

    sys.A = A;
    sys.B = B;
end


%    x = msspoly('x',n);
%    u = msspoly('u',m);
%    xdot = msspoly(sys.f(0,x,u));

%    AB = double(subs(diff(xdot,[x;u]),[x;u],[x0;u0]));
%    A = AB(1:n,1:n);
%    B = AB(1:n,n+(1:m));


