% Time Invariant Domain of Attraction Estimate
%
% rho = ti_poly_doa(x,f,V,Lmonom)
% 
% Input:
%    x  -- an n-by-1 of free msspoly (the state variables)
%    f  -- an n-by-1 msspoly of differential equations (\dot x)
%    V  -- a positive definite candidate Lyapunov function.
%    Lmonom -- monomials for the LaGrange multiplier.
%
% Returns:
%    rho -- a conservative quantity such that { x | V(x) < rho}
%           is in the domain of attraction of the origin.
%
% Inspired by Parillo '05.
%
function rho = ti_poly_roa_sosgp(x0,f,S,K)


%     rho = ti_poly_roa(xbar,f(x,u0-K*xbar),);
    pvar x0 x1 x2 x3 x4 x5 u t;
    x = [x0 x1,x2,x3,x4,x5]';
    xbar = x;
    x = xbar+x0;
    u0 = zeros(3,1);
    xdot = f(x,u0-K*xbar);
    V = xbar'*S*xbar;
    if size(x,2) ~= 1, error('x must be a column'); end
    if size(xdot,2) ~= 1, error('f must be a column'); end
    if size(x,1) ~= size(xdot,1),
        error('x and xdot not matching sizes'); 
    end
    if size(V) ~= [1 1], error('V must be 1-by-1'); end
    
    Vdot = jacobian(V,x)*xdot;
    prog = sosprogram(x);
    [prog,s]=sossosvar(prog,x');
    prog = sosdecvar(prog,s);
%     s = monomials(x,0:7);
    prog = sosdecvar(prog,s);
    
    % Initialize the sum of squares program
    prog = sosprogram(x);
    % estimate g (the uperbound of lyap function)
    % Set options
    opts = gsosoptions;
    opts.minobj = -100000; % minmum value of obj
    opts.maxobj = 0;   % maxmum value of obj
    
    % Form bilinear constraints.  The variable t:=-gamma is used to convert
    % the maximization of gamma into a minimization of t.
    
    [prog,s]=sossosvar(prog,x');
    prog = sosdecvar(prog,s);
    % =============================================
    % Declare decision variable gam too
    sosc = polyconstr;
    sosc(1) = s>=0;
%     sosc(2) = -1 <= -(state'*state) + s*(V+t);
    sosc(2) = -s*Vdot <=  (x'*x)*(V+t);

    [info,dopt,sossol] = gsosopt(sosc,x,t,opts);
    s = subs(s,dopt);
    rho = -info.tbnds(2);
    

end