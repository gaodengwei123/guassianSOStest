function V = regionOfAttraction(sys,varargin)
% Estimates the region of attraction.

ok_sedumi = checkDependency('sedumi');
ok_mosek = checkDependency('mosek');
% checkDependency('yalmip');
if ~ok_sedumi && ~ok_mosek
    error('You need either MOSEK or SeDuMi installed to use this function.');
end

%% get Lyapunov candidate
num_x = sys.getNumStates;
num_u = sys.getNumInput;
t = msspoly('t');
x = msspoly('x',num_x);
u = msspoly('u',num_u);

% lyapunov equation
V = varargin{1};

% ploynominal or not
try
    f= sys.dynamics(t,x,u);
catch
    f = TaylorApproxSys(sys,sys.xG,sys.uG,3);
end
% default is origin(if not trans sys to V frame)
if ~all(V.x0==0)
    f = TransFormSystem(f,V);
end
%% zero all inputs
if (sys.getNumInput>0)
    f = subs(f,u,zeros(sys.getNumInput,1));
end

%% handle options
if (nargin>2) options=varargin{2};
else options=struct(); end
if (~isfield(options,'method'))
    options.method={'levelset'};
elseif (~iscell(options.method))
    options.method={options.method};
end
if (~isfield(options,'degV')) options.degV = 4; end
if (~isfield(options,'max_iterations')) options.max_iterations=10; end
if (~isfield(options,'converged_tol')) options.converged_tol=.01; end
if (~isfield(options,'optimize')) options.optimize=true; end
if (~isfield(options,'clean_tol')) options.clean_tol = 1e-6; end % tolerance for cleaning small terms
if (~isfield(options,'numSamples')) options.numSamples = 10^(num_x+1); end

if (~isfield(options,'degL1'))
    options.degL1 = options.degV-1 + deg(f,x);  % just a guess
end
if (~isfield(options,'degL2'))
    options.degL2 = options.degL1;
end

% Choose sdp solver
if ok_mosek
    options.solver = @spot_mosek;
else
    options.solver = @spot_sedumi;
end

for i=1:length(options.method)
    %% compute level set
    switch (lower(options.method{i}))
        case 'bilinear'
            V = bilinear(V,f,options);
        case 'levelset'
            V = levelSetMethod(V,f,options);
        case 'stochasticsearch'
            V = stochasticsearch(V,f,options);
%         case 'parsearch'
%             V = parsearch(V,f,options);
        otherwise
            error(['don''t know method: ', options.method]);
    end
end

end

function [T,Vbal,fbal,S,A] = balance(x,V,f,S,A)
if (nargin<4 || isempty(S))
    S=.5*doubleSafe(subs(diff(diff(V,x)',x),x,0*x));  % extract Hessian
end
if (nargin<5 || isempty(A))
    A = doubleSafe(subs(diff(f,x),x,0*x));
end

[T,D] = balanceQuadForm(S,(S*A+A'*S));

if (nargout>1)
    Vbal=subs(V,x,T*x);
    if (nargout>2)
        fbal=inv(T)*subs(f,x,T*x);
    end
end
end

%% for the bilinear search
function V = bilinear(V0,f,options)

x = V0.x;
num_x = length(x);
V=V0;

[T,V0bal,fbal,S0,A] = balance(x,V0.getPoly,f);
rho = 1;

L1monom = monomials(x,0:options.degL1);
L2monom = monomials(x,0:options.degL2);
Vmonom = monomials(x,0:options.degV);

vol=0;
Vpoly = V.getPoly;
for iter=1:options.max_iterations
    last_vol = vol;
    
    % balance on every iteration (since V and Vdot are changing):
    [T,Vbal,fbal]=balance(x,Vpoly,f,S0/rho,A);
    V0bal=subs(V0.getPoly,x,T*x);
    
    [L1,sigma1] = findL1(x,fbal,Vbal,L1monom,options);
    L2 = findL2(x,Vbal,V0bal,rho,L2monom,options);
    [Vbal,rho] = optimizeV(x,fbal,L1,L2,V0bal,sigma1,Vmonom,options);
    vol = rho;
    Vpoly = subs(Vbal,x,inv(T)*x);
    
    V = V0.getPoly/rho;
    
    % check for convergence
    if ((vol - last_vol) < options.converged_tol*last_vol)
        break;
    end
end
xfun0 = getLevelSet(x,Vpoly,options);
plot(xfun0(1,:),xfun0(2,:),'b','lineWidth',2)
end

function [L1,sigma1] = findL1(x,f,V,Lxmonom,options)
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% construct multipliers for Vdot
[prog,L1] = prog.newFreePoly(Lxmonom);

% construct Vdot
Vdot = clean(diff(V,x)*f);

% construct slack var
[prog,sigma1] = prog.newPos(1);

% setup SOS constraints
prog = prog.withSOS(-Vdot + L1*(V - 1) - sigma1*V);
prog = prog.withSOS(L1);

% run SeDuMi/MOSEK and check output
solver = options.solver;
options = spot_sdp_default_options();
sol = prog.minimize(-sigma1,solver,options);

if sol.status == spotsolstatus.STATUS_SOLVER_ERROR
    error('The solver threw an internal error.');
end
if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end

L1 = sol.eval(L1);
sigma1 = sol.eval(sigma1);
end

function L2 = findL2(x,V,V0,rho,Lxmonom,options)
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% construct multipliers
[prog,L2] = prog.newFreePoly(Lxmonom);

[prog,slack] = prog.newPos(1);

prog = prog.withSOS(-(V-1) + L2*(V0-rho));
prog = prog.withSOS(L2);

solver = options.solver;
options = spot_sdp_default_options();
% options.verbose = 1;
sol = prog.minimize(slack,solver,options);% keyboard;

if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end

L2 = sol.eval(L2);
end

function [V,rho]=optimizeV(x,f,L1,L2,V0,sigma1,Vxmonom,options)
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% construct V
[prog,V] = prog.newFreePoly(Vxmonom);
Vdot = diff(V,x)*f;

% construct rho
[prog,rho] = prog.newPos(1);

% setup SOS constraints
prog = prog.withSOS(-Vdot + L1*(V - 1) - sigma1*V/2);
prog = prog.withSOS(-(V-1) + L2*(V0 - rho));
prog = prog.withSOS(V);

% run SeDuMi/MOSEK and check output
solver = options.solver;
options = spot_sdp_default_options();
options.verbose=1;
sol = prog.minimize(-rho,solver,options);

if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded.');
end

V = sol.eval(V);
rho = double(sol.eval(rho));
end


%% Pablo's method (jointly convex in rho and lagrange multipliers)
function V = levelSetMethod(V0,f,options)

x = V0.x;
[T,V,f] = balance(x,V0.getPoly,f);

%% compute Vdot
Vdot = diff(V,x)*f;

% check Hessian Vdot at origin, to make sure it's negative def.
H=.5*doubleSafe(subs(diff(diff(Vdot,x)',x),x,0*x));  % extract Hessian
if (~isPositiveDefinite(-H)) error('Vdot must be negative definite at the origin'); end

prog = spotsosprog;
prog = prog.withIndeterminate(x);
Lmonom = monomials(x,0:options.degL1);
%  Lmonom = hermite_basis(monomials(x,0:options.degL1));

[prog,rho] = prog.newFree(1);

[prog,L] = prog.newFreePoly(Lmonom);

prog = prog.withSOS((x'*x)^floor((options.degL1 + deg(Vdot)-deg(V))/2)*(V - rho) +  L*Vdot);

solver = options.solver;
options = spot_sdp_default_options();
options.verbose = 1;
sol = prog.minimize(-rho,solver,options);

if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end

rho = doubleSafe(sol.eval(rho));
if (rho<=0) error('optimization failed'); end

V = V/rho;

%% undo balancing
V = subs(V,x,inv(T)*x);

%   V = V_function(t,x,S,[]);


end

%% Stochastic search
function V = stochasticsearch(V0,f,options)
% w = msspoly('w',1);
x = V0.x;
V=V0;
% f = f+[0;w];
w=[];
rho = 1;

L1monom = monomials([x;w],0:2);
Lw0 = monomials([x;w],0:2);

vol=0;
Vpoly = V0.getPoly;
for iter=1:options.max_iterations
    last_vol = vol;
    
    Vbal = Vpoly;
    fbal = f;
    [L1,L2] = findLw(x,w,fbal,Vbal,rho,L1monom,Lw0,options);
%     uf = findLrhow(x,w,fbal,Vbal,L1,L2,Lw0,options);
    uf = 0;
    [Vbal,rho] = optimizeVw(x,w,fbal,L1,L2,Lw0,uf,V0,options);
    vol = rho;
%     Vpoly = subs(Vbal,x,inv(T)*x);
        Vpoly = Vbal;
        V = V0.getPoly/rho;
    % check for convergence
    if ((vol - last_vol) < options.converged_tol*last_vol)
        break;
    end
end
xfun0 = getLevelSet(x,Vpoly,options);
plot(xfun0(1,:),xfun0(2,:),'b','lineWidth',2)
end

function [L1,L2] = findLw(x,w,f,V,rho,Lxmonom,Lwmonom,options)
prog = spotsosprog;
prog = prog.withIndeterminate([x;w]);

% construct multipliers
[prog,L1] = prog.newFreePoly(Lxmonom);
% [prog,Lw] = prog.newFreePoly(Lwmonom);
% [prog,Lw2] = prog.newFreePoly(Lwmonom);
prog = prog.withSOS(L1);
% prog = prog.withSOS(Lw);
% prog = prog.withSOS(Lw2);
% construct u
% um = monomials(x,1:1);
% [prog,uf] = prog.newFreePoly(um);
% f= f+[0;uf];
% construct Vdot
V = clean(V,options.clean_tol);
Vdot = clean(diff(V,x)*f,options.clean_tol);
% Vdot0 = subs(Vdot,w,0*w);
% [prog,gamma] = prog.newPos(1);
[prog,gamma0] = prog.newPos(1);

% Vdot0 = subs(Vdot,w,0*w);
% Lxmonom = subs(Lxmonom,w,0*w);
% [prog,L10] = prog.newFreePoly(Lxmonom);
[prog,L2] = prog.newFreePoly(Lxmonom);


% prog = prog.withSOS(L10);
prog = prog.withSOS(L2);

prog = prog.withSOS(-gamma0*(x'*x)^2- Vdot + L1*(V-rho) +L2*(x'*x-5));
% prog = prog.withSOS(- Vdot );+Lw*(w-0.001)+Lw2*(-0.001-w )




% run SeDuMi/MOSEK and check output   L1*(V-rho) + 
solver = options.solver;
pars = spot_sdp_default_options();

[prog,dummy] = prog.newPos(1);
sol = prog.minimize(dummy,solver,pars);

if sol.status == spotsolstatus.STATUS_SOLVER_ERROR
    error('The solver threw an internal error.');
end
if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end

L1 = sol.eval(L1);
L2 = sol.eval(L2);
% Lw = sol.eval(Lw);

end

function uf= findLrhow(x,w,f,V,L1,L2,Lwmonom,options)
% min rho
% fix L,V,search for rho,Lw

prog = spotsosprog;
prog = prog.withIndeterminate([x;w]);

% construct u
um = monomials(x,1:1);
[prog,uf] = prog.newFreePoly(um);
f= f+[0;uf];

Vdot = clean(diff(V,x)*f);

[prog,Lw] = prog.newFreePoly(Lwmonom);
prog = prog.withSOS(Lw);

[prog,rho] = prog.newPos(1);
% [prog,gamma] = prog.newPos(1);

prog = prog.withSOS(- Vdot +L1*(V-rho) +L2*(x'*x-5)+Lw*(w'*w-1));

solver = options.solver;
options = spot_sdp_default_options();

sol = prog.minimize(-rho,solver,options);

if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end

uf = sol.eval(uf);
end

function [V,rho]=optimizeVw(x,w,f,L,L2,Lwmonom,uf,V0,options)

% fix L
% search for: V rho

prog = spotsosprog;
prog = prog.withIndeterminate([x;w]);

% Create Phi and new V
[prog,Phi] = prog.newPSD(length(x));
V = x'*Phi*x;
% [prog,V] = prog.newFreePoly(Vxmonom);
f = f+[0;uf];
Vdot = diff(V,x)*f;

V = clean(V,options.clean_tol);
Vdot = clean(Vdot,options.clean_tol);

% construct uncertainty
[prog,Lw] = prog.newFreePoly(Lwmonom);
[prog,Lw2] = prog.newFreePoly(Lwmonom);
prog = prog.withSOS(Lw);
prog = prog.withSOS(Lw2);
% construct rho
[prog,rho] = prog.newPos(1);

% [prog,gamma] = prog.newPos(1);

% setup SOS constraints
prog = prog.withEqs(trace(Phi) - trace(V0.Spp));
prog = prog.withSOS( - Vdot + L*(V-rho)+L2*(x'*x-5));
prog = prog.withSOS(V);

% run SeDuMi/MOSEK and check output
solver = options.solver;
options = spot_sdp_default_options();
sol = prog.minimize(-rho,solver,options);

if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end

V = sol.eval(V);
rho = double(sol.eval(rho));
end

function V = maxROAfeedback()
% First step: Fix rho and V, search for L and u
[L1,uf] = findLU(V,rho,x,f,u,options);
plot(iter,rho,'ro');
title(['iteration ' num2str(iter)])
drawnow;
hold on

% Second step: Fix V and L, search for u and rho
[uf,rho] = findURho(V,x,u,f,L1,options);
plot(iter+0.33,rho,'ro');
title(['iteration ' num2str(iter)])
drawnow;

% Third step: Fix u and L, search for V and rho
[V,rho,Phi] = findVRho(x,u,uf,f,L1,S0,options);
plot(iter+0.66,rho,'ro');
title(['iteration ' num2str(iter)])
drawnow;




end






function y=doubleSafe(x)
y=double(x);
if (~isa(y,'double')) error('double failed'); end
end
