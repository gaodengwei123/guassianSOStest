function [c,V] = maxROAFeedbackcontrol(sys,K,V,options)

% The approach searches for a polynomial controller and a Lyapunov function
% and attempts to maximize the size of the ROA. 


ok_sedumi = checkDependency('sedumi');
ok_mosek = checkDependency('mosek');
if ~ok_sedumi && ~ok_mosek
    error('You need either MOSEK or SeDuMi installed to use this function.');
end
typecheck(V,'V_function'); 

% Default options
if (nargin<4) options = struct(); end
if (~isfield(options,'controller_deg')) options.controller_deg = 1; end % Degree of polynomial controller to search for
if (~isfield(options,'max_iterations')) options.max_iterations = 10; end % Maximum number of iterations (3 steps per iteration)
if (~isfield(options,'converged_tol')) options.converged_tol = 1e-3; end % Tolerance for checking convergence
if (~isfield(options,'rho0')) options.rho0 = 0.01; end % Initial "guessed" rho 
if (~isfield(options,'clean_tol')) options.clean_tol = 1e-6; end % tolerance for cleaning small terms
if (~isfield(options,'backoff_percent')) options.backoff_percent = 5; end % 5 percent backing off
if (~isfield(options,'degL1')) options.degL1 = options.controller_deg + 1; end % Chosen to do degree matching
if (~isfield(options,'degLu')) options.degLu = options.controller_deg - 1; end 
if (~isfield(options,'degLu1')) options.degLu1 = 2; end 
if (~isfield(options,'degLu2')) options.degLu2 = 2; end
if (~isfield(options,'degLup')) options.degLup = 2; end
if (~isfield(options,'degLum')) options.degLum = 2; end

%% get Lyapunov candidate
num_x = sys.getNumStates;
num_u = sys.getNumInput;
t = msspoly('t');
u = msspoly('u',num_u);
x = msspoly('x',num_x);

% test the system must be polynominal
try
    f= sys.dynamics(t,x,u);
catch
    f = TaylorApproxSys(sys,sys.xG,sys.uG,3);
end
% default is origin
if ~isempty(V.x0)||~all(V.x0==0)
    f = TransFormSystem(f,V);
end


%% Choose sdp solver (perfer mosek)
if ok_mosek
    options.solver = @spot_mosek;
else
    options.solver = @spot_sedumi;
end

% Get V 
S0 = V.Spp;
x0 = V.x0;
V = V.getPoly;

% Compute fmax and fmin if saturations
if ~all(isinf([sys.umax;sys.umin]))
  if getNumInputs(sys)>1
    error('saturations only implemented for a single input (so far)');
  end
  if any(isinf([sys.umax;sys.umin]))
    error('input saturations must be all inf or not inf (for now)');
  end
  f_umax = subs(f,u,sys.umax);
  f_umin = subs(f,u,sys.umin);
  saturations = true;
else
  saturations = false;
end


% Initialize u, Phi and rho
uf = -K*x;   % control input
rho = options.rho0;

rho_last = rho;
for iter=1:options.max_iterations
  
  if saturations
       % First step: Fix rho and V, search for L 
       [L1,Lu1,Lu2,Lp,Lup,Lm,Lum] = findLwithSats(V,rho,x,u,f,f_umax,f_umin,sys.umax,sys.umin,uf,options);
       plot(iter,rho,'ro');
       title(['iteration ' num2str(iter)])
       drawnow;
       hold on

      % Second step: Fix V and L, search for u and rho
      [uf,rho] = findURhowithSats(V,x,u,f,f_umax,f_umin,sys.umax,sys.umin,L1,Lu1,Lu2,Lp,Lup,Lm,Lum,options);
      plot(iter+0.33,rho,'ro');
      title(['iteration ' num2str(iter)])
      drawnow;

      % Third step: Fix u and L, search for V and rho
      [V,rho,Phi] = findVRhowithSats(x,u,uf,f,f_umax,f_umin,sys.umax,sys.umin,L1,Lp,Lm,S0,options);
      plot(iter+0.66,rho,'ro');
      title(['iteration ' num2str(iter)])
      drawnow;
   else
       % First step: Fix rho and V, search for L and u
       [L1,uf] = findLU(V,rho,x,f,u,options);
%        plot(iter,rho,'ro');
%        title(['iteration ' num2str(iter)])
%        drawnow;
%        hold on

      % Second step: Fix V and L, search for u and rho
      [uf,rho] = findURho(V,x,u,f,L1,options);
%       plot(iter+0.33,rho,'ro');
%       title(['iteration ' num2str(iter)])
%       drawnow;

      % Third step: Fix u and L, search for V and rho
      [V,rho,Phi] = findVRho(x,u,uf,f,L1,S0,options);
%       plot(iter+0.66,rho,'ro');
%       title(['iteration ' num2str(iter)])
%       drawnow;
   end

  % check for convergence
  if (double(rho - rho_last) < options.converged_tol*rho_last)  % see if it's converged
   break;
  end
  
  rho_last = double(rho);
  
end

% Optimized Lyapunov function
V = V_function(t,x,double(Phi),[]);
% c = double(diff(uf,x));
c = uf;
end

function [L1f,Lu1f,Lu2f,Lpf,Lupf,Lmf,Lumf] = findLwithSats(V,rho,x,u,f,f_umax,f_umin,umax,umin,uopt,options)

    disp('Step 1: Searching for multipliers...')
        
    % SOS program 
    prog = mssprog;
    
    % Compute Vdot  
    Vdot = diff(V,x)*subs(f,u,uopt);
    Vdotp = diff(V,x)*f_umax;
    Vdotm = diff(V,x)*f_umin; 
       
    % Clean stuff
    V = clean(V,options.clean_tol);
    Vdot = clean(Vdot,options.clean_tol);
    Vdotp = clean(Vdotp,options.clean_tol);
    Vdotm = clean(Vdotm,options.clean_tol);
    
    % Declare multipliers
    L1m = monomials(x,options.degL1:options.degL1);
    [prog,l1] = new(prog,length(L1m),'free');
    L1 = l1'*L1m;
    prog.sos = L1;
       
    if options.degLu1 == 0 % This is a bit silly to have to special case
        [prog,Ls] = new(prog,4,'pos');
        Lu1 = Ls(1);
        Lu2 = Ls(2);
        Lup = Ls(3);
        Lum = Ls(4);
    else
        Lu1m = monomials(x,0:options.degLu1);
        [prog,lu1] = new(prog,length(Lu1m),'free');
        Lu1 = lu1'*Lu1m;

        Lu2m = monomials(x,0:options.degLu2);
        [prog,lu2] = new(prog,length(Lu2m),'free');
        Lu2 = lu2'*Lu2m;


        Lupm = monomials(x,0:options.degLup);
        [prog,lup] = new(prog,length(Lupm),'free');
        Lup = lup'*Lupm;

        Lumm = monomials(x,0:options.degLum);
        [prog,lum] = new(prog,length(Lumm),'free');
        Lum = lum'*Lumm;
        
        prog.sos = Lu1;
        prog.sos = Lu2;
        prog.sos = Lup;
        prog.sos = Lum; 
    end
    
    Lpm = monomials(x,options.degL1:options.degL1);
    [prog,lp] = new(prog,length(Lpm),'free');
    Lp = lp'*Lpm;
    prog.sos = Lp;
    
    Lmm = monomials(x,options.degL1:options.degL1);
    [prog,lm] = new(prog,length(Lmm),'free');
    Lm = lm'*Lmm;
    prog.sos = Lm;
    
    % Create gammas
    [prog,gamma] = new(prog,1,'pos');
    [prog,gammap] = new(prog,1,'pos');
    [prog,gammam] = new(prog,1,'pos'); 
    
    % Declare SOS conditions
    prog.sos =  -gamma*(x'*x)^(deg(Vdot)/2) - Vdot + L1*(V-rho) + Lu1*(uopt - umax) + Lu2*(umin - uopt);
    prog.sos =  -gammap*(x'*x)^(deg(Vdotp)/2) - Vdotp + Lp*(V-rho) + Lup*(umax - uopt);
    prog.sos =  -gammam*(x'*x)^(deg(Vdotm)/2) - Vdotm + Lm*(V-rho) + Lum*(uopt - umin);
              
   % keyboard;
    % Solve SOS program  
    pars = struct();
    pars.fid = 0;
    [prog,dummy] = new(prog,1,'pos');
    % keyboard;
    [prog,info] = sedumi(prog,dummy,1,pars,0) 
    
    % Optimized multipliers
    L1f = prog(L1);
    Lu1f = prog(Lu1);
    Lu2f = prog(Lu2);
    Lpf = prog(Lp);
    Lupf = prog(Lup);
    Lmf = prog(Lm);
    Lumf = prog(Lum);

end

function [L1f,uf] = findLU(V,rho,x,f,u,options)

    disp('Step 1: Searching for multipliers and controller...')
    
    % SOS program 
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    % Create u
    um = monomials(x,1:options.controller_deg);
    
    uf = [];
    for j = 1:length(u)
        [prog,un] = prog.newFreePoly(um);
        uf = [uf;un];
    end

    % Compute Vdot  
    Vdot = diff(V,x)*subs(f,u,uf);
       
    % Clean stuff
    V = clean(V,options.clean_tol);
    Vdot = clean(Vdot,options.clean_tol);
    
    % Declare multipliers
    L1m = monomials(x,options.degL1:options.degL1); % Make it homogeneous
    [prog,L1] = prog.newFreePoly(L1m);
    prog = prog.withSOS(L1);
    
    % Create gammas
    [prog,gamma] = prog.newPos(1);

    prog = prog.withSOS(-gamma*(x'*x)^(deg(Vdot,x)/2) - Vdot + L1*(V-rho));

    % Solve SOS program  
    solver = options.solver;
    pars = spot_sdp_default_options();
    pars.verbose = 0;
    sol = prog.minimize(0,solver,pars);  %¡¡no optimaize obj
    
    if sol.status == spotsolstatus.STATUS_SOLVER_ERROR
        error('The solver threw an internal error.');
    end
    if ~sol.isPrimalFeasible
        error('Problem looks primal infeasible');
    end
    if ~sol.isDualFeasible
        error('Problem looks dual infeasible. It is probably unbounded. ');
    end
    
    % Optimized multipliers
    L1f = sol.eval(L1);
% % SOS program 
%     prog = mssprog;
% %     w = msspoly('w',2);
%     % Create u
%     um = monomials(x,1:options.controller_deg);
%     lu = [];
%     for j = 1:length(u)
%         [prog,luj] = new(prog,length(um),'free');
%         lu = [lu;luj'];
%     end
%     uf = lu*um;
%     
%     % Compute Vdot  
% %     f = f+[0;0;w];
%     Vdot = diff(V,x)*subs(f,u,uf);
%        
%     % Clean stuff
%     V = clean(V,options.clean_tol);
%     Vdot = clean(Vdot,options.clean_tol);
%     
%     % Declare multipliers
%     L1m = monomials(x,options.degL1:options.degL1); % Make it homogeneous
%     [prog,l1] = new(prog,length(L1m),'free');
%     L1 = l1'*L1m;
%     prog.sos = L1;
%     
%     % Create gammas
%     [prog,gamma] = new(prog,1,'pos');
%     
%     % Declare SOS conditions
%     
% %     Lw0 = monomials([x;w],0:2); % Make it homogeneous  +Lw*(w'*w-0.01)
% %     [prog,lw1] = new(prog,length(Lw0),'free');
% %     Lw = lw1'*Lw0;
% %     prog.sos = Lw;
% %     [prog,lw2] = new(prog,length(Lw0),'free');
% %     Lw2 = lw2'*Lw0;
% %     prog.sos = Lw2;
%     
% %     +Lw*(w-0.0001)+Lw2*(-0.0001-w)
%     prog.sos =  -gamma*(x'*x)^(deg(Vdot,x)/2) - Vdot + L1*(V-rho);
% 
%     % Solve SOS program  
%     pars = struct();
%     pars.fid = 0;
%     [prog,dummy] = new(prog,1,'free');
%     [prog,info] = sedumi(prog,dummy,1,pars,0) 
%     
%     % Optimized multipliers
%     L1f = prog(L1);
end


function [uf,rho] = findURhowithSats(V,x,u,f,f_umax,f_umin,umax,umin,L1,Lu1,Lu2,Lp,Lup,Lm,Lum,options)

    disp('Step 2: Searching for controller and rho...')

    % SOS program 
    prog = mssprog;
    
    % Variables we need later
    [prog,bslack] = new(prog,1,'pos');
    [prog,dummy] = new(prog,1,'pos');

    % Declare rho
    [prog,rho] = new(prog,1,'pos');

    % Create u
    um = monomials(x,1:options.controller_deg);
    [prog,lu] = new(prog,length(um),'free');
    uf = lu'*um;
    
    % Compute Vdot  
    Vdot = diff(V,x)*subs(f,u,uf);
    Vdotp = diff(V,x)*f_umax;
    Vdotm = diff(V,x)*f_umin;
    
    % Clean stuff
    V = clean(V,options.clean_tol);
    Vdot = clean(Vdot,options.clean_tol);
    Vdotp = clean(Vdotp,options.clean_tol);
    Vdotm = clean(Vdotm,options.clean_tol);
    
    % Declare SOS conditions
    prog.sos = -Vdot + L1*(V-rho) + Lu1*(uf - umax) + Lu2*(umin - uf);
    prog.sos = -Vdotp + Lp*(V-rho) + Lup*(umax - uf);
    prog.sos = -Vdotm + Lm*(V-rho) + Lum*(uf - umin);
             
    % Solve SOS program
    pars.fid = 0;
    [prog,info] = sedumi(prog,-rho,1,pars,1)   
    
    if info.numerr == 2 || info.pinf == 1 || info.dinf == 1
      error('No solution found');
    end
        
    disp('Backing off now...')
        
    % Back off on objective
    rhof = double(prog(rho));
       
    prog.eq = bslack - (-(1-options.backoff_percent/100)*rhof + rho);
   
    [prog,info] = sedumi(prog,dummy,1,pars,0) 
    
    if info.numerr == 2 || info.pinf == 1 || info.dinf == 1
        keyboard;
    end
    
    rho = double(prog(rho));
    uf = prog(uf);
%     Lu1f = prog(Lu1);
%     Lu2f = prog(Lu2);
%     ufcoeffs = double(prog(lu));
        
end

function [uf,rho] = findURho(V,x,u,f,L1,options)

    disp('Step 2: Searching for controller and rho...')

    % SOS program 
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    
    % Some variables needed later
    [prog,bslack] = prog.newPos(1);

    % Declare rho
    [prog,rho] = prog.newPos(1);

    % Create u
    um = monomials(x,1:options.controller_deg);
    uf = [];
    for j = 1:length(u)
        [prog,un] = prog.newFreePoly(um);
        uf = [uf;un];
    end
    
    % Compute Vdot  
    Vdot = diff(V,x)*subs(f,u,uf);
    
    % Clean stuff
    V = clean(V,options.clean_tol);
    Vdot = clean(Vdot,options.clean_tol);
    
    % Declare SOS conditions
    prog = prog.withSOS(-Vdot + L1*(V-rho));
             
    % Solve SOS program

    solver = options.solver;
    pars = spot_sdp_default_options();
    pars.verbose = 1; % obj function
    
    sol = prog.minimize(-rho,solver,pars); 

    disp('Backing off now...')
        
    % Back off on objective
    rhof = sol.eval(rho);
        
    prog = prog.withEqs(bslack - (-(1-options.backoff_percent/100)*rhof + rho));
   
    pars.verbose = 0;
    sol = prog.minimize(0,solver,pars);
    
    if sol.status == spotsolstatus.STATUS_SOLVER_ERROR
        error('The solver threw an internal error.');
    end
    if ~sol.isPrimalFeasible
        error('Problem looks primal infeasible');
    end
    
    if ~sol.isDualFeasible
        error('Problem looks dual infeasible. It is probably unbounded. ');
    end

    rho = sol.eval(rho);
    uf = sol.eval(uf);
% % SOS program 
%     prog = mssprog;
%     
%     % Some variables needed later
%     [prog,bslack] = new(prog,1,'pos');
%     [prog,dummy] = new(prog,1,'pos');
% 
%     % Declare rho
%     [prog,rho] = new(prog,1,'pos');
% 
%     % Create u
%     um = monomials(x,1:options.controller_deg);
%     lu = [];
%     for j = 1:length(u)
%         [prog,luj] = new(prog,length(um),'free');
%         lu = [lu;luj'];
%     end
%     uf = lu*um;
%     
%     % Compute Vdot  
%     Vdot = diff(V,x)*subs(f,u,uf);
%     
%     % Clean stuff
%     V = clean(V,options.clean_tol);
%     Vdot = clean(Vdot,options.clean_tol);
%     
%     % Declare SOS conditions
%     prog.sos = -Vdot + L1*(V-rho);
%              
%     % Solve SOS program
%     pars.fid = 0;
%     [prog,info] = sedumi(prog,-rho,1,pars,1)   
%     
%     if info.numerr == 2 || info.pinf == 1 || info.dinf == 1
%         keyboard;
%     end
%         
%     disp('Backing off now...')
%         
%     % Back off on objective
%     rhof = double(prog(rho));
%         
%     prog.eq = bslack - (-(1-options.backoff_percent/100)*rhof + rho);
%    
%     [prog,info] = sedumi(prog,dummy,1,pars,0) 
%     
%     if info.numerr == 2 || info.pinf == 1 || info.dinf == 1
%         keyboard;
%     end
%     
%     rho = double(prog(rho));
%     uf = prog(uf);
end


function [V,rho,Phi] = findVRhowithSats(x,u,uf,f,f_umax,f_umin,umax,umin,L1,Lp,Lm,S0,options)

    disp('Step 3: Searching for V and rho...')

    % SOS program 
    prog = mssprog;
    
    % Some variables needed later
    [prog,bslack] = new(prog,1,'pos');
    [prog,dummy] = new(prog,1,'pos');
    
    % Create rho
    [prog,rho] = new(prog,1,'pos');
    
    % Create Phi and new V
    [prog,Phi] = new(prog,length(x),'psd');
    V = x'*Phi*x;
    
    % Normalization constraint
    prog.eq = trace(Phi) - trace(S0);
    
    % Compute Vdot  
    Vdot = diff(V,x)*subs(f,u,uf);
    Vdotp = diff(V,x)*f_umax;
    Vdotm = diff(V,x)*f_umin;
       
    % Clean stuff
    V = clean(V,options.clean_tol);
    Vdot = clean(Vdot,options.clean_tol);
    Vdotp = clean(Vdotp,options.clean_tol);
    Vdotm = clean(Vdotm,options.clean_tol);
    L1 = clean(L1,options.clean_tol);
    Lp = clean(Lp,options.clean_tol);
    Lm = clean(Lm,options.clean_tol);
    
    % Declare multipliers
    if options.degLu1 == 0 % This is a bit silly to have to special case
        [prog,Ls] = new(prog,4,'pos');
        Lu1 = Ls(1);
        Lu2 = Ls(2);
        Lup = Ls(3);
        Lum = Ls(4);
    else
        Lu1m = monomials(x,0:options.degLu1);
        [prog,lu1] = new(prog,length(Lu1m),'free');
        Lu1 = lu1'*Lu1m;

        Lu2m = monomials(x,0:options.degLu2);
        [prog,lu2] = new(prog,length(Lu2m),'free');
        Lu2 = lu2'*Lu2m;

        Lupm = monomials(x,0:options.degLup);
        [prog,lup] = new(prog,length(Lupm),'free');
        Lup = lup'*Lupm;

        Lumm = monomials(x,0:options.degLum);
        [prog,lum] = new(prog,length(Lumm),'free');
        Lum = lum'*Lumm;
        
        prog.sos = Lu1;
        prog.sos = Lu2;
        prog.sos = Lup;
        prog.sos = Lum;
    end
    
    % Declare SOS conditions
    prog.sos = -Vdot + L1*(V-rho) + Lu1*(uf - umax) + Lu2*(umin - uf);
    prog.sos = -Vdotp + Lp*(V-rho) + Lup*(umax - uf);
    prog.sos = -Vdotm + Lm*(V-rho) + Lum*(uf - umin);
            
    % Solve SOS program  
    pars = struct();
    pars.fid = 0;
    [prog,info] = sedumi(prog,-rho,1,pars,1) 
            
    disp('Backing off now...')
       
    % Back off on objective
    rhof = double(prog(rho));
        
    prog.eq = bslack - (-(1 - options.backoff_percent/100)*rhof + rho);
   
    [prog,info] = sedumi(prog,dummy,1,pars,0) 
    if info.numerr == 2 || info.pinf == 1 || info.dinf == 1
        keyboard;
    end
     
    rho = double(prog(rho));
    V = x'*double(prog(Phi))*x;
    Phi = double(prog(Phi));

end

function [V,rho,Phi] = findVRho(x,u,uf,f,L1,S0,options)

    disp('Step 3: Searching for V and rho...')

    % SOS program 
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    
    % Some variables needed later
    [prog,bslack] =  prog.newPos(1);
    
    % Create rho
    [prog,rho] = prog.newPos(1);
    
    % Create Phi and new V
    [prog,Phi] = prog.newPSD(length(x));
    V = x'*Phi*x;
    
    % Normalization constraint
     prog = prog.withEqs(trace(Phi) - trace(S0));
    
    % Compute Vdot  
    Vdot = diff(V,x)*subs(f,u,uf);
       
    % Clean stuff
    V = clean(V,options.clean_tol);
    Vdot = clean(Vdot,options.clean_tol);
    
    % Declare SOS conditions
    prog = prog.withSOS(-Vdot + L1*(V-rho));
            
    % Solve SOS program  
    solver = options.solver;
    pars = spot_sdp_default_options();
    pars.verbose = 1;
    sol = prog.minimize(-rho,solver,pars);
    
    disp('Backing off now...')
       
    % Back off on objective
    rhof = sol.eval(rho);
    
    prog = prog.withEqs(bslack - (-(1 - options.backoff_percent/100)*rhof + rho));
    pars.verbose = 0;
    sol = prog.minimize(0,solver,pars);
    if sol.status == spotsolstatus.STATUS_SOLVER_ERROR
        error('The solver threw an internal error.');
    end
    if ~sol.isPrimalFeasible
        error('Problem looks primal infeasible');
    end
    
    if ~sol.isDualFeasible
        error('Problem looks dual infeasible. It is probably unbounded. ');
    end
    rho = sol.eval(rho);
    Phi = sol.eval(Phi);
    V = x'*Phi*x;
%     % SOS program 
%     prog = mssprog;
%     
%     % Some variables needed later
%     [prog,bslack] = new(prog,1,'pos');
%     [prog,dummy] = new(prog,1,'pos');
%     
%     % Create rho
%     [prog,rho] = new(prog,1,'pos');
%     
%     % Create Phi and new V
%     [prog,Phi] = new(prog,length(x),'psd');
%     V = x'*Phi*x;
%     
%     % Normalization constraint
%      prog.eq = trace(Phi) - trace(S0);
%     
%     % Compute Vdot  
%     Vdot = diff(V,x)*subs(f,u,uf);
%        
%     % Clean stuff
%     V = clean(V,options.clean_tol);
%     Vdot = clean(Vdot,options.clean_tol);
%     
%     % Declare SOS conditions
%     prog.sos = -Vdot + L1*(V-rho);
%             
%     % Solve SOS program  
%     pars = struct();
%     pars.fid = 0;
%     [prog,info] = sedumi(prog,-rho,1,pars,1) 
%         
%     disp('Backing off now...')
%        
%     % Back off on objective
%     rhof = double(prog(rho));
%         
%     prog.eq = bslack - (-(1 - options.backoff_percent/100)*rhof + rho);
%    
%     [prog,info] = sedumi(prog,dummy,1,pars,0) 
%     if info.numerr == 2 || info.pinf == 1 || info.dinf == 1
%         keyboard;
%     end
%      
%     rho = double(prog(rho));
%     V = x'*double(prog(Phi))*x;
%     Phi = double(prog(Phi));
end




