% /*! @FunnelImprovebyLagrange.m
% *************************************************************************
% <PRE>
% file.name       : FunnelImprovebyLagrange.m
% related files   :
% function&ablity : Calculate Invariant_Funnels(reference:Control Design along Trajectories with Sums of Squares Programming)
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/4/26     1.00                     build
% 2016/7/20     1.00                     rebuild form drake
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
function [c,V,VS1] = FunnelImprovebyLagrange(sys,V,c0,tv,ts,xx,u0,options)
% (sys,V{i},feedbackU{i},tv,taus(i),xx,u0,options)
checkDependency('mosek');
checkDependency('sedumi');

t = msspoly('t');
u = msspoly('u',options.m);
x = msspoly('x',options.n);
% sdpvar t
% u = sdpvar(options.m,1);
% x = sdpvar(options.n,1);
%% start Time-Varying Design
% Compute fmax and fmin if saturations
uf = c0;
f = sys.f(ts,x,u);
% u0 = subs(uf,x,x0);% current state
% if uf>sys.umax||uf<sys.umin
if ~all(isinf([sys.umax;sys.umin]))
    f_umax = subs(f,u,sys.umax);
    f_umin = subs(f,u,sys.umin);
    saturations = true;
else
    saturations = false;
end
    saturations = false;
% Initialize u, Phi and rho

S0 = V.gettime(ts).S;
rho = options.rho0_tau;
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
       [L1,~] = findLU(V,rho,x,f,u,options);
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

  % check for convergence
  if ((rho - rho_last) < options.converged_tol*rho_last)  % see if it's converged
   break;
  end
  
  rho_last = rho;
  
end

% Optimized Lyapunov function
% V = QuadraticLyapunovFunction(Vframe,double(Phi),zeros(length(x),1),0);
VS1 = double(Phi);
% c = SpotPolynomialSystem(Vframe,[],OrigUFrame,[],[],[],uf,[]);
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
    prog = mssprog;
    
    % Create u
    um = monomials(x,1:options.controller_deg);
    lu = [];
    for j = 1:length(u)
        [prog,luj] = new(prog,length(um),'free');
        lu = [lu;luj'];
    end
    uf = lu*um;
    
    % Compute Vdot  
    Vdot = diff(V,x)*subs(f,u,uf);
       
    % Clean stuff
    V = clean(V,options.clean_tol);
    Vdot = clean(Vdot,options.clean_tol);
    
    % Declare multipliers
    L1m = monomials(x,options.degL1:options.degL1); % Make it homogeneous
    [prog,l1] = new(prog,length(L1m),'free');
    L1 = l1'*L1m;
    prog.sos = L1;
    
    % Create gammas
    [prog,gamma] = new(prog,1,'pos');
    
    % Declare SOS conditions
    prog.sos =  -gamma*(x'*x)^(deg(Vdot,x)/2) - Vdot + L1*(V-rho);
            
    % Solve SOS program  
%     pars = struct();
%     pars.fid = 0;
%     [prog,dummy] = new(prog,1,'free');
%     [prog,info] = sedumi(prog,dummy,1,pars,0) 
%     
%     % Optimized multipliers
%     L1f = prog(L1);
    pars = spot_sdp_default_options();
    pars.verbose = 0;
    sol = prog.minimize(0,@spot_mosek,pars);
    
    % Optimized multipliers
    L1f{k} = sol.eval(L1);  

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
    prog = mssprog;
    
    % Some variables needed later
    [prog,bslack] = new(prog,1,'pos');
    [prog,dummy] = new(prog,1,'pos');

    % Declare rho
    [prog,rho] = new(prog,1,'pos');

    % Create u
    um = monomials(x,1:options.controller_deg);
    lu = [];
    for j = 1:length(u)
        [prog,luj] = new(prog,length(um),'free');
        lu = [lu;luj'];
    end
    uf = lu*um;
    
    % Compute Vdot  
    Vdot = diff(V,x)*subs(f,u,uf);
    
    % Clean stuff
    V = clean(V,options.clean_tol);
    Vdot = clean(Vdot,options.clean_tol);
    
    % Declare SOS conditions
    prog.sos = -Vdot + L1*(V-rho);
             
    % Solve SOS program
    pars.fid = 0;
    [prog,info] = sedumi(prog,-rho,1,pars,1)   
    
    if info.numerr == 2 || info.pinf == 1 || info.dinf == 1
        keyboard;
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
       
    % Clean stuff
    V = clean(V,options.clean_tol);
    Vdot = clean(Vdot,options.clean_tol);
    
    % Declare SOS conditions
    prog.sos = -Vdot + L1*(V-rho);
            
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
%% add some part with V
function [Vmin,b] = minimumV(x,V)
    if (deg(V,x)>2)
        prog = mssprog;
        [prog,slack] = new(prog,1,'free');
        prog.sos = slack + V;
        [prog,info] = sedumi(prog,slack,0);
        Vmin = -doubleSafe(prog(slack));
    else
        H = doubleSafe(0.5*diff(diff(V,x)',x));
        b = -0.5*(H\doubleSafe(subs(diff(V,x),x,0*x)'));
        Vmin = subs(V,x,b);
    end
end

%%  Initialize rho
function rho = rhoLineSearch(Vtraj0,Vy,utraj,ts,forig_u,Phi,dts,options,u,ui,x,psys)

    N = length(ts)-1;
    disp('Step 0: Initializing rho with bisection search...')

    rho = zeros(N+1,1);
    rho(1) = options.rho0;

    for k = 1:N
        rhonow = rho(k);
        rhomin = 0.5*rhonow;
        rhomax = 10*rhonow;
        rho(k+1) = fzero(@(rhonext) checkRho(Vtraj0,Vy,rhonext,utraj,ts,forig_u,Phi,dts,options,u,ui,x,k,rhonow,psys),[rhomin rhomax],optimset('TolX',1e-5));
        rho(k+1) = 1.01*rho(k+1) % 1.01 To ensure feasibility
        
    end

end