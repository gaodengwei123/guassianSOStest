% /*! @myReachabilityFunnel.m
% *************************************************************************
% <PRE>
% file.name       : myReachabilityFunnel.m
% related files   :
% function&ablity : Calculate Invariant_Funnels(Funnel Libraries for Real-Time Robust Feedback Motion Planning)
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/4/26     1.00                     build
% 2016/7/26     1.00                     rebuild form drake
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
function [Vtraj0,rho] = myReachabilityFunnel(sys,Vtraj0,ts,options)
dbstop if error

checkDependency('mosek');   % this toolbox have dowloaded in 2016.7.27, and have one year lincese.[javaaddpath('/full/path/to/mosekmatlab.jar')]
checkDependency('sedumi');  % this is free

% Get the necessary variables
N = length(ts);
x = Vtraj0.x;
u= sys.p_u;
xdot = cell(1,N);
xdotOrig = cell(1,N);
V = cell(1,N);
Vdot = cell(1,N);
[Tayxdot, TayxdotOrig]= TaylorExpansion(sys,3); % taylor expantion in order

num_u = sys.getNumInput;

if options.saturations && (num_u > 1)
    error('Sorry, I cannot handle actuator saturations for systems with more than one actuator.')
end

Vmin = zeros(N-1,1);
%% evaluate dynamics and Vtraj at every ts once (for efficiency/clarity)
for i=1:N
    V{i} = Vtraj0.getPoly(ts(i));
    xdot{i} = Tayxdot(ts(i));
    xdotOrig{i} = TayxdotOrig(ts(i));
    if (sys.getNumInput>0)   % zero all inputs
        xdot{i} = subs(xdot{i},sys.p_u,zeros(sys.getNumInput,1));
    end
    
    dVdt=Vtraj0.getPolyTimeDeriv(ts(i));
    Vdot{i}=diff(V{i},x)*xdot{i} + dVdt;
    
    % balancing
    if (strcmp(options.lyap_parameterization,'rho'))
        S1=.5*doubleSafe(subs(diff(diff(V{i},x)',x),x,0*x));
        S2=.5*doubleSafe(subs(diff(diff(Vdot{i},x)',x),x,0*x));
        [T,D] = balanceQuadForm(S1,S2);
%         Ts{i}=T;
        V{i}=clean(subss(V{i},x,T*x));
        Vdot{i}=clean(subss(Vdot{i},x,T*x));
%         f{i}=clean(subss(xdot{i},x,T*x));
    end
    %       Vmin(i) = minimumV(x,V{i});
    i
    
    % form new rho finding
    ui{i} = sys.FunContr.eval(ts(i)) + sys.K(ts(i))*x;
    if i > 1
        Phi{i} = 0.5*eye(length(x));
    else
        Phi{i} = zeros(length(x),length(x));
    end
    Vy{i} = Vtraj0.getPoly(ts(i)) + x'*Phi{i}*x;
    Vmin(i) =  minimumV(x,Vy{i});
end
%%
% update the Vdot and V
% Vtraj0.Vdot = Vdot;
%% Initialize rho with "tube" constraints
if (~isfield(options,'degL1'))
    options.degL1 = 4;%deg(Vdot{1},x);  % just a guess
end

options.converged_tol = -Inf;

%% Initialize rho
dts = diff(ts);
% rhof = 1;   % todo: handle the more general case and get it from containment
% rho = flipud(rhof*exp(-options.rho0_tau*(ts-ts(1))/(ts(end)-ts(1))))+max(Vmin);
% rhodot = diff(rho)./dts;
rho = rhoLineSearch(Vtraj0,Vy,sys.FunTraj,ts,xdotOrig,Phi,dts,options,u,ui,x,Tayxdot);
% load rho
rhof = 1;
rhodot = diff(rho)./dts;
% check accuracy by sampling
for i=1:N-1
    m(i)=sampleCheck(x,V{i},Vdot{i},rho(i),rhodot(i));
end
if (max(m)>0)
  figure(4);clf;fnplt(foh(ts,rho')); 
  figure(5);clf;plot(ts(1:end-1),m); drawnow;
  error('infeasible rho. increase options.rho0_tau');
end

% perform bilinear search the actual verification here
rhointegral = 0;
%%
for iter=1:options.max_iterations
    last_rhointegral = rhointegral;
    L=findMultipliers(x,V,Vdot,rho,rhodot,options);
    
    [rho,rhointegral]=optimizeRho(x,V,Vdot,L,dts,Vmin,rhof,options);
    
    rhodot = diff(rho)./dts;
    % plot current rho
    if (options.plot_rho)
        rhopp=foh(ts,rho');
        figure; fnplt(rhopp); title(['iteration ',num2str(iter)]); drawnow;
    end
    save datarho
    % check for convergence
    if ((rhointegral - last_rhointegral) < options.converged_tol*last_rhointegral)  % see if it's converged
        break;
    end
end
save data20170122
% load('data20170121.mat')



end

% fix lagrange multipliers, optimize rho
function [rho,rhointegral]=optimizeRho(x,V,Vdot,L,dts,Vmin,rhof,options)
N = length(V)-1;
prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog,rho] = prog.newPos(N,1);
rho = [rho;rhof]+Vmin;
rhointegral=0;
for i=1:N
    rhodot(i,1) = (rho(i+1)-rho(i))/dts(i);
    rhointegral = rhointegral+rho(i)*dts(i)+.5*rhodot(i)*dts(i)^2;
    if (options.stability)
        Vdot{i} = Vdot{i}*rho(i)-V{i}*rhodot(i);
        prog = prog.withSOS(-Vdot{i} + L{i}*(V{i}-rho(i)));
    else
        prog = prog.withSOS(-(Vdot{i}-rhodot(i)+L{i}*(V{i}-rho(i))));
    end
end
pars = spot_sdp_default_options();
% pars.verbose = 1;
sol = prog.minimize(-rhointegral,@spot_mosek,pars);
if strcmp(sol.status,'STATUS_PRIMAL_INFEASIBLE')
    load chirp
    sound(y,Fs)
end
rho = double(sol.eval(rho));
rhointegral = double(sol.eval(rhointegral));

%% ========================
% N = length(V)-1;
% 
% prog = mssprog;
% [prog,rho] = new(prog,N,'pos');
% rho = [rho;rhof]+Vmin;
% rhointegral=0;
% for i=1:N
%     rhodot(i,1) = (rho(i+1)-rho(i))/dts(i);
%     rhointegral = rhointegral+rho(i)*dts(i)+.5*rhodot(i)*dts(i)^2;
%     
%     if (options.stability)
%         Vdot{i} = Vdot{i}*rho(i)-V{i}*rhodot(i);
%         prog.sos = -Vdot{i} + L{i}*(V{i}-rho(i));
%     else
%         prog.sos = -(Vdot{i}-rhodot(i)+L{i}*(V{i}-rho(i)));
%     end
% end
% %  prog = sedumi(prog,ones(1,N)*slack-rhointegral,1,struct());
% [prog,info] = sedumi(prog,-rhointegral,0); %1,struct());
% if (info.numerr>1)
%     keyboard;
%     error('sedumi failed due to numerical errors');
% end
% rho = doubleSafe(prog(rho));
% rhointegral = doubleSafe(prog(rhointegral));
end

% fix rho, optimize lagrange multipliers
function L=findMultipliers(x,V,Vdot,rho,rhodot,options)
% note: compute L for each sample point
disp('Step 1: Searching for multipliers.')
N = length(V)-1;
%% mousk comptete
parfor i=1:N
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    % Declare multipliers
    L1m = monomials(x,0:options.degL1);
    [prog,l1] = prog.newFree(length(L1m));
    L1 = l1'*L1m;
    % Create gammas
    [prog,gamma] = prog.newFree(1);
    % Declare SOS conditions
    if (options.stability)
        Vdot{i} = Vdot{i}*rho(i)-V{i}*rhodot(i);  % intentionally skip 1/rho^2 term
        prog = prog.withSOS(gamma*V{i} - Vdot{i} + L1*(V{i}-rho(i)));% + L2*(.1-V{i});
        prog = prog.withSOS(L1);
    else
        prog = prog.withSOS(gamma-(Vdot{i}-rhodot(i) + L1*(V{i}-rho(i))));
%         prog = prog.withSOS(-gamma*(x'*x)^(deg(Vdot{i},x)/2) - Vdot{i} + rhodot(i) + L1*(V{i}-rho(i)))
    end
    % Solve SOS program
    pars = spot_sdp_default_options();
%     pars.verbose = 1;
    sol = prog.minimize(gamma,@spot_mosek,pars);
    L{i} = sol.eval(L1);
%     slack{i}=double(sol.eval(gamma));
end

%% =========================
% N=10;
% parfor i=1:N
%     prog = mssprog;
%     Lxmonom = monomials(x,0:options.degL1);
%     [prog,l] = new(prog,length(Lxmonom),'free');
%     L1 = l'*Lxmonom;
%     
%     [prog,gamma] = new(prog,1,'free');
%     if (options.stability)
%         
%         %      [prog,l2] = new(prog,length(Lxmonom),'free');
%         %      L2 = l2'*Lxmonom;
%         
%         Vdot{i} = Vdot{i}*rho(i)-V{i}*rhodot(i);  % intentionally skip 1/rho^2 term
%         prog.sos = gamma*V{i} - Vdot{i} + L1*(V{i}-rho(i));% + L2*(.1-V{i});
%         prog.sos = L1;
%         %      prog.sos = L2;
%     else
%         prog.sos = gamma-(Vdot{i}-rhodot(i) + L1*(V{i}-rho(i)));
%     end
%     
%     [prog,info{i}] = sedumi(prog,gamma,0);
%     
%     if (info{i}.pinf==0 && info{i}.dinf==0)
%         slack{i}=double(prog(gamma));
%         L{i} = prog(L1);
%     end
% end
%slack

% for i=fliplr(1:N)
%     if (slack{i}>1e-4 || info{i}.pinf~=0 || info{i}.dinf~=0)
%         if (length(x)~=2)
%             dims=[2;4]; d=ones(length(x),1); d(dims)=0; d=logical(d);
%             [m,b]=minimumV(x,V{i});
%             V{i}=subs(V{i},x(d),b(d));
%             Vdot{i}=subs(Vdot{i},x(d),b(d));
%             x=x(dims);
%         end
%         figure(2); clf; [minVmRho,junk]=plotPoly(x,V{i}-rho(i));
%         if (options.stability)
%             figure(1); clf; [junk,maxVdot]=plotPoly(x,Vdot{i});
%             fprintf(1,'segment %d of %d, slack=%f,\n sampledmax(d/dt(V/rho))=%f, sampledmin(V-rho)=%f\n info:\n',i,N,maxVdot,full(minVmRho),slack{i});
%         else
%             figure(1); clf; [junk,maxVdotmRhodot]=plotPoly(x,Vdot{i}-rhodot(i));
%             fprintf(1,'segment %d of %d, slack=%f,\n sampledmax(Vdot-rhodot)=%f, sampledmin(V-rho)=%f\n info:\n',i,N,full(maxVdotmRhodot),full(minVmRho),slack{i});
%         end
%         info{i}
%         error('rho is infeasible');
%     end
% end
end


function [mi,ma]=plotPoly(x,P,rho)
if(nargin<3) rho=0; end
[X1,X2]=ndgrid(-2:.1:2,-2:.1:2);
Ps=reshape(doubleSafe(msubs(P,x,[X1(:)';X2(:)'])),size(X1));
mi=min(min(Ps));
ma=max(max(Ps));
surf(X1,X2,Ps); colorbar;
view(0,90);
hold on;
[c,h]=contour3(X1,X2,Ps,[rho,rho]);
set(h,'EdgeColor',[1 1 1],'LineWidth',4);
end

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

function m=sampleCheck(x,V,Vdot,rho,rhodot)
if (deg(V,x)>2) error('only checks quadratics'); end

n=length(x);
K=100;
X = randn(n,K);
X = X./repmat(sqrt(sum(X.^2,1)),n,1);

H = doubleSafe(0.5*diff(diff(V,x)',x));
b = -0.5*(H\doubleSafe(subs(diff(V,x),x,0*x)'));

try
    X = repmat(b,1,K) + (H/(doubleSafe(rho-subs(V,x,b))+eps))^(-1/2)*X;
catch
    keyboard;
end
m=max(doubleSafe(msubs(Vdot,x,X))) - rhodot;
if (m>0)
    warning('found a positive Vdot');
end
end


function y=doubleSafe(x)
y=double(x);
if (~isa(y,'double')) error('double failed'); end
end

%% add form sampledFiniteTimeReachabilityFunnel to find a feasible rho

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
    rho(k+1) = 1.01*rho(k+1); % 1.01 To ensure feasibility
    
end

end

function gamma = checkRho(Vtraj0,Vy,rhonext,utraj,ts,forig_u,Phi,dts,options,u,ui,x,k,rho,psys)
% Compute rhodot
rhodot = (rhonext - rho)/dts(k);

prog = spotsosprog;
prog = prog.withIndeterminate(x);

Phidotk = (Phi{k+1} - Phi{k})/(ts(k+1) - ts(k));
V0k = Vtraj0.getPoly(ts(k));
V0dotk = Vtraj0.getPolyTimeDeriv(ts(k));

% Compute Vdot
Vdoty = diff(V0k,x)*subss(forig_u{k},u,ui{k}) + V0dotk + 2*x'*Phi{k}*subss(forig_u{k},u,ui{k}) + x'*Phidotk*x;

% Clean stuff
V = clean(Vy{k},options.clean_tol);
Vdoty = clean(Vdoty,options.clean_tol);

% Declare multipliers
L1m = monomials(x,0:options.degL1);
[prog,l1] = prog.newFree(length(L1m));
L1 = l1'*L1m;

% Create gammas
[prog,gamma] = prog.newFree(1);

% Declare SOS conditions
prog = prog.withSOS(-gamma*(x'*x)^(deg(Vdoty,x)/2) - Vdoty + rhodot + L1*(V-rho));

% Solve SOS program
pars = spot_sdp_default_options();
pars.verbose = 1;
sol = prog.minimize(-gamma,@spot_mosek,pars);

if strcmp(sol.info.solverInfo.itr.prosta,'PRIMAL_INFEASIBLE')
    gamma = -1.0;
else
    
    gamma = double(sol.eval(gamma));
    
    if strcmp(sol.info.solverInfo.itr.prosta,'UNKNOWN')
        gamma = -1.0;
    end
end
end