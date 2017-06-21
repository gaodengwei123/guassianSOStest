% /*! @SampleReachabilityFunnel.m
% *************************************************************************
% <PRE>
% file.name       : SampleReachabilityFunnel.m
% related files   :
% function&ablity : Calculate Invariant_Funnels by sample based method(Funnel Libraries for Real-Time Robust Feedback Motion Planning)
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/11/12    1.00                     build
% 2016/12/06    2.00                     add Taylor Expansion in the system
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
function [Vtraj0,rho] = SampleReachabilityFunnel(sys,Vtraj0,ts,options)
N = length(ts);
x = Vtraj0.x;
rho = zeros(N,1);
rhof = 1;
xdot = cell(1,N);
Tayxdot = TaylorExpansion(sys,3); % taylor expantion in order
V = cell(1,N);
Vdot = cell(1,N);
for i = 1:N
    V{i} = Vtraj0.getPoly(ts(i));
    xdot{i} = Tayxdot(ts(i));
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
        V{i}=clean(subss(V{i},x,T*x));
        Vdot{i}=clean(subss(Vdot{i},x,T*x));
        %         f{i}=clean(subss(xdot{i},x,T*x));
    end
end

parfor i = 1:N
    rho(i) = SampleBasedRho(x,V{i},Vdot{i});
    i
end

for i = 1:N
    rho(i) = rho(i)/rho(end)*rhof;
end
% Vtraj0 = Vtraj0.updateV(foh(ts,1./rho'));
if (options.plot_rho)
    rhopp=foh(ts,rho');
    figure(10);
    fnplt(rhopp);
    drawnow;
end
end

function rho = SampleBasedRho(x,V,dV)
SmpleNum = 5000;
c = 1000000;
vx = 10-20*rand(4,SmpleNum);
for i = 1:SmpleNum
    if double(subs(dV,x,vx(:,i)))>=0&&double(subs(V,x,vx(:,i)))<c
        c = double(subs(V,x,vx(:,i)));
    end
end
rho = c;
end

function y=doubleSafe(x)
y=double(x);
if (~isa(y,'double')) error('double failed'); end
end
