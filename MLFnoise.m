function Ell = MLFnoise(sys,Vtraj)
ts = sys.breaks;
N = length(ts);
x = Vtraj.x;
xdot = cell(1,N);
Tayxdot = TaylorExpansion(sys,3); % taylor expantion in order
V = cell(1,N);
Vdot = cell(1,N);
for i = 1:N
    V{i} = Vtraj.getPoly(ts(i));
    xdot{i} = Tayxdot(ts(i));
    if (sys.getNumInput>0)   % zero all inputs
        xdot{i} = subs(xdot{i},sys.p_u,zeros(sys.getNumInput,1));
    end
    dVdt=Vtraj.getPolyTimeDeriv(ts(i));
    Vdot{i}=diff(V{i},x)*xdot{i} + dVdt;
    
    % balancing
    
    S1=.5*doubleSafe(subs(diff(diff(V{i},x)',x),x,0*x));
    S2=.5*doubleSafe(subs(diff(diff(Vdot{i},x)',x),x,0*x));
    [T,D] = balanceQuadForm(S1,S2);
    V{i}=clean(subss(V{i},x,T*x));
    Vdot{i}=clean(subss(Vdot{i},x,T*x));
    
end
xs = sym('x',[1 4],'real');
xs = xs(:)

for i = 1:length(ts)
    V{i} = msspoly2sym(x,xs,V{i});
    Vdot{i} = msspoly2sym(x,xs,Vdot{i});
    xdot{i} = msspoly2sym(x,xs,xdot{i});
    
    % error ellipsoid
    Q = diag([0.1,0.1,0.1,0.1]);
    ErrorEll{i} = xs'*Q*xs;
    ErrorElldot{i} = (2*Q*xs)'*xdot{i};
    
    VN{i} = V{i}/(1-ErrorEll{i});
    VNdot{i} = Vdot{i}/(1-ErrorEll{i})+V{i}*ErrorElldot{i}/(1-ErrorEll{i})^2;
    
end

for i = 1:N
    rho(i) = SampleBasedRho(xs,VN{i},VNdot{i});
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


