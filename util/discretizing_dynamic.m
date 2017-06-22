% transform the dynamic system to discret
% build by dengwei 2017
function g = discretizing_dynamic(sys,x0,u0,temp)

k1 = sys.dynamics(t,x0,u0);
k2 = sys.dynamics(t,x0+k1/2,u0);
k3 = sys.dynamics(t,x0+k2/2,u0);
k4 = sys.dynamics(t,x0+k3,u0);
g = x0+temp/6*(k1+2*k2+2*k3+k4);

end




