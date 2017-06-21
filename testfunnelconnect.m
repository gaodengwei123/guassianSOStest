clear all
clc
close all
%%
x = msspoly('x',2);
% x0 = [2;3];
% P1 = [1 0;0 1];
% P2 = [4 1;1 3];
% V1 = (x-x0)'*P1*(x-x0);
% V2 = (x-x0)'*P2*(x-x0);
% x0 = [22.6000;38.4000];%;0.4266;0
x0 = [0;0];
x1 = x(1);x2=x(2);
V1 = (1925.6)+(-43.839)*x1+(0.9699)*x1^2+(-74.488)*x2+(0.9699)*x2^2+(-8.4585e-14)*x2*x1;
V2 = (3045.3)+(-68.921)*x1+(1.5259)*x1^2+(-118.05)*x2+(1.5374)*x2^2+(-0.0012602)*x2*x1;
P1 = double(0.5*diff(diff(V1,x)',x));
P2 = double(0.5*diff(diff(V2,x)',x));
V1 = (x-x0)'*P1*(x-x0);
V2 = (x-x0)'*P2*(x-x0);

options.x0 = x0;
options.plotdims=[1 2];
xfun1 = getLevelSet(x,V1,options);
xfun2 = getLevelSet(x,V2,options);
figure
hold on
plot(xfun1(1,:),xfun1(2,:),'b','LineWidth',1);
plot(xfun2(1,:),xfun2(2,:),'r','LineWidth',1);
tic
prog = spotsosprog;
prog = prog.withIndeterminate(x);
Lmonom = monomials(x,0:2);
[prog,rho] = prog.newFree(1);
[prog,L] = prog.newSOSPoly(Lmonom);
prog = prog.withSOS(L*(1 - V2)-(rho-V1));% V1 is in V2
option = spot_sdp_default_options();
sol = prog.minimize(-rho,@spot_mosek,option);
rho = double(sol.eval(rho));
toc
V1 = (x-x0)'*(P1/rho)*(x-x0);
xfun3 = getLevelSet(x,V1,options);
plot(xfun3(1,:),xfun3(2,:),'k','LineWidth',1);
