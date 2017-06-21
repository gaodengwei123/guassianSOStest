clear all
close all
clc
dbstop if error
%%
Pendsys = Pendulumdyn();
% Pendsys2 = doublePendulumdyn();

A = @(x)([0  1; -9.8*cos(x(1)) -0.1]);
B = [0;1];
Q = diag([1,1]);R=1;
range = [-pi -2.5*pi pi 2.5*pi]; 
%% dynmaic plot
% 
% figure
% hold on
% axHandle = gca;
% xlabel(axHandle,'\theta');
% ylabel(axHandle,'$\dot{\theta}$','Interpreter','latex')
% [x1,x2] = meshgrid(-4:0.2:4,-8:0.2:8);
% u = x2;
% v = -0.1*x2 - 9.8*sin(x1) - 19.6010*x1 - 6.1635*x2 + 19.6010*pi;
% sumuv = sqrt(u.^2+v.^2);
% u = u./sumuv;
% v = v./sumuv;
% quiver(x1,x2,u,v,0.5)


%% ROA calculation
% [K,S,~] = lqr(A([pi;0]),B,Q,R);
% % K = [19.6010    6.1635];
% Pend = feedbackSys(Pendsys,K,B);
% 
% options.max_iterations = 100;
% options.degV = 6;
% options.converged_tol = 1e-2;
% options.method={'levelset'};
% 
% x0 = [pi;0];
% t = msspoly('t');
% u = msspoly('u',1);
% x = msspoly('x',2);
% % x-[1;1]
% V0 = V_function(t,x,S,[]);
% V0.x0 = x0;
% V = regionOfAttraction(Pend,V0,options);
% % plot lyapunov function
% options.x0 = x0;
% xfun0 = getLevelSet(x,V,options);
% plot(xfun0(1,:),xfun0(2,:),'b','lineWidth',2)

%% plot the trajectory
% odefun = @(t,x)Pendsys.dynamics(t,x,-K*(x-[pi;0]));
% sol = ode45(odefun,[0,2],[0;0]);
% xnew = deval(sol,0:0.01:2);
% plot(xnew(1,:),xnew(2,:),'r')

%% ===========  lqr rrt star=======
sys.A = A;
sys.B = @(x)[0;1];
sys.f = Pendsys;
sys.Q = Q;
sys.R = R;
sys.uppercontrol = 3;
Num = 5000;

xT = [pi;0];
x0 = [0;0];
% plot(xT(1),xT(2),'o');

G = LQR_trees_star.LQR_trees_search(sys,Num,range,xT,x0);


