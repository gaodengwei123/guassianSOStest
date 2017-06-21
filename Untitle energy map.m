% drow energy map
clear
close all
clc
A = @(x)([0  1; -9.8*cos(x(1)) -0.1]);
eng = @(x1,x2)(4.9*(1-cos(x1))+1/8*x2.^2);

B = @(x)[0;1];
Q = diag([1,1]);
R=1;
figure
hold on

range = [-2*pi -2*pi 2*pi 2*pi];
Num = 100;
Randomp = haltonset(2);
Randomp = scramble(Randomp,'RR2');
X0 = net(Randomp,Num+1);

axHandle = gca;
xlabel(axHandle,'\theta');
ylabel(axHandle,'$\dot{\theta}$','Interpreter','latex')
[x1,x2] = meshgrid(range(1):0.5:range(3),range(2):0.5:range(4));
u = x2;
v = -0.1*x2 - 9.8*sin(x1);
sumuv = 0.1;%sqrt(u.^2+v.^2);
u = u./sumuv;
v = v./sumuv;
colorA = colormap(jet);
maxeng = max(max(eng(x1,x2)));
engindex = @(x1,x2)(ceil(eng(x1,x2)/maxeng*64));
% quiver3(x1,x2,eng(x1,x2),u,v,zeros(size(u)),0.5)%,'Color',colorA(engindex(x1,x2),:)

%% energy map

mesh(x1,x2,eng(x1,x2))

load('traj2017516', 'traj')
plot3(traj.state(:,1),traj.state(:,2),eng(traj.state(:,1),traj.state(:,2))+0.1,'r','LineWidth',2);
%%
K = 20;
plot(pi,0,'o')
for i = 1:Num
    xrand = [range(1)+(range(3)-range(1))*X0(i,1);range(2)+(range(4)-range(2))*X0(i,2)];

    Ax = A(xrand);
    Bx = B(xrand);
    [~,S0,~] = lqr(Ax,Bx,Q,R);
    th=linspace(0,2*pi,K);
    X = [sin(th);cos(th)];
    S = S0/.1;
    y = repmat(xrand,1,K) + S^(-1/2)*X;
    z = eng(y(1,:),y(2,:));
    plot3(y(1,:),y(2,:),z,'r')
    

    [~,S1,~] = lqr(Ax,Bx,Q,R*50);
    th=linspace(0,2*pi,K);
    X = [sin(th);cos(th)];
    S = S1/.1;
    y = repmat(xrand,1,K) + S^(-1/2)*X;
    z = eng(y(1,:),y(2,:));
    plot3(y(1,:),y(2,:),z,'r')
    
    
end
axis equal
grid on











