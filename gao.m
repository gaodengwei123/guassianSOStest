clear all
close all
clc
dbstop if error
figrue
hold on
axHandle = gca;
xlabel(axHandle,'x1');
ylabel(axHandle,'x2');
[x1,x2] = meshgrid(-5:0.2:5,-5:0.2:5);
u = -x2;
v = -x2.*(1-x1.^2)+x1;
sumuv = sqrt(u.^2+v.^2);
u = u./sumuv;
v = v./sumuv;
quiver(x1,x2,u,v,0.5)