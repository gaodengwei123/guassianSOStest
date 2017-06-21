clear
clc
figure(1)
dbstop if error
syms x1 x2 real
f = @(x1,x2)([-x2
    x1+x2*(x1^2-1)]);
f1 = -x2;
f2 = x1-x2+x1^2*x2;
% f = @(x1,x2)([x2
%               -sin(x1)-0.5*x2]);
% f1 = x2;
% f2 = -sin(x1)-0.5*x2;
% f = @(x1,x2)([-2*x2+x1*x2
%               -x2+x1*x2]);
% f1 = -2*x2+x1*x2;
% f2 = -x2+x1*x2;
% f = @(x1,x2)([x1
%               -x1+1/3*x1^3-x2]);
% f1 = x1;
% f2 = -x1+1/3*x1^3-x2;
% f = @(x1,x2,x3)([-x1+x2*x3^2
%                  -x2+x1*x2
%                  -x3]);
% f1 = -x1+x2*x3^2;
% f2 =-x2+x1*x2;
% f3 =-x3;

temp = 0.1;
Lmax = 20;
tmax = 10;
r0 = 2;
derx = 0.05;
R = [];
lx1 = -3:derx:3;
lx2 = -3:derx:3;
[sx1,sx2]=meshgrid(lx1, lx2);
num1 = size(sx1,1);
num2 = size(sx1,2);
tic
for i =1:num1*num2
    t = 0;
    L = 0;
    xnext = [sx1(i) sx2(i)]';
    while t<tmax&&xnext'*xnext>2
        t = t+temp;
        xnext = xnext+f(xnext(1),xnext(2))*temp;
        L = L+sqrt(f(xnext(1),xnext(2))'*f(xnext(1),xnext(2)))*temp;
        if L>Lmax
            L = 2*Lmax;
            break;
        end
    end
    
    R = [R;L];
end
toc
sx3 = reshape(R,num1,num2);
contour(sx1,sx2, sx3)
figure
surf(sx1,sx2, sx3); % 画出立体曲面图
colorbar





