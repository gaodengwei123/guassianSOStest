% /*! @lqrbaseddistance.m
% *************************************************************************
% <PRE>
% file.name       : lqrbaseddistance.m
% related files   :     
% function&ablity : 
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         : 
% --------------------------------------------------------------------------------
% record of modify : 
% date          version     name         content 
% 2016/04/25    1.00                     build
% </PRE>
% ********************************************************************************
% 
% * right(c)
% ¡¡
% *************************************************************************
% input:
%        
% output:
% *************************************************************************
function [minindex, mindistance] = lqrbaseddistance(nodes, stateno, point, speed1,speed2)

mindistance = -1;
minindex = -1;


for l=1:stateno
   xbar(1:2) = nodes(l,1:2) - point;
   xbar(3:4) = nodes(l,5:6) - [speed1 speed2];
   % We have a polynomial dynamical system:
   f = @(x,u) [x(3);  x(1)+(x(3)+u(1)/2); x(4);  x(2)+(x(4)+u(2)/2)];
   
   u0 = [0;0];

   % Pick Q and R
   Q = eye(4); R = 1;
   R = [1 0; 0 1];
   if size(xbar,2) ~= 1, xbar = xbar'; end
   if size(xbar,2) ~= 1, error('x0 must be a column.'); end
   if size(u0,2) ~= 1, error('u0 must be a column.'); end
   
   n = length(xbar); m = length(u0);
   
   if size(Q) ~= [n n], error('Q must be n-by-n.'); end
   if size(R) ~= [m m], error('R must be m-by-m.'); end
    
   x = msspoly('x',n);
   u = msspoly('u',m);
   xdot = msspoly(f(x,u));
% 
%    if norm(double(subs(xdot,[x;u],[xbar;u0])),1) > 1e-16
%        error('(x0,u0) is not an equilibrium');
%    end
    
   AB = double(subs(diff(xdot,[x;u]),[x;u],[xbar;u0]));
   A = AB(1:n,1:n);
   B = AB(1:n,n+(1:m));
   [K,S] = lqr(A,B,Q,R);
   syms r;
   d(:,:,1) = exp(A*100) * nodes(l,1) + int(exp(A*(100-r)), 0, 100); 
   d(:,:,2) = exp(A*100) * nodes(l,2) + int(exp(A*(100-r)), 0, 100);
   mind = -1;
   d1 = d(:,:,1);
   d2 = d(:,:,2);
   for k=1:1
      current1(:,:) = 1/2 * double(d1' * S *d1);
      current2(:,:) = 1/2 * double(d2' * S *d2);
      c1 = norm(current1) / 2 +k;
      c2 = norm(current2) / 2 +k;
      c=norm([c1 c2]);
      if mind == -1 || c < mind
          mind = c;
      end
   end
   if mindistance == -1 || mind < mindistance
      minindex = l;
      mindistance = mind;
   end
%    stateno2 = 
%    for k = 1:
   
end