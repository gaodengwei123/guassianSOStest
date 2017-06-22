% /*! @ellipsoid2Polyhedron.m
% *************************************************************************
% <PRE>
% file.name       : ellipsoid2Polyhedron.m
% related files   :
% function&ablity : transprot the ellipsoid to Polyhedron
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/8/25     1.00                     build
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

function [Poly, Rho] = ellipsoid2Polyhedron(origin0,x,ellipse)

options.inclusion = 'slice';
options.color=.7*[1 1 1];
options.tol = .01;

options.plotdims = [1;2];
options.x0 = origin0;
Poly1 = getLevelSet(x,ellipse,options);
num = size(Poly1,2);
Rho1 = Poly1-repmat(options.x0,1,num);
Rho = zeros(1,num);
Poly = zeros(2,num);        % 2 is the display dimension
for i = 1:num
    Rho(i) = norm(Rho1(:,i));
    Poly(:,i) = Poly1(1:2,i)/Rho(i);
end

