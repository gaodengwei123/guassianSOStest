% /*! @plotEllipsoid.m
% *************************************************************************
% <PRE>
% file.name       : plotEllipsoid.m
% related files   :
% function&ablity : plot Ellipsoid
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/1/29     1.00        dengwei      build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:
% control;
%
% output:
% *************************************************************************
function plotEllipsoid(x)
% must 3 dim plot
[n,N] = size(x);
if n~=3
    error('Ellipsoid must be 3 dim')
end
clr = [1 0 0]; % red color
Alpha = 0.4;   % shadow
chll = convhulln(x');
patch('Vertices', x', 'Faces', chll, ...
    'FaceVertexCData', repmat(clr,N,1), 'FaceColor', 'flat', ...
    'FaceAlpha', Alpha);
shading interp;
lighting phong;
material('metal');
view(3);

end