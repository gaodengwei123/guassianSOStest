% /*! @plotROA.m
% *************************************************************************
% <PRE>
% file.name       : plotROA.m
% related files   :
% function&ablity :
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/6/15     1.00                     build
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
function h = plotROA(obj,options)
if (nargin<2) options=struct(); end
if ~isfield(options,'plotdims') options.plotdims = [1;2]; end
if ~isfield(options,'x0') options.x0 = zeros(length(obj.x),1); end
if (~isfield(options,'inclusion')) options.inclusion = 'slice'; end
if (~isfield(options,'color')) options.color=.7*[1 1 1]; end
if (~isfield(options,'tol')) options.tol = .01; end

hold on;
view(0,90);
% TODO: Here we split between projection and slice.
if strcmp(options.inclusion,'slice')
    x = getLevelSet(obj.x,obj.getPoly,options);
elseif strcmp(options.inclusion,'projection')
    x=getProjection(obj.x,obj.getPoly,options.x0,options.plotdims,options);
else
    error(['Unknown inclusion method: ' options.inclusion]);
end
h=fill3(x(1,:),x(2,:),zeros(1,size(x,2)),options.color,'LineStyle','-','LineWidth',2);

end