% /*! @plot_myFunnel.m
% *************************************************************************
% <PRE>
% file.name       : plot_myFunnel.m
% related files   :
% function&ablity :
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/8/1      1.00                     build
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
function h = plot_myFunnel(sys,obj,options)
if (nargin<2) options=struct(); end
if ~isfield(options,'plotdims') options.plotdims = [1;2]; end
if ~isfield(options,'x0') options.x0 = zeros(length(obj.x),1); end
if (~isfield(options,'inclusion')) options.inclusion = 'slice'; end
% if (~isfield(options,'color')) options.color=.7*[1 1 1]; end
if (~isfield(options,'color')) options.color= [0.0710    0.3060    0.8750]; end%[0 0 1]
if (~isfield(options,'tol')) options.tol = .01; end
hold on;
view(0,90);

h=[];
ts = obj.getbreak;
x0 = options.x0;

for i=length(ts)-1:-1:1
    if strcmp(options.inclusion,'slice')
        options.x0 = x0.eval(ts(i));
        xfun0 = getLevelSet(obj.x,obj.getPoly(ts(i)),options);
        options.x0 = x0.eval(ts(i+1));
        xfun1 = getLevelSet(obj.x,obj.getPoly(ts(i+1)),options);
    elseif strcmp(options.inclusion,'projection')
        options.x0 = x0.eval(ts(i));
        xfun0 = getProjection(obj.x,obj.getPoly(ts(i)),options);
        options.x0 = x0.eval(ts(i+1));
        xfun1 = getProjection(obj.x,obj.getPoly(ts(i+1)),options);
    else
        error(['Unknown inclusion method: ' options.inclusion]);
    end
    if (i==length(ts)-1)
        % draw level-set at the end
        %         h1 = plot3(xfun1(1,:),xfun1(2,:),repmat(.1,1,size(xfun1,2)),'Color',[1 1 1]);
    end
    % i
    x = [xfun0,xfun1(:,end:-1:1)];
    k = convhull(x(1,:),x(2,:));
    AAA = colormap(gray);
    AAA(40:end,:) =[];
    AAA(1:10,:) =[];
    j = fix(i/(length(ts)-1)*length(AAA(:,1)));
    h=[h;fill3(x(1,k),x(2,k),repmat(0,1,length(k)),AAA(j,:),'LineStyle','none')];
    
    % plot convhulls
    %     h=[h;plot3(x(1,k),x(2,k),repmat(-.1,1,length(k)),'Color',[0 0 0])];
    
    
    %     h=[h;plot3(xfun0(1,:),xfun0(2,:),repmat(.1,1,size(xfun0,2)),'Color',[.5 .5 .5])];
end
% set(h,'alpha',0.5);
% draw level-set at the beginning
% h=[h;plot3(xfun0(1,:),xfun0(2,:),repmat(.1,1,size(xfun0,2)),'Color',[1 1 1])];
% plot trajectory
h0=plot(sys.PlotObj.solutionPlot.state(: ,1), sys.PlotObj.solutionPlot.state(: ,2));
% set(h0,'Color',[1 0 0]);
set(h0,'Color',[0 1 0]);

% h = [h;h1];
alpha(h,.5)
% global hPdiff
% hPdiff;
uistack(h,'bottom')
% if ellipsoid use plotEllipsoid.m
end





