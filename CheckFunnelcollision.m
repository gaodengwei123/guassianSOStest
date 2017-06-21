% /*! @CheckFunnelcollision.m
% *************************************************************************
% <PRE>
% file.name       : CheckFunnelcollision.m
% related files   :
% function&ablity :
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/03/6     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output:
% *************************************************************************
function [Xadd,flag]= CheckFunnelcollision(x_new, field)
% shrink funnel with obs
flag = 0;
vert=[];
for i = 1:field.number_of_obstacles
    obstacle = field.obstacles{i};
    xv = obstacle(1,:)';
    yv = obstacle(2,:)';
    state = x_new.Ftuple.traj.eval(x_new.Ftuple.breaks);
    xq = state(1,:)';
    yq = state(2,:)';
    [in, on] = inpolygon(xq,yq,xv,yv);
    flag = sum(in)+sum(on);
    if flag  % trajectory is cllision return empty
        Xadd = [];
        return
    end
    trajxy = [xq yq];
    [~,d] = dsearchn(trajxy,obstacle');
    Index = find(d<2); % 2 is a funnel radius guess
    if ~isempty(Index)
        vert = [vert, obstacle(:,Index)];
    end
end


% shrink the Funnel based on obstacles
N = length(x_new.Ftuple.breaks);
ts = x_new.Ftuple.breaks;
V = x_new.Ftuple.V;
Traj = x_new.Ftuple.traj;
rho(1:N)=1;
% iterate over t
if ~isempty(vert)
    for i = fliplr(1:N-1)
        rho(i) = rho(i+1);
        x0 = Traj.eval(ts(i));
        x = [vert;repmat(x0(3:end),1,size(vert,2))];
        Vvert = [];
        for k = 1:size(vert,2)
            Vvert = [Vvert,V.eval(ts(i),x(:,k))];
        end
        if (min(Vvert)<rho(i))
            rho(i) = min(Vvert);
        end
    end
end
Xadd = x_new;
% funnel is too small we stop it here to decrease the computation
if rho(1)<0.1
    Xadd = [];
    flag = 1;
    return
end

if ~all(rho==1)   % shirnk the funnel
    Xadd.Ftuple.V = V.updateV(foh(V.getbreak,1./rho));
end
end








