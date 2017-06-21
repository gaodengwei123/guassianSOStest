function [output, gpopsHistory] = GaoMap_gpopsSolveCar(sys,Time,point)

global CONSTANS
CONSTANS = sys.INPUTS;
CONSTANS.PRMflag = 1;

% xmin = min(point(:,1));
% xmax = max(point(:,1));
% ymin = min(point(:,2));
% ymax = max(point(:,2));
phimin = -pi;
phimax = pi;
dphimin = -pi;
dphimax = pi;

param_min = [];
param_max = [];
path_min = [];
path_max = [];

event_min = [];
event_max = [];
duration_min = [];
duration_max = [];
iphaseSize = size(point,1)-1;
for iphase = 1:iphaseSize

    x0 = point(iphase,:);
    xT = point(iphase+1,:);
    CONSTANS.state_goal_phase(iphase,:) = xT;
    t0 = Time(iphase);
    tf = Time(iphase+1)+iphaseSize;    % time upper bound must be extend to ensure time cost in interval
    
    limits(iphase).time.min = [t0 tf];
    limits(iphase).time.max = [t0 tf];
    limits(iphase).control.min = -pi;
    limits(iphase).control.max = pi;
    
    limits(iphase).state.min (1 ,:) = [ x0(1) x0(1)-10 xT(1)];
    limits(iphase).state.max (1 ,:) = [ x0(1) xT(1)+10 xT(1)];
    limits(iphase).state.min (2 ,:) = [ x0(2) x0(2)-10 xT(2)];
    limits(iphase).state.max (2 ,:) = [ x0(2) xT(2)+10 xT(2)];
    limits(iphase).state.min (3 ,:) = [ x0(3) phimin xT(3)];
    limits(iphase).state.max (3 ,:) = [ x0(3) phimax xT(3)];
    limits(iphase).state.min (4 ,:) = [ x0(4) dphimin xT(4)];
    limits(iphase).state.max (4 ,:) = [ x0(4) dphimax xT(4)];
    
    limits(iphase).parameter.min = param_min ;
    limits(iphase).parameter.max = param_max ;
    limits(iphase).path.min = path_min ;
    limits(iphase).path.max = path_max ;
    limits(iphase).event .min = event_min ;
    limits(iphase).event.max = event_max ;
    limits(iphase).duration.min = duration_min ;
    limits(iphase).duration.max = duration_max ;
    
    guess(iphase).time = [t0; tf];
    guess(iphase).state(: ,1) = [x0(1); xT(1)];
    guess(iphase).state(: ,2) = [x0(2); xT(2)];
    guess(iphase).state(: ,3) = [x0(3); xT(3)];
    guess(iphase).state(: ,4) = [x0(4); xT(4)];
    guess(iphase).control = [0;0];
    guess(iphase).parameter = [];
    
%     if iphase<iphaseSize-1;
%         ipair = iphase; % pair of phases to link
%         linkages(ipair).left.phase = iphase;
%         linkages(ipair).right.phase = iphase+1;
%         linkages(ipair).min = [-0.1; -0.1; -0.1; -0.1];
%         linkages(ipair).max = [0.1; 0.1; 0.1; 0.1];
%     end
end

setup.name = 'CarControl-Problem';
setup.funcs.cost = 'Car_controlCost_Phases';
setup.funcs.dae = 'Car_controlDae';
% setup.funcs.link = 'CarLink';
setup.limits = limits ;
setup.guess = guess ;
setup.derivatives = 'finite-difference';
setup.direction = 'increasing';
setup.linkages = [];%linkages;
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-3;
setup.mesh.maxiterations = 10;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;
[output, gpopsHistory] = gpops(setup);


end




