% /*! @GMPLdatabase.m
% *************************************************************************
% <PRE>
% file.name       : GMPLdatabase.m
% related files   :
% function&ablity : Guessian estimate
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/08/8     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output:

% *************************************************************************
function V = GMPLdatabase(sys,V)
state = sys.FunTraj;
%% GPtrainng
x = [];
y = [];
% for jj=1:4
%     switch jj
%         case 1
%             load('pathdata1.mat','xx','Vtrim')
%         case 2
%             load('pathdata2.mat','xx','Vtrim')
%         case 3
%             load('pathdata3.mat','xx','Vtrim')
%         case 4
%             load('pathdata4.mat','xx','Vtrim')
%     end
%             state = xx;
%             V = Vtrim;
%%
n = size(state.eval(0),1);
ts = V.getbreak;
% ts = linspace(ts(1),ts(end),1000);
ts=ts(:);
N = length(ts);

for i=1:N
    stateK(:,i) = state.eval(ts(i)); 
    Vbase(:,i) = reshape(ppval(V.Spp,ts(i)),n^2,1);
end
% database time
% for i = 1:100
%     testSample = i*1;
%     testTime(i) = V.getbreak(testSample);
%     Vbase(:,i) = reshape(V.S{testSample},n^2,1);
%     real_rho(i) = full(V.gettime(statest(:,testSample),testSample).Vvalue);
%     %% correct rho
%     if real_rho(i)<0
%         fprintf(1,'rho = %d at i = %d\n',real_rho(i),i);
%         real_rho(i)=-real_rho(i);
%     end
% %%
% stateK(:,i) = state(testTime(i));
% end
%% ========guessian estimate============
load('datafunnel0122.mat','rho')
% fitting rho is just test=========
Fitrho = spline(V.getbreak,rho);
rho = ppval(Fitrho,ts);
%============================
sys = sys.timecalculation(100);% calcluate the time points
taus = sys.breaks;
sys = sys.tv_poly_linearize; % linearize the sys
A = sys.A;
B = sys.B;
Qf = 1*sys.INPUTS.Q(0); % this is just a initial guess of the target ROA
[tv,Ss] = tv_lqr_riccati(sys.Maxinterval,A,B,sys.INPUTS.Q,sys.INPUTS.R,Qf);
Spp = spline(tv,Ss);
V = V_function(sys.p_t,sys.p_x,Spp,taus);

for i=1:length(taus)
    Vtest(:,i) = reshape(ppval(V.Spp,taus(i)),n^2,1);
end
%%
x = [x,stateK];           % state
y = [y,double(rho)];      % real data base
t = [state.eval(taus)];   % estimate

% x = [x,Vbase];           % state
% y = [y,double(rho)];     % real data base
% t = [Vtest];             % estimate
%%
GMPLprogama
new_rho = MI_Cand_gp;
figure
hold on
plot(taus,new_rho,'ko')
plot(ts,rho,'r')
% plot(new_rho,'o')
plot(taus,new_rho+10*sqrt(gp_v),'b--')
plot(taus,new_rho-10*sqrt(gp_v),'b--')
hullupper = new_rho+10*sqrt(gp_v);
hulldown = new_rho-10*sqrt(gp_v);
fill([taus;flip(taus)],[hullupper;flip(hulldown)],'b','FaceAlpha',0.2,'EdgeColor','none');

xlabel('sample points');
ylabel('\rho'); 

% =============end==================
for i = 1:size(new_rho,1)
    shrinkRho = new_rho./real_rho';
end
foh(taus,1./rho')
V = V.updateV(shrinkRho);
end
