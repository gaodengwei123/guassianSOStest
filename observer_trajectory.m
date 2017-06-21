clear
clc
close all
load('datavehicle20170123.mat')
figure
sys = sysploy;
field = ObstacleField();
field.range = [0 0 50 50];
field = field.addsensors([20 25],8);
show(field);
sys = sys.setfield(field);
xtraj = sys.FunTraj;
utraj = sys.FunContr;
taus = sys.breaks;

INPUTS = sys.INPUTS;
% integral the state with feedback
n = sys.getNumStates;
m = sys.getNumInput;
obsevern = length(sys.INPUTS.noise_v);
N = sys.timeNum;

% dynamic integral
x0 = xtraj.eval(xtraj.pp.breaks);
plot(x0(1,:),x0(2,:),'r','lineWidth',4)

rho = 1*[3 3 0.2 0.2]';

th = 0:pi/50:2*pi;
xunit = rho(1) * cos(th) + x0(1);
yunit = rho(2) * sin(th) + x0(2);
h = plot(xunit, yunit);

INPUTS.noise_w = 10*INPUTS.noise_w;
INPUTS.Qk = diag(INPUTS.noise_w.^2);
sys.INPUTS.noise_w = INPUTS.noise_w;
sys.INPUTS.Qk = INPUTS.Qk;

plotDynNum = 1;    % Monte Carlo simulation
SampleNum = 100;
ts = linspace(taus(1),taus(end),SampleNum);  % set sample time

for i = 1:plotDynNum
    
    InitialState = xtraj.eval(0)+rho.*(2*rand(n,1)-1);
    
    FeedbackU = @(t,x)(utraj.eval(t)+sys.K(t)*(x-xtraj.eval(t)));
    odefun = @(t,x)sys.dynamics(t,x,FeedbackU(t,x));

    state(:,1) = InitialState;      % estimate
    realstate(:,1) = InitialState;  % real
    
    for j = 1:SampleNum-1
        
        Tkf = ts(j+1)-ts(j);
        [~,XS] = ode45(odefun,[ts(j), ts(j+1)],state(:,j));
        state(:,j+1) = XS(end,:)';   
        % observer trajectory
        observer = GaoObserver(odefun,sys.INPUTS.Rk);
        [z,H,Obs] = observer.measure(ts(j),realstate(:,j),FeedbackU(ts(j),state(:,j)),INPUTS.noise_v.*(randn(obsevern,1)));
        if Obs
            z = state(1:3,j+1)+INPUTS.noise_v.*(randn(obsevern,1));
            [L, P] = GaoObserver.filter(obj, A, Qt, C, R, Pt_1, Tkf, Obs);
            realstate(:,j+1) = realstate(:,j) + L*(z - H*state(:,j+1));
        else
            realstate(:,j+1) = state(:,j) + 
        end
        
        % real trajectory
        odefun = @(t,x)Stochsys.dynamics(t,x,FeedbackU(t,state(:,j+1)), INPUTS.noise_w.*(randn(n,1)));
        [t,XS] = ode45(odefun,[ts(j), ts(j+1)],realstate(:,j));
        realstate(:,j+1) = XS(end,:)';
        
    end
    plot(state(1,:),state(2,:),'k')
    plot(realstate(1,:),realstate(2,:),'g')
    drawnow
end


% Tkf = ts(2)-ts(1);
% Pt_1 = INPUTS.Qk;
% pdim =  [1 0 0 0;0 1 0 0]';
% for i=1:100
%     A = sys.A(ts(i));
%     Qt = Stochsys.Wc;
%     % kalman filter
%     In = eye(size(A));
%     At = In +Tkf*A +Tkf^2/2*A^2 +Tkf^3/6*A^3 +Tkf^4/24*A^4 +Tkf^5/120*A^5;
%     Qt = Qt*Tkf+(At*Qt+(At*Qt)')*Tkf^2/2+(At*(At*Qt+(At*Qt)')+At*(At*Qt+(At*Qt)')')*Tkf^3*6;
%     Ptbar = At*Pt_1*At'+Qt;
%     errorEllipose(i) = error_ellipse(state(1:2,i),Ptbar(1:2,1:2));
%     plot(errorEllipose(i),'y')
% end