% GMPL
clear all
clc
load('par_state.mat')
load('control_gain.mat')
load('g.mat')
load('par_state2.mat')
load('control_gain2.mat')
load('g2.mat')
% Initialize Gaussian Process
tic
t = par_state;        % estimate
likfunc = @likGauss;
hyp.lik = log(0.0); % input is signal noise
inffunc = @infExact;
meanfunc = @meanConst; 
covfunc = {@covMaternard,3};
hyp.cov = log([1 1 0.01 1 1 0.01 1]');
x = par_state2(:,1:50);            % state
y = g2(:,1:50);                     % real data base
hyp.mean = mean(y');
hyp = minimize(hyp,'gp', -100, inffunc, meanfunc, covfunc, likfunc, x', y');    % par training
[MI_Cand_gp, gp_v, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x', y', t');
toc
figure
plot(fmu(1:30),'r')
hold on
plot(g)
plot(fmu(1:30)+sqrt(fs2(1:30)),'y')
plot(fmu(1:30)-sqrt(fs2(1:30)),'y')