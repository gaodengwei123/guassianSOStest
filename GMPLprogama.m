
% Initialize Gaussian Process
tic
likfunc = @likGauss;
hyp.lik = log(0.0); % input is signal noise
inffunc = @infExact;
meanfunc = @meanConst; 
covfunc = {@covMaternard,3};
hyp.cov = log(ones(size(x,1)+1,1));
hyp.mean = mean(y');
hyp = minimize(hyp,'gp', -100, inffunc, meanfunc, covfunc, likfunc, x', y);    % par training
[MI_Cand_gp, gp_v, fmu, fs2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x', y, t');
toc
