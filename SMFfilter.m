function [Pt,Lambda] = SMFfilter(A, BK, Qt, C, R, Pt_1, Lambda0, Tkf, Obs)
% kalman filter map
n = size(A);
In = eye(n);
At = In +Tkf*A +Tkf^2/2*A^2 +Tkf^3/6*A^3 +Tkf^4/24*A^4 +Tkf^5/120*A^5;
AK = At+BK*Tkf;
Qt = Qt*Tkf+(At*Qt+(At*Qt)')*Tkf^2/2+(At*(At*Qt+(At*Qt)')+At*(At*Qt+(At*Qt)')')*Tkf^3*6;

beta = sqrt(trace(Qt))/(sqrt(trace(Qt))+sqrt(trace(At*Pt_1*At')));
% Pbar = At*Pt_1/(1-beta)*At'+Qt/beta;
[U,D] = UDFactor(Pt_1);
[UQ,DQ] = UDFactor(Qt);
S = [At*U UQ];
Dtilde = blkdiag(D/(1-beta), DQ/beta);
Pbar = S*Dtilde*S';
%     [Ubar,Dbar] = UDFactor(Pbar);
%     s = Dbar*Ubar'*C'/(1-rho);
%     c = C*Ubar*s+R/rho;
if Obs == 0 % no measurments  only prediect
    Pt = Pbar;
    Pt = chol(Pt)'*chol(Pt);
else % prediecting
    rho = sqrt(max(svd(R)))/(sqrt(max(svd(C*Pbar*C')))+sqrt(max(svd(R))));
    Wk = C*Pbar/(1-rho)*C'+R/rho;
    Phat = Pbar/(1-rho)-Pbar/(1-rho)*C'*(Wk)^-1*C*Pbar/(1-rho);
    delta = 1;
    Pt = delta*Phat;
    Pt = chol(Pt)'*chol(Pt);
%     [Uhat,Dhat] = UDFactor(Phat);
%     s = Dhat*Uhat'*C'/(1-rho);
%     c = C*Uhat*s+R/rho;
%     Pt = Uhat*(Dhat/(1-rho)-1/c*s*s')*Uhat';
end

% SMF filter
QKt = eps*eye(n)-BK*Tkf*Pt;
QKt = chol(QKt)'*chol(QKt);
beta1 = sqrt(trace(Qt))/(sqrt(trace(QKt))+sqrt(trace(Qt)));
Qmao = QKt/(1-beta1)+Qt/beta1;

Lnew = AK*Lambda0*AK';
beta2 = sqrt(trace(Qmao))/(sqrt(trace(Lnew))+sqrt(trace(Qmao)));
[Ut,Dt] = UDFactor(Lambda0);
[UmQ,DmQ] = UDFactor(Qmao);
S2 = [AK*Ut UmQ];
Dtilde = blkdiag(Dt/(1-beta2), DmQ/beta2);
Lambda = double(S2*Dtilde*S2');
Lambda = chol(Lambda)'*chol(Lambda);

% can not be smaller than zeros
if min(eig(Lambda))<0||min(eig(Pt))<0
    error('must passive')
end

end










