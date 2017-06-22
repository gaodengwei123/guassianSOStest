function [Pt,Lambda] = kfilter2(A, BK, Qt, C, R, Pt_1, Lambda0, Tkf, Obs)
% kalman filter map
In = eye(size(A));
At = In +Tkf*A +Tkf^2/2*A^2 +Tkf^3/6*A^3 +Tkf^4/24*A^4 +Tkf^5/120*A^5;
AK = At+BK*Tkf;
Qt = Qt*Tkf+(At*Qt+(At*Qt)')*Tkf^2/2+(At*(At*Qt+(At*Qt)')+At*(At*Qt+(At*Qt)')')*Tkf^3*6;

Ptbar = At*Pt_1*At'+Qt;
if Obs == 0 % no measurments  only prediecting
    Pt = Ptbar;
    Lt = zeros(size(C'));
else % prediecting
    St = C*Ptbar*C'+R;
    Lt = Ptbar*C'*(St^-1);
    Pt = Ptbar - Lt*C*Ptbar;
end
Pt = chol(Pt)'*chol(Pt);
Lambda = AK*Lambda0*AK'+Lt*C*Ptbar;
Lambda = chol(Lambda)'*chol(Lambda);
% Qk = (-BK*Tkf)*Pt*(-BK*Tkf)'+Qt;
% Lambda = AK*Lambda0*AK'+Qk;
end










