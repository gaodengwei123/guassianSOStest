function [Pk,Kk,Fikk_1] = kfilter(Ft, Xk_1, Qt, Hk, Zk, Rk, Pk_1, Tkf, n, Obs)
% kalman filter
% Ft��state matrix��
% G:  system noise matrix
% Qt��noise variance matrix
% Hk��measurement matrix��
% Rk��measurement variance matrix��
% X��state��
% P��is the predicted state covariance��
% Tkf��period��
% n:dimension.
    In = eye(n);
    Fikk_1 = In +Tkf*Ft +Tkf^2/2*Ft^2 +Tkf^3/6*Ft^3 +Tkf^4/24*Ft^4 +Tkf^5/120*Ft^5;
    Qk_1 = Qt*Tkf+(Ft*Qt+(Ft*Qt)')*Tkf^2/2+(Ft*(Ft*Qt+(Ft*Qt)')+Ft*(Ft*Qt+(Ft*Qt)')')*Tkf^3*6;
    
    if Obs == 0 % no measurments  only prediecting
        Pkk_1=Fikk_1*Pk_1*Fikk_1'+Qk_1;
        Pk=Pkk_1;
        Kk = zeros(size(Hk'));
%         Xkk_1=Fikk_1*Xk_1;
%         Xk=Xkk_1;

    else % prediecting
        Pkk_1=Fikk_1*Pk_1*Fikk_1'+Qk_1;
        Kk=Pkk_1*Hk'*(Hk*Pkk_1*Hk'+Rk)^-1;
        Pk=(In-Kk*Hk)*Pkk_1*(In-Kk*Hk)'+Kk*Rk*Kk';

%         Xkk_1=Fikk_1*Xk_1;
%         Xk=Xkk_1+Kk*(Zk-Hk*Xkk_1);
    end
end