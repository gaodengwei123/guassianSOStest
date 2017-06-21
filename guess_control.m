% /*! @guess_control.m
% *************************************************************************
% <PRE>
% file.name       : guess_control.m
% related files   :
% function&ablity : guess initial control for videoRay-IV
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2015/10/29    1.00                     build
% 2016/8/31     2.00                     extand to GPOPS initial guess
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input :
% output:control:
% *************************************************************************
function control_output = guess_control(x0,x1,INPUTS)
m = INPUTS.m;
Izz = INPUTS.Izz;

% set control parameter
Kp = 10*diag([1 1 0.1]);       % proportional gain
Kd = 0.1*diag([1 1 0.1]);     % differential gain

% transform matrix B->I
psi = x0(3);
J = [cos(psi)        -sin(psi)     0
    sin(psi)         cos(psi)      0
    0                0             1];

e = x0(1:3)-x1(1:3);
de = J*(x0(4:6)-x1(4:6));

Contr = -Kp*e-Kd*de;
% psid = atan2(Contr(2),Contr(1));

% control force and torque trans to body frame
control = J'*[Contr(1:2)
    Contr(3)];

% saturation
upForce = INPUTS.controlSat*INPUTS.m/2;
if INPUTS.underactuation == 1  % underactuation control system
    if abs(control(3)) > INPUTS.momentSat;
        control(3)= sign(control(3))*INPUTS.momentSat;
    end
    F1 = 0.5*(control(1)+control(3)/INPUTS.rG);
    F2 = 0.5*(control(1)-control(3)/INPUTS.rG);

    if abs(F1) > upForce
        F1= sign(F1)*upForce;
        F2 = F1-control(3)/INPUTS.rG;
    end
    if abs(F2) > upForce
        F2 = sign(F2)*upForce;
        F1 = F2+control(3)/INPUTS.rG;
    end
    control_output = [sign(control(1))*(F1+F2)/m
        0
        (F1-F2)*INPUTS.rG];
else
    if abs(control(1)) > upForce;
        control(1) = sign(control(1))*INPUTS.controlSat;
    end
    if abs(control(2)) > upForce;
        control(2) = sign(control(2))*INPUTS.controlSat;
    end
    if abs(control(3)) > INPUTS.momentSat;
        control(3) = sign(control(3))*INPUTS.momentSat;
    end
    control_output = [control(1:2)/m;control(3)];
end



