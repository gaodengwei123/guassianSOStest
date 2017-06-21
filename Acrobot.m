classdef Acrobot<GaoSystem
    properties
        % parameters from Spong95 (except inertias are now relative to the
        % joints)
        % axis)
        l1 = 1; l2 = 2;
        m1 = 1; m2 = 1;
        g = 9.81;
        b1=.1;  b2=.1;
        %    b1=0; b2=0;
        lc1 = .5; lc2 = 1;
        Ic1 = .083;  Ic2 = .33;
        
        xG
        uG
    end
    methods
        function obj = Acrobot()
            mark = 0;
            obj = obj@GaoSystem(mark,4,1);
            obj.xG = [pi;0;0;0];
            obj.uG = 0;
        end
        
        function [xdot,df]= dynamics(obj,t,x,u)
            q = x(1:2);
            v = x(3:4);
            [H,C,B] = manipulatorDynamics(obj,q,v);
            
            Hinv = inv(H);
            if (obj.getNumInput>0) tau=B*u - C; else tau=-C; end
            % tau = tau + computeConstraintForce(obj,q,v,H,tau,Hinv);
            
            vdot = Hinv*tau;
            % note that I used to do this (instead of calling inv(H)):
            %   vdot = H\tau
            % but I already have and use Hinv, so use it again here
            
            xdot = [vToqdot(obj,q)*v; vdot];
            if (nargout>1)
                df = dynamicsGradients(obj, t, x, u);
            end
            % quiver(x(1),x(2),0,u,0.5,'y')
        end
        
        % systme output
        function y = output(obj,t,x,u)
            y = x;
            nX = obj.getNumStates;
            nU = obj.getNumInput;
            [~,df] = geval(@obj.dynamics,t0,x0,u0);
            obj.A = df(:,1+(1:nX));
            obj.B = df(:,nX+1+(1:nU));
        end
        
         function A = dynamicA(obj, t, x, u)
            [~,df,~,~] = obj.dynamics(t,x,u);
            nX = obj.getNumStates;
            A = full(df(:,1+(1:nX)));
        end
        function B = dynamicB(obj, t, x, u)
            [~,df,~,~] = obj.dynamics(t,x,u);
            nX = obj.getNumStates;
            nU = obj.getNumInput;
            B = full(df(:,nX+1+(1:nU)));
        end
    end
    
    methods (Access = private)
        function [H,C,B] = manipulatorDynamics(obj,q,qd)
            % keep it readable:
            m1=obj.m1; m2=obj.m2; l1=obj.l1; g=obj.g; lc1=obj.lc1; lc2=obj.lc2; b1=obj.b1; b2=obj.b2;
            I1 = obj.Ic1 + obj.m1*obj.lc1^2; I2 = obj.Ic2 + obj.m2*obj.lc2^2;
            m2l1lc2 = m2*l1*lc2;  % occurs often!
            
            c = cos(q(1:2,:));  s = sin(q(1:2,:));  s12 = sin(q(1,:)+q(2,:));
            
            h12 = I2 + m2l1lc2*c(2);
            H = [ I1 + I2 + m2*l1^2 + 2*m2l1lc2*c(2), h12; h12, I2 ];
            
            C = [ -2*m2l1lc2*s(2)*qd(2), -m2l1lc2*s(2)*qd(2); m2l1lc2*s(2)*qd(1), 0 ];
            G = g*[ m1*lc1*s(1) + m2*(l1*s(1)+lc2*s12); m2*lc2*s12 ];
            
            % accumate total C and add a damping term:
            C = C*qd + G + [b1;b2].*qd;
            
            B = [0; 1];
        end
        
        function [VqInv,dVqInv] = vToqdot(obj, q)
            % defines the linear map qdot = Vqinv * v
            % default relationship is that v = qdot
            VqInv = eye(length(q));
            dVqInv = zeros(numel(VqInv), 2);
        end
        
        function df = dynamicsGradients(a1, a2, a3, a4)
            % Symbol table:
            Ic1=a1.Ic1;
            Ic2=a1.Ic2;
            a3_1=a3(1);
            a3_2=a3(2);
            a3_3=a3(3);
            a3_4=a3(4);
            a4_1=a4(1);
            b1=a1.b1;
            b2=a1.b2;
            g=a1.g;
            l1=a1.l1;
            lc1=a1.lc1;
            lc2=a1.lc2;
            m1=a1.m1;
            m2=a1.m2;
            % Compute Gradients:
            df = sparse(4,6);
            df(1,4) = 1;
            df(2,5) = 1;
            df(3,2) = -(g*l1*(lc2^2*m2^2*cos(a3_1) + 2*Ic2*m2*cos(a3_1) - lc2^2*m2^2*cos(a3_1 + 2*a3_2)) + g*lc1*(2*Ic2*m1*cos(a3_1) + 2*lc2^2*m1*m2*cos(a3_1)))/(2*Ic1*Ic2 + 2*l1^2*lc2^2*m2^2 + 2*Ic2*l1^2*m2 + 2*Ic2*lc1^2*m1 + 2*Ic1*lc2^2*m2 + 2*lc1^2*lc2^2*m1*m2 - 2*l1^2*lc2^2*m2^2*cos(a3_2)^2);
            df(3,3) = (lc2*m2*(Ic2 + lc2^2*m2)*(a3_4^2*l1*cos(a3_2) - g*cos(a3_1 + a3_2) + 2*a3_3*a3_4*l1*cos(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2) + (lc2*m2*(g*cos(a3_1 + a3_2) + a3_3^2*l1*cos(a3_2))*(Ic2 + lc2^2*m2 + l1*lc2*m2*cos(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2) - (l1*lc2*m2*sin(a3_2)*(a3_4*b2 - a4_1 + g*lc2*m2*sin(a3_1 + a3_2) + a3_3^2*l1*lc2*m2*sin(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2) + (2*l1^2*lc2^2*m2^2*cos(a3_2)*sin(a3_2)*(Ic2 + lc2^2*m2)*(a3_3*b1 + g*lc2*m2*sin(a3_1 + a3_2) + g*l1*m2*sin(a3_1) + g*lc1*m1*sin(a3_1) - a3_4^2*l1*lc2*m2*sin(a3_2) - 2*a3_3*a3_4*l1*lc2*m2*sin(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2)^2 - (2*l1^2*lc2^2*m2^2*cos(a3_2)*sin(a3_2)*(Ic2 + lc2^2*m2 + l1*lc2*m2*cos(a3_2))*(a3_4*b2 - a4_1 + g*lc2*m2*sin(a3_1 + a3_2) + a3_3^2*l1*lc2*m2*sin(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2)^2;
            df(3,4) = -(Ic2*b1 - sin(a3_2)*(2*a3_3*l1*lc2^3*m2^2 + 2*a3_4*l1*lc2^3*m2^2 + 2*Ic2*a3_3*l1*lc2*m2 + 2*Ic2*a3_4*l1*lc2*m2) + b1*lc2^2*m2 - a3_3*l1^2*lc2^2*m2^2*sin(2*a3_2))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + l1^2*lc2^2*m2^2*(sin(a3_2)^2 - 1) + lc1^2*lc2^2*m1*m2);
            df(3,5) = (Ic2*b2 + l1*(sin(a3_2)*(2*a3_3*lc2^3*m2^2 + 2*a3_4*lc2^3*m2^2 + 2*Ic2*a3_3*lc2*m2 + 2*Ic2*a3_4*lc2*m2) - b2*lc2*m2*(2*sin(a3_2/2)^2 - 1)) + b2*lc2^2*m2)/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + l1^2*lc2^2*m2^2*(sin(a3_2)^2 - 1) + lc1^2*lc2^2*m1*m2);
            df(3,6) = -(Ic2 + m2*(lc2^2 + l1*lc2*cos(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2);
            df(4,2) = (g*(l1*lc2^2*m2^2*cos(a3_1) - Ic1*lc2*m2*cos(a3_1 + a3_2) + Ic2*l1*m2*cos(a3_1) + Ic2*lc1*m1*cos(a3_1) - l1^2*lc2*m2^2*cos(a3_1 + a3_2) + lc1*lc2^2*m1*m2*cos(a3_1) - lc1^2*lc2*m1*m2*cos(a3_1 + a3_2)) + g*cos(a3_2)*(l1^2*lc2*m2^2*cos(a3_1) - l1*lc2^2*m2^2*cos(a3_1 + a3_2) + l1*lc1*lc2*m1*m2*cos(a3_1)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2);
            df(4,3) = (2*l1*lc2*m2*sin(a3_2)*(a3_4*b2 - a4_1 + g*lc2*m2*sin(a3_1 + a3_2) + a3_3^2*l1*lc2*m2*sin(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2) - ((Ic2 + lc2^2*m2 + l1*lc2*m2*cos(a3_2))*(a3_4^2*l1*lc2*m2*cos(a3_2) - g*lc2*m2*cos(a3_1 + a3_2) + 2*a3_3*a3_4*l1*lc2*m2*cos(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2) - (l1*lc2*m2*sin(a3_2)*(a3_3*b1 + g*(m2*(lc2*sin(a3_1 + a3_2) + l1*sin(a3_1)) + lc1*m1*sin(a3_1)) - a3_4^2*l1*lc2*m2*sin(a3_2) - 2*a3_3*a3_4*l1*lc2*m2*sin(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2) - ((g*lc2*m2*cos(a3_1 + a3_2) + a3_3^2*l1*lc2*m2*cos(a3_2))*(Ic1 + Ic2 + l1^2*m2 + lc1^2*m1 + lc2^2*m2 + 2*l1*lc2*m2*cos(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2) + (2*l1^2*lc2^2*m2^2*cos(a3_2)*sin(a3_2)*(a3_4*b2 - a4_1 + g*lc2*m2*sin(a3_1 + a3_2) + a3_3^2*l1*lc2*m2*sin(a3_2))*(Ic1 + Ic2 + l1^2*m2 + lc1^2*m1 + lc2^2*m2 + 2*l1*lc2*m2*cos(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2)^2 - (2*l1^2*lc2^2*m2^2*cos(a3_2)*sin(a3_2)*(Ic2 + lc2^2*m2 + l1*lc2*m2*cos(a3_2))*(a3_3*b1 + g*(m2*(lc2*sin(a3_1 + a3_2) + l1*sin(a3_1)) + lc1*m1*sin(a3_1)) - a3_4^2*l1*lc2*m2*sin(a3_2) - 2*a3_3*a3_4*l1*lc2*m2*sin(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2)^2;
            df(4,4) = ((b1 - 2*a3_4*l1*lc2*m2*sin(a3_2))*(Ic2 + lc2^2*m2 + l1*lc2*m2*cos(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2) - (2*a3_3*l1*lc2*m2*sin(a3_2)*(Ic1 + Ic2 + l1^2*m2 + lc1^2*m1 + lc2^2*m2 + 2*l1*lc2*m2*cos(a3_2)))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2);
            df(4,5) = -(Ic1*b2 + Ic2*b2 + cos(a3_2)*(2*b2*l1*lc2*m2 + 2*a3_3*l1^2*lc2^2*m2^2*sin(a3_2) + 2*a3_4*l1^2*lc2^2*m2^2*sin(a3_2)) + b2*l1^2*m2 + b2*lc1^2*m1 + b2*lc2^2*m2 + 2*a3_3*l1*lc2^3*m2^2*sin(a3_2) + 2*a3_4*l1*lc2^3*m2^2*sin(a3_2) + 2*Ic2*a3_3*l1*lc2*m2*sin(a3_2) + 2*Ic2*a3_4*l1*lc2*m2*sin(a3_2))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2);
            df(4,6) = (Ic1 + Ic2 + lc1^2*m1 + m2*(l1^2 + lc2^2) + 2*l1*lc2*m2*cos(a3_2))/(Ic1*Ic2 + l1^2*lc2^2*m2^2 + Ic2*l1^2*m2 + Ic2*lc1^2*m1 + Ic1*lc2^2*m2 + lc1^2*lc2^2*m1*m2 - l1^2*lc2^2*m2^2*cos(a3_2)^2);

        end
    end
    
end