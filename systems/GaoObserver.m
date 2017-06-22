classdef GaoObserver
    
    
    properties
        forward_model
        V
        P
    end
    
    methods
        function obj = GaoObserver(dynamic,output_covariance)
            obj.forward_model = dynamic;
            obj.V = chol(output_covariance+eps*eye(size(output_covariance,1)));
            obj.P = [];
        end
        
        function [xcdot, P] = dynamics(obj,t,x,u)
            % Continuous time observer update
            yhat = output(obj,t,x,u);
            % forecast
            xcdot = dynamics(obj.forward_model,t,x,u);
            % update
            [L, P]= filter(obj, A, Qt, H, R, P, Tkf, Obs);
            xcdot = xcdot + L*(y-yhat);
            
        end
               
        function [Lt, Pt]= filter(obj, A, Qt, C, R, Pt_1, Tkf, Obs)
            % kalman filter
            In = eye(size(A));
            At = In +Tkf*A +Tkf^2/2*A^2 +Tkf^3/6*A^3 +Tkf^4/24*A^4 +Tkf^5/120*A^5;
            Qt = Qt*Tkf+(At*Qt+(At*Qt)')*Tkf^2/2+(At*(At*Qt+(At*Qt)')+At*(At*Qt+(At*Qt)')')*Tkf^3*6;
            Ptbar = At*Pt_1*At'+Qt;
            if Obs
                St = C*Ptbar*C'+R;
                Lt = Ptbar*C'*(St^-1);
                Pt = Ptbar - Lt*C*Ptbar;
            else
                Pt = Ptbar;
                Lt = [];
            end
        end
        
        function [z,H,Obs] = measure(obj,t,x,u,v)
            [z,H,Obs] = measurement(obj.sys,t,x,u);
            if length(v)~=length(z)
                error('observer must match noise')
            end
            z = z + obj.V*v; % only one sensor in simuliation
        end
        
        function obj = mimoCascade(obj,Stochsys)
            % combine sys1 and sys2: output of sys1 is the input of sys2 
            % x1dot = A1x1 + B1u + xcdot01
            % x2dot = A2x2 + B2(C1x1+D1u+y01) + xcdot02
            % y = C2x2 + D2(C1x1+D1u + y01) + y02
            u = (utraj.eval(t)+sys.K(t)*(x-xtraj.eval(t)));
            Stochsys.sys.dynamic(obj,t,x,u);
            obj.dynamics = dynamics(obj,t,x,u);
        end
    end
    
end



if length(varargin) > 1
    % there is sensor here
    v = varargin{2};  % measurement error
    [z,H,Obs] = Observer(obj,t,x,u,v);
    % add an filter here
    
    % LQG control
    A = obj.sys.A(t);
    Qt = obj.Wc;
    R = obj.V;
    
    % EKF
    [L, P]= filter(obj, A, Qt, H, R, P, Tkf, Obs);
    if Obs % get measurement
        x = x + L*(z-H*x);
        u = obj.sys.FunContr.eval(t)+obj.sys.K(t)*(x-obj.sys.FunTraj.eval(t));
        xcdot = dynamics(obj.sys,t,x,u)+L*(z-H*x);
    else
        % forecast
        x = Observer(obj,t,x,u,obj.sys.INPUTS.noise_w);
        u = obj.sys.FunContr.eval(t)+obj.sys.K(t)*(x-obj.sys.FunTraj.eval(t));
        xcdot = dynamics(obj.sys,t,x,u) + obj.Wc*w;
    end
else