classdef StochasticSystem<GaoSystem
    % Implements
    %   xcdot = Ac*x + Bc*u + xcdot0
    %   y = C*x + D*u + y0
    properties
        sys
        Wc
        V
        P
        t0
    end
    
    
    methods
        function obj = StochasticSystem(sys,dynamics_covariance,output_covariance)
            obj = obj@GaoSystem(sys.mark, sys.getNumStates, sys.getNumInput);
            obj.sys = sys;   % dynamic system of robot
            obj.Wc = chol(dynamics_covariance+eps*eye(sys.getNumStates));
            obj.V = chol(output_covariance+eps*eye(size(output_covariance,1)));
        end
        
        function xcdot = dynamics(obj,t,x,u,w)
            % forecast without obsever: u0+K(x-x0)
            xcdot = dynamics(obj.sys,t,x,u)+ obj.Wc*w;
            if Obs == NearSensor(obj.INPUTS.field,x(1:2))
                u = (obj.sys.FunContr.eval(t)+obj.sys.K(t)*(x-obj.sys.FunTraj.eval(t)));
            end
        end
       
    end
end

