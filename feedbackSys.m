classdef feedbackSys<GaoSystem
    properties
        sys
        xG
        uG
    end
    methods
        function obj = feedbackSys(sys,K,B)
            obj = obj@GaoSystem(sys.mark,sys.getNumStates,sys.getNumInput);
            obj.K = K;
            obj.B = B;
            obj.sys = sys;
            try
                obj.xG = sys.xG;
                obj.uG = sys.uG;
            catch
                obj.xG = zeros(sys.getNumStates);
                obj.uG = zeros(sys.getNumInput);
            end
        end
        function [xdot,df]= dynamics(obj,t,x,~)
            if isa(obj.K,'function_handle')
                u = -obj.K(x(1),x(2));
            else
                u = -obj.K*x;
            end
            xdot = obj.sys.dynamics(t,x,u);

            if (nargout>1)
                [~,df0] = obj.sys.dynamics(t,x,0);
                df = df0-diag(obj.Gain);
            end
        end
        
        % systme output
        function y = output(obj,t,x,u)
            y = x;
        end
 
    end
end