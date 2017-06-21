classdef StochasticVandpol<GaoSystem
    properties
        w
        v
    end
    methods
        function obj = StochasticVandpol()
            mark = 0;
            obj = obj@GaoSystem(mark,2,1);
            w = msspoly('w',2);
            v = msspoly('v',2);
            obj.w = w;
            obj.v = v;
        end
        function [xdot,df]= dynamics(obj,t,x,u)
            w = obj.w;
            xdot = [-x(2)+w(1);
                -x(2)*(1-x(1)^2)+x(1)+u+w(2)];
            if (nargout>1)
                df = [0 -1; 1+2*x(1)*x(2) -1+x(1)^2];
            end
        end
        
        % systme output
        function y = output(obj,t,x,u)
            y = x;
        end
        
    end
end