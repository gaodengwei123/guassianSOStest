classdef doublePendulumdyn<GaoSystem
    properties
        b = 0.1;
        g=9.8;
    end
    
    
    methods
        function obj = doublePendulumdyn()
            mark = 0;
            obj = obj@GaoSystem(mark,2,1);
        end
        function [xdot,df,d2f,d3f]= dynamics(obj,t,x,u)
            xdot = [x(2);
                    u-obj.b*x(2)-obj.g*sin(x(1))];
            if (nargout>1)
                df = [0 1; -obj.g*cos(x(1)) -obj.b];
            end
            if (nargout>2)
                d2f = [0 0; obj.g*sin(x(1)) 0];
            end
            if (nargout>3)
                d3f = [0 0; obj.g*cos(x(1)) 0];
            end
        end
        
        % systme output
        function y = output(obj,t,x,u)
            y = x;
        end
        

    end
end