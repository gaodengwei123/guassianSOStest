classdef polynomialSystemDynamic<Gaosystem
    properties
        Ac
        Bc
        xcdot0
        Ad
        Bd
        xdn0
        C
        D
        y0
    end
    methods
        function g = newdyn(t,sys,frame,ctf,dtf)
            x = sys.getStateFrame.getPoly; z=frame.getPoly;
            p_t = msspoly('t',1);
            [d,dc]=output(t,[],sys.p_x);
            dcdt = dc(:,1); 
            dcdx = dc(:,2:end);
            g = dcdt + dcdx*sys.getPolyDynamics(t);
            g = subss(g,x,d);
        end
        
        function [y,dy] = output(obj,t,x,u)
            if (isTI(obj))
                y=zeros(obj.num_y,1);
                if (obj.num_x) y=y+obj.C*x; end
                if (obj.num_u) y=y+obj.D*u; end
                if ~isempty(obj.y0) y=y+obj.y0; end
                if (nargout>1)
                    dy = [0*y,obj.C,obj.D];
                end
            else
                y=zeros(obj.num_y,1);
                if (obj.num_x)
                    y=y+eval(obj.C,t)*x;
                end
                if (obj.num_u)
                    y=y+eval(obj.D,t)*u;
                end
                if ~isempty(obj.y0)
                    y=y+eval(obj.y0,t);
                end
                if (nargout>1)
                    dy = zeros(obj.num_y,1+obj.num_x+obj.num_u);
                    if (obj.num_x)
                        dy(:,1:1+obj.num_x)=dy(:,1:1+obj.num_x) + [deriv(obj.C,t)*x, eval(obj.C,t)];
                    end
                    if (obj.num_u)
                        dy(:,[1,2+obj.num_x:end])=dy(:,[1,2+obj.num_x:end])+[deriv(obj.D,t)*u,eval(obj.D,t)];
                    end
                    if ~isempty(obj.y0)
                        dy(:,1)=dy(:,1)+deriv(obj.y0,t);
                    end
                end
            end
        end
        
    end
end