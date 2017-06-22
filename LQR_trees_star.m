classdef LQR_trees_star
    %% ==============dynmaic system======================
    properties
        A;B;Q;R;sys;
        Xx;% msspoly
        limits; % control boundary
        upperCost = 1;
    end
    %% start method
    methods 
        function obj = LQR_trees_star(sys)
            obj.A = sys.A;
            obj.B = sys.B;
            obj.Q = sys.Q;
            obj.R = sys.R;
            obj.sys = sys.f;
            obj.limits = sys.uppercontrol;  % control less than the boundary will be accept
            obj.Xx =  msspoly('x',size(sys.Q,1));
        end
    end
    methods (Static=true)
        function G = LQR_trees_search(sys,Num,range,xT,x0)
            obj = LQR_trees_star(sys);
            Randomp = haltonset(2);
            Randomp = scramble(Randomp,'RR2');
            X0 = net(Randomp,Num+1);
            %% initila V and E======
            Ax = obj.A(x0);
            Bx = obj.B(x0);
            [~,S,~] = lqr(Ax,Bx,obj.Q,obj.R);
            V = LQR_Nodes(x0,S);
            % plot ROA of xT
%             K = 20;
%             th=linspace(0,2*pi,K);
%             X = [sin(th);cos(th)];
%             S = S/0.1;
%             y = repmat([pi;0],1,K) + S^(-1/2)*X;
%             plot(y(1,:),y(2,:),'r')
            drawnow
            %% ===================
            for i = 2:Num
                i
                xrand = [range(1)+(range(3)-range(1))*X0(i,1);range(2)+(range(4)-range(2))*X0(i,2)];
                if mod(i,100) == 0
                    xrand = [pi;0];
                end
                [xnearest,xRand,flag,K] = LQRNearest(obj,V,xrand);
                if flag==1   % over control limits
                    continue
                end
                xnew = LQRSteer(obj,xnearest,xRand,K);
%                 plot(xnew.node(1),xnew.node(2),'o')
                if collisionFree(obj,xnew)
                    Xnear = LQRNear(obj,V,xnew);
                    if isempty(Xnear)
                        Xnear = xnearest;
                    end
                    xnew = chooseparnet(obj,Xnear,xnew);
                    %                     plot(xnew);
                    % Append to Tree
                    V = [V xnew];
                    
                    V = rewire(obj,V,Xnear,xnew);
                end
                if norm(xnew - xT)<0.0  % is near to target
                    break
                end
            end
            G = V;   % find a link x0 and xT
            plotpath(V(end))
        end
    end
    %%
    methods(Access = private)
        function [xnearest,xrand,flag,K]= LQRNearest(obj,V,xrand)
            % V: vectix, xrand: random state
            % linearize sys at xrand
            Ax = obj.A(xrand);
            Bx = obj.B(xrand);
            [K,S,~] = lqr(Ax,Bx,obj.Q,obj.R);
%             S = eye(length(xrand));
            % find a near V
            Num = length(V);
            X = repmat(xrand,1,Num);
            Estate = [V(:).node]-X;
            Sa = Estate'*S;
            LQRdist = sum(Sa.*Estate',2);
            [~,Index]= min(LQRdist);
            xnearest = V(Index);
            if abs(K*(xnearest - xrand))>obj.limits
                flag =1;
                return
            end
%             xrand = LQR_Nodes(xrand,S);
            flag = 0;
        end
        
        function Xnear = LQRNear(obj,V,xnew)
            Num = length(V);
            X = repmat(xnew.node,1,Num);
            Estate = V-X;
            Sa = Estate'*xnew.S;
            LQRdist = sum(Sa.*Estate',2);
%             gamma   = 1;  d   = 2;
%             ner = gamma*( log(Num+1)/Num )^(1/d);

            Xnear = V(LQRdist<=0.1);
        end
        
        function xnew = LQRSteer(obj,xnear,xrand,K)
            % find a node xnew to xnear with a certain engergy
            %             A0 = obj.A(xrand);
            %             B0 = obj.B(xrand);
            %             E = Ax-A0;
            %             D1 = 0.9*(-S0*Ax-Ax'*S0);
            %             Ac = (A0-B0*(obj.R)^-1*B0'*S0)';
            %             FF = -S0*Ax-Ax'*S0-D1;
            %             S1 = lyap(Ac,FF);
            %             S=S0+S1;
            
            %             S = xnear.S;
            %             S0 = S/obj.upperCost;
            %             xrand = rand(xnear.dim);
            %             xR = ones(size(xrand))-2*rand(size(xrand));
            %             X = xR/norm(xR);
            %             xnew = xnear + S0^(-1/2)*X;
            odefun = @(t,x)obj.sys.dynamics(t,x,-K*(x-xrand));
            sol = ode45(odefun,[0, 0.1],[xnear.node]);
            xnew = deval(sol,0.1);
%             xdot = obj.sys.dynamics(0,xnear.node,-K*(xnear-xrand));
%             xnew = xnear+xdot*0.1;
            Ax = obj.A(xnew);
            Bx = obj.B(xnew);
            [~,Snew,~] = lqr(Ax,Bx,obj.Q,obj.R);
            
            CostNow = Cost(xnear)+(xnear-xnew)'*Snew*(xnear-xnew);
            xnew = LQR_Nodes(xnew,Snew,CostNow);
            %             xnew = addparent(xnew,xnear);
        end
        
        function xnew = chooseparnet(obj,Xnear,xnew)
            % find a best parent for traj
            
            Num = length(Xnear);
            X = repmat(xnew,1,Num);
            Estate = Xnear-X;
            Sa = Estate'*xnew.S;
            LQRdist = sum(Sa.*Estate',2);
            [CostNow, index]= min(Cost(Xnear) + LQRdist');
            xnear = Xnear(index);
            xnew = addparent(xnew,xnear,CostNow);
        end
        
        function flag = collisionFree(obj,traj)
            % check the collision
            flag = 1;
        end
        
        function V = rewire(obj,V,Xnear,xnew)
            for i=1:length(Xnear)   % number of xnew neighbor
                xnear = Xnear(i);
                %                 if xnear.parents~=xnew  % not rewire parents
                Estate = xnew - xnear;
                LQRdist = Estate'*xnear.S*Estate;
                
                CostNow = Cost(xnew) + LQRdist;
                if  CostNow <Cost(xnear)
                    V = V.update(xnear,xnew,CostNow);
                end
                %                 end
            end
        end
        
        
        %%
        
    end
end