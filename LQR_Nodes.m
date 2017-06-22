% /*! @LQR_Nodes.m
% *************************************************************************
% <PRE>
% file.name       : LQR_Nodes.m
% related files   :
% function&ablity : solve with LQR_trees_star.m
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/5/4      1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:
%
% output:
% *************************************************************************

classdef LQR_Nodes
    
    properties
        node
        dim
        parents
        energy = 0.1; % default
        S
        cost = 0;% default
        pow = 1;% default: not in the tree
        h;      %¡¡plot handle
    end
    %% start method
    methods
        function obj = LQR_Nodes(a,varargin)
            if isa(a,'LQR_Nodes')
                obj.node = a.node;
                obj.dim = size(a.node);
            else
                obj.node = a;
                obj.dim = size(a);
            end
            if nargin > 1
                obj.S = varargin{1};
            else
                obj.S = [];
            end
            if nargin > 2
                obj.cost = varargin{2};
            end
            if nargin > 3
                obj.energy = varargin{3};
            end
        end
        
        function cost = Cost(obj)
            cost = [obj(:).cost];
        end
        
        function obj = addparent(obj,a,varargin)
            %             a: aprent; b:cost
            obj.parents = a;
            %             update cost
            if nargin>2
                obj.cost = varargin{1};
            end
            delete(obj.h);
            colorA = colormap(jet);
            i = 64-floor(obj.cost/100*64);
            if i<1, i=1; end
            obj.h = line([obj.node(1) a.node(1)],[obj.node(2) a.node(2)],'Color',colorA(i,:));
            drawnow
        end
        
        function obj = update(obj,x,x0,Cost)
            % find x in the tree
            obj(x.pow) = obj(x.pow).addparent(x0,Cost);
        end
        
        function h = plot(a)
            h = plot(a.node(1),a.node(2),'o');
            
        end
        
        function h = plotpath(a)
            x = a;
            while ~isempty(x.parents)
                h = line([x.node(1),x.parents.node(1)],[x.node(2),x.parents.node(2)],'Color',[1 0 0],'LineWidth',2);
                x = x.parents;
            end
            
        end
        
        function c = minus(a,b)
            if isa(a,'LQR_Nodes')
                c = [a(:).node]+uminus(b);
            else
                c = a+uminus(b);
            end
        end
        
        function c = uminus(a)
                c = -[a(:).node];
        end
        
        function c = plus(a,b)
            if isa(b,'LQR_Nodes')
                c = [a(:).node]+[b(:).node];
            else
                c = [a(:).node]+b;
            end
        end
        
        function c = ne(a,b)
            if isempty(a)
                c = 1;
            else
                c = all((a.node-b.node)~=0);
            end
        end
        
        function [x,y]=size(p,n)
            if nargin == 1
                if nargout == 2
                    x = p.dim(1);
                    y = p.dim(2);
                else
                    x = p.dim;
                end
            else
                x = p.dim(n);
            end
        end
        
        function x=vertcat(varargin)
            if ~isa(varargin{2},'LQR_Nodes')
                error('cannot coabine number')
            end
            x = cat(1,varargin{:});
            x(end).pow = length(x);
        end
        
        function x=horzcat(varargin)
            if ~isa(varargin{2},'LQR_Nodes')
                error('cannot coabine number')
            end
            x = cat(2,varargin{:});
            x(end).pow = length(x);
        end
        
        function disp(obj)
            N = length(obj);
            if N==0
                return
            end
            for i=1:obj(1).dim(1)
                for j=1:obj(1).dim(2),
                    for k = 1:N
                        Str = num2str(obj(k).node);
                        fprintf('%s  ',Str(i,:));
                    end
                end
                fprintf('\n');
            end
        end
        
        %         function q=subsref(p,s)
        %             switch char(s.type),
        %                 case '()',
        %                     q = p.node(s.subs{:});
        %                 otherwise
        %                     error('option not supported')
        %             end
        %         end
        
        
    end
end
