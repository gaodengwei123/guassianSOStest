classdef (InferiorClasses = {?ConstantTrajectory})polyniminalTrajectory
    
    properties
        pp
        tspan
        dim
    end
    
    methods
        function obj = polyniminalTrajectory(ppform)
            if isnumeric(ppform)
                ppform = ConstantTrajectory(ppform);
            end
            if isa(ppform,'ConstantTrajectory')
                ppform = mkpp([-inf,inf],ppform.pt,ppform.dim);
            end
            obj.pp = ppform;
            obj.tspan = [min(obj.pp.breaks) max(obj.pp.breaks)];
            obj.dim = obj.pp.dim;
        end
        

        
        function y = eval(obj,t)
            t=max(min(t,obj.tspan(end)),obj.tspan(1));
            y = ppval(obj.pp,t);  % still benefits from being safe (e.g. for supporting TaylorVar)
        end
        
        function dtraj = fnder(obj,order)
            if nargin<2
                order = 1;
            end
            % this requires the curve-fitting toolbox, so i'm implementing it myself below
            %      dtraj = polyniminalTrajectory(fnder(obj.pp,order));
            
            % first handle order=1
            [b,c,l,k,d] = unmkpp(obj.pp);
            if (order>k-1)  % handle the case of too-high order
                dtraj = polyniminalTrajectory(mkpp(b,0*c(:,1),d));
                return;
            end
            for i=1:k-1
                cnew(:,i) = (k-i)*c(:,i);
            end
            dtraj = polyniminalTrajectory(mkpp(b,cnew,d));
            
            if (order>1)
                dtraj = fnder(dtraj,order-1);
            end
        end
        
        function df = deriv(obj,t)
            [b,c,~,k,d] = unmkpp(obj.pp);
            if (k==1)  % handle the case of too-high order
                ppform = mkpp(b,0*c(:,1),d);
            else
                for i=1:k-1
                    cnew(:,i) = (k-i)*c(:,i);
                end
                ppform = mkpp(b,cnew,d);
            end
            df = ppval(ppform,t);
        end
        
        function nobj = uminus(obj)
            obj.pp.coefs = -obj.pp.coefs;
            nobj = polyniminalTrajectory(obj.pp);
        end
        
        function tf = eq(a,b)
            % only implement the trivial case of pptraj=scalar const
            % (which is what I need right now)
            if isscalar(a)
                tmp=b;b=a;a=tmp;
            end
            if isscalar(b)
                % first check if it's a constant
                if any(a.pp.coefs(:,1:end-1)~=0), tf=false; return; end
                tf = all(a.pp.coefs(:,end)==b);
            else
                error('not implemented yet');
            end
        end
        
        function tf = ne(a,b)
            tf = ~eq(a,b);
        end
        
        function t = getBreaks(obj)
            t = obj.pp.breaks;
        end
        
        function c = plus(a,b)
            if ~isequal(size(a),size(b))
                error('must be the same size');  % should support scalars, too (but don't yet)
            end
            if any(size(a)==0)  % handle the empty case
                c = ConstantTrajectory(zeros(size(a)));
                return;
            end
            if isa(a,'ConstantTrajectory') a=double(a); end
            if isa(b,'ConstantTrajectory') b=double(b); end
            
            if isnumeric(a)  % then only b is a polyniminalTrajectory
                [breaks,coefs,l,k,d] = unmkpp(b.pp);
                if length(d)<2, d=[d 1]; elseif length(d)>2, error('plus is not defined for ND arrays'); end
                coefs = reshape(coefs,[d,l,k]);
                for i=1:l,
                    coefs(:,:,i,end)=a+coefs(:,:,i,end);
                end
                c=polyniminalTrajectory(mkpp(breaks,coefs,[size(a,1) d(2)]));
                return;
            elseif isnumeric(b) % then only a is a polyniminalTrajectory
                [breaks,coefs,l,k,d] = unmkpp(a.pp);
                if length(d)<2, d=[d 1]; elseif length(d)>2, error('plus is not defined for ND arrays'); end
                coefs = reshape(coefs,[d,l,k]);
                for i=1:l,
                    coefs(:,:,i,end)=coefs(:,:,i,end)+b;
                end
                c=polyniminalTrajectory(mkpp(breaks,coefs,[d(1) size(b,2)]));
                return;
            end
            
            [abreaks,acoefs,al,ak,ad] = unmkpp(a.pp);
            [bbreaks,bcoefs,bl,bk,bd] = unmkpp(b.pp);
            
            if ~isequal(abreaks,bbreaks)
                breaks = unique([abreaks,bbreaks]);
                a = refine(a,breaks);
                b = refine(b,breaks);
                
                [abreaks,acoefs,al,ak,ad] = unmkpp(a.pp);
                [bbreaks,bcoefs,bl,bk,bd] = unmkpp(b.pp);
            end
            if (ak>=bk)
                coefs=acoefs; coefs(:,end-bk+1:end)=coefs(:,end-bk+1:end)+bcoefs;
            else
                coefs=bcoefs; coefs(:,end-ak+1:end)=coefs(:,end-ak+1:end)+acoefs;
            end
            
            c = polyniminalTrajectory(mkpp(abreaks,coefs,ad));
        end
        
        function a = inv(a)
            if a.pp.order==1
                [breaks,coefs,l,k,d] = unmkpp(a.pp);
                coefs = reshape(coefs,[d,l]);  % k should be 1
                for i=1:l
                    coefs(:,:,i) = inv(coefs(:,:,i));
                end
                a = polyniminalTrajectory(mkpp(breaks,coefs,d));
            else
                a = inv@Trajectory(a);
            end
        end
        
        function c = minus(a,b)
            c = plus(a,uminus(b));
        end
        
        function c = mtimes(a,b)
            if any([size(a,1) size(b,2)]==0)  % handle the empty case
                c = sparse(zeros(size(a,1),size(b,2)));
                return;
            end
            if isa(a,'ConstantTrajectory') a=double(a); end
            if isa(b,'ConstantTrajectory') b=double(b); end
            
            if isnumeric(a)  % then only b is a polyniminalTrajectory
                [breaks,coefs,l,k,d] = unmkpp(b.pp);
                if length(d)<2, d=[d 1]; elseif length(d)>2, error('mtimes is not defined for ND arrays'); end
                if isscalar(a), cd = d; elseif isscalar(b), cd = size(a); else cd = [size(a,1),d(2)]; end
                coefs = full(reshape(coefs,[d,l,k])); a=full(a);
                for i=1:l, for j=1:k,
                        c(:,:,i,j)=a*coefs(:,:,i,j);
                    end, end
                c=polyniminalTrajectory(mkpp(breaks,c,cd));
                return;
            elseif isnumeric(b) % then only a is a polyniminalTrajectory
                [breaks,coefs,l,k,d] = unmkpp(a.pp);
                if length(d)<2, d=[d 1]; elseif length(d)>2, error('mtimes is not defined for ND arrays'); end
                if isscalar(a), cd = d; elseif isscalar(b), cd = size(a); else cd = [size(a,1),d(2)]; end
                coefs = full(reshape(coefs,[d,l,k])); b=full(b);
                for i=1:l, for j=1:k,
                        c(:,:,i,j)=coefs(:,:,i,j)*b;
                    end, end
                c=polyniminalTrajectory(mkpp(breaks,c,[d(1) size(b,2)]));
                return;
            end
            
            [abreaks,acoefs,al,ak,ad] = unmkpp(a.pp);
            [bbreaks,bcoefs,bl,bk,bd] = unmkpp(b.pp);
            
            if ~isequal(abreaks,bbreaks)
                breaks = unique([abreaks,bbreaks]);
                a = refine(a,breaks);
                b = refine(b,breaks);
                
                [abreaks,acoefs,al,ak,ad] = unmkpp(a.pp);
                [bbreaks,bcoefs,bl,bk,bd] = unmkpp(b.pp);
            end
            
            if (length(ad)<2) ad=[ad 1];
            elseif (length(ad)>2) error('mtimes not defined for ND arrays'); end
            if (length(bd)<2) bd=[bd 1];
            elseif (length(bd)>2) error('mtimes not defined for ND arrays'); end
            
            acoefs = reshape(acoefs,[ad,al,ak]);
            bcoefs = reshape(bcoefs,[bd,bl,bk]);
            
            cbreaks = abreaks; % also bbreaks, by our assumption above
            if isscalar(a)
                cd = bd;
            elseif isscalar(b)
                cd = ad;
            else
                cd = [ad(1) bd(2)];
            end
            cl = al;  % also bl, by our assumption that abreaks==bbreaks
            ck = ak+bk-1;
            
            ccoefs = zeros([cd,cl,ck]);
            for l=1:cl
                for j=1:ak  % note: could probably vectorize at least the inner loops
                    for k=1:bk
                        ccoefs(:,:,l,ck-(ak-j)-(bk-k))=ccoefs(:,:,l,ck-(ak-j)-(bk-k)) + acoefs(:,:,l,j)*bcoefs(:,:,l,k);
                    end
                end
            end
            c = polyniminalTrajectory(mkpp(cbreaks,ccoefs,cd));
        end
        
        function obj = refine(obj,newbreaks)
            obj = polyniminalTrajectory(pprfn(obj.pp,newbreaks));
        end
        function s = size(obj,dim)
            s=obj.dim;
            if (length(s)==1) s=[s,1]; end
            if (nargin>1) s=s(dim); end
        end
        
        function varargout = subsref(a,s)
            if (length(s)==1 && strcmp(s(1).type,'()'))
                [breaks,coefs,l,k,d] = unmkpp(a.pp);
                coefs = reshape(coefs,[d,l,k]);
                if (length(s.subs)==1 && length(d)>1)
                    subs = cell(1,length(d));
                    [subs{:}] = ind2sub(d,s.subs{:});
                    s.subs = subs;
                end
                s.subs = {s.subs{:},':',':'};
                coefs = subsref(coefs,s);
                d=size(subsref(a.eval(a.tspan(1)),s));
                if numel(d)==2 && d(2)==1, d = d(1); end  % column vectors are a special case that's handled differently by the spline class
                varargout{1} = PPTrajectory(mkpp(breaks,coefs,d));
            else % use builtin
                varargout=cell(1,max(nargout,1));
                [varargout{:}] = builtin('subsref',a,s);
            end
        end
        function a = subsasgn(a,s,b)
            if (length(s)==1 && strcmp(s(1).type,'()'))
                if isempty(a) % handle the special case
                    [breaks,coefs,l,k,d] = unmkpp(b.pp);
                    e=[];
                    d_extended = [d,l,k];
                    coefs = reshape(coefs,d_extended);
                    s.subs = {s.subs{:},':',':'};
                    e = subsasgn(e,s,coefs);
                    a = PPTrajectory(mkpp(breaks,e,d_extended(1:end-2)));
                    return;
                end
                if isnumeric(a) % then b must be a PPTrajectory
                    breaks = b.getBreaks();
                    a = PPTrajectory(zoh(breaks,repmat(a,[1+0*size(a),length(breaks)])));
                elseif isa(a,'ConstantTrajectory')
                    breaks = b.getBreaks();
                    a = PPTrajectory(zoh(breaks,repmat(a.pt,[1+0*size(a),length(breaks)])));
                end
                typecheck(a,'PPTrajectory');  % i believe this is the only way this method would get called
                [breaks,coefs,l,k,d] = unmkpp(a.pp);
                if isnumeric(b)
                    b = PPTrajectory(zoh(breaks,repmat(b,[1+d*0,length(breaks)])));
                elseif isa(b,'ConstantTrajectory')
                    b = PPTrajectory(zoh(breaks,repmat(b.pt,[1+d*0,length(breaks)])));
                end
                typecheck(b,'PPTrajectory');
                [breaks2,coefs2,l2,k2,d2] = unmkpp(b.pp);
                if ~isequal(breaks,breaks2)
                    a = subsasgn@Trajectory(a,s,b);
                end
                if (k<k2)
                    coefs =  [zeros(prod(d)*l,k2-k),coefs];
                    k=k2;
                elseif (k2<k)
                    coefs2 = [zeros(prod(d2)*l2,k-k2),coefs2];  % pad with zeros
                    k2=k;
                end
                coefs = reshape(coefs,[d,l,k]);
                coefs2 = reshape(coefs2,[d2,l2,k2]);
                s.subs = {s.subs{:},':',':'};
                coefs = subsasgn(coefs,s,coefs2);
                d = size(coefs); d=d(1:end-2);
                a = PPTrajectory(mkpp(breaks,coefs,d));
            else
                a = subsasgn@Trajectory(a,s,b);
            end
        end
    end
end
