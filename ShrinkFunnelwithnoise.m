function Vtrim =  ShrinkFunnelwithnoise(V,bound)
% shrink the Funnel based on obstacles
ts = V.getbreak;
rho(length(ts))=1;
vert=[];

p0 = [1 1 1 1];
for i = 1:length(bound)
    vertbound = bound{i};
%     x0 = V.x0.eval(ts(i));
%     F=@(p,x)p(1)*x(:,1).^2+p(2)*x(:,1).*x(:,2)+p(3)*x(:,2).^2+p(4)*x(:,1)+p(5)*x(:,2)+p(6);
%     F=@(p,x)(p(1)*(x0(1) - x(:,1)) + p(3)*(x0(2) - x(:,2)))*(x0(1) - x(:,1)) + (p(2)*(x0(1) - x(:,1)) + p(4)*(x0(2) - x(:,2)))*(x0(2) - x(:,2));
%     pr = nlinfit(vertbound',zeros(size(vertbound,2),1),F,p0);
%     Vfun = F(pr,V.x');
%     H = double(0.5*diff(diff(Vfun,V.x)',V.x));
    
    vertbound = vertbound(:,randperm(size(vertbound,2),10));
    vert = [vert,[vertbound(1,:);vertbound(2,:)]];

end

% iterate over t
for i = fliplr(1:length(ts)-1)
    rho(i) = rho(i+1);
    x0 = V.x0.eval(ts(i));
    x = [vert;repmat(x0(3:end),1,size(vert,2))];
    Vvert = [];
    for k = 1:length(vert)
        Vvert = [Vvert,V.eval(ts(i),x(:,k))];
    end
    if (min(Vvert)<rho(i))
        rho(i) = min(Vvert);
    end
end
% shrink the last funnel
x0 = V.x0.eval(ts(end));
x = [vert;repmat(x0(3:end),1,size(vert,2))];
Vvert = [];
for k = 1:length(vert)
    Vvert = [Vvert,V.eval(ts(end),x(:,k))];
end
if (min(Vvert)<1)
    rho(end) = min(Vvert);
end
% todo: still need to update this

Vtrim = V.updateV(foh(ts,1./rho));
end