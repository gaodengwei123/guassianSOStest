function Vtrim =  ShrinkFunnelwithObstacles(sys,V)
% shrink the Funnel based on obstacles
ts = V.getbreak;
rho(length(ts))=1;
vert=[];
% transform the QR obstacle to Polygon
for i = 1:sys.getNumObstacles
    obstacle = obstaclePolygon(sys.INPUTS.obstacle(i,:),sys.INPUTS.obstacleRadius);
    vert = [vert,[obstacle(1,:);obstacle(2,:)]];
end

% iterate over t
for i = fliplr(1:length(ts)-1)
    rho(i) = rho(i+1);
    x0 = sys.FunTraj.eval(ts(i));
    x = [vert;repmat(x0(3:end),1,size(vert,2))];
    Vvert = [];
    for k = 1:length(vert)
        Vvert = [Vvert,V.eval(ts(i),x(:,k))];
    end
    if (min(Vvert)<rho(i))
        rho(i) = min(Vvert);
    end
end

% todo: still need to update this

Vtrim = V.updateV(foh(ts,1./rho));
end