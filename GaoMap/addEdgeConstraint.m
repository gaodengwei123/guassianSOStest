function obj = addEdgeConstraint(obj)
X = obj.X;
vertex = obj.vertex;
% find nodes are closest to 4 vertex
for j = 1:4
    dist = inf;
    for i = 1:size(X,1)
        distnow = norm(X(i,:)-vertex(j,:));
        if dist > distnow
            dist = distnow;
            vertexnear = i;
        end
    end
    X = [X;[X(vertexnear,1),vertex(j,2)];[vertex(j,1),X(vertexnear,2)]];
end

% add nodes in edge
%  --------------3-------------
% |  2                      3   |
% |                             |
% 2                             4
% |                             |
% |  1                      4   |
%  --------------1-------------

dist(1) = min(abs(X(1:end-8,2)-vertex(1,2)));
dist(2) = min(abs(X(1:end-8,1)-vertex(1,1)));
dist(3) = min(abs(X(1:end-8,2)-vertex(3,2)));
dist(4) = min(abs(X(1:end-8,1)-vertex(3,1)));

X1(1) = X(end-1,1); X2(1) = X(end-7,1);
X1(2) = X(end-4,2); X2(2) = X(end-6,2);
X1(3) = X(end-3,1); X2(3) = X(end-5,1);
X1(4) = X(end-2,2); X2(4) = X(end,2);
for i = 1:4
    n = 1;
    distmin =inf;
    while dist(i)<distmin
        distmin = min(X1(i)-X2(i))/n;
        if distmin>dist(i)
            n = n+1;
        else
            break
        end
    end
    if i ==1||i==3
        X = [X;[linspace(X1(i),X2(i),n+1)' repmat(vertex(i,2),n+1,1)]];
    else
        X = [X;[repmat(vertex(i,1),n+1,1) linspace(X1(i),X2(i),n+1)']];
    end
end
X = unique(X,'rows');
obj.X = X;

obj.vertex = vertex;
end