function GaoGraph = creat_voronoi_graph(Map,varargin)
X = Map.X;
range = Map.range;

figure
hold on
rectangle('Position',range)

[vxx,vy] = voronoi(X(:,1),X(:,2));
voronoi(X(:,1),X(:,2));

% remove link which is intesect with obstacles
if nargin>2
    Obsline = varargin{2};
    intrsct = [];       % index of line intesect
    for j = 1:size(Obsline,1)
        for i = 1:size(vxx,2)
            x=[Obsline(j,1), Obsline(j,3),vxx(1,i),vxx(2,i)];
            y=[Obsline(j,2), Obsline(j,4),vy(1,i),vy(2,i)];
            
            dt1=det([1,1,1;x(1),x(2),x(3);y(1),y(2),y(3)])*det([1,1,1;x(1),x(2),x(4);y(1),y(2),y(4)]);
            dt2=det([1,1,1;x(1),x(3),x(4);y(1),y(3),y(4)])*det([1,1,1;x(2),x(3),x(4);y(2),y(3),y(4)]);
            if(dt1<=0 && dt2<=0)
                intrsct = [intrsct,i];         %If lines intesect
            end
        end
    end
    vxx(:,intrsct) = [];
    vy(:,intrsct) = [];
end

%% remove nodes which are close
if nargin>1
    distRemove = varargin{1};
    i = 1;
    while i ~= size(vxx,2)
        distx = (vxx(1,i)-vxx(2,i))^2;
        disty = (vy(1,i)-vy(2,i))^2;
        dist = sqrt(distx+disty);
        if dist<distRemove
%             replace the node combine
            conindexx = find(vxx==vxx(2,i));
            conindexy = find(vy==vy(2,i));
            conindex = intersect(conindexx,conindexy);
            vxx(conindex) = vxx(1,i);
            vy(conindex) = vy(1,i);
%             remove the too close link
            vxx(:,i) = [];
            vy(:,i) = [];
            i=i-1;
        end
        i=i+1;
    end
end

V = [vxx(:),vy(:)];
%% remove the infinte node
removeidex=[];
for i = 1:2:size(V,1)-1
%     if any(abs(V(i,:))>6)||any(abs(V(i+1,:))>6)
    if V(i,1)<range(1)||V(i,2)<range(2)||V(i,1)>range(1)+range(3)||V(i,2)>range(2)+range(4)||...
            V(i+1,1)<range(1)||V(i+1,2)<range(2)||V(i+1,1)>range(1)+range(3)||V(i+1,2)>range(2)+range(4)
        removeidex = [removeidex,i,i+1];
    end
end
V(removeidex,:)=[];
%%
vxx = reshape(V(:,1),[2,size(V,1)/2]);
vy = reshape(V(:,2),[2,size(V,1)/2]);

tV = unique(V,'rows');% all the nodes

nume = size(vxx,2);
vx = [vxx; NaN(1,nume)];
vx = vx(:);
vy = [vy; NaN(1,nume)];
vy = vy(:);
Vxy = [vx,vy];  % all the link breaks with NAN
% plot feasible link
line(vx,vy,'linewidth',2)
plot(tV(:,1),tV(:,2),'y*')
NodeNum = size(tV,1);

% find which is the same points in Vxy and return the index
% tV is the vertex
% vx is the graph
for i = 1:NodeNum
    Tag = [];
    weights = [];
    tindex = intersect(find(vx==tV(i,1)),find(vy==tV(i,2)));
    node = Vxy(tindex(1),:);
    Node.node = node;   % read the value of each node
    Node.Tag = i;       % mark the index of each node
    for j = 1:size(tindex,1)
        if rem(tindex(j),3)==1
            linkNode = Vxy(tindex(j)+1,:);
            node = [node;linkNode];
        else
            linkNode = Vxy(tindex(j)-1,:);
            node = [node;linkNode];
        end
        Tag(j) = intersect(find(tV(:,1)==linkNode(1)),find(tV(:,2)==linkNode(2)));
        weights(j) = norm(Node.node - linkNode);
    end
    Node.link = node(2:end,:); % the next node form this vertex
    Node.linkTag = Tag;
    Node.linkNum = size(Tag,2);
    Node.linkweights = weights;
    Nodes(i) = Node;
    plot(Node.node(1),Node.node(2),'y*')
    text(Node.node(1),Node.node(2),num2str(i))
end
 
GaoGraph.Nodes = Nodes;
GaoGraph.vertex = tV;
GaoGraph = GaoNodes();






