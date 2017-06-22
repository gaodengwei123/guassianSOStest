function node = NearestNode(Nodes,state)
if max(size(state))~=2 
    error('node must be a 2-D position')
end
dist = inf;
for i = 1:length(Nodes)
    dist0 = norm(Nodes(i).node-state(:)');
    if dist0<dist
        dist = dist0;
        node = i;
    end
end

end