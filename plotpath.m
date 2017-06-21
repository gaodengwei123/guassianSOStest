function FindPath = plotpath(Vertex, fatherInedx,GaoNodes, target)

Ftuples = [];
index = Vertex(fatherInedx).path;
for i = 1:length(index)
    Ftuples = [Ftuples;Vertex(index(i)).Ftuple];
   h = plot(Vertex(index(i)).Ftuple);
end
LI =  target(2)-Vertex(fatherInedx).pos(2) + 5;
Ftuple = ShitFunnel(GaoNodes.FunnleLibrary(LI),GaoNodes.trajLibrary(LI), Vertex(fatherInedx));
h= plot(Ftuple);
FindPath.parent = index;
FindPath.Ftuple = [Ftuples;Ftuple];
end