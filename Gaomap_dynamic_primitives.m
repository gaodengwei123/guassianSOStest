function output = Gaomap_dynamic_primitives(sys,Path,pathsort,Nodes,start,target)
global CONSTANS
CONSTANS = sys.INPUTS;
for k = 1:length(Path)
    index = pathsort(k);
    pathNode = [start';Nodes(Path{index}.path,:);target'];
    plot(pathNode(:,1),pathNode(:,2),'linewidth',3,'color','r')
    [guessTime, guesspoint, guessControl ] = OptGuessEnhancer(sys, sys.INPUTS, pathNode);
    hold on
    [output, ~] = gpopsSolveCar(guessTime,guesspoint,guessControl,guesspoint(1,:),guesspoint(end,:));
    plot(output.solutionPlot.state(:,1),output.solutionPlot.state(:,2))
end

end