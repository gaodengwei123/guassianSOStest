function output = Gaomap_dynamic_planning_Phases(sys,Path,pathsort,Nodes,start,target)
global CONSTANS
CONSTANS = sys.INPUTS;

for k = 1:length(Path)
    index = pathsort(k);
    pathNode = [start';Nodes(Path{index}.path,:);target'];
    plot(pathNode(:,1),pathNode(:,2),'linewidth',3,'color','r')
    [PhasesTime, Phasespoint] = OptGuessPhases(sys, pathNode);
    hold on
 
    [output, ~] = GaoMap_gpopsSolveCar(sys,PhasesTime,Phasespoint);

    for j = 1:6
        plot(output.solutionPlot(j).state(:,1),output.solutionPlot(j).state(:,2))
    end
end






end