function prmSHOW = PRMnominal(sys,x0)
% node = sys.PlotObj.prmSimple.getEdges();

% plot(node(1,:),node(2,:),'*')
start = sys.INPUTS.state_initial;
termnel = sys.INPUTS.state_goal;
[sys,guessTime,guessstate,guessControl] = Optguesspoint(sys,start(1:2),termnel(1:2));

hold on
for i = 1:size(guessstate,1)-1
    x0 = guessstate(i,:)';
    xT = guessstate(i+1,:)';
    [output{i}, ~] = PRM_gpopsSolveCar(sys,x0,xT);
    plot(output{i}.solutionPlot.state(:,1),output{i}.solutionPlot.state(:,2),'r')
    
end