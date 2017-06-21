function MA = trajectory_Libiry()
sys = CarDynamic();

figure
hold on
speed = sys.INPUTS.speed;
x0 = [0 0 0 0];
length = 5;
dertX = 0.5;
for i = -5:5
    xT = [length,dertX*i,0,0];
    point = [x0;xT];
    Time = [0,length/speed];
    output = trajLibrary_gpopsSolveCar(sys, Time, point);
    plot(output.solutionPlot.state(:,1),output.solutionPlot.state(:,2))
    drawnow
end
MA




end