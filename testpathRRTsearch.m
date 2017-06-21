clear
close all
dbstop if error
clc

GaoGraph = GaoNodes();
sys = CarDynamic();

% set initial state

% GaoGraph = GaoGraph.RRT_trajLibrary(sys,5);
% save RRTdata2017Traj
% load('RRTdata2017Traj.mat')
% GaoGraph = GaoGraph.RRT_FunnelLibrary(sys);
% save RRTdata2017funnel
load('RRTdata2017funnel.mat')

vertex = [0 0;0 50;50 50;50 0];
OBSpos = [15 30;30 25;40 10;50 20];
field = ObstacleField();
% field = field.GenerateObstacles(vertex,OBSpos,size(OBSpos,1));
%% random obs
field.number_of_obstacles = 4;
field.obstacles{1}=[13.5793	32.6698
11.9557	31.0304
15.4982	28.4543
18.0074	29.9766
16.5314	32.7869
13.5793	32.6698]';

field.obstacles{2} = [28.7823	27.6347
30.9963	28.103
32.6199	23.185
28.1919	23.0679
28.7823	27.6347]';

field.obstacles{3} = [52.1033	22.1311
48.1181	19.7892
51.0701	17.4473
53.1365	19.6721
52.1033	22.1311]';

field.obstacles{4} = [37.1956	11.1241
40.4428	8.07963
43.5424	10.5386
37.1956	11.1241]';


range = [0 0 80 50];
resolution = 0.1;
field = field.buildQRMap(range, resolution);
field = field.addsensors(OBSpos,8);
show(field);
% GaoPath = GaoPath();
GaoPath = GaoPath.setfield(field,[],[]);
GaoPath.Nodes = [];    % Node is form BRM so here it must be []
% figure
% hold on
% for i=1:9
%     plotFunnelLibary(sys,GaoGraph.FunnleLibrary(i).V,GaoGraph.FunnleLibrary(i).rho,GaoGraph.trajLibrary(i).time,GaoGraph.trajLibrary(i).state, 1);
% end
sys.INPUTS.state_goal =  [0;40;0;0];
sys.INPUTS.state_initial = [80;40;0;0];
sys = sys.setfield(field);
method = 'FunnelRRT';
start = sys.INPUTS.state_goal(1:2)';
target = sys.INPUTS.state_initial(1:2)';
% add the possible trajectories region
hold on
% Trajx = [start(1);target(1)];
% Trajyup = [start(2);start(2)+4*target(1)/5];
% Trajydown = [start(2);start(2)-4*target(1)/5];
% line(Trajx,Trajyup);
% line(Trajx,Trajydown);

AvoidRegin = [0 0;50 0;0 35;0 0];
AvoidRegin1 = [0 45;15 50;0 50;0 40];
fill(AvoidRegin(:,1),AvoidRegin(:,2),[0.5 0.5 0.5])
fill(AvoidRegin1(:,1),AvoidRegin1(:,2),[0.5 0.5 0.5])



plot(target(1),target(2),'r*');
GaoGraph = FeedbackGain(GaoGraph,sys);  % add the feedback gain into funnel

sys.INPUTS.Qk = 2*sys.INPUTS.Qk;
GaoPath = GaoPath.findpath(GaoGraph, sys, start, target, method);







