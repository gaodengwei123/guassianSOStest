clear
close all
dbstop if error
clc
%% test
X = [-2.3 -0.7;-1.5 1.3;0 -0.5;0.8 1.2;2 -1.5];
range = [-3.5,-3,7,6];
Map = GaoMap(X,range,[0 -3]);

%% vedioRay
% X = [3 0.5;3 1;3 1.5;3 2;3 2.5;3 3;3 3.5;3 4;3 4.5;3 5;3 5.5;3 6;
%      6 6;6 5.5;6 5;6 4.5;6 4;
%      7 9.5;7 9;7 8.5;7 8];
% eX= [0 2; 0 4;0 6;0 8; 0 10;2 0;4 0;6 0;8 0;10 0;10 2; 10 4;10 6;10 8; 10 10;2 10;4 10;6 10;8 10;10 10];
% range = [0,0,10,10];
% Map = GaoMap(X,range);
% Obsline = [3 6 3 0;6 6 6 4;7 8 7 10]; % [x1 y1 x2 y2;...]
%% car
% X = [10 30; 30 40; 25 20];
% range = [0,0,50,50];
% Map = GaoMap(X,range);
% add constraint in edge when use voronoi to build nodes
% Map = addEdgeConstraint(Map);
field = ObstacleField();
% field = field.GenerateObstacles(getvertex(Map),[],[]);
% simpleMap = buildQRMap(field, range, 0.1);

load('BRMmap.mat'); % from image_to_QRmap
field.Map = BRMmap;
field.range = [0 0 50 50];
field.resolution = 0.1;
field = field.addsensors([12 11;16 12;20 15;23 18;33 18;33 11;37 25;40 32],2);
% ;30 22
% figure
% hold on
field = field.Map2Obstacle();
show(field);

boundary = [1.0514	13.2009
8.52804	17.1729
8.76168	40.6542
12.2664	39.2523
14.8364	44.5093
25.1168	44.1589
30.4907	42.6402
35.2804	38.2009
43.1075	37.3832
46.729	35.3972
48.3645	28.972
44.0421	26.8692
45.0935	23.014
37.0327	20.9112
37.5	19.0421
35.0467	18.2243
35.7477	14.9533
37.8505	11.9159
36.3318	11.3318
36.6822	10.1636
32.5935	7.71028
28.6215	15.1869
12.2664	5.60748
9.11215	4.6729
5.02336	4.78972
2.80374	7.59346
1.0514	13.2009];
plot(boundary(:,1),boundary(:,2),'k','linewidth',3)
GaoPath = GaoPath();
GaoPath = GaoPath.setfield(field,35,6);
GaoGraph = GaoNodes();
GaoGraph = GaoGraph.getEdges(GaoPath);
sys = CarDynamic();
sys.INPUTS.state_goal =  [25.9;39.9;0;0];
sys.INPUTS.state_initial = [8.5;8.5;pi/2;0];

% start = 1;
% target = 17;
% GaoGraph = GaoGraph.Calculate_trajLibrary(sys, GaoPath, start, target);
% save data2017Traj
% load('data2017Traj.mat')
% GaoGraph = GaoGraph.Calculate_FunnelLibrary(sys);
% save data2017funnel
load('data2017funnel.mat')
method = 'FunnelAstar';

sys = sys.setfield(field);
GaoPath = GaoPath.findpath(GaoGraph, sys, start, target, method);

% if node is short than this distant they will be merged
% distRemove = 1;
% Graph= creat_voronoi_graph(Map,distRemove);
% GaoGraph= creat_PRM_graph(Map,distRemove);
% Depth-First-Search
% startNode = [0;0];
% targetNode = [50;40];
% MA = trajectory_primitives(GaoGraph,startNode,targetNode);
% MA = trajectory_Libiry(sys);
% Path = DFS_gao_node(Map,Nodes,NearestNode(Nodes,startNode),NearestNode(Nodes,targetNode));
% Path = BFS_gao_node(Map,Nodes,NearestNode(Nodes,startNode),NearestNode(Nodes,targetNode));

% sort the path with cost
% [~,pathsort] = sort(cellfun(@(Path) Path.cost,Path));

Nodes = vec2mat([GaoGraph.Graph(:).node],2);
Vec = vec2mat([GaoGraph.trajLibrary(:).mark],3);
for i = 1:length(GaoPath.path)-1
    if i==1
        Trajvec = [0,GaoPath.path(i:i+1)];
    else
        Trajvec = [GaoPath.path(i-1:i+1)];
    end
    index(i) = find(ismember(Vec,Trajvec,'rows'),1);
    
end
% intersect
% Traj = vec2mat([GaoGraph.trajLibrary(:).state(:,1)],2);

%     index = pathsort(i);
%     index = i;
% a = plot(Nodes(GaoPath.path,1),Nodes(GaoPath.path,2),'linewidth',1,'color','r');
% drawnow
% Frame(i) = getframe(gcf);
% delete(a)
ShrinkRho = 1;
% index = [2     9    13    20    28    38    48    55]; %shrotestpath
% index = [2     6   167   165   163   161   140   133    63    77];
for j = length(index):-1:1
    % plot trajectory
%     plot(GaoGraph.trajLibrary(index(j)).state(:,1),GaoGraph.trajLibrary(index(j)).state(:,2),'linewidth',2,'color','r');
    % plot funnel
    [Pini,rhoini] = plotFunnelLibary(sys,GaoGraph.FunnleLibrary(index(j)),GaoGraph.trajLibrary(index(j)).time,GaoGraph.trajLibrary(index(j)).state, ShrinkRho,j);
%     if j>1
%         ShrinkRho = funnel_composition(sys,GaoGraph.FunnleLibrary(index(j-1)).V,GaoGraph.trajLibrary(index(j-1)).time,GaoGraph.trajLibrary(index(j-1)).state, Pini, rhoini);
%     end
    if j==1
        
    end
end

% writegif('pathtest1',Frame,1);
% sys = CarDynamic();
% % prmSHOW = Gaomap_dynamic_planning(sys,Path,pathsort,Link,startNode,targetNode);
% prmSHOW = Gaomap_dynamic_planning_Phases(sys,Path,pathsort,Link,startNode,targetNode);























%% plot node in gragh
% s = [];
% t = [];
% for i = 1:length(Nodes)
%     s = [s, repmat(Nodes{i}.Tag,1,Nodes{i}.linkNum)];
%     t = [t, Nodes{i}.linkTag];
% end
% G = digraph(s,t);