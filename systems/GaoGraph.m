function [GaoGraph] = GaoGraph(varargin)

X = varargin{1};
if nargin < 2
    range = [];
end
if nargin < 2
    range = [];
end

Map = GaoMap(X,range);
GaoGraph.Map = Map;
GaoGraph.Nodes = [];
GaoGraph.Graph = [];
GaoGraph = class(GaoGraph, 'GaoGraph');
end