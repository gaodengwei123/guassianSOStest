function [Map] = GaoMap(varargin)
X = varargin{1};
if isempty(varargin{2})
    % add range
    extand = 0.1*[(max(X(:,1))-min(X(:,1))) max(X(:,2))-min(X(:,2))];
    range = [min(X(:,1))-extand(1), min(X(:,2))-extand(2),max(X(:,1))-min(X(:,1))+2*extand(1), max(X(:,2))-min(X(:,2))+2*extand(2)];
else
    range = varargin{2};
end
% find the 4 vertex for range in map
vertex = [range(1:2);range(1:2)+[0, range(4)];range(1:2)+[range(3), range(4)];range(1:2)+[range(3), 0]];

% creat class
Map.X = X;              % obstacles in  map
Map.range = range;      % range of map
Map.vertex = vertex;        % road vertex

Map = class(Map, 'GaoMap');

return;


