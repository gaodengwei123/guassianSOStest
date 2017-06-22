function Path = DFS_gao(G,start,target)

if isa(G ,'digraph')
    flagG = 1;
elseif isa(G ,'graph')
    flagG = 0;
else
    error('G must be graph or digraph')
end

Edge = table2array(G.Edges);
m = max(Edge(:));                % nodes number
A = compresstable2matrix(Edge);  % make graph as matrix

top = 1;                         % top of stack
k = 1;
stack{top}.node = start;         % add first node to stack
stack{top}.path = start;

while top ~= 0                  % judge stack is empty or not
    pre_len = length(stack);    % length of stack before next serach
    i = stack{top}.node;        % get top in stack
    for j = 1:m
        if A(i,j)==1 && isempty(find(stack{top}.path==j,1))    % if node is link and do not in stack
            top=top+1;                              % extend
            stack{top}.node = j;                    % add new node in stack
            if ~isfield(stack{top},'path')
                stack{top}.path = stack{top-1}.path;
            end
            stack{top}.path = [stack{top}.path, j];
            break;
        end
    end
    
    if length(stack)==pre_len   % if length of stack is not grow out of stack
        stack{top} = [];
        top=top-1;
    end
    if top~=0
        if stack{top}.node == target
            Path{k} = stack{top}.path;
            k = k+1;
            % out of stack
            stack(top) = [];
            top=top-1;
            continue
        end
    end
end

% A=compresstable2matrix(re);
% figure;
% netplot(A,1)

end

function A=compresstable2matrix(b)
n=size(b,1);
m=max(b(:));
A=zeros(m,m);

for i=1:n
    A(b(i,1),b(i,2))=1;
    A(b(i,2),b(i,1))=1;
end

end