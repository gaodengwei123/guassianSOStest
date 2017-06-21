function BestFather = BestPredecessor(Tree,IndexTarget)

fatherInedx = find(IndexTarget==1);  % the index of target's father in the tree
% find all the possible fathers
Cost = [Tree(fatherInedx(:)).cost];
[~,index] = min(Cost);
BestFather = fatherInedx(index);
% for i = 1:length(fatherInedx)
%     x_Predecessors = Tree(fatherInedx(i));      % find vetrex in the Tree
%     [~,Inedx]= min([x_Predecessors(:).cost]);   % find the minimum cost
%     
%     BestFather = Tree(fatherInedx(Inedx));
% 
% end
end