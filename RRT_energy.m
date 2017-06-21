function [path, path_length] = RRT_energy(start_state, goal_region, obstacles)
%RRT Summary of this function goes here
%   Detailed explanation goes here
%% Initializing 
epsilon=2;
x=start_state;
reach=0;
idx=0;
%% loop for reaching goal
while reach==0
    x_rand=100 * rand(1,2);
    %% extend
    % find nearest neighbor
    [idx_near,~]=knnsearch(x, x_rand);
    x_near=x(idx_near,:);
    % grow tree to get new x
    vector=x_rand-x_near;
    vector_len=sqrt(sum((vector).^2));
    x_new=x_near+vector./vector_len.*epsilon;
    % new node collision check
    collision=collision_check_point(x_new(1), x_new(2), obstacles);
    if collision==1
       continue;
    end
    % add to list
    idx=cat(1,idx,idx_near); 
    x=cat(1,x,x_new);
    plot(x_new(1), x_new(2), 'c*');
    drawnow;
    %% reach goal check
    collision=collision_check_point(x_new(1), x_new(2), goal_region);
    if collision==1
       [x_path,~,x_idx]=intersect(x_new,x,'rows');
       path=x_path;
       while 1
          x_path=x(idx(x_idx),:);
          [x_path,~,x_idx]=intersect(x_path,x,'rows');
          path=cat(1,path,x_path);
           if x_idx==1
               break;
           end
       end
       path_length=(size(path,1)-1)*epsilon;
       reach=1;
    end
    
end
end

