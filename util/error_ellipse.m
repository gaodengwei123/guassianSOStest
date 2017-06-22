% /*! @error_ellipse.m
% *************************************************************************
% <PRE>
% file.name       : error_ellipse.m
% related files   :
% function&ablity :
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2017/2/23     3.00                     error_ellipse
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input :

% output:
% *************************************************************************
function E = error_ellipse(avg,covariance)

[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, ~] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
%     smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
%     smallest_eigenvec = eigenvec(1,:);
end
if ~isreal(smallest_eigenval)||~isreal(largest_eigenval)
   error('covaniance is wrong') 
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;  % sqrt(5.991)
phi = angle;
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

D = diag([1/a^2 1/b^2]);
P = R'*D*R;
% make symmetric
Q = P^-1;
Q = chol(Q);
Q = Q'*Q;
E = ellipsoid(avg,Q);
% plot(EE,'r')
end
