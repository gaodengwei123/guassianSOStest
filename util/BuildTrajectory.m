% build spline trajectory with matrix
% build by dengwei 2016
function ppform = BuildTrajectory(array)
[a,b] = size(array);
if a>=b
    warning('the spline array must be row direction');
    array = array';
end
Num = max(a,b);
ppform = spline(Num,array);
warning('donot know the breaks');
end