% transform system to spline form
% build by dengwei 2017
% row and column have some worning here, be careful!!!!
function ploy = handle2ploy(t,a)
N = size(t,1);
[s1,s2] = size(a(0));
if s1>=s2
    b = zeros(s1,N);
    for i = 1:N
        b(:,i) = a(t(i));
    end
    ploy = spline(t,b);
else
    b = zeros(N,s2);
    for i = 1:N
        b(i,:) = a(t(i));
    end
    ploy = spline(t',b');
    [breaks,coefs,~,~,~] = unmkpp(ploy);
    ploy = mkpp(breaks,coefs,[s1 s2]);
end
end