function V = extractV_Function(Vpoly,x,x0)

t = msspoly('t');

if (deg(Vpoly,x)>2) error('not quadratic'); end

S = double(.5*subs(diff(diff(Vpoly,x)',x),x,0*x));

V = V_function(t,x,S,[]);
V.x0 = x0;
end