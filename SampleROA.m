function SampleROA()
x = msspoly('x',2);

Q = eye(2);
R = eye(2);
sys = vandpol();
[f, df]= sys.dynamics(0,x,0);
SmpleNum = 100;
Randomp = haltonset(2);
Randomp = scramble(Randomp,'RR2');
X0 = net(Randomp,SmpleNum+1);
sumV = 0;
for i = 1:SmpleNum
    x0 = X0(i,:)';
    A = double(subs(df,x,x0));
    S=lyap(A',Q); 
    V = (x-x0)'*S*(x-x0);
    sumV = sumV+V;
end

prog = spotsosprog;
prog = prog.withIndeterminate(x);

[prog,p] = prog.newFreePoly(monomials(x,0:4));
[prog,q] = prog.newFreePoly(monomials(x,0:4));
[prog,rho] = prog.newFree(1);

prog = prog.withSOS(rho*q-p);

solver = @spot_mosek;
options = spot_sdp_default_options();
sol = prog.minimize(-rho,solver,options);



dV =  [diff(V,x(1)) diff(V,x(2))]*f;



end

% function rho = SampleBasedRho(x,V,dV)
% SmpleNum = 1000;
% Randomp = haltonset(2);
% Randomp = scramble(Randomp,'RR2');
% X0 = net(Randomp,SmpleNum+1);
% c = 1000000;
% vx = 5-10*rand(2,SmpleNum);
% for i = 1:SmpleNum
%     x_rand = 5-10*[X0(rj,1) X0(rj,2)];
%     if double(subs(dV,x,vx(:,i)))>=0&&double(subs(V,x,vx(:,i)))<c
%         c = double(subs(V,x,vx(:,i)));
%     end
% end
% rho = c;
% end