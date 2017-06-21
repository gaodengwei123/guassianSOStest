prog = mssprog;

l=1; %linear 
[prog,rhocoeff] = new(prog,(N-1)*(l+1),'free');
% Build the segmentaus
rhos = reshape(rhocoeff,(N-1),l+1)*monomials(t,0:l);

dtaus = diff(taus);

% Necessary relationships between segmentaus, rho(T), and 0
[prog,po] = new(prog,N-1,'pos');

prog.eq = (1 - subs(rhos(N-1),t,dtaus(end)));
prog.eq = po(N-1) - subs(rhos(1),t,0);
for i = 1:N-2
    prog.eq = po(i) - subs(rhos(i),t,dtaus(i));
    prog.eq = (subs(rhos(i+1),t,0) - subs(rhos(i),t,dtaus(i)));
end

% Build rhos and rhodotaus
rhoss = repmat(rhos,1,k);
rhodotauss = repmat(rhos,1,k);
for i = 1:N-1
    tt = linspace(intervals{i}(1),intervals{i}(2),k);
    rhoss(i,:) = msubs(rhos(i),t,tt);
    rhodotauss(i,:) = msubs(diff(rhos(i),t),t,tt);
end

Lis = repmat(Ls{1},(N-1),k);
for i = 1:M
    Lis(i) = Ls{i}';
end

nsd =  Vdots - rhodotauss+Lis.*(Vs-rhoss);
% if (options.stability)
%     Vdot{i} = Vdot{i}*rho(i)-V{i}*rhodot(i);
%     prog.sos = -Vdot{i} + L{i}*(V{i}-rho(i));
% else
%     prog.sos = -(Vdot{i}-rhodot(i)+L{i}*(V{i}-rho(i)));
% end
prog.sos = -nsd;

obj = 0;
rhoint = reshape(rhocoeff,(N-1),l+1)*diag(1./(1:l+1))*(t*monomials(t,0:l));
for i = 1:N-1
    obj = obj + subs(rhoint(i),t,dtaus(i));
end

[prog,info] = sedumi(prog,-obj,1,struct());

rho = cell(N-1,1);
for i = 1:N-1
    rho{i} = prog(rhos(i));
end

rhopp = sim_ps2pp(t,prog(rhos'),taus');




