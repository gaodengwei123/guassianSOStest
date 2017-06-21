for i = 1:N-1
    intervals{i} = [taus(i) taus(i+1)]-taus(i);
    Pii = reshape(Pi(:,i),n,n);
    xbari = (x-x0i(:,i));
    xdot{i} = f0(t,x,ui(:,i)-K(taus(i))*xbari);
    V{i} = xbari'*reshape(Pi(:,i),n,n)*xbari;
    rho{i} = rhoi(:,i);
end

k=2;
Vs = repmat(V{1},N-1,k);
rhos = repmat(V{1},N-1,k);
Vdots = repmat(V{1},N-1,k);
rhodots = repmat(V{1},N-1,k);

for i = 1:N-1
    tt = linspace(intervals{i}(1),intervals{i}(2),k);
    tts(i,:) = tt;
    Vs(i,:) = msubs(V{i},t,tt);
    Vdots(i,:) = msubs(diff(V{i},t)+diff(V{i},x)*xdot{i},t,tt);
    rhos(i,:)  = msubs(rho{i},t,tt);
    rhodots(i,:) = msubs(diff(rho{i},t),t,tt);
end

M = length(Vs(:));
for i =1:M
    warning('off','MATLAB:nearlySingularMatrix');
    [gammas{i},Ls{i}] = tv_poly_rho_L_samp(t,x,tts(i),rhos(i),rhodots(i),Vs(i),Vdots(i));
    if mod(i,10) == 0
        fprintf([ int2str(i) '/' int2str(M) '\n']);
    end
end
if ~all([gammas{:}] < 0)
    if iter==1
        error('Initial rho(t) infeasible, increase c.');
    end
end
figure(3)
plot(ts,ppval(rhopp,ts),'r'); drawnow;
hold on
improve_rho   %¡¡improve the rho of ROA

plot(taus,ppval(rhopp,taus),'b'); drawnow;
save rhopp
rho0 = ppval(rhopp,taus);
for i=1:length(rho0)
    S0(:,:,i) = 1.01*Ps0(:,:,i)/rho0(i);
end