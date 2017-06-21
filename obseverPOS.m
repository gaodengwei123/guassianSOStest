function Obs = obseverPOS(sys,t)
checkstate = xtraj.eval(t);
obstacle = sys.INPUTS.obstacle;
for j=1:size(obstacle,1)
    if norm(checkstate(1:2,i)'-obstacle(j,:))<15
        Obs = 1;
        %             plot(checkstate(1,i),checkstate(2,i),'*')
        break;
    else
        Obs=0;
    end
end
end