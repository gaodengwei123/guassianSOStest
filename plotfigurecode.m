% plot all the picture what we want
function plotfigurecode(sys,options)
if (~isfield(options,'plotPRM')), options.plotPRM =0; end
if (~isfield(options,'plotMap')), options.plotMap =0; end
 
if options.plotMap == 1
    figure(1)
    hold on
    show(sys.PlotObj.prmSimple.Map)
    xlabel('x(m)');
    ylabel('y(m)');
    title('')
end
if options.plotPRM == 1
    figure(2)
    show(sys.PlotObj.prmSimple)
end

% plot(sys.PlotObj.path(:,1),sys.PlotObj.path(:,2),'g','LineWidth',2)
% plot(sys.PlotObj.path(:,1),sys.PlotObj.path(:,2),'go','MarkerFaceColor','g')
% plot(sys.PlotObj.solutionPlot.state(: ,1), sys.PlotObj.solutionPlot.state(: ,2),'r');
% load('optimalpath2.mat')
% hold on
% plot(sys.PlotObj.path(:,1),sys.PlotObj.path(:,2),'g','LineWidth',2)
% plot(sys.PlotObj.path(:,1),sys.PlotObj.path(:,2),'go','MarkerFaceColor','g')
% plot(sys.PlotObj.solutionPlot.state(: ,1), sys.PlotObj.solutionPlot.state(: ,2),'r');
% load('optimalpath3.mat')
% hold on
% plot(sys.PlotObj.path(:,1),sys.PlotObj.path(:,2),'g','LineWidth',2)
% plot(sys.PlotObj.path(:,1),sys.PlotObj.path(:,2),'go','MarkerFaceColor','g')
% plot(sys.PlotObj.solutionPlot.state(: ,1), sys.PlotObj.solutionPlot.state(: ,2),'r');
% load('optimalpath4.mat')
% hold on
% plot(sys.PlotObj.path(:,1),sys.PlotObj.path(:,2),'g','LineWidth',2)
% plot(sys.PlotObj.path(:,1),sys.PlotObj.path(:,2),'go','MarkerFaceColor','g')
% plot(sys.PlotObj.solutionPlot.state(: ,1), sys.PlotObj.solutionPlot.state(: ,2),'r');
% %%
% load('pathdata1.mat')
% V = V.shrinkV(repmat(1.025,V.timeNum,1));
% options.plotdims=[1 2];
% options.inclusion='slice';
% plot_myFunnel(V,options);
% load('pathdata2.mat')
% V = V.shrinkV(repmat(1.025,V.timeNum,1));
% options.plotdims=[1 2];
% options.inclusion='slice';
% plot_myFunnel(V,options);
% load('pathdata3.mat')
% V = V.shrinkV(repmat(1.025,V.timeNum,1));
% options.plotdims=[1 2];
% options.inclusion='slice';
% plot_myFunnel(V,options);
% load('pathdata4.mat')
% V = V.shrinkV(repmat(1.025,V.timeNum,1));
% options.plotdims=[1 2];
% options.inclusion='slice';
% plot_myFunnel(V,options);
% 
% axis([-1,11,-1,11])
% scale = 0.1;
% quiver(sys.PlotObj.solutionPlot.state(1:10:end,1),sys.PlotObj.solutionPlot.state(1:10:end,2),cos(sys.PlotObj.solutionPlot.state(1:10:end,3)),sin(sys.PlotObj.solutionPlot.state(1:10:end,3)),scale )
% 
% %% videoRay model special
% figure(2)
% subplot(3,1,1)
% plot(sys.PlotObj.solutionPlot.time,sys.PlotObj.solutionPlot.state(: ,4))
% title('velocity of x')
% subplot(3,1,2)
% plot(sys.PlotObj.solutionPlot.time,sys.PlotObj.solutionPlot.state(: ,5))
% title('velocity of y')
% subplot(3,1,3)
% plot(sys.PlotObj.solutionPlot.time,sys.PlotObj.solutionPlot.state(: ,3))
% title('psi')
% %%
% 
% plot(sys.PlotObj.state(1,:),sys.PlotObj.state(2,:),'LineWidth',4)
% plot(sys.PlotObj.stateOL(1,:),sys.PlotObj.stateOL(2,:),'k')
% for i = 1:sys.PlotObj.Max_t
%     plot(sys.PlotObj.Eerror(i),'y')
%     plot(sys.PlotObj.Esetdiff(i),'r')
%     drawnow;
%     frame(i)=getframe(gcf);
% end
% writegif('test1.gif',frame,0.1);
% axis equal
% title('Optimized local control')
% xlabel('x- position')
% ylabel('y- position')
% 
% figure(3)
% plot(sys.PlotObj.solutionPlot.time, sys.PlotObj.solutionPlot.control(: ,1),'r')
% hold on
% plot(sys.PlotObj.solutionPlot.time, sys.PlotObj.solutionPlot.control(: ,2),'b')
% plot(sys.PlotObj.solutionPlot.time, sys.PlotObj.solutionPlot.control(: ,3),'k')
% legend('u','v','psi')
% ylabel('acceleration(m/s^2) (rad/s^2)');
% xlabel('t(s)');

end

