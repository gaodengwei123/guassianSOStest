% /*! @OptGuessEnhancer.m
% *************************************************************************
% <PRE>
% file.name       : OptGuessEnhancer.m
% related files   :
% function&ablity : to guess the state for guess
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/06/17     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output:

% *************************************************************************
function [nuTime , nuGuess , nuControl ] = OptGuessEnhancer(sys, INPUTS, guessPoints )
% --------------------------------------
% BEGIN : function OptGuessEnhancer .m
%
% This script takes the main waypoints provided in the courseOptShell
% script and added in additional points between them to provide a more
% accurate initial guess to the GPOPS algorthim .
%
% --------------------------------------
switch (sys.mark)
    case 1  % car
        N = 1; % Number of steps to add between each waypoint
        speed = INPUTS.speed;
        K = size(guessPoints,1)-1;
        nuGuess = [guessPoints(1,1:2),atan2(guessPoints(2 ,2)-guessPoints(1 ,2), guessPoints(2 ,1)-guessPoints(1 ,1)),0];
        for i = 1:K
            dist(i) = norm(guessPoints(i+1,1:2)-guessPoints(i,1:2));
            deltaTime(i) = dist(i)/speed;
            h = atan2(guessPoints(i+1,2) -guessPoints(i ,2) ,guessPoints(i+1 ,1) -guessPoints(i ,1));
            deltaDist = dist(i)/N;
            for j = 1:N
                nuGuess = [nuGuess; nuGuess(end,1)+cos(h)*deltaDist,nuGuess(end,2)+sin(h)*deltaDist,h,0];
            end
        end
        totalTime = sum(deltaTime(1 ,:));
        delTime = totalTime/(N*K);
        
        nuTime(1) = 0;
        nuControl(:,1) = zeros(sys.getNumInput,1);
        for i = 1:N*K
            nuTime(i+1) = nuTime(end)+ delTime ;
            nuControl(:,i+1) = zeros(sys.getNumInput,1);
        end
        nuTime = nuTime';
        nuControl = nuControl';
        
    case 2  % VideoRay2D
%         [nuTime , nuGuess , nuControl ] = guessVideoRayState(INPUTS,guessPoints);
        N = 1; % Number of steps to add between each waypoint
        speed = INPUTS.speed;
        K = size(guessPoints,1) -1;
        h(1) = pi/3;
        nuGuess = [guessPoints(1,1:2),h(1),0, 0, 0];
        nuTime(1) = 0;
        
        for i = 1:K
            dist(i) = norm(guessPoints(i+1,1:2)-guessPoints(i,1:2)); % distance between two PRM guess points
            deltaDist = dist(i)/N;                                   % distance between two real guess points
            h(i+1) = atan2(guessPoints(i+1,2) - guessPoints(i,2), guessPoints(i+1,1) -guessPoints(i,1));
            if (i == 1||i==K)&&(dist(i)<3||nuTime(i)<1)   % donnot satisfy condition
                deltaTime(i) = sqrt(2*dist(i)/INPUTS.controlSat);
            else
                deltaTime(i) = dist(i)/speed;
            end

            for j = 1:N
                if i==K&&j==N
                    nuGuess = [nuGuess; nuGuess(end,1)+cos(h(i+1))*deltaDist, nuGuess(end,2)+sin(h(i+1))*deltaDist, 0, 0, 0.0, 0];
                else
                    nuGuess = [nuGuess; nuGuess(end,1)+cos(h(i+1))*deltaDist, nuGuess(end,2)+sin(h(i+1))*deltaDist, h(i+1), speed, 0.0, (h(i+1)-h(i))/deltaTime(i)];
                end
                nuTime((i-1)*N+j+1) = nuTime(end)+ deltaTime(i);
                nuControl(:,(i-1)*N+j) = guess_control(nuGuess(i,:)',nuGuess(i+1,:)',INPUTS);
            end
        end
        nuControl(:,N*K+1) = zeros(3,1);
        nuTime = nuTime';
        nuControl = nuControl';
end

end






