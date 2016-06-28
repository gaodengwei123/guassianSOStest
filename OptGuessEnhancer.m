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
function [nuTime , nuGuess , nuControl ] = OptGuessEnhancer(INPUTS, guessPoints )
% --------------------------------------
% BEGIN : function OptGuessEnhancer .m
%
% This script takes the main waypoints provided in the courseOptShell
% script and added in additional points between them to provide a more
% accurate initial guess to the GPOPS algorthim .
%
% --------------------------------------
N = 1; % Number of steps to add between each waypoint
speed = INPUTS.speed;
K = size(guessPoints,1) -1;
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

nuTime = 0;
nuControl = 0;
for i = 1:N*K
    nuTime(i+1) = nuTime(end)+ delTime ;
    nuControl(i+1) = 0;
end
nuTime = nuTime';
nuControl = nuControl';