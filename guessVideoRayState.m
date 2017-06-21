% /*! @guessVideoRayState.m
% *************************************************************************
% <PRE>
% file.name       : guessVideoRayState.m
% related files   :
% function&ablity : to guess the state for VideoRay
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/09/9     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output:

% *************************************************************************

function [nuTime , nuGuess , nuControl ] = guessVideoRayState(INPUTS,guessPoints)
speed = INPUTS.speed;
K = size(guessPoints,1) -1;
h(1) = pi/3;  % intial angle
nuGuess = [guessPoints(1,1:2),h(1),0, 0, 0]; % intial state
nuTime(1) = 0;

for i = 1:K
    N=1;
    % guess angle
    h(i+1) = atan2(guessPoints(i+1,2) - guessPoints(i,2), guessPoints(i+1,1) -guessPoints(i,1));
    % distance between two PRM guess points
    dist(i) = norm(guessPoints(i+1,1:2)-guessPoints(i,1:2)); 
    % distance between two real guess points
    deltaDist = dist(i)/N;                                   
    while deltaDist>3
        N = N+1;
        deltaDist = dist(i)/N;
    end
    if (i ==1||i ==K)&&deltaDist<3
        deltaTime = (h(i+1)-h(i))/INPUTS.momentSat+sqrt(2*deltaDist/INPUTS.controlSat);
    else
        deltaTime = (h(i+1)-h(i))/INPUTS.momentSat+deltaDist/speed;
    end

    
    nuGuess = [nuGuess; nuGuess(end,1)+cos(h(i+1))*deltaDist, nuGuess(end,2)+sin(h(i+1))*deltaDist, 0, 0, 0.0, 0];
    
    
    

end
nuControl(:,N*K+1) = zeros(3,1);
nuTime = nuTime';
nuControl = nuControl';
end

function my_model



end