%% ==============ground Vehicle model==============================
clear all
close all
% close all
clc
dbstop if error
global Gaussiangao
Gaussiangao = 1;

% show the equation in local_controlDae
%% ================================================================
% 1; 2VideoRay2D: is underwater vehicle£»3 Sapcecraft: is the spacecraft tracking;4 VideoRay3D
% -------------------------------------------------------------------------
dynamic_mark = 1;
switch (dynamic_mark)
    case 1  %car: is ground vehicle
        sys = CarDynamic();
    case 2  %VideoRay2D: is underwater vehicle
        sys = VideoRay2DDynamic();
    case 3  %Sapcecraft: is the spacecraft tracking
        sys = SapcecraftDynamic();
    case 4  %VideoRay3D: is underwater vehicle
        sys = VideoRay3DDynamic();
end

%% inital state points
OL_state(:,1) = sys.INPUTS.state_initial;           % start point
U_int(:,1) = zeros(sys.getNumInput,1);              % initialize control law

current_ell = ellipsoid(OL_state(:,1),diag(sys.INPUTS.noise_w.^2));
currentq = OL_state(:,1);

%% ===================================trajectory===================

if 1
    %% GPOPS to solve the open-loop trajectory and control
    [sys,OL_state,OL_control,OL_time] = optimal_local_programming(sys, OL_state(:,1), current_ell);
    sys = sys.updateNominal(OL_time,OL_state,OL_control);   % set the sys trajectory to ploynominal form
    %% LQR initialized to calculate the ROA in every nonminal points
    [Funnel_V, sysploy] = Invariant_Funnels(sys);
%     save datavehicle20170123
else
    %     load('videoRayTrajectory.mat')
    %     load('carTrajectory.mat')
    %     load('test1.mat')
    %     load('test2.mat')
    load('datavehicle20170123.mat')
end
% nominal in PRM node
% prmSHOW = PRMnominal(sys);

% figure
% show(sys.PlotObj.prmSimple)
% v1 = sys.PlotObj.prmSimple.getEdges();
%     Funnel_V = GMPLdatabase(sysploy,Funnel_V);
%     figure
%     options.plotdims=[1 2];
%     options.inclusion='slice';
%     options.x0 = sysploy.FunTraj;
%     VFrame = Funnel_V.inFrame(sysploy.FunTraj);
%     plot_myFunnel(sysploy,VFrame,options);
%     current_ell = MLFnoise(sysploy,Funnel_V);
current_ell = Pontryagin_set_difference(sysploy,Funnel_V);


%% ===================end =======================================