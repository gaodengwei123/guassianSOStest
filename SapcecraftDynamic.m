% /*! @SapcecraftDynamic.m
% *************************************************************************
% <PRE>
% file.name       : SapcecraftDynamic.m
% related files   :
% function&ablity : build system
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         : not finished
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/11/10     1.00       gaodengwei   build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:
%
% output:
% *************************************************************************
classdef SapcecraftDynamic<GaoSystem
    methods
        function obj = SapcecraftDynamic()
            mark = 3;                  % 1 is the Car 2-D model
            obj = obj@GaoSystem(mark,6,3);
            % parameters(INPUTS)
            obj = SetParameter(obj);   % set parameters for the vechicle
            obj = Dyanmic_f(obj);      % build the dynamic function
            obj.getNumObstacles = size(obj.INPUTS.obstacle,1);   % Number of Obstacle
        end
        
        function obj = SetParameter(obj)
            deg = pi/180;              % rad->degree
            INPUTS.state_initial = zeros(obj.getNumStates,1);%initilal points
            % obstacles location
            INPUTS.obstacle = [-600 40 60];
            INPUTS.obstacle_velocity=[0;0;0];
            
            INPUTS.obstacleRadius = 4;              % obstacles radius
            INPUTS.obstacleRadius_multiplier = 1;
            INPUTS.state_goal = [0;0;0;0;0;0];      % terminal point
            
            INPUTS.miu = 3.986004405e14;                   % gravitational constant
            INPUTS.Mt = 1000;                              % target mass/kg
            % =======================target imformation============
            INPUTS.target_a = 673e3+6371e3;                % target semi-major axis
            INPUTS.target_e = 0.003;                       % target orbit eccentricity
            INPUTS.target_RA = 277.5*deg;                  % Right Ascension of Ascending Node\目标航天器初始升交点赤经°
            INPUTS.target_incl = 98.05*deg;                % orbit inclination\轨道倾角
            INPUTS.target_Omega = 0*deg;                   % argument of perigee\近地点幅角
            INPUTS.target_TA = 0*deg;                      % true anomaly\真近点角
            INPUTS.target_h = sqrt(INPUTS.miu*INPUTS.target_a*(1-INPUTS.target_e^2));% angular momentum\角动量
            % target location and vecolcity in inertial frame
            [initRt_i,initVt_i] = sv_from_coe(INPUTS.target_h,INPUTS.target_e,INPUTS.target_RA,INPUTS.target_incl,INPUTS.target_Omega,INPUTS.target_TA,INPUTS.miu);
            INPUTS.initRt_i = initRt_i;
            INPUTS.initVt_i = initVt_i;
            % =================chaser imformation==================
            INPUTS.Mc = 1400;                          % kg/追踪航天器的质量
            INPUTS.Jc = diag([954,1030,605]);          % 追踪航天器的转动惯量
            INPUTS.Qc0=[0.3062;0.8777;0.2637;0.2575];  % 追踪航天器的初始姿态
            INPUTS.Qc0 = Qc0/norm(Qc0);
            INPUTS.Wc0 = [0.000;0.000;0.000];          % 追踪航天器的初始姿态角速度
            INPUTS.initRc_i = initRt_i + Ci_o'*initState(1:3);
            INPUTS.initVc_i = initVt_i + initState(4:6) + cross(Wt,initState(1:3));
            obj.INPUTS = INPUTS;
        end
        
        function obj = Dyanmic_f(obj)
            error('this dynamic have not been accomplished')

        end
    end
    
end






















