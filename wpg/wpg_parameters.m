classdef wpg_parameters<handle
    %%% walking pattern generator parameters class
    properties
        robot                   % preset robot chosen
        type_traj               %preset trajectory chosen
        firstSS                 % preset firstSS chosen
        poly_degree             % degree of ZMP polynomials
        nbpolyssp               % number of polynomial for ZMP in ssp
        nbpolydsp               % number of polynomial for ZMP in dsp
        nbpolypi                % number of polynomial for starting dsp
        nbpolypf                % number of polynomial for stopping dsp
        nbpolyzmp1              % number of polynomial for ZMP1 in dsp
        backtoankle             % distance from back to ankle of foot
        fronttoankle            % distance from  front to ankle of foot
        exttoankle              % distance from exterior to ankle of foot
        inttoankle              % distance from interior to ankle of foot
        sole_margin             % sole margin to keep stability
        xankmax                 % stepping forward max
        xankmin                 % stepping forward min (if negative, it means stepping backward max)
        yankmin                 % width min between ankles
        yankmax                 % width max between ankles
        nbstep                  % number of foot step
        pankinit_firstSS        % initial position of the first SS foot
        pankinit_firstinair     % initial position of non first SS foot
        pankfin_lastSS          % position of the last SS foot
        pankfin_lastinair       % position of the non last SS foot
        xpsa_zmpinit            % initial position, speed and acceleration of ZMP on x axis
        xpsa_zmpfin             % final position, speed and acceleration of ZMP on x axis
        ypsa_zmpinit            % initial position, speed and acceleration of ZMP on y axis
        ypsa_zmpfin             % final position, speed and acceleration of ZMP on y axis
        xpcominit               % initial position of COM on x axis
        xpcomfin                % final position of COM on x axis
        ypcominit               % initial position of COM on y axis
        ypcomfin                % final position of COM on y axis
        xscominit               % initial speed of COM on x axis
        xscomfin                % final speed of COM on x axis
        yscominit               % initial speed of COM on y axis
        yscomfin                % final speed of COM on y axis
        xpsa_zmp1init           % initial position, speed and acceleration of ZMP1 on x axis
        xpsa_zmp1fin            % final position, speed and acceleration of ZMP1 on x axis
        ypsa_zmp1init           % initial position, speed and acceleration of ZMP1 on y axis
        ypsa_zmp1fin            % final position, speed and acceleration of ZMP1 on y axis
        z                       % COM altitude in m
        g                       % gravity constant en m/s²
        Tc                      % time constant
        w                       % inverse time constant
        m                       % robot mass
        mg                      % robot weight
        ha                      % ankle height in m
        he                      % max ankle height in air in m
        frequency               % discretization frequency
        tss                     % time duration of SSP in ms
        tss_divide              % divided duration of SSP in ms
        tds                     % time duration of DSP in ms
        tds_divide              % divided duration of DSP in ms
        tpi                     % time duration of starting phase in ms
        tpi_divide              % divided duration of starting phase in ms
        tpf                     % time duration of stopping phase in ms
        tpf_divide              % divided duration of stopping phase in ms
        tpassage                % time border of phases
        tpassage_ghost          % time border of phases with ghost phases
        nbphases                % number of phases
        discretization          % number of discrete point of each phases. !Warning!, during phase one, time=0 is not counted
        type_phase              % type of each phase, 0 means DSP, 1 means SSP after landing, 2 means SSP before take off
        nbpointdiscret          % number of discrete points
        dt_type_phase           % discretization of type_phase
        nbcontrolpointzmp       % number of optimization parameter linked to points ABCD
        nbcontrolpointzmp1      % number of optimization parameters linked to point B' of ZMP1 aka number of double support phases +2 (because ZMP1 is cut during starting and stopping phases), +6 because initial and final DSP are cut in two
        nbparamank              % number of ankle position
        nbparamtotal            % number of optimization parameters in one direction
        nbpankle                % number of foot positions on the floor
        step_number_pankle_fixed% fixed ankle positions
        psi_zmp                 % fixed ankle angles used for ZMP in SSP
        psi_zmp1                % fixed ankle angles used for ZMP1 in DSP
        psi_zmp2                % fixed ankle angles used for ZMP2 in DSP
        psi_pstep               % fixed ankle angles used for pstep
        lambda                  % weight of fcom
        mu                      % weight of ankle torques
        epsilon                 % weight of ZMP1&2 acceleration
        psa_abcdDSP             % all optimized parameters
        psa_abcd                % ABCD boundary conditions optimized parameters
        pstep                   % step positions
        pabcd                   % via-point positions
        e                       % real sole width without constraint
        pos_ankle_loc           % ankle position respect to sole frame
        rightorleft             % foot (left or right)
    end
    methods
        function obj = wpg_parameters(robot,type_traj,firstSS,frequency,lambda,epsilon,e,rightorleft,poly_degree,nbpolyssp,nbpolydsp,nbpolypi,nbpolypf)
            %% %%%choose the trajectory%%%
            obj.type_traj=type_traj;
            %1 : forward walking 0.5m, 4step
            %2 : forward walking 1m, 10step
            %3 : forward walking 3m, 10step
            %4 : forward walking 1.5m, 9step
            %5 : forward walking 1m, 10 steps, diagonal
            %6 : forward walking 1m, 10 steps, lateral
            %7 : forward walking 2m, 10 steps, 'human like'
            %8 : forward+diagonal+stretch, 30step
% %             %1 : forward walking 1m, 10step
% %             %2 : forward+diagonal+stretch, 30step
% %             %3 : forward walking 1.5m, 9step
% %             %4 : forward walking 3m, 10 steps, human like
% %             %5 : arc circle walking
% %             %6 : 1 step walking
% %             %7 : lateral walking
% %             %8 : arc circle little
            %%%%%%%%%%%
            
            obj.g = 9.81; % gravity constant en m/s²
            
            % sole width
            obj.e = e; %real sole width
            
            %%%define the number of foot step%%%
			obj.nbstep = 4;
            
            switch (obj.type_traj)
                case 1
                    obj.nbstep = 4;
                case 2
                    obj.nbstep = 10;
                case 3
                    obj.nbstep = 10;
                case 4
                    obj.nbstep = 9;
                case 5
                    obj.nbstep = 10;
                case 6
                    obj.nbstep = 10;
                case 7
                    obj.nbstep = 10;
                case 8
                    obj.nbstep = 30;
                case 9
                    obj.nbstep = 4;
                    
                otherwise
                    msg = 'choose a preset trajectory';
                    errormsg = [msg];
                    error(errormsg,[])
            end
            
            % foot (left or right)
            obj.rightorleft = rightorleft; %+1 for right foot and -1 for left foot
            
            % robot
            obj.robot_config(robot);
			
			obj.firstSS = firstSS;
			
            %%% weight of optimization criterion %%%
            obj.lambda = lambda; %weight of fcom
            obj.mu = 1-obj.lambda; %weight of ankle torques
            obj.epsilon = epsilon; %weight of ZMP1&2 acceleration
            
            % polynomial parameters
%             obj.poly_degree=6;
            obj.poly_degree=poly_degree;
            obj.nbpolyssp=nbpolyssp;
            obj.nbpolydsp=nbpolydsp;
            obj.nbpolypi=nbpolypi;
            obj.nbpolypf=nbpolypf;
            obj.nbpolyzmp1=obj.nbpolydsp;
			
			%%% discretization frequency %%%
            obj.frequency = frequency;
			
            % time configuration
            obj.time_config();
            
            % sole margin to keep stability
            obj.sole_margin = 0.02; 
        end
        function obj=robot_config(obj,robot)
            %%%choose the robot%%%
            if isnan(robot) || isempty(robot) || ~(any(robot==[1 2]))
                msg = 'choose a robot: \n';
                msg1 = '1: hrp2 \n';
                msg2 = '2: hrp4';
                errormsg = [msg msg0 msg1 msg2];
                error(errormsg,[])
            end           
            obj.robot = robot;
			
            switch(obj.robot)
                case 1
					%%% foot dimension %%%
					obj.backtoankle = 0.1; % from back to the ankle
                    obj.fronttoankle = 0.13; % from front to the ankle
                    obj.exttoankle = 0.075; % from exterior to the ankle
                    obj.inttoankle = 0.055; % from interior to the ankle
                    obj.ha = 0.12; % ankle height in m
                    obj.he = 0.2; % max ankle height during movement in m 
                    if obj.rightorleft==1
                        obj.pos_ankle_loc = [obj.backtoankle;obj.exttoankle-((obj.exttoankle+obj.inttoankle)/2);obj.e+obj.ha]; % ankle position in sole frame
                    else
                        obj.pos_ankle_loc = [obj.backtoankle;((obj.exttoankle+obj.inttoankle)/2)-obj.exttoankle;obj.e+obj.ha]; % ankle position in sole frame                        
                    end
                    % COM altitude in m
					obj.z = 0.808511;
					% robot mass in kg
					obj.m = 80;
					%%%foot step positions constraint%%%
                    obj.xankmax = 0.4; %stepping forward max
                    obj.xankmin = -0.4; %stepping forward min (if negative, it means stepping backward max)
                    obj.yankmin = obj.inttoankle+0.03; %width min between ankles
                    obj.yankmax = obj.inttoankle+0.4; %width max between ankles
					%%% define initial - final positions of right and left foot%%%	
                    obj.pankinit_firstSS = [0.0095 -obj.rightorleft*0.095]; %initial position of left foot. 
                    obj.pankinit_firstinair = [0.0095 obj.rightorleft*0.095]; %right foot is supposed to be symmetrical by x axis to the left foot with initial pose
                    switch (obj.type_traj)
                        case 1
                            obj.pankfin_lastSS = [0.5095 (obj.rightorleft)^obj.nbstep*0.095];
                            obj.pankfin_lastinair = [0.5095 -(obj.rightorleft)^obj.nbstep*0.095];
                        case 2
                            obj.pankfin_lastSS = [1.0095 (obj.rightorleft)^obj.nbstep*0.095];
                            obj.pankfin_lastinair = [1.0095 -(obj.rightorleft)^obj.nbstep*0.095];
                        case 3
                            obj.pankfin_lastSS = [2.5095 (obj.rightorleft)^obj.nbstep*0.095];
                            obj.pankfin_lastinair = [2.5095 -(obj.rightorleft)^obj.nbstep*0.095];
                        case 4
                            obj.pankfin_lastSS = [1.5095 (obj.rightorleft)^obj.nbstep*0.095];
                            obj.pankfin_lastinair = [1.5095 -(obj.rightorleft)^obj.nbstep*0.095];
                        case 5
                            obj.pankfin_lastSS = [1.0095 (obj.rightorleft)^obj.nbstep*0.095-1];
                            obj.pankfin_lastinair = [1.0095 -(obj.rightorleft)^obj.nbstep*0.095-1];
                        case 6
                            obj.pankfin_lastSS = [0.0095 (obj.rightorleft)^obj.nbstep*0.095-1];
                            obj.pankfin_lastinair = [0.0095 -(obj.rightorleft)^obj.nbstep*0.095-1];
                        case 7
                            obj.pankfin_lastSS = [1.0095 (obj.rightorleft)^obj.nbstep*0.095];
                            obj.pankfin_lastinair = [1.0095 -(obj.rightorleft)^obj.nbstep*0.095];
                        case 8
                            obj.pankfin_lastSS = [1.0095 (obj.rightorleft)^obj.nbstep*0.095-0.2];
                            obj.pankfin_lastinair = [1.0095 -(obj.rightorleft)^obj.nbstep*0.095-0.2];
                        case 9
                            obj.pankfin_lastSS = [0.5095 (obj.rightorleft)^obj.nbstep*0.095-0.2];
                            obj.pankfin_lastinair = [0.5095 -(obj.rightorleft)^obj.nbstep*0.095-0.2];

                        otherwise
                            msg = 'choose a preset trajectory';
                            errormsg = [msg];
                            error(errormsg,[])
                    end
                case 2
					%%% foot dimension %%%
                    obj.backtoankle = 0.098; %from back to the ankle
                    obj.fronttoankle = 0.128; %from front to the ankle
                    obj.exttoankle = 0.076; %from exterior to the ankle
                    obj.inttoankle = 0.054; %from interior to the ankle
					obj.ha = 0.093; %ankle height in m
                    obj.he = obj.ha+0.08; %max ankle height during movement in m
                    if obj.rightorleft==1
                        obj.pos_ankle_loc = [obj.backtoankle;obj.exttoankle-((obj.exttoankle+obj.inttoankle)/2);obj.e+obj.ha]; % ankle position in sole frame
                    else
                        obj.pos_ankle_loc = [obj.backtoankle;((obj.exttoankle+obj.inttoankle)/2)-obj.exttoankle;obj.e+obj.ha]; % ankle position in sole frame                        
                    end
                    %COM altitude in m
					obj.z = 0.781739;	
					%robot mass in kg
					obj.m = 43;
					%%%foot step positions constraint%%%	
                    obj.xankmax = 0.4; %stepping forward max
                    obj.xankmin = -0.4; %stepping forward min (if negative, it means stepping backward max)
                    obj.yankmin = obj.inttoankle+0.056; %width min between ankles
                    obj.yankmax = obj.inttoankle+0.4; %width max between ankles	
					%%% define initial - final positions of right and left foot%%%
                    obj.pankinit_firstSS = [0.0095 -obj.rightorleft*0.0815817]; %initial position of left foot. 
                    obj.pankinit_firstinair = [0.0095 obj.rightorleft*0.0815817]; %right foot is supposed to be symmetrical by x axis to the left foot with initial pose
                    switch (obj.type_traj)
                        case 1
                            obj.pankfin_lastSS = [0.5095 (obj.rightorleft)^obj.nbstep*0.0815817];
                            obj.pankfin_lastinair = [0.5095 -(obj.rightorleft)^obj.nbstep*0.0815817];
                        case 2
                            obj.pankfin_lastSS = [0.5095 (obj.rightorleft)^obj.nbstep*0.0815817];
                            obj.pankfin_lastinair = [0.5095 -(obj.rightorleft)^obj.nbstep*0.0815817];
                        case 3
                            obj.pankfin_lastSS = [2.5095 (obj.rightorleft)^obj.nbstep*0.0815817];
                            obj.pankfin_lastinair = [2.5095 -(obj.rightorleft)^obj.nbstep*0.0815817];
                        case 4
                            obj.pankfin_lastSS = [1.5095 (obj.rightorleft)^obj.nbstep*0.0815817];
                            obj.pankfin_lastinair = [1.5095 -(obj.rightorleft)^obj.nbstep*0.0815817];
                        case 5
                            obj.pankfin_lastSS = [1.0095 (obj.rightorleft)^obj.nbstep*0.0815817-1];
                            obj.pankfin_lastinair = [1.0095 -(obj.rightorleft)^obj.nbstep*0.0815817-1];
                        case 6
                            obj.pankfin_lastSS = [0.0095 (obj.rightorleft)^obj.nbstep*0.0815817-1];
                            obj.pankfin_lastinair = [0.0095 -(obj.rightorleft)^obj.nbstep*0.0815817-1];
                        case 7
                            obj.pankfin_lastSS = [2.0095 (obj.rightorleft)^obj.nbstep*0.0815817];
                            obj.pankfin_lastinair = [2.0095 -(obj.rightorleft)^obj.nbstep*0.0815817];
                        case 8
                            obj.pankfin_lastSS = [1.0095 (obj.rightorleft)^obj.nbstep*0.0815817-0.2];
                            obj.pankfin_lastinair = [1.0095 -(obj.rightorleft)^obj.nbstep*0.0815817-0.2];
                        case 9
                            obj.pankfin_lastSS = [0.5095 (obj.rightorleft)^obj.nbstep*0.0815817];
                            obj.pankfin_lastinair = [0.5095 -(obj.rightorleft)^obj.nbstep*0.0815817];
                        otherwise
                            msg = 'choose a preset trajectory';
                            errormsg = [msg];
                            error(errormsg,[])
                    end
            end
            obj.Tc = sqrt(obj.z/obj.g); %time constant
            obj.w = 1/obj.Tc;
			obj.mg = obj.m*obj.g; %robot weight in N

			%%% initial and final zmp %%%
            obj.xpsa_zmpinit = [(obj.pankinit_firstSS(1)+obj.pankinit_firstinair(1))/2+0.0095;0;0];    
			obj.xpsa_zmpfin = [(obj.pankfin_lastSS(1)+obj.pankfin_lastinair(1))/2+0.0095;0;0];
            obj.ypsa_zmpinit = [(obj.pankinit_firstSS(2)+obj.pankinit_firstinair(2))/2;0;0];           
			obj.ypsa_zmpfin = [(obj.pankfin_lastSS(2)+obj.pankfin_lastinair(2))/2;0;0];
			%%%initial and final com position%%%
            obj.xpcominit = (obj.pankinit_firstSS(1)+obj.pankinit_firstinair(1))/2+0.0095;
			obj.xpcomfin = (obj.pankfin_lastSS(1)+obj.pankfin_lastinair(1))/2+0.0095;
            obj.ypcominit = (obj.pankinit_firstSS(2)+obj.pankinit_firstinair(2))/2;
			obj.ypcomfin = (obj.pankfin_lastSS(2)+obj.pankfin_lastinair(2))/2;
            %%%initial and final com velocity%%%
            obj.xscominit = 0;
			obj.xscomfin = 0;
            obj.yscominit = 0;
			obj.yscomfin = 0;
			%initial and final position of zmp1
			obj.xpsa_zmp1init = obj.xpsa_zmpinit;
			obj.ypsa_zmp1init = [obj.pankinit_firstinair(2);0;0];
			obj.xpsa_zmp1fin = obj.xpsa_zmpfin;
			obj.ypsa_zmp1fin = [obj.pankfin_lastSS(2);0;0];		
        end
		
        function obj=time_config(obj)
            %%% time on passage points
            obj.tss = .90; %time duration of single support phases in ms %WARNING:  * frequency / nbpolyssp = positive integer
            obj.tss_divide = obj.tss/obj.nbpolyssp;
            obj.tds = .400; %time duration of double support phases in ms %WARNING:  * frequency / nbpolydsp = positive integer
%             obj.tds = 1.2;
            obj.tds_divide = obj.tds/obj.nbpolydsp;
            obj.tpi = .800; %time duration of starting phase in ms %WARNING:  * frequency / nbpolypi = positive integer
            obj.tpi_divide = obj.tpi/obj.nbpolypi;
            obj.tpf = .800; %time duration of stopping phase in ms %WARNING:  * frequency / nbpolypf = positive integer
            obj.tpf_divide = obj.tpf/obj.nbpolypf;

            %%% computation of phase time duration without starting and stopping%%%
            tpassagessp=[0];
            for i=1:obj.nbpolyssp
                tpassagessp(i+1)=tpassagessp(i)+obj.tss_divide;
            end
            tpassagedsp=[0];
            for i=1:obj.nbpolydsp
                tpassagedsp(i+1)=tpassagedsp(i)+obj.tds_divide;
            end
            tpassagepi=[0];
            for i=1:obj.nbpolypi
                tpassagepi(i+1)=tpassagepi(i)+obj.tpi_divide;
            end
            tpassagepf=[0];
            for i=1:obj.nbpolypf
                tpassagepf(i+1)=tpassagepf(i)+obj.tpf_divide;
            end
            
            obj.tpassage = [0];
            for i=1:obj.nbstep
                obj.tpassage=[obj.tpassage obj.tpassage(end)+tpassagessp(2:end)];
                
                obj.tpassage=[obj.tpassage obj.tpassage(end)+tpassagedsp(2:end)];
            end
            obj.tpassage=[obj.tpassage obj.tpassage(end)+tpassagessp(2:end)];
            
            %%% add starting and stopping duration%%%
            obj.tpassage = [tpassagepi obj.tpassage(2:end)+tpassagepi(end) obj.tpassage(end)+tpassagepi(end)+tpassagepf(2:end)];
            
            %%% add starting and stopping ghost control point time %%%
            tpassagepi_=[0];
            for i=1:obj.poly_degree
                tpassagepi_(i+1)=tpassagepi_(i)-obj.tpi;
            end
            tpassagepf_=[0];
            for i=1:obj.poly_degree
                tpassagepf_(i+1)=tpassagepf_(i)+obj.tpf;
            end
            obj.tpassage_ghost = [tpassagepi_(end:-1:2) obj.tpassage obj.tpassage(end)+tpassagepf_(2:end)];
            %%%%%%%
            
            %%% number of phases in the whole walk%%%
            obj.nbphases = length(obj.tpassage)-1;
			
			obj.compute_dicretization_type_phase();            
        end

        function obj=compute_dicretization_type_phase(obj)
            %%computation of the number of points during each phase and the type of phase:%%%
            %0 means DSP
            %!=0 means SSP
            obj.discretization = zeros(1,obj.nbphases);
            for i=1:obj.nbphases
                obj.discretization(i) = size(obj.tpassage(i):1/obj.frequency:obj.tpassage(i+1)-1/obj.frequency,2);
            end
            
            obj.type_phase = zeros(1,obj.nbpolypi);
            obj.type_phase = [obj.type_phase repmat([1:obj.nbpolyssp zeros(1,obj.nbpolydsp)],1,obj.nbstep)];
            obj.type_phase = [obj.type_phase 1:obj.nbpolyssp zeros(1,obj.nbpolypf)];
            %%%%%%%
            
            %%%number of discretization points%%%
            obj.nbpointdiscret = sum(obj.discretization)+1;
            %%%%%%%

            %%%discretize the type of phase%%%
            obj.dt_type_phase = obj.compute_dt_type_phase(obj.type_phase,obj.discretization,obj.nbphases,obj.nbpointdiscret);
            %%%%%%%

            %%%constants%%%
            %number of optimization parameter linked to points ABCD
            obj.nbcontrolpointzmp = (obj.nbstep+1)*(obj.nbpolyssp)+(obj.nbstep)*(obj.nbpolydsp)+obj.nbpolypi+obj.nbpolypf+obj.poly_degree;
            %number of optimization parameters linked to point B' of ZMP1 are the number of double support phases +2 (because ZMP1 is cut during starting and stopping phases)
            obj.nbcontrolpointzmp1 = (obj.nbstep)*(obj.nbpolydsp+obj.poly_degree)+2*obj.poly_degree+obj.nbpolypi+obj.nbpolypf;
            %number of ankle position we are looking at
            obj.nbparamank = obj.nbstep-1;
            %number of foot positions on the floor
            obj.nbpankle = obj.nbstep+3;
            %number of optimization parameters in one direction
            obj.nbparamtotal = obj.nbcontrolpointzmp + obj.nbcontrolpointzmp1 + obj.nbparamank +4;%+4 init and final COM conditions

            %%%fixed ankle position%%%
            obj.step_number_pankle_fixed = [];
            switch (obj.type_traj)
                case 1
                    obj.step_number_pankle_fixed = [1 obj.pankinit_firstinair;
                                                    2 obj.pankinit_firstSS;
                                                    3 0.125 -0.0815817;
                                                    4 0.250 0.0815817;
                                                    5 0.375 -0.0815817;
                                                    obj.nbpankle-1 obj.pankfin_lastSS;
                                                    obj.nbpankle obj.pankfin_lastinair];
                case 8
                     obj.step_number_pankle_fixed=[1 obj.pankinit_firstinair;
                                                   2 obj.pankinit_firstSS;
                                                   8 0.5 -obj.firstSS*0.095;
                                                   19 2 -obj.firstSS*0.75;
                                                   25 2 -obj.firstSS*0.2
                                                   obj.nbpankle-1 obj.pankfin_lastSS;
                                                   obj.nbpankle obj.pankfin_lastinair]; 
                otherwise
                    obj.step_number_pankle_fixed = [1 obj.pankinit_firstinair;
                                                    2 obj.pankinit_firstSS;
                                                    obj.nbpankle-1 obj.pankfin_lastSS;
                                                    obj.nbpankle obj.pankfin_lastinair]; 
            end
			
            %%%fixed ankle angular position%%%
            obj.psi_zmp=zeros(1,length(obj.discretization)+1);
            obj.psi_zmp1=zeros(1,length(obj.discretization)+1);
            obj.psi_zmp2=zeros(1,length(obj.discretization)+1);
            obj.psi_pstep=zeros(1,obj.nbpankle);
            switch (obj.type_traj)
                case 2
%                     obj.psi_zmp(31:36)=pi/2;
%                     obj.psi_zmp1(37:44)=pi/2;
%                     obj.psi_zmp2(23:30)=pi/2;
%                     obj.psi_pstep(4)=pi/2;
                case 7
                    for i=3:obj.nbpankle-2
                        obj.psi_pstep(i)=obj.firstSS^(i)*pi/12;
                        
                        n=obj.nbpolypi+(i-2)*(obj.nbpolyssp+obj.nbpolydsp);
                        tssp=ones(1,obj.nbpolyssp)*obj.firstSS^(i)*pi/12;
                        tdsp=ones(1,obj.nbpolydsp)*obj.firstSS^(i)*pi/12;
                        obj.psi_zmp(n+1:n+obj.nbpolyssp)=tssp;
                        obj.psi_zmp1(n+obj.nbpolyssp+1:n+obj.nbpolyssp+obj.nbpolydsp)=tdsp;
                        obj.psi_zmp2(n-obj.nbpolydsp+1:n)=tdsp;
                    end
                case 8
                    for i=6:17
                        obj.psi_pstep(i)=-obj.firstSS*(pi/4-abs((i-6-6)*(pi/4/(17-6)*2)));

                        n=obj.nbpolypi+(i-2)*(obj.nbpolyssp+obj.nbpolydsp);
                        tssp=-obj.firstSS*ones(1,obj.nbpolyssp)*(pi/4-abs((i-6-6)*(pi/4/(17-6)*2)));
                        tdsp=-obj.firstSS*ones(1,obj.nbpolydsp)*(pi/4-abs((i-6-6)*(pi/4/(17-6)*2)));
                        obj.psi_zmp(n+1:n+obj.nbpolyssp)=tssp;
                        obj.psi_zmp1(n+obj.nbpolyssp+1:n+obj.nbpolyssp+obj.nbpolydsp)=tdsp;
                        obj.psi_zmp2(n-obj.nbpolydsp+1:n)=tdsp;
                    end
                    for i=23:obj.nbstep-1
                        obj.psi_pstep(i)=-obj.firstSS*(pi/4-abs((i-23-3)*(pi/4/(obj.nbstep-1-23)*2)));

                        n=obj.nbpolypi+(i-2)*(obj.nbpolyssp+obj.nbpolydsp);
                        tssp=-obj.firstSS*ones(1,obj.nbpolyssp)*(pi/4-abs((i-23-3)*(pi/4/(obj.nbstep-1-23)*2)));
                        tdsp=-obj.firstSS*ones(1,obj.nbpolydsp)*(pi/4-abs((i-23-3)*(pi/4/(obj.nbstep-1-23)*2)));
                        obj.psi_zmp(n+1:n+obj.nbpolyssp)=tssp;
                        obj.psi_zmp1(n+obj.nbpolyssp+1:n+obj.nbpolyssp+obj.nbpolydsp)=tdsp;
                        obj.psi_zmp2(n-obj.nbpolydsp+1:n)=tdsp;
                    end
            end
        end

    end
    methods (Static)
        function dt_type_phase=compute_dt_type_phase(type_phase,discretization,nbphases,nbpointdiscret)
            dt_type_phase=zeros(nbpointdiscret,size(type_phase,1));
            nb=0;
            for j=1:nbphases
                dt_type_phase(nb+1:nb+discretization(:,j),:)=type_phase(:,j)';
                nb=nb+discretization(:,j);
            end

            dt_type_phase(nbpointdiscret,1)=type_phase(:,end);
        end
    end
end