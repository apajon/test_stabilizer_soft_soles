classdef wpg_trajectories<handle

    properties
        xpzmp       % zmp positions on x axis
        ypzmp       % zmp positions on y axis
        xszmp       % zmp speed on x axis
        yszmp       % zmp speed on y axis
        xazmp       % zmp acceleration on x axis
        yazmp       % zmp acceleration on y axis
        xpzmp1      % zmp1 positions on x axis
        ypzmp1      % zmp1 positions on y axis
        xpzmp2      % zmp2 positions on x axis
        ypzmp2      % zmp2 positions on y axis
        xpcom       % com position on x axis
        ypcom       % com position on y axis
        xpankle_l   % left ankle positions on x axis
        ypankle_l   % left ankle positions on y axis
        zpankle_l   % left ankle positions on z axis
        xpankle_r   % right ankle positions on x axis
        ypankle_r   % right ankle positions on y axis
        zpankle_r   % right ankle positions on z axis
        rpsi_l      % left ankle orientation around z axis
        rphi_l      % left ankle orientation around y axis
        rtheta_l    % left ankle orientation around x axis
        rpsi_r      % right ankle orientation around z axis
        rphi_r      % right ankle orientation around y axis
        rtheta_r    % right ankle orientation around x axis
    end
    methods
        function obj=wpg_trajectories()
        end
        
        function obj=drawing(obj,walking_param,zmp,display,figure_number)
            %draw ZMP, COM, ZMP1, ZMP2, foot step position in walking_param
            %display :
            %0 : don't display trajectories
            %1 : display trajectories
            % load('trajectories.mat')
            walking_param.nbparamtotal=walking_param.nbcontrolpointzmp+walking_param.nbpankle+walking_param.nbcontrolpointzmp1;

            walking_param.psa_abcd=[walking_param.psa_abcdDSP(1:walking_param.nbcontrolpointzmp);walking_param.psa_abcdDSP(walking_param.nbparamtotal+1:walking_param.nbparamtotal+walking_param.nbcontrolpointzmp)];
%             walking_param.pstep=[walking_param.psa_abcdDSP(walking_param.nbcontrolpointzmp+1:walking_param.nbcontrolpointzmp+walking_param.nbpankle) walking_param.psa_abcdDSP(walking_param.nbparamtotal+walking_param.nbcontrolpointzmp+1:walking_param.nbparamtotal+walking_param.nbcontrolpointzmp+walking_param.nbpankle)];
            walking_param.pstep=[walking_param.psa_abcdDSP(walking_param.nbcontrolpointzmp+1:walking_param.nbcontrolpointzmp+walking_param.nbpankle) walking_param.psa_abcdDSP(walking_param.nbparamtotal+walking_param.nbcontrolpointzmp+1:walking_param.nbparamtotal+walking_param.nbcontrolpointzmp+walking_param.nbpankle)];
%             walking_param.pstep=[walking_param.pankinit_firstinair;
%                 walking_param.pankinit_firstSS;
%                 walking_param.pstep;
%                 walking_param.pankfin_lastSS;
%                 walking_param.pankfin_lastinair];
            if display
                figure_num=figure_number;
                %%%Clear figure(3) to draw y(x) coordinates%%%
                figure(figure_num);
                clf;
                set(gca,'fontsize',14)
                axis equal;
                title('ZMP projected on the ground')
                xlabel('x(m)')
                ylabel('y(m)')
                hold on;
                plot(walking_param.pstep(:,1),walking_param.pstep(:,2),'+')
                hold off
            end
            %%%using to draw zmp and com with new ZMP generator%%%

            obj.xpzmp=zmp.A_xzmp*walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)+zmp.B_xzmp;
            obj.ypzmp=zmp.A_yzmp*walking_param.psa_abcd(length(walking_param.psa_abcd)/2+1:end)+zmp.B_yzmp;
            obj.xszmp=zmp.A_xzmp_spd*walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)+zmp.B_xzmp_spd;
            obj.yszmp=zmp.A_yzmp_spd*walking_param.psa_abcd(length(walking_param.psa_abcd)/2+1:end)+zmp.B_yzmp_spd;
            obj.xazmp=zmp.A_xzmp_acc*walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)+zmp.B_xzmp_acc;
            obj.yazmp=zmp.A_yzmp_acc*walking_param.psa_abcd(length(walking_param.psa_abcd)/2+1:end)+zmp.B_yzmp_acc;



            obj.xpzmp1=zmp.A_xzmp1*walking_param.psa_abcdDSP(1:walking_param.nbparamtotal)+zmp.B_xzmp1;
            obj.ypzmp1=zmp.A_yzmp1*walking_param.psa_abcdDSP(walking_param.nbparamtotal+1:walking_param.nbparamtotal+walking_param.nbparamtotal)+zmp.B_yzmp1;


            obj.xpzmp2=zmp.A_xzmp2*walking_param.psa_abcdDSP(1:walking_param.nbparamtotal)+zmp.B_xzmp2;
            obj.ypzmp2=zmp.A_yzmp2*walking_param.psa_abcdDSP(walking_param.nbparamtotal+1:walking_param.nbparamtotal+walking_param.nbparamtotal)+zmp.B_yzmp2;



            % xpzmp=trj(:,1);
            % ypzmp=trj(:,2);

            % %%%drawing xzmp(t)%%%
            % figure(6)
            % clf 
            % axis auto
            % title('x(t) of ZMP and COM (optimization on COM work)')
            % xlabel('t(s)')
            % ylabel('x(m)')
            % hold on
            % plot(xtimezmp,xpzmp)
            % hold off
            % %%%%%%%
            % 
            % %%%drawing yzmp(t)%%%
            % figure(7)
            % clf 
            % axis auto
            % title('y(t) of ZMP and COM (optimization on COM work)')
            % xlabel('t[s]')
            % ylabel('y[m]')
            % hold on
            % plot(ytimezmp,ypzmp)
            % hold off
            % %%%%%%%

            if display
                %%%drawing ZMP positions%%%
                figure(figure_num);
                hold on;
                plot(obj.xpzmp,obj.ypzmp,'-b','LineWidth',2)
%                 plot(obj.xpzmp(any(walking_param.dt_type_phase==0,2)),obj.ypzmp(any(walking_param.dt_type_phase==0,2)),'*g','LineWidth',2);
                % plot(xpzmp1,ypzmp1,'-*r')
                % plot(xpzmp2,ypzmp2,'-*m')
                plot(obj.xpzmp1,obj.ypzmp1,'*r')
                plot(obj.xpzmp2,obj.ypzmp2,'*m')
                legend('viapoints','ZMP')
                hold off
                %%%%%%%%%%%
            end

            %%%COM computation%%%
            obj.xpcom=zmp.A_xcom*walking_param.psa_abcd(1:end/2)+zmp.B_xcom;
            obj.ypcom=zmp.A_ycom*walking_param.psa_abcd(end/2+1:end)+zmp.B_ycom;
            %%%%%%%%%%%

            if display
                %%%drawing COM positions%%%
                figure(figure_num)
                hold on
                plot(obj.xpcom,obj.ypcom,'green','LineWidth',2);
                % hleg = legend('ankle position','Via points','ZMP','COM','Location','NorthEast');
                hleg = legend('ankle position','ZMP','ZMP1','ZMP2','COM','Location','EastOutside');
                % Make the text of the legend italic and color it brown
                set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
                hold off
                %%%%%%%%%%%
                % figure(6)
                % hold on
                % plot(time,xpcom,'green');
                % hleg = legend('ZMP','COM','Location','Northwest');
                % % Make the text of the legend italic and color it brown
                % set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
                % hold off
                % figure(7)
                % hold on
                % plot(time,ypcom,'green');
                % hleg = legend('ZMP','COM','Location','Northwest');
                % % Make the text of the legend italic and color it brown
                % set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
                % hold off
                % %%%%%%
                %%%%%%
            end

            if display
                %%%drawing foot steps%%%
                figure(figure_num)
                hold on;
                XY=drawing_rectangle_rotate(walking_param.pstep,[walking_param.psi_pstep 0],walking_param.backtoankle,walking_param.fronttoankle,walking_param.exttoankle,walking_param.inttoankle,walking_param.firstSS);
                for i=1:length(walking_param.pstep)
                    plot(XY(i,1:5),XY(i,6:10),'-k','LineWidth',2)
                end
                
                XY=drawing_rectangle_rotate(walking_param.pstep,[walking_param.psi_pstep 0],walking_param.backtoankle-walking_param.sole_margin,walking_param.fronttoankle-walking_param.sole_margin,walking_param.exttoankle-walking_param.sole_margin,walking_param.inttoankle-walking_param.sole_margin,walking_param.firstSS);
                for i=1:length(walking_param.pstep)
%                 for i=6:7
                    plot(XY(i,1:5),XY(i,6:10),':k','LineWidth',2)
                end
                hold off
            end
        end
    end
end