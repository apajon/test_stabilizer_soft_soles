classdef wpg_qp_problem_nl<handle

    properties
        H
        f
        H_nl
        f_nl
        G
        A
        b
        Aeq
        beq
        C
    end
    methods
        function obj=wpg_qp_problem_nl(walking_param,gradients)
            
            obj.qp_generating_SSP(walking_param,gradients);

            %%
            %%% modifying some number of parameters to extract one foot
            % step optimization. never look at it
            %use to adapt the WPg for a 1 step optimization
            walking_param.nbparamBb = (walking_param.nbstep+2+1)*3;
            walking_param.nbparamBb = walking_param.nbparamBb-15-3;
            walking_param.nbparamABCD = (walking_param.nbstep*3+2+1)*3-6;
            Cssp = [10:walking_param.nbparamABCD-9-9 walking_param.nbparamABCD+1:walking_param.nbparamABCD+4];
            T = sum(walking_param.discretization(1:4))+2:sum(walking_param.discretization(1:7))+1;
            
            %%
            %%% compute the matrix of time discretization %%%
            M_t=obj.compute_M_t(walking_param.discretization,walking_param.frequency,walking_param.nbphases,walking_param.nbpointdiscret);
            M_t=M_t(T,:);
            
            obj.compute_ZMP(walking_param,gradients,M_t,Cssp);
            obj.compute_COM(walking_param,gradients,M_t,Cssp,T);
            obj.compute_cost_force(walking_param,gradients);        
            obj.compute_cost_torque_ssp(walking_param,gradients,T);
            obj.compute_cons_zmp(walking_param,gradients,T);
            
            %%
            obj.qp_generating_DSP(walking_param,gradients);
            
            %%
            %%% modifying some number of parameters to extract one foot
            % step optimization. never look at it
            %use to adapt the WPG for a 1 step optimization
            walking_param.nbparamABCD=(walking_param.nbstep*3+2+1)*3-6;
            Cssp=[10:walking_param.nbparamABCD-9-9 walking_param.nbparamABCD+1:walking_param.nbparamABCD+4];
            walking_param.nbparamABCD=walking_param.nbparamABCD-18-9;
            walking_param.nbparamBb=(walking_param.nbstep+2+1)*3;
            L=sum(walking_param.discretization([1 3]))+2:sum(walking_param.discretization([1 3 6]))+1;
            Cdsp=[Cssp Cssp(end)+9+1:Cssp(end)+9+walking_param.nbparamank walking_param.nbparamtotal+4-walking_param.nbparamBb+10:walking_param.nbparamtotal+4-6-3];
            walking_param.nbparamBb=walking_param.nbparamBb-15-3;
            T=sum(walking_param.discretization(1:4))+2:sum(walking_param.discretization(1:7))+1;
            
            %%
            %where to cut the beginning and ending of zmp1 in percentage
            cutting=1/2;
            %%% compute the matrix of time discretization for ZMP1%%%
            %change from M_t because durig the initial and final phases
            %ZMP1 is cut into 2 polynomials
            M_t1=obj.compute_M_t1(walking_param.discretization,walking_param.frequency,cutting);
            M_t1=M_t1(any(walking_param.dt_type_phase==0,2),:);%suppress SSP discretization steps because ZMP1&2 are only defined in DSP
            M_t1=M_t1(L,:);
            
            %matrix of force repartition in DSP
            gradients.Ac=obj.compute_force_repartition(walking_param.tpassage,cutting);
            gradients.k_diag=diag(M_t1*gradients.Ac);
            
            obj.compute_ZMP1(walking_param,gradients,M_t1,cutting,Cdsp);
            obj.compute_cost_torque_zmp1(walking_param,gradients,T,L);
            obj.compute_cost_acc_zmp1(walking_param,gradients,cutting,L,Cdsp);
            obj.compute_ZMP2(walking_param,gradients,T);
            obj.compute_cost_torque_zmp2(walking_param,gradients,T,L);
            obj.compute_cost_acc_zmp2(walking_param,gradients,cutting,Cssp,T,L,Cdsp);
            obj.compute_cons_zmp12(walking_param,gradients,L);
            
            %%
            obj.A=gradients.AconstraintDSP;
            obj.b=gradients.BconstraintDSP;
           
            obj.qp_generating_problem_matrix(walking_param,gradients);
            
            walking_param.nbparamtotal=walking_param.nbparamABCD+walking_param.nbparamank+walking_param.nbparamBb;

            obj.C=walking_param.lambda*(gradients.xBfcom'*gradients.xBfcom+gradients.yBfcom'*gradients.yBfcom) ...
                +(1-walking_param.lambda)*(gradients.xBtorquessp'*gradients.xBtorquessp+gradients.yBtorquessp'*gradients.yBtorquessp) ...
                +walking_param.epsilon*(gradients.xBazmp1'*gradients.xBazmp1+gradients.yBazmp1'*gradients.yBazmp1+gradients.xBazmp2'*gradients.xBazmp2+gradients.yBazmp2'*gradients.yBazmp2) ...
                +0.5*(1-walking_param.lambda)*(gradients.xBtorquedsp1'*gradients.xBtorquedsp1+gradients.yBtorquedsp1'*gradients.yBtorquedsp1+gradients.xBtorquedsp2'*gradients.xBtorquedsp2 ...
                +gradients.yBtorquedsp2'*gradients.yBtorquedsp2);

            obj.compute_H_G_C(walking_param,gradients);
            obj.compute_A_b(walking_param,gradients);
            obj.compute_Aeq_beq(walking_param,gradients);
            
            
            %%
            walking_param.nbparamABCD = walking_param.nbparamABCD+4;
        end
        
        function obj=qp_generating_SSP(obj,walking_param,gradients)
            %% % modifying some number of parameters to extract one foot
            % step optimization. never look at it
            walking_param.nbparamBb = (walking_param.nbstep+2+1)*3;
            walking_param.nbparamBb = walking_param.nbparamBb-15-3;
            walking_param.nbparamABCD = (walking_param.nbstep*3+2+1)*3-6;
            Cssp = [10:walking_param.nbparamABCD-9-9 walking_param.nbparamABCD+1:walking_param.nbparamABCD+4];
            T = sum(walking_param.discretization(1:4))+2:sum(walking_param.discretization(1:7))+1;
            
            %% %%% compute the matrix of time discretization
            dt=obj.compute_coeff_dt_matrix(walking_param.discretization,walking_param.frequency,walking_param.nbphases,walking_param.nbpointdiscret);
            dt=dt(T,:);

            %% %matrix of sinh(wt) and cosh(wt)
            scwt = obj.compute_coeff_swt_cwt_matrix(walking_param.discretization,walking_param.w,walking_param.frequency,walking_param.nbphases);
            scwt=scwt(T,9:end-8-6);

            %% %coeff gradient computation of zmp_poly
            [gradients.xAzmp_gradient gradients.xBzmp_gradient]=obj.compute_coeff_zmp_gradient(walking_param.tpassage,walking_param.xpsa_zmpinit,walking_param.xpsa_zmpfin,walking_param.nbphases,walking_param.nbparamABCD);
            [gradients.yAzmp_gradient gradients.yBzmp_gradient]=obj.compute_coeff_zmp_gradient(walking_param.tpassage,walking_param.ypsa_zmpinit,walking_param.ypsa_zmpfin,walking_param.nbphases,walking_param.nbparamABCD);
            %% %A coeff of COM poly
            % A_gradient=com_morisawa_A_gradient(length(walking_param.tpassage)-1,walking_param.w);
            A_gradient=obj.com_morisawa_A_gradient(length(walking_param.tpassage)-1-8-3,walking_param.w);

            %% %VW gradient
            [xAVW xBVW]=obj.com_morisawa_VW_gradient(gradients.xAzmp_gradient(25:end-24-18,Cssp),gradients.xBzmp_gradient(25:end-24-18),A_gradient,walking_param.tpassage(5:end-7),walking_param.w);
            [yAVW yBVW]=obj.com_morisawa_VW_gradient(gradients.yAzmp_gradient(25:end-24-18,Cssp),gradients.yBzmp_gradient(25:end-24-18),A_gradient,walking_param.tpassage(5:end-7),walking_param.w);

            %% %position of ZMP and COM in x-axis
            [gradients.xApzmp gradients.xBpzmp] = obj.compute_traj_discrete(gradients.xAzmp_gradient(:,Cssp),gradients.xBzmp_gradient,dt);
            [gradients.xApcom gradients.xBpcom] = obj.pcom_generator_morisawa_gradient(gradients.xAzmp_gradient(25:end-24-18,Cssp),gradients.xBzmp_gradient(25:end-24-18),xAVW,xBVW,A_gradient,dt(:,25:end-24-18),scwt);
            %% %position of ZMP and COM in y-axis
            [gradients.yApzmp gradients.yBpzmp] = obj.compute_traj_discrete(gradients.yAzmp_gradient(:,Cssp),gradients.yBzmp_gradient,dt);
            [gradients.yApcom gradients.yBpcom] = obj.pcom_generator_morisawa_gradient(gradients.yAzmp_gradient(25:end-24-18,Cssp),gradients.yBzmp_gradient(25:end-24-18),yAVW,yBVW,A_gradient,dt(:,25:end-24-18),scwt);

            %% %force on com in x and y axis
            [gradients.xAfcom gradients.xBfcom] = obj.fcom_gradient(gradients.xApzmp,gradients.xBpzmp,gradients.xApcom,gradients.xBpcom,walking_param.mg,walking_param.z);
            [gradients.yAfcom gradients.yBfcom] = obj.fcom_gradient(gradients.yApzmp,gradients.yBpzmp,gradients.yApcom,gradients.yBpcom,walking_param.mg,walking_param.z);

            %% %force2 on com
            [gradients.xAfcom2 gradients.xBfcom2] = obj.compute_quad_matrix(gradients.xAfcom,gradients.xBfcom);
            [gradients.yAfcom2 gradients.yBfcom2] = obj.compute_quad_matrix(gradients.yAfcom,gradients.yBfcom);

            %% %torque in ankle
            [xAtorquessp gradients.xBtorquessp] = obj.torque_ankle_SSP_gradient(gradients.xApzmp,gradients.xBpzmp,gradients.xAfcom,gradients.xBfcom,walking_param.pankinit_firstSS(1),walking_param.pankfin_lastSS(1),walking_param.discretization,walking_param.mg,walking_param.ha,T);
            [yAtorquessp gradients.yBtorquessp] = obj.torque_ankle_SSP_gradient(gradients.yApzmp,gradients.yBpzmp,gradients.yAfcom,gradients.yBfcom,walking_param.pankinit_firstSS(2),walking_param.pankfin_lastSS(2),walking_param.discretization,walking_param.mg,walking_param.ha,T);

            %% %torque 2 in ankle
            [gradients.xAtorque2ssp gradients.xBtorque2ssp] = obj.compute_quad_matrix(xAtorquessp,gradients.xBtorquessp);
            [gradients.yAtorque2ssp gradients.yBtorque2ssp] = obj.compute_quad_matrix(yAtorquessp,gradients.yBtorquessp);

            %% %constraint ineq in SSP
            [xApankle xBpankle]=obj.torque_ankle_positions_SSP(walking_param.pankinit_firstSS(1),walking_param.pankfin_lastSS(1),walking_param.discretization);
            [yApankle yBpankle]=obj.torque_ankle_positions_SSP(walking_param.pankinit_firstSS(2),walking_param.pankfin_lastSS(2),walking_param.discretization);
            xApankle=xApankle(T,:);
            xBpankle=xBpankle(T,:);
            yApankle=yApankle(T,:);
            yBpankle=yBpankle(T,:);

            [gradients.Apzmpconstraintssp gradients.Bpzmpconstraintssp]=obj.zmp_constraint_stability_SSP_xy(gradients.xApzmp,gradients.xBpzmp,gradients.yApzmp,gradients.yBpzmp,xApankle,xBpankle,yApankle,yBpankle,walking_param.discretization(5:end-7),walking_param.backtoankle,walking_param.fronttoankle,walking_param.exttoankle,walking_param.inttoankle,walking_param.sole_margin,walking_param.psi,walking_param.type_phase(5:end-7),walking_param.nbparamBb,walking_param.firstSS);
            %% %Constraint eq for com speed initial and final
            %Warning : A*P+B=0 => put -B in constraint in optimized function
            [gradients.xAscomeq xBscomeq] = obj.scom_eqcontraint(gradients.xAzmp_gradient(25:end-24-18,Cssp),gradients.xBzmp_gradient(25:end-24-18),xAVW,xBVW,A_gradient,walking_param.tpassage(5:end-7),3,walking_param.w);
            [gradients.yAscomeq yBscomeq] = obj.scom_eqcontraint(gradients.yAzmp_gradient(25:end-24-18,Cssp),gradients.yBzmp_gradient(25:end-24-18),yAVW,yBVW,A_gradient,walking_param.tpassage(5:end-7),3,walking_param.w);
            %% %%%concatenation%%%
            gradients.Bscomeq=[xBscomeq;yBscomeq];
            %%%%%%%
        end
        
        function obj=compute_ZMP(obj,walking_param,gradients,M_t,Cssp)
            %The ZMP is computed from : zmp=A_zmp*X+B_zmp
            %A_zmp=M_t*A_Ca
            %B_zmp=M_t*B_Ca
            %A_Ca*X+B_Ca gives the zmp polynomial coefficients
            %In the case of this simplified version of the code, B_zmp and
            %B_Ca are null. Thus, A_Ca is equivalent to Ca.
            %% %compute matrix Ca
            %Ca is the matrix which multiplied by zmp via-points boundary
            %conditions gives a vector with the zmp polynomial coefficients
            [gradients.A_xCa gradients.B_xCa]=obj.compute_zmp_Ca(walking_param.tpassage,walking_param.xpsa_zmpinit,walking_param.xpsa_zmpfin,walking_param.nbphases,walking_param.nbparamABCD);
            [gradients.A_yCa gradients.B_yCa]=obj.compute_zmp_Ca(walking_param.tpassage,walking_param.ypsa_zmpinit,walking_param.ypsa_zmpfin,walking_param.nbphases,walking_param.nbparamABCD);
            %% %compute matrix A_zmp and B_zmp
            [gradients.A_xzmp gradients.B_xzmp] = obj.compute_traj_discrete(gradients.A_xCa(:,Cssp),gradients.B_xCa,M_t);
            [gradients.A_yzmp gradients.B_yzmp] = obj.compute_traj_discrete(gradients.A_yCa(:,Cssp),gradients.B_yCa,M_t);
        end
        function obj=compute_COM(obj,walking_param,gradients,M_t,Cssp,T)
            %the COM is computed from : com=A_com*X+B_com
            %com=M_cs*y+M_t*l
            %com=M_cs*(A_y*X+B_y)+M_t*(A_l*X+B_l)
            %A_y=G^-1*[A_A.A_Ca N 0]
            %B_y=G^-1*[A_A.B_Ca]
            %A_l=[A_A.A_Ca 0]
            %B_l=[A_A.B_Ca]
            %In the case of this simplified version of the code, B_Ca is
            %null. Thus, B_y and B_l are null and A_y and A_l are
            %equivalent to y and l.
            %% %compute matrix G composed with sinh(wt) and cosh(wt)
            M_cs = obj.compute_com_M_cs(walking_param.discretization,walking_param.w,walking_param.frequency,walking_param.nbphases);
            M_cs=M_cs(T,9:end-8-6);

            %% %A coeff of COM poly
            % A_gradient=com_morisawa_A_gradient(length(walking_param.tpassage)-1,walking_param.w);
            A_A=obj.compute_com_A_A(length(walking_param.tpassage)-1-8-3,walking_param.w);
            %% %VW gradient
            [A_xy B_xy]=obj.compute_com_y(gradients.A_xCa(25:end-24-18,Cssp),gradients.B_xCa(25:end-24-18),A_A,walking_param.tpassage(5:end-7),walking_param.w);
            [A_yy B_yy]=obj.compute_com_y(gradients.A_yCa(25:end-24-18,Cssp),gradients.B_yCa(25:end-24-18),A_A,walking_param.tpassage(5:end-7),walking_param.w);
            %% %Compute matrix A_com and B_com
            [gradients.A_xcom gradients.B_xcom]=obj.compute_com(gradients.A_xCa(25:end-24-18,Cssp),gradients.B_xCa(25:end-24-18),A_xy,B_xy,A_A,M_t(:,25:end-24-18),M_cs);
            [gradients.A_ycom gradients.B_ycom]=obj.compute_com(gradients.A_yCa(25:end-24-18,Cssp),gradients.B_yCa(25:end-24-18),A_yy,B_yy,A_A,M_t(:,25:end-24-18),M_cs);
            
            %% %Constraint eq for com speed initial and final
            %These constraints matrix are defined here because they are a
            %part of COM definition
            %Warning : A*P+B=0 => put -B in constraint in optimized function
            [gradients.A_xcons_scom gradients.B_xcons_scom] = obj.cons_scom_initfin(gradients.A_xCa(25:end-24-18,Cssp),gradients.B_xCa(25:end-24-18),A_xy,B_xy,A_A,walking_param.tpassage(5:end-7),3,walking_param.w);
            [gradients.A_ycons_scom gradients.B_ycons_scom] = obj.cons_scom_initfin(gradients.A_yCa(25:end-24-18,Cssp),gradients.B_yCa(25:end-24-18),A_yy,B_yy,A_A,walking_param.tpassage(5:end-7),3,walking_param.w);
%             %% %%%concatenation%%%
            gradients.B_cons_scom=[gradients.B_xcons_scom;gradients.B_ycons_scom];
            gradients.A_cons_scom=[gradients.A_xcons_scom                 zeros(size(gradients.A_xcons_scom,1),walking_param.nbparamank) zeros(size(gradients.A_xcons_scom,1),walking_param.nbparamBb) zeros(size(gradients.A_ycons_scom)) zeros(size(gradients.A_xcons_scom,1),walking_param.nbparamank) zeros(size(gradients.A_xcons_scom,1),walking_param.nbparamBb);
                                  zeros(size(gradients.A_xcons_scom))    zeros(size(gradients.A_ycons_scom,1),walking_param.nbparamank) zeros(size(gradients.A_ycons_scom,1),walking_param.nbparamBb) gradients.A_ycons_scom              zeros(size(gradients.A_ycons_scom,1),walking_param.nbparamank) zeros(size(gradients.A_ycons_scom,1),walking_param.nbparamBb)];
        end
        function obj=compute_cost_force(obj,walking_param,gradients)
            %% %compute force on com in x and y axis
            %fcom=A_fcom*X+B_fcom=m*com_acceleration=m*com_a
            %com_a=g/z_com*(com-zmp)=g/z_com*[(A_com-A_zmp)*X+B_com-B_zmp]
            [gradients.A_xfcom gradients.B_xfcom]=obj.compute_force(gradients.A_xzmp,gradients.B_xzmp,gradients.A_xcom,gradients.B_xcom,walking_param.mg,walking_param.z);
            [gradients.A_yfcom gradients.B_yfcom]=obj.compute_force(gradients.A_yzmp,gradients.B_yzmp,gradients.A_ycom,gradients.B_ycom,walking_param.mg,walking_param.z);

            %% %compute the cost function criteria linked to the force on coml
            %compute the square of the force on com
            %fcom²=X^T*A_fcom2*X+2*B_fcom2*X+C_fcom2
            [gradients.A_xfcom2 gradients.B_xfcom2 gradients.C_xfcom2] = obj.compute_quad_matrix(gradients.A_xfcom,gradients.B_xfcom);
            [gradients.A_yfcom2 gradients.B_yfcom2 gradients.C_yfcom2] = obj.compute_quad_matrix(gradients.A_yfcom,gradients.B_yfcom);

        end
        function obj=compute_cost_torque_ssp(obj,walking_param,gradients,T)
            %% %compute torques in ankle during SSP in x and y axis
            
            %% %compute the ankle position for the torque computation
            %an_ssp=A_an_ssp*X+B_an_ssp
            [A_xan_ssp B_xan_ssp]=obj.compute_ankle_positions_SSP(walking_param.pankinit_firstSS(1),walking_param.pankfin_lastSS(1),walking_param.discretization);
            [A_yan_ssp B_yan_ssp]=obj.compute_ankle_positions_SSP(walking_param.pankinit_firstSS(2),walking_param.pankfin_lastSS(2),walking_param.discretization);
            A_xan_ssp=A_xan_ssp(T,:);
            B_xan_ssp=B_xan_ssp(T,:);
            A_yan_ssp=A_yan_ssp(T,:);
            B_yan_ssp=B_yan_ssp(T,:);
            
            %% %torque in ankle
            %compute the torques in ankle in x and y axis
            %?_SSP=A_t_SSP*X+B_t_SSP
            %Warning : torques in x-axis are function of optimization
            %parameters in y-axis and reverse for the torques in y-axis.
            [gradients.A_yt_ssp gradients.B_yt_ssp] = obj.compute_torque_ssp(gradients.A_xzmp,gradients.B_xzmp,gradients.A_xfcom,gradients.B_xfcom,A_xan_ssp,B_xan_ssp,walking_param.mg,walking_param.ha);
            [gradients.A_xt_ssp gradients.B_xt_ssp] = obj.compute_torque_ssp(gradients.A_yzmp,gradients.B_yzmp,gradients.A_yfcom,gradients.B_yfcom,A_yan_ssp,B_yan_ssp,walking_param.mg,walking_param.ha);

            %% %torque 2 in ankle
            %% %compute the cost function criteria linked to the torques in ankle in SSP
            %compute the square of the force on com
            %T_SSP²=X^T*A_t2_SSP*X+2*B_t2_SSP*X+C_t?2_SSP
            %Warning : torques in x-axis are function of optimization
            %parameters in y-axis and reverse for the torques in y-axis.
            [gradients.A_yt2_ssp gradients.B_yt2_ssp gradients.C_yt2_ssp] = obj.compute_quad_matrix(gradients.A_yt_ssp,gradients.B_yt_ssp);
            [gradients.A_xt2_ssp gradients.B_xt2_ssp gradients.C_xt2_ssp] = obj.compute_quad_matrix(gradients.A_xt_ssp,gradients.B_xt_ssp);
        end 
        function obj=compute_cons_zmp(obj,walking_param,gradients,T)
            %% %compute the inequality constraint equations in SSP
            %% %compute the ankle position for the torque computation used in the stability constraint
            %an_ssp=A_an_ssp*X+B_an_ssp
            [A_xan_ssp B_xan_ssp]=obj.compute_ankle_positions_SSP(walking_param.pankinit_firstSS(1),walking_param.pankfin_lastSS(1),walking_param.discretization);
            [A_yan_ssp B_yan_ssp]=obj.compute_ankle_positions_SSP(walking_param.pankinit_firstSS(2),walking_param.pankfin_lastSS(2),walking_param.discretization);
            A_xan_ssp=A_xan_ssp(T,:);
            B_xan_ssp=B_xan_ssp(T,:);
            A_yan_ssp=A_yan_ssp(T,:);
            B_yan_ssp=B_yan_ssp(T,:);
            %% %compute stability constraint which keep the ZPM under the support foot
            %Hypo: the foot is a rectangle of length l and width w
            %To keep the stability, the ZMP must stay in the rectangle.
            %O is the foot step ankle position and Z the ZMP. OZ is the
            %vector from O to Z.
            %If the the foot is oriented with an angle a, we have 4
            %conditions to verify at each time step:
            %cos(a)*OZ<=fronttoankle
            %cos(-a)*OZ>=-backtoankle
            %if the support foot is the left foot:
            %sin(a)*OZ<=exttoankle
            %sint(-a)*OZ>=-inttoankle
            %if the suppport foot is the right foot:
            %sin(a)*OZ<=inttoankle
            %sin(-a)*OZ>=-exttoankle
            %
            %we must verify :
            %A_cons_zmp_stability*X <= B_cons_zmp_stability
            [gradients.A_cons_zmp_stability gradients.B_cons_zmp_stability]=obj.zmp_constraint_stability_SSP_xy(gradients.A_xzmp,gradients.B_xzmp,gradients.A_yzmp,gradients.B_yzmp,A_xan_ssp,B_xan_ssp,A_yan_ssp,B_yan_ssp,walking_param.discretization(5:end-7),walking_param.backtoankle,walking_param.fronttoankle,walking_param.exttoankle,walking_param.inttoankle,walking_param.sole_margin,walking_param.psi,walking_param.type_phase(5:end-7),walking_param.nbparamBb,walking_param.firstSS);
        end
        
        function obj=qp_generating_DSP(obj,walking_param,gradients)

            walking_param.nbparamABCD=(walking_param.nbstep*3+2+1)*3-6;
            Cssp=[10:walking_param.nbparamABCD-9-9 walking_param.nbparamABCD+1:walking_param.nbparamABCD+4];
            walking_param.nbparamABCD=walking_param.nbparamABCD-18-9;
            walking_param.nbparamBb=(walking_param.nbstep+2+1)*3;
            L=sum(walking_param.discretization([1 3]))+2:sum(walking_param.discretization([1 3 6]))+1;
            Cdsp=[Cssp Cssp(end)+9+1:Cssp(end)+9+walking_param.nbparamank walking_param.nbparamtotal+4-walking_param.nbparamBb+10:walking_param.nbparamtotal+4-6-3];
            walking_param.nbparamBb=walking_param.nbparamBb-15-3;
            T=sum(walking_param.discretization(1:4))+2:sum(walking_param.discretization(1:7))+1;

            %%%gestion du DSP%%%
            %%%gestion de zmp1%%%
            %where to cut the beginning and ending of zmp1 in %
            cutting=1/2;
            %compute the matrix of time discretization for zmp1
            dt1=obj.compute_coeff_dt1_matrix(walking_param.discretization,walking_param.frequency,cutting);
            dt1=dt1(any(walking_param.dt_type_phase==0,2),:);
            dt1=dt1(L,:);
            %coeff gradient computation of zmp1_poly
            [xAzmp1_coeff_gradient xBzmp1_coeff_gradient]=obj.compute_coeff_zmp1_gradient(walking_param.tpassage,walking_param.xpsa_zmp1init,walking_param.xpsa_zmp1fin,cutting);
            [yAzmp1_coeff_gradient yBzmp1_coeff_gradient]=obj.compute_coeff_zmp1_gradient(walking_param.tpassage,walking_param.ypsa_zmp1init,walking_param.ypsa_zmp1fin,cutting);
            %position of ZMP1
            [gradients.xApzmp1 gradients.xBpzmp1]=obj.compute_traj_discrete(xAzmp1_coeff_gradient(:,Cdsp),xBzmp1_coeff_gradient,dt1);
            [gradients.yApzmp1 gradients.yBpzmp1]=obj.compute_traj_discrete(yAzmp1_coeff_gradient(:,Cdsp),yBzmp1_coeff_gradient,dt1);
            %matrix of force repartition in DSP
            [Ac]=obj.compute_repartition_matrix(walking_param.tpassage,cutting);
            Afb=diag(dt1*Ac);
            %% %ankle position uses for torques by zmp1
            [xApankle1 xBpankle1]=obj.torque_ankle_positions_DSP_1(walking_param.pankinit_firstinair(1),walking_param.pankinit_firstSS(1),walking_param.pankfin_lastSS(1),walking_param.discretization);
            [yApankle1 yBpankle1]=obj.torque_ankle_positions_DSP_1(walking_param.pankinit_firstinair(2),walking_param.pankinit_firstSS(2),walking_param.pankfin_lastSS(2),walking_param.discretization);
    %         yBpankle1(1:walking_param.discretization(1)+1)=-yBpankle1(1:walking_param.discretization(1)+1);%les pieds sont parallèles au debut, La création des pankle commence par le pied gauche or zmp1 est sous le pieds droit
            xApankle1=xApankle1(any(walking_param.dt_type_phase==0,2),:);
            xBpankle1=xBpankle1(any(walking_param.dt_type_phase==0,2),:);
            yApankle1=yApankle1(any(walking_param.dt_type_phase==0,2),:);
            yBpankle1=yBpankle1(any(walking_param.dt_type_phase==0,2),:);

            xApankle1=xApankle1(L,:);
            xBpankle1=xBpankle1(L,:);
            yApankle1=yApankle1(L,:);
            yBpankle1=yBpankle1(L,:);
            %% %torques by zmp1
            xAfcom_=gradients.xAfcom(any(walking_param.dt_type_phase(T,:)==0,2),:);
            xBfcom_=gradients.xBfcom(any(walking_param.dt_type_phase(T,:)==0,2),:);
            yAfcom_=gradients.yAfcom(any(walking_param.dt_type_phase(T,:)==0,2),:);
            yBfcom_=gradients.yBfcom(any(walking_param.dt_type_phase(T,:)==0,2),:);

            [xAtorquedsp1 gradients.xBtorquedsp1] = obj.torque_ankle_DSP_gradient_1(gradients.xApzmp1,gradients.xBpzmp1,xAfcom_,xBfcom_,xApankle1,xBpankle1,walking_param.mg,walking_param.ha,Afb);
            [yAtorquedsp1 gradients.yBtorquedsp1] = obj.torque_ankle_DSP_gradient_1(gradients.yApzmp1,gradients.yBpzmp1,yAfcom_,yBfcom_,yApankle1,yBpankle1,walking_param.mg,walking_param.ha,Afb);

            [gradients.xAtorque2_1 gradients.xBtorque2_1] = obj.compute_quad_matrix(xAtorquedsp1,gradients.xBtorquedsp1);
            [gradients.yAtorque2_1 gradients.yBtorque2_1] = obj.compute_quad_matrix(yAtorquedsp1,gradients.yBtorquedsp1);

            %% %%%gestion de zmp2%%%
            %% %position of ZMP2
            xApzmp_=gradients.xApzmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            xBpzmp_=gradients.xBpzmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            yApzmp_=gradients.yApzmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            yBpzmp_=gradients.yBpzmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            [gradients.xApzmp2 gradients.xBpzmp2]=obj.compute_pzmp2(xApzmp_,xBpzmp_,gradients.xApzmp1,gradients.xBpzmp1,Afb,walking_param.mg);
            [gradients.yApzmp2 gradients.yBpzmp2]=obj.compute_pzmp2(yApzmp_,yBpzmp_,gradients.yApzmp1,gradients.yBpzmp1,Afb,walking_param.mg);
            %% %ankle position uses for torques by zmp2
            [xApankle2 xBpankle2]=obj.torque_ankle_positions_DSP_2(walking_param.pankinit_firstSS(1),walking_param.pankfin_lastSS(1),walking_param.pankfin_lastinair(1),walking_param.discretization);
            [yApankle2 yBpankle2]=obj.torque_ankle_positions_DSP_2(walking_param.pankinit_firstSS(2),walking_param.pankfin_lastSS(2),walking_param.pankfin_lastinair(2),walking_param.discretization);
    %         yBpankle2(end-discretization(end):end)=-yBpankle2(end-discretization(end):end);
    %         yBpankle2(end-walking_param.discretization(end)+1:end)=yBpankle2(end-walking_param.discretization(end)+1:end)-(-1)^(walking_param.nbstep)*0.095*2;
            xApankle2=xApankle2(any(walking_param.dt_type_phase==0,2),:);
            xBpankle2=xBpankle2(any(walking_param.dt_type_phase==0,2),:);
            yApankle2=yApankle2(any(walking_param.dt_type_phase==0,2),:);
            yBpankle2=yBpankle2(any(walking_param.dt_type_phase==0,2),:);

            xApankle2=xApankle2(L,:);
            xBpankle2=xBpankle2(L,:);
            yApankle2=yApankle2(L,:);
            yBpankle2=yBpankle2(L,:);
            %% %torques by zmp2
            [xAtorquedsp2 gradients.xBtorquedsp2] = obj.torque_ankle_DSP_gradient_2(gradients.xApzmp2,gradients.xBpzmp2,xAfcom_,xBfcom_,xApankle2,xBpankle2,walking_param.mg,walking_param.ha,Afb);
            [yAtorquedsp2 gradients.yBtorquedsp2] = obj.torque_ankle_DSP_gradient_2(gradients.yApzmp2,gradients.yBpzmp2,yAfcom_,yBfcom_,yApankle2,yBpankle2,walking_param.mg,walking_param.ha,Afb);

            [gradients.xAtorque2_2 gradients.xBtorque2_2] = obj.compute_quad_matrix(xAtorquedsp2,gradients.xBtorquedsp2);
            [gradients.yAtorque2_2 gradients.yBtorque2_2] = obj.compute_quad_matrix(yAtorquedsp2,gradients.yBtorquedsp2);

            %% %%%cost function acceleration zmp1 & 2%%%
            dddt1=obj.compute_coeff_dddt1_matrix(walking_param.discretization,walking_param.frequency,cutting);
            dddt1=dddt1(any(walking_param.dt_type_phase==0,2),:);
            dddt1=dddt1(L,:);
            %% %acceleration of zmp1
            [xAazmp1 gradients.xBazmp1]=obj.compute_traj_discrete(xAzmp1_coeff_gradient(:,Cdsp),xBzmp1_coeff_gradient,dddt1);
            [yAazmp1 gradients.yBazmp1]=obj.compute_traj_discrete(yAzmp1_coeff_gradient(:,Cdsp),yBzmp1_coeff_gradient,dddt1);
            %% %acceleration2 of zmp1
            [gradients.xAa2zmp1 gradients.xBa2zmp1]=obj.compute_quad_matrix(xAazmp1,gradients.xBazmp1);
            [gradients.yAa2zmp1 gradients.yBa2zmp1]=obj.compute_quad_matrix(yAazmp1,gradients.yBazmp1);

            %% %speed of zmp1
            ddt1=obj.compute_coeff_ddt1_matrix(walking_param.discretization,walking_param.frequency,cutting);
            ddt1=ddt1(any(walking_param.dt_type_phase==0,2),:);
            ddt1=ddt1(L,:);
            [xAszmp1 xBszmp1]=obj.compute_traj_discrete(xAzmp1_coeff_gradient(:,Cdsp),xBzmp1_coeff_gradient,ddt1);
            [yAszmp1 yBszmp1]=obj.compute_traj_discrete(yAzmp1_coeff_gradient(:,Cdsp),yBzmp1_coeff_gradient,ddt1);
            %% %speed of zmp
            ddt=obj.compute_coeff_ddt_matrix(walking_param.discretization,walking_param.frequency);
            ddt=ddt(T,:);
            [xAszmp xBszmp] = obj.compute_traj_discrete(gradients.xAzmp_gradient(:,Cssp),gradients.xBzmp_gradient,ddt);
            [yAszmp yBszmp] = obj.compute_traj_discrete(gradients.yAzmp_gradient(:,Cssp),gradients.yBzmp_gradient,ddt);

            xAszmp_=xAszmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            xBszmp_=xBszmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            yAszmp_=yAszmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            yBszmp_=yBszmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            %% %acceleration of zmp
            dddt=obj.compute_coeff_dddt_matrix(walking_param.discretization,walking_param.frequency);
            dddt=dddt(T,:);
            [xAazmp xBazmp] = obj.compute_traj_discrete(gradients.xAzmp_gradient(:,Cssp),gradients.xBzmp_gradient,dddt);
            [yAazmp yBazmp] = obj.compute_traj_discrete(gradients.yAzmp_gradient(:,Cssp),gradients.yBzmp_gradient,dddt);

            xAazmp_=xAazmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            xBazmp_=xBazmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            yAazmp_=yAazmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            yBazmp_=yBazmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            %% %acceleration of zmp2
            dAfb=diag(ddt1*Ac);
            ddAfb=diag(dddt1*Ac);
            
            [xAazmp2 gradients.xBazmp2]=obj.compute_azmp2_gradient(xApzmp_,xBpzmp_,gradients.xApzmp1,gradients.xBpzmp1,xAszmp_,xBszmp_,xAszmp1,xBszmp1,xAazmp_,xBazmp_,xAazmp1,gradients.xBazmp1,Afb,dAfb,ddAfb);
            [yAazmp2 gradients.yBazmp2]=obj.compute_azmp2_gradient(yApzmp_,yBpzmp_,gradients.yApzmp1,gradients.yBpzmp1,yAszmp_,yBszmp_,yAszmp1,yBszmp1,yAazmp_,yBazmp_,yAazmp1,gradients.yBazmp1,Afb,dAfb,ddAfb);
            %acceleration2 of zmp1
            [gradients.xAa2zmp2 gradients.xBa2zmp2]=obj.compute_quad_matrix(xAazmp2,gradients.xBazmp2);
            [gradients.yAa2zmp2 gradients.yBa2zmp2]=obj.compute_quad_matrix(yAazmp2,gradients.yBazmp2);

            %% %%%constraint stability zmp1&2
            [Aconstraint1 Bconstraint1]=obj.zmp_constraint_stability_DSP_xy_1(gradients.xApzmp1,gradients.xBpzmp1,gradients.yApzmp1,gradients.yBpzmp1,xApankle1,xBpankle1,yApankle1,yBpankle1,walking_param.discretization(5:end-7),walking_param.backtoankle,walking_param.fronttoankle,walking_param.exttoankle,walking_param.inttoankle,walking_param.sole_margin,walking_param.psi,walking_param.type_phase(5:end-7),walking_param.nbparamABCD,walking_param.nbparamBb,walking_param.firstSS);
            [Aconstraint2 Bconstraint2]=obj.zmp_constraint_stability_DSP_xy_2(gradients.xApzmp2,gradients.xBpzmp2,gradients.yApzmp2,gradients.yBpzmp2,xApankle2,xBpankle2,yApankle2,yBpankle2,walking_param.discretization(5:end-7),walking_param.backtoankle,walking_param.fronttoankle,walking_param.exttoankle,walking_param.inttoankle,walking_param.sole_margin,walking_param.psi,walking_param.type_phase(5:end-7),walking_param.nbparamABCD,walking_param.nbparamBb,walking_param.firstSS);
            %% %%%contrainte de non chevauchement des pieds%%%
            [Apankconst Bpankconst]=obj.zmp_constraint_ankle_pos(walking_param.pankinit_firstSS,walking_param.pankfin_lastSS,walking_param.backtoankle,walking_param.fronttoankle,walking_param.exttoankle,walking_param.inttoankle,walking_param.xankmax,walking_param.xankmin,walking_param.yankmax,walking_param.yankmin,walking_param.psi,walking_param.nbparamABCD,walking_param.nbparamank,walking_param.nbparamBb);

            gradients.AconstraintDSP=[gradients.Apzmpconstraintssp;    Aconstraint1;   Aconstraint2;   Apankconst];
            gradients.BconstraintDSP=[gradients.Bpzmpconstraintssp;    Bconstraint1;   Bconstraint2;   Bpankconst];


            %adaptation de Aeq avec DSP
            gradients.AscomeqDSP=[gradients.xAscomeq                 zeros(size(gradients.xAscomeq,1),walking_param.nbparamank) zeros(size(gradients.xAscomeq,1),walking_param.nbparamBb) zeros(size(gradients.yAscomeq)) zeros(size(gradients.xAscomeq,1),walking_param.nbparamank) zeros(size(gradients.xAscomeq,1),walking_param.nbparamBb);
                                  zeros(size(gradients.xAscomeq))    zeros(size(gradients.yAscomeq,1),walking_param.nbparamank) zeros(size(gradients.yAscomeq,1),walking_param.nbparamBb) gradients.yAscomeq              zeros(size(gradients.yAscomeq,1),walking_param.nbparamank) zeros(size(gradients.yAscomeq,1),walking_param.nbparamBb)];
        end
        
        function obj=compute_ZMP1(obj,walking_param,gradients,M_t1,cutting,Cdsp)
            %coeff gradient computation of zmp1_poly
            [gradients.A_xCa1 gradients.B_xCa1]=obj.compute_zmp1_Ca1(walking_param.tpassage,walking_param.xpsa_zmp1init,walking_param.xpsa_zmp1fin,cutting);
            [gradients.A_yCa1 gradients.B_yCa1]=obj.compute_zmp1_Ca1(walking_param.tpassage,walking_param.ypsa_zmp1init,walking_param.ypsa_zmp1fin,cutting);
            %position of ZMP1
            [gradients.A_xzmp1 gradients.B_xzmp1]=obj.compute_traj_discrete(gradients.A_xCa1(:,Cdsp),gradients.B_xCa1,M_t1);
            [gradients.A_yzmp1 gradients.B_yzmp1]=obj.compute_traj_discrete(gradients.A_yCa1(:,Cdsp),gradients.B_yCa1,M_t1);
        end
        function obj=compute_cost_torque_zmp1(obj,walking_param,gradients,T,L)
            %% %ankle position uses for torques by zmp1
            [A_xan_dsp B_xan_dsp]=obj.compute_ankle_positions_zmp1(walking_param.pankinit_firstinair(1),walking_param.pankinit_firstSS(1),walking_param.pankfin_lastSS(1),walking_param.discretization);
            [A_yan_dsp B_yan_dsp]=obj.compute_ankle_positions_zmp1(walking_param.pankinit_firstinair(2),walking_param.pankinit_firstSS(2),walking_param.pankfin_lastSS(2),walking_param.discretization);
            %yBpankle1(1:walking_param.discretization(1)+1)=-yBpankle1(1:walking_param.discretization(1)+1);%les pieds sont parallèles au debut, La création des pankle commence par le pied gauche or zmp1 est sous le pieds droit
            
            %suppress SSP discretization steps because ZMP1 is only defined in DSP
            A_xan_dsp=A_xan_dsp(any(walking_param.dt_type_phase==0,2),:);
            B_xan_dsp=B_xan_dsp(any(walking_param.dt_type_phase==0,2),:);
            A_yan_dsp=A_yan_dsp(any(walking_param.dt_type_phase==0,2),:);
            B_yan_dsp=B_yan_dsp(any(walking_param.dt_type_phase==0,2),:);

            %suppress value to keep only one foot step
            A_xan_dsp=A_xan_dsp(L,:);
            B_xan_dsp=B_xan_dsp(L,:);
            A_yan_dsp=A_yan_dsp(L,:);
            B_yan_dsp=B_yan_dsp(L,:);
            %% %torques by zmp1
            %suppress SSP discretization steps because ZMP1 is only defined in DSP
            A_xfcom_=gradients.A_xfcom(any(walking_param.dt_type_phase(T,:)==0,2),:);
            B_xfcom_=gradients.B_xfcom(any(walking_param.dt_type_phase(T,:)==0,2),:);
            A_yfcom_=gradients.A_yfcom(any(walking_param.dt_type_phase(T,:)==0,2),:);
            B_yfcom_=gradients.B_yfcom(any(walking_param.dt_type_phase(T,:)==0,2),:);
            %Warning : torques in x-axis are function of optimization
            %parameters in y-axis and reverse for the torques in y-axis.
            [gradients.A_yt_zmp1 gradients.B_yt_zmp1] = obj.compute_torque_zmp1(gradients.A_xzmp1,gradients.B_xzmp1,A_xfcom_,B_xfcom_,A_xan_dsp,B_xan_dsp,walking_param.mg,walking_param.ha,gradients.k_diag);
            [gradients.A_xt_zmp1 gradients.B_xt_zmp1] = obj.compute_torque_zmp1(gradients.A_yzmp1,gradients.B_yzmp1,A_yfcom_,B_yfcom_,A_yan_dsp,B_yan_dsp,walking_param.mg,walking_param.ha,gradients.k_diag);
            %Warning : torques in x-axis are function of optimization
            %parameters in y-axis and reverse for the torques in y-axis.
            [gradients.A_yt2_zmp1 gradients.B_yt2_zmp1 gradients.C_yt2_zmp1] = obj.compute_quad_matrix(gradients.A_yt_zmp1,gradients.B_yt_zmp1);
            [gradients.A_xt2_zmp1 gradients.B_xt2_zmp1 gradients.C_xt2_zmp1] = obj.compute_quad_matrix(gradients.A_xt_zmp1,gradients.B_xt_zmp1);
        end
        function obj=compute_cost_acc_zmp1(obj,walking_param,gradients,cutting,L,Cdsp)
            %% %%%cost function acceleration zmp1 & 2%%%
            %compute the double derivative of M_t1
            M_t1_dd=obj.compute_M_t1_dd(walking_param.discretization,walking_param.frequency,cutting);
            M_t1_dd=M_t1_dd(any(walking_param.dt_type_phase==0,2),:);
            M_t1_dd=M_t1_dd(L,:);
            %% %acceleration of zmp1
            [gradients.A_xzmp1_acc gradients.B_xzmp1_acc]=obj.compute_traj_discrete(gradients.A_xCa1(:,Cdsp),gradients.B_xCa1,M_t1_dd);
            [gradients.A_yzmp1_acc gradients.B_yzmp1_acc]=obj.compute_traj_discrete(gradients.A_yCa1(:,Cdsp),gradients.B_yCa1,M_t1_dd);
            %% %acceleration2 of zmp1
            [gradients.A_xzmp1_acc2 gradients.B_xzmp1_acc2 gradients.C_xzmp1_acc2]=obj.compute_quad_matrix(gradients.A_xzmp1_acc,gradients.B_xzmp1_acc);
            [gradients.A_yzmp1_acc2 gradients.B_yzmp1_acc2 gradients.C_yzmp1_acc2]=obj.compute_quad_matrix(gradients.A_yzmp1_acc,gradients.B_yzmp1_acc);
        end
        function obj=compute_ZMP2(obj,walking_param,gradients,T)
            %% %position of ZMP2
            A_xzmp_=gradients.A_xzmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            B_xzmp_=gradients.B_xzmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            A_yzmp_=gradients.A_yzmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            B_yzmp_=gradients.B_yzmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            [gradients.A_xzmp2 gradients.B_xzmp2]=obj.compute_zmp2(A_xzmp_,B_xzmp_,gradients.A_xzmp1,gradients.B_xzmp1,gradients.k_diag,walking_param.mg);
            [gradients.A_yzmp2 gradients.B_yzmp2]=obj.compute_zmp2(A_yzmp_,B_yzmp_,gradients.A_yzmp1,gradients.B_yzmp1,gradients.k_diag,walking_param.mg);
        end
        function obj=compute_cost_torque_zmp2(obj,walking_param,gradients,T,L)
            %% %ankle position uses for torques by zmp2
            [A_xan_dsp B_xan_dsp]=obj.compute_ankle_positions_zmp2(walking_param.pankinit_firstSS(1),walking_param.pankfin_lastSS(1),walking_param.pankfin_lastinair(1),walking_param.discretization);
            [A_yan_dsp B_yan_dsp]=obj.compute_ankle_positions_zmp2(walking_param.pankinit_firstSS(2),walking_param.pankfin_lastSS(2),walking_param.pankfin_lastinair(2),walking_param.discretization);
    %         yBpankle2(end-discretization(end):end)=-yBpankle2(end-discretization(end):end);
    %         yBpankle2(end-walking_param.discretization(end)+1:end)=yBpankle2(end-walking_param.discretization(end)+1:end)-(-1)^(walking_param.nbstep)*0.095*2;
            A_xan_dsp=A_xan_dsp(any(walking_param.dt_type_phase==0,2),:);
            B_xan_dsp=B_xan_dsp(any(walking_param.dt_type_phase==0,2),:);
            A_yan_dsp=A_yan_dsp(any(walking_param.dt_type_phase==0,2),:);
            B_yan_dsp=B_yan_dsp(any(walking_param.dt_type_phase==0,2),:);

            A_xan_dsp=A_xan_dsp(L,:);
            B_xan_dsp=B_xan_dsp(L,:);
            A_yan_dsp=A_yan_dsp(L,:);
            B_yan_dsp=B_yan_dsp(L,:);
            %% %torques by zmp1
            %suppress SSP discretization steps because ZMP1 is only defined in DSP
            A_xfcom_=gradients.A_xfcom(any(walking_param.dt_type_phase(T,:)==0,2),:);
            B_xfcom_=gradients.B_xfcom(any(walking_param.dt_type_phase(T,:)==0,2),:);
            A_yfcom_=gradients.A_yfcom(any(walking_param.dt_type_phase(T,:)==0,2),:);
            B_yfcom_=gradients.B_yfcom(any(walking_param.dt_type_phase(T,:)==0,2),:);
            %% %torques by zmp2
            %Warning : torques in x-axis are function of optimization
            %parameters in y-axis and reverse for the torques in y-axis.
            [gradients.A_yt_zmp2 gradients.B_yt_zmp2] = obj.compute_torque_zmp2(gradients.A_xzmp2,gradients.B_xzmp2,A_xfcom_,B_xfcom_,A_xan_dsp,B_xan_dsp,walking_param.mg,walking_param.ha,gradients.k_diag);
            [gradients.A_xt_zmp2 gradients.B_xt_zmp2] = obj.compute_torque_zmp2(gradients.A_yzmp2,gradients.B_yzmp2,A_yfcom_,B_yfcom_,A_yan_dsp,B_yan_dsp,walking_param.mg,walking_param.ha,gradients.k_diag);
            %Warning : torques in x-axis are function of optimization
            %parameters in y-axis and reverse for the torques in y-axis.
            [gradients.A_yt2_zmp2 gradients.B_yt2_zmp2 gradients.C_yt2_zmp2] = obj.compute_quad_matrix(gradients.A_yt_zmp2,gradients.B_yt_zmp2);
            [gradients.A_xt2_zmp2 gradients.B_xt2_zmp2 gradients.C_xt2_zmp2] = obj.compute_quad_matrix(gradients.A_xt_zmp2,gradients.B_xt_zmp2);
        end
        function obj=compute_cost_acc_zmp2(obj,walking_param,gradients,cutting,Cssp,T,L,Cdsp)
            %compute the double derivative of M_t1
            M_t1_dd=obj.compute_M_t1_dd(walking_param.discretization,walking_param.frequency,cutting);
            M_t1_dd=M_t1_dd(any(walking_param.dt_type_phase==0,2),:);
            M_t1_dd=M_t1_dd(L,:);
            %% %speed of zmp1
            M_t1_d=obj.compute_coeff_ddt1_matrix(walking_param.discretization,walking_param.frequency,cutting);
            M_t1_d=M_t1_d(any(walking_param.dt_type_phase==0,2),:);
            M_t1_d=M_t1_d(L,:);
            [xAszmp1 xBszmp1]=obj.compute_traj_discrete(gradients.A_xCa1(:,Cdsp),gradients.B_xCa1,M_t1_d);
            [yAszmp1 yBszmp1]=obj.compute_traj_discrete(gradients.A_yCa1(:,Cdsp),gradients.B_yCa1,M_t1_d);
            
            %% %position of zmp
            A_xzmp_=gradients.A_xzmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            B_xzmp_=gradients.B_xzmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            A_yzmp_=gradients.A_yzmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            B_yzmp_=gradients.B_yzmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            %% %speed of zmp
            M_t_d=obj.compute_coeff_ddt_matrix(walking_param.discretization,walking_param.frequency);
            M_t_d=M_t_d(T,:);
            [xAszmp xBszmp] = obj.compute_traj_discrete(gradients.A_xCa(:,Cssp),gradients.B_xCa,M_t_d);
            [yAszmp yBszmp] = obj.compute_traj_discrete(gradients.A_yCa(:,Cssp),gradients.B_yCa,M_t_d);

            xAszmp_=xAszmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            xBszmp_=xBszmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            yAszmp_=yAszmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            yBszmp_=yBszmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            %% %acceleration of zmp
            M_t_dd=obj.compute_coeff_dddt_matrix(walking_param.discretization,walking_param.frequency);
            M_t_dd=M_t_dd(T,:);
            [xAazmp xBazmp] = obj.compute_traj_discrete(gradients.A_xCa(:,Cssp),gradients.B_xCa,M_t_dd);
            [yAazmp yBazmp] = obj.compute_traj_discrete(gradients.A_yCa(:,Cssp),gradients.B_yCa,M_t_dd);

            xAazmp_=xAazmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            xBazmp_=xBazmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            yAazmp_=yAazmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            yBazmp_=yBazmp(any(walking_param.dt_type_phase(T,:)==0,2),:);
            %% %acceleration of zmp2
            k_diag_d=diag(M_t1_d*gradients.Ac);
            k_diag_dd=diag(M_t1_dd*gradients.Ac);
            
            [gradients.A_xzmp2_acc gradients.B_xzmp2_acc]=obj.compute_zmp2_acc(A_xzmp_,B_xzmp_,gradients.A_xzmp1,gradients.B_xzmp1,xAszmp_,xBszmp_,xAszmp1,xBszmp1,xAazmp_,xBazmp_,gradients.A_xzmp1_acc,gradients.B_xzmp1_acc,gradients.k_diag,k_diag_d,k_diag_dd);
            [gradients.A_yzmp2_acc gradients.B_yzmp2_acc]=obj.compute_zmp2_acc(A_yzmp_,B_yzmp_,gradients.A_yzmp1,gradients.B_yzmp1,yAszmp_,yBszmp_,yAszmp1,yBszmp1,yAazmp_,yBazmp_,gradients.A_yzmp1_acc,gradients.B_yzmp1_acc,gradients.k_diag,k_diag_d,k_diag_dd);
            %acceleration2 of zmp1
            [gradients.A_xzmp2_acc2 gradients.B_xzmp2_acc2 gradients.C_xzmp2_acc2]=obj.compute_quad_matrix(gradients.A_xzmp2_acc,gradients.B_xzmp2_acc);
            [gradients.A_yzmp2_acc2 gradients.B_yzmp2_acc2 gradients.C_yzmp2_acc2]=obj.compute_quad_matrix(gradients.A_yzmp2_acc,gradients.B_yzmp2_acc);
        end
        function obj=compute_cons_zmp12(obj,walking_param,gradients,L)
            %% %ankle position uses for torques by zmp1
            [A_xan_dsp B_xan_dsp]=obj.compute_ankle_positions_zmp1(walking_param.pankinit_firstinair(1),walking_param.pankinit_firstSS(1),walking_param.pankfin_lastSS(1),walking_param.discretization);
            [A_yan_dsp B_yan_dsp]=obj.compute_ankle_positions_zmp1(walking_param.pankinit_firstinair(2),walking_param.pankinit_firstSS(2),walking_param.pankfin_lastSS(2),walking_param.discretization);
            %yBpankle1(1:walking_param.discretization(1)+1)=-yBpankle1(1:walking_param.discretization(1)+1);%les pieds sont parallèles au debut, La création des pankle commence par le pied gauche or zmp1 est sous le pieds droit
            
            %suppress SSP discretization steps because ZMP1 is only defined in DSP
            A_xan_dsp=A_xan_dsp(any(walking_param.dt_type_phase==0,2),:);
            B_xan_dsp=B_xan_dsp(any(walking_param.dt_type_phase==0,2),:);
            A_yan_dsp=A_yan_dsp(any(walking_param.dt_type_phase==0,2),:);
            B_yan_dsp=B_yan_dsp(any(walking_param.dt_type_phase==0,2),:);

            %suppress value to keep only one foot step
            A_xan_dsp1=A_xan_dsp(L,:);
            B_xan_dsp1=B_xan_dsp(L,:);
            A_yan_dsp1=A_yan_dsp(L,:);
            B_yan_dsp1=B_yan_dsp(L,:);
            %% %ankle position uses for torques by zmp2
            [A_xan_dsp B_xan_dsp]=obj.compute_ankle_positions_zmp2(walking_param.pankinit_firstSS(1),walking_param.pankfin_lastSS(1),walking_param.pankfin_lastinair(1),walking_param.discretization);
            [A_yan_dsp B_yan_dsp]=obj.compute_ankle_positions_zmp2(walking_param.pankinit_firstSS(2),walking_param.pankfin_lastSS(2),walking_param.pankfin_lastinair(2),walking_param.discretization);
    %         yBpankle2(end-discretization(end):end)=-yBpankle2(end-discretization(end):end);
    %         yBpankle2(end-walking_param.discretization(end)+1:end)=yBpankle2(end-walking_param.discretization(end)+1:end)-(-1)^(walking_param.nbstep)*0.095*2;
            A_xan_dsp=A_xan_dsp(any(walking_param.dt_type_phase==0,2),:);
            B_xan_dsp=B_xan_dsp(any(walking_param.dt_type_phase==0,2),:);
            A_yan_dsp=A_yan_dsp(any(walking_param.dt_type_phase==0,2),:);
            B_yan_dsp=B_yan_dsp(any(walking_param.dt_type_phase==0,2),:);

            A_xan_dsp2=A_xan_dsp(L,:);
            B_xan_dsp2=B_xan_dsp(L,:);
            A_yan_dsp2=A_yan_dsp(L,:);
            B_yan_dsp2=B_yan_dsp(L,:);
            %% %%%constraint stability zmp1&2
            [gradients.A_cons_zmp1_stab gradients.B_cons_zmp1_stab]=obj.cons_xy_zmp1_stability(gradients.A_xzmp1,gradients.B_xzmp1,gradients.A_yzmp1,gradients.B_yzmp1,A_xan_dsp1,B_xan_dsp1,A_yan_dsp1,B_yan_dsp1,walking_param.discretization(5:end-7),walking_param.backtoankle,walking_param.fronttoankle,walking_param.exttoankle,walking_param.inttoankle,walking_param.sole_margin,walking_param.psi,walking_param.type_phase(5:end-7),walking_param.nbparamABCD,walking_param.nbparamBb,walking_param.firstSS);
            [gradients.A_cons_zmp2_stab gradients.B_cons_zmp2_stab]=obj.cons_xy_zmp2_stability(gradients.A_xzmp2,gradients.B_xzmp2,gradients.A_yzmp2,gradients.B_yzmp2,A_xan_dsp2,B_xan_dsp2,A_yan_dsp2,B_yan_dsp2,walking_param.discretization(5:end-7),walking_param.backtoankle,walking_param.fronttoankle,walking_param.exttoankle,walking_param.inttoankle,walking_param.sole_margin,walking_param.psi,walking_param.type_phase(5:end-7),walking_param.nbparamABCD,walking_param.nbparamBb,walking_param.firstSS);
        end
        
        function obj=qp_generating_problem_matrix(obj,walking_param,gradients)
            [xAviapoint_nl xBviapoint_nl] = obj.cost_viapoint_gradient_nl(gradients.xAfcom2,gradients.xBfcom2,gradients.xAtorque2ssp,gradients.xBtorque2ssp,walking_param.lambda,walking_param.mu);
            [yAviapoint_nl yBviapoint_nl] = obj.cost_viapoint_gradient_nl(gradients.yAfcom2,gradients.yBfcom2,gradients.yAtorque2ssp,gradients.yBtorque2ssp,walking_param.lambda,walking_param.mu);
            
            %%% cost taking into account torque in DSP and azmp1 %%%
            [xAviapointDSP_nl xBviapointDSP_nl]=obj.cost_viapoint_gradient_withDSP_nl(xAviapoint_nl,xBviapoint_nl,gradients.xAtorque2_1,gradients.xBtorque2_1,gradients.xAtorque2_2,gradients.xBtorque2_2,gradients.xAa2zmp1+gradients.xAa2zmp2,gradients.xBa2zmp1+gradients.xBa2zmp2,walking_param.epsilon,0.5*walking_param.mu);
            [yAviapointDSP_nl yBviapointDSP_nl]=obj.cost_viapoint_gradient_withDSP_nl(yAviapoint_nl,yBviapoint_nl,gradients.yAtorque2_1,gradients.yBtorque2_1,gradients.yAtorque2_2,gradients.yBtorque2_2,gradients.yAa2zmp1+gradients.yAa2zmp2,gradients.yBa2zmp1+gradients.yBa2zmp2,walking_param.epsilon,0.5*walking_param.mu);

            %%%position, speed and acceleration gradient of viapoints%%%
            [xAviapoint xBviapoint] = obj.cost_viapoint_gradient(gradients.xAfcom2,gradients.xBfcom2,gradients.xAtorque2ssp,gradients.xBtorque2ssp,walking_param.lambda,walking_param.mu);
            [yAviapoint yBviapoint] = obj.cost_viapoint_gradient(gradients.yAfcom2,gradients.yBfcom2,gradients.yAtorque2ssp,gradients.yBtorque2ssp,walking_param.lambda,walking_param.mu);

            %%% cost taking into account torque in DSP and azmp1 %%%
            [xAviapointDSP xBviapointDSP]=obj.cost_viapoint_gradient_withDSP_azmp1(xAviapoint,xBviapoint,gradients.xAtorque2_1,gradients.xBtorque2_1,gradients.xAtorque2_2,gradients.xBtorque2_2,gradients.xAa2zmp1+gradients.xAa2zmp2,gradients.xBa2zmp1+gradients.xBa2zmp2,walking_param.epsilon,0.5*walking_param.mu);
            [yAviapointDSP yBviapointDSP]=obj.cost_viapoint_gradient_withDSP_azmp1(yAviapoint,yBviapoint,gradients.yAtorque2_1,gradients.yBtorque2_1,gradients.yAtorque2_2,gradients.yBtorque2_2,gradients.yAa2zmp1+gradients.yAa2zmp2,gradients.yBa2zmp1+gradients.yBa2zmp2,walking_param.epsilon,0.5*walking_param.mu);

            %%%viapointcost with DSP%%%
            obj.H=[xAviapointDSP zeros(size(yAviapointDSP));zeros(size(xAviapointDSP)) yAviapointDSP];
            obj.f=[xBviapointDSP yBviapointDSP]';

            obj.H_nl=[xAviapointDSP_nl zeros(size(yAviapointDSP_nl));zeros(size(xAviapointDSP_nl)) yAviapointDSP_nl];
            obj.f_nl=[xBviapointDSP_nl yBviapointDSP_nl]';            
            
            AscomeqDSP_path=gradients.AscomeqDSP;
            Bscomeq_path=gradients.Bscomeq;

            %%%%Contrainte sur des ankle positions%%%
            [Aeq beq]=obj.pankle_fixed_path(AscomeqDSP_path,Bscomeq_path,walking_param.nbparamABCD,walking_param.nbparamank,walking_param.nbparamBb,walking_param.step_number_pankle_fixed);


            %% 1 4 7
            xpA=[[-1 0 0] zeros(1,6) [1 0 0] zeros(1,4+walking_param.nbparamank) zeros(1,walking_param.nbparamBb) zeros(1,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb)];
            Aeq=[Aeq;xpA];
            beq=[beq;-0.125];
            ypA=[zeros(1,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb) [1 0 0] zeros(1,6) [1 0 0] zeros(1,4+walking_param.nbparamank) zeros(1,walking_param.nbparamBb)];
            Aeq=[Aeq;ypA];
            beq=[beq;0];

            xsA=[[0 -1 0] zeros(1,6) [0 1 0] zeros(1,4+walking_param.nbparamank) zeros(1,walking_param.nbparamBb) zeros(1,4+walking_param.nbparamABCD+walking_param.nbparamank+walking_param.nbparamBb)];
            Aeq=[Aeq;xsA];
            beq=[beq;0];
            ysA=[zeros(1,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb) [0 1 0] zeros(1,6) [0 1 0] zeros(1,4+walking_param.nbparamank) zeros(1,walking_param.nbparamBb) ];
            Aeq=[Aeq;ysA];
            beq=[beq;0];

            xaA=[[0 0 -1] zeros(1,6) [0 0 1] zeros(1,4+walking_param.nbparamank) zeros(1,walking_param.nbparamBb) zeros(1,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb)];
            Aeq=[Aeq;xaA];
            beq=[beq;0];
            yaA=[zeros(1,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb) [0 0 1] zeros(1,6) [0 0 1] zeros(1,4+walking_param.nbparamank) zeros(1,walking_param.nbparamBb) ];
            Aeq=[Aeq;yaA];
            beq=[beq;0];
 
            xpA=[zeros(2,walking_param.nbparamABCD) [-1 1 0 0;0 0 -1 1] zeros(2,walking_param.nbparamank) zeros(2,walking_param.nbparamBb) zeros(2,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb)];
            Aeq=[Aeq;xpA];
            beq=[beq;-0.125;0];
            ypA=[zeros(2,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb) zeros(2,walking_param.nbparamABCD) [1 1 0 0;0 0 1 1] zeros(2,walking_param.nbparamank) zeros(2,walking_param.nbparamBb)];
            Aeq=[Aeq;ypA];
            beq=[beq;0;0];

            obj.Aeq=Aeq;
            obj.beq=beq;
        end
        function obj=compute_H_G_C_nl(obj,walking_param,gradients)
            %compute the coefficient of the cost function
            %cost=1/2*X^T*H*X+G*X+1/2*C
            %cost_quad=1/2*X^T*H*X+f'*X=1/2*X^T*H*X+G'*X
            %% Compute H
            xH_ssp=walking_param.lambda*[gradients.A_xfcom2 zeros(size(gradients.A_xfcom2,1),size(gradients.A_yt2_ssp,2)-size(gradients.A_xfcom2,2)); ...
                zeros(size(gradients.A_yt2_ssp,1)-size(gradients.A_xfcom2,1),size(gradients.A_xfcom2,2)) ... 
                zeros(size(gradients.A_yt2_ssp,1)-size(gradients.A_xfcom2,1),size(gradients.A_yt2_ssp,2)-size(gradients.A_xfcom2,2))] ...
                +walking_param.mu*gradients.A_yt2_ssp;
            xH=[xH_ssp zeros(size(xH_ssp,1),size(gradients.A_yt2_zmp1,2)-size(xH_ssp,2)); ...
                zeros(size(gradients.A_yt2_zmp1,1)-size(xH_ssp,1),size(gradients.A_yt2_zmp1,2))] ...
                +0.5*walking_param.mu*(gradients.A_yt2_zmp1+gradients.A_yt2_zmp2) ...
                +walking_param.epsilon*(gradients.A_xzmp1_acc2+gradients.A_xzmp2_acc2);
            
            yH_ssp=walking_param.lambda*[gradients.A_yfcom2 zeros(size(gradients.A_yfcom2,1),size(gradients.A_xt2_ssp,2)-size(gradients.A_yfcom2,2)); ...
                zeros(size(gradients.A_xt2_ssp,1)-size(gradients.A_yfcom2,1),size(gradients.A_yfcom2,2)) ...
                zeros(size(gradients.A_xt2_ssp,1)-size(gradients.A_yfcom2,1),size(gradients.A_xt2_ssp,2)-size(gradients.A_yfcom2,2))] ...
                +walking_param.mu*gradients.A_xt2_ssp;
            yH=[yH_ssp zeros(size(yH_ssp,1),size(gradients.A_xt2_zmp1,2)-size(yH_ssp,2)); ...
                zeros(size(gradients.A_xt2_zmp1,1)-size(yH_ssp,1),size(gradients.A_xt2_zmp1,2))] ...
                +0.5*walking_param.mu*(gradients.A_xt2_zmp1+gradients.A_xt2_zmp2) ...
                +walking_param.epsilon*(gradients.A_yzmp1_acc2+gradients.A_yzmp2_acc2);
            
            obj.H=[xH zeros(size(yH));zeros(size(xH)) yH];
            
            %% Compute G
            xG_ssp=walking_param.lambda*[gradients.B_xfcom2 zeros(1,size(gradients.B_yt2_ssp,2)-size(gradients.B_xfcom2,2))] ...
                +walking_param.mu*gradients.B_yt2_ssp;
            xG=[xG_ssp zeros(1,size(gradients.B_yt2_zmp1,2)-size(xG_ssp,2))] ...
                +0.5*walking_param.mu*(gradients.B_yt2_zmp1+gradients.B_yt2_zmp2) ...
                +walking_param.epsilon*(gradients.B_xzmp1_acc2+gradients.B_xzmp2_acc2);
            
            yG_ssp=walking_param.lambda*[gradients.B_yfcom2 zeros(1,size(gradients.B_xt2_ssp,2)-size(gradients.B_yfcom2,2))] ...
                +walking_param.mu*gradients.B_xt2_ssp;
            yG=[yG_ssp zeros(1,size(gradients.B_xt2_zmp1,2)-size(yG_ssp,2))] ...
                +0.5*walking_param.mu*(gradients.B_xt2_zmp1+gradients.B_xt2_zmp2) ...
                +walking_param.epsilon*(gradients.B_yzmp1_acc2+gradients.B_yzmp2_acc2);
            
            obj.f=[xG yG]';
            obj.G=[xG yG];
            
            %% Compute C
            xC=walking_param.lambda*(gradients.C_xfcom2) ...
                +(1-walking_param.lambda)*(gradients.C_yt2_ssp) ...
                +0.5*(1-walking_param.lambda)*(gradients.C_yt2_zmp1+gradients.C_yt2_zmp2) ...    
                +walking_param.epsilon*(gradients.C_xzmp1_acc2+gradients.C_xzmp2_acc2);
            
            yC=walking_param.lambda*(gradients.C_yfcom2) ...
                +(1-walking_param.lambda)*(gradients.C_xt2_ssp) ...
                +0.5*(1-walking_param.lambda)*(gradients.C_xt2_zmp1+gradients.C_xt2_zmp2) ...    
                +walking_param.epsilon*(gradients.C_yzmp1_acc2+gradients.C_yzmp2_acc2);
            
            obj.C=xC+yC;
        end
        function obj=compute_H_G_C(obj,walking_param,gradients)
            %compute the coefficient of the cost function
            %cost=1/2*X^T*H*X+G*X+1/2*C
            %cost_quad=1/2*X^T*H*X+f'*X=1/2*X^T*H*X+G'*X
            %% Compute H
            xH_ssp=walking_param.lambda*[gradients.A_xfcom2 zeros(size(gradients.A_xfcom2,1),size(gradients.A_yt2_ssp,2)-size(gradients.A_xfcom2,2)); ...
                zeros(size(gradients.A_yt2_ssp,1)-size(gradients.A_xfcom2,1),size(gradients.A_xfcom2,2)) ... 
                zeros(size(gradients.A_yt2_ssp,1)-size(gradients.A_xfcom2,1),size(gradients.A_yt2_ssp,2)-size(gradients.A_xfcom2,2))] ...
                +walking_param.mu*gradients.A_yt2_ssp;
            xH=[xH_ssp zeros(size(xH_ssp,1),size(gradients.A_yt2_zmp1,2)-size(xH_ssp,2)); ...
                zeros(size(gradients.A_yt2_zmp1,1)-size(xH_ssp,1),size(gradients.A_yt2_zmp1,2))] ...
                +0.5*walking_param.mu*(gradients.A_yt2_zmp1+gradients.A_yt2_zmp2) ...
                +walking_param.epsilon*(gradients.A_xzmp1_acc2+gradients.A_xzmp2_acc2);
            
            yH_ssp=walking_param.lambda*[gradients.A_yfcom2 zeros(size(gradients.A_yfcom2,1),size(gradients.A_xt2_ssp,2)-size(gradients.A_yfcom2,2)); ...
                zeros(size(gradients.A_xt2_ssp,1)-size(gradients.A_yfcom2,1),size(gradients.A_yfcom2,2)) ...
                zeros(size(gradients.A_xt2_ssp,1)-size(gradients.A_yfcom2,1),size(gradients.A_xt2_ssp,2)-size(gradients.A_yfcom2,2))] ...
                +walking_param.mu*gradients.A_xt2_ssp;
            yH=[yH_ssp zeros(size(yH_ssp,1),size(gradients.A_xt2_zmp1,2)-size(yH_ssp,2)); ...
                zeros(size(gradients.A_xt2_zmp1,1)-size(yH_ssp,1),size(gradients.A_xt2_zmp1,2))] ...
                +0.5*walking_param.mu*(gradients.A_xt2_zmp1+gradients.A_xt2_zmp2) ...
                +walking_param.epsilon*(gradients.A_yzmp1_acc2+gradients.A_yzmp2_acc2);
            
            obj.H=[xH zeros(size(yH));zeros(size(xH)) yH];
            
            %% Compute G
            xG_ssp=walking_param.lambda*[gradients.B_xfcom2 zeros(1,size(gradients.B_yt2_ssp,2)-size(gradients.B_xfcom2,2))] ...
                +walking_param.mu*gradients.B_yt2_ssp;
            xG=[xG_ssp zeros(1,size(gradients.B_yt2_zmp1,2)-size(xG_ssp,2))] ...
                +0.5*walking_param.mu*(gradients.B_yt2_zmp1+gradients.B_yt2_zmp2) ...
                +walking_param.epsilon*(gradients.B_xzmp1_acc2+gradients.B_xzmp2_acc2);
            
            yG_ssp=walking_param.lambda*[gradients.B_yfcom2 zeros(1,size(gradients.B_xt2_ssp,2)-size(gradients.B_yfcom2,2))] ...
                +walking_param.mu*gradients.B_xt2_ssp;
            yG=[yG_ssp zeros(1,size(gradients.B_xt2_zmp1,2)-size(yG_ssp,2))] ...
                +0.5*walking_param.mu*(gradients.B_xt2_zmp1+gradients.B_xt2_zmp2) ...
                +walking_param.epsilon*(gradients.B_yzmp1_acc2+gradients.B_yzmp2_acc2);
            
            obj.f=[xG yG]';
            obj.G=[xG yG];
            
            %% Compute C
            xC=walking_param.lambda*(gradients.C_xfcom2) ...
                +(1-walking_param.lambda)*(gradients.C_yt2_ssp) ...
                +0.5*(1-walking_param.lambda)*(gradients.C_yt2_zmp1+gradients.C_yt2_zmp2) ...    
                +walking_param.epsilon*(gradients.C_xzmp1_acc2+gradients.C_xzmp2_acc2);
            
            yC=walking_param.lambda*(gradients.C_yfcom2) ...
                +(1-walking_param.lambda)*(gradients.C_xt2_ssp) ...
                +0.5*(1-walking_param.lambda)*(gradients.C_xt2_zmp1+gradients.C_xt2_zmp2) ...    
                +walking_param.epsilon*(gradients.C_yzmp1_acc2+gradients.C_yzmp2_acc2);
            
            obj.C=xC+yC;
        end
        function obj=compute_A_b(obj,walking_param,gradients)
            %% %%%contrainte de non chevauchement des pieds%%%
            [gradients.A_cons_stretch gradients.B_cons_stretch]=obj.cons_stretching(walking_param.pankinit_firstSS,walking_param.pankfin_lastSS,walking_param.backtoankle,walking_param.fronttoankle,walking_param.exttoankle,walking_param.inttoankle,walking_param.xankmax,walking_param.xankmin,walking_param.yankmax,walking_param.yankmin,walking_param.psi,walking_param.nbparamABCD,walking_param.nbparamank,walking_param.nbparamBb);
            
            obj.A=[gradients.A_cons_zmp_stability;    gradients.A_cons_zmp1_stab;   gradients.A_cons_zmp2_stab;   gradients.A_cons_stretch];
            obj.b=[gradients.B_cons_zmp_stability;    gradients.B_cons_zmp1_stab;   gradients.B_cons_zmp2_stab;   gradients.B_cons_stretch];
        end
        function obj=compute_Aeq_beq(obj,walking_param,gradients)
            %%%%Contrainte sur des ankle positions%%%
            [gradients.A_cons_ankle gradients.B_cons_ankle]=obj.cons_ankle_fixed_path(walking_param.nbparamABCD,walking_param.nbparamank,walking_param.nbparamBb,walking_param.step_number_pankle_fixed);
            
            Aeq=[gradients.A_cons_scom;gradients.A_cons_ankle];
            beq=[gradients.B_cons_scom;gradients.B_cons_ankle];
            
            %% Add constraint to obtain a cyclical walk
            %force a distance in x of 0.125m between the initial and final
            %via-point
            %force the equality in x in speed and acceleration between the
            %initial and final via-point
            %force a symetry in y in position, speed and acceleration
            %between the initial and final via-point
            %1 4 7
            xpA=[[-1 0 0] zeros(1,6) [1 0 0] zeros(1,4+walking_param.nbparamank) zeros(1,walking_param.nbparamBb) zeros(1,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb)];
            Aeq=[Aeq;xpA];
            beq=[beq;-0.125];
            ypA=[zeros(1,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb) [1 0 0] zeros(1,6) [1 0 0] zeros(1,4+walking_param.nbparamank) zeros(1,walking_param.nbparamBb)];
            Aeq=[Aeq;ypA];
            beq=[beq;0];

            xsA=[[0 -1 0] zeros(1,6) [0 1 0] zeros(1,4+walking_param.nbparamank) zeros(1,walking_param.nbparamBb) zeros(1,4+walking_param.nbparamABCD+walking_param.nbparamank+walking_param.nbparamBb)];
            Aeq=[Aeq;xsA];
            beq=[beq;0];
            ysA=[zeros(1,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb) [0 1 0] zeros(1,6) [0 1 0] zeros(1,4+walking_param.nbparamank) zeros(1,walking_param.nbparamBb) ];
            Aeq=[Aeq;ysA];
            beq=[beq;0];

            xaA=[[0 0 -1] zeros(1,6) [0 0 1] zeros(1,4+walking_param.nbparamank) zeros(1,walking_param.nbparamBb) zeros(1,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb)];
            Aeq=[Aeq;xaA];
            beq=[beq;0];
            yaA=[zeros(1,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb) [0 0 1] zeros(1,6) [0 0 1] zeros(1,4+walking_param.nbparamank) zeros(1,walking_param.nbparamBb) ];
            Aeq=[Aeq;yaA];
            beq=[beq;0];
            
            %force a distance in x of 0.125m between the initial and final COM
            %force an equal speed in x between the initial and final COM
            %force a symetry in position and speed in y between the initial
            %and final COM
            xpA=[zeros(2,walking_param.nbparamABCD) [-1 1 0 0;0 0 -1 1] zeros(2,walking_param.nbparamank) zeros(2,walking_param.nbparamBb) zeros(2,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb)];
            Aeq=[Aeq;xpA];
            beq=[beq;-0.125;0];
            ypA=[zeros(2,walking_param.nbparamABCD+4+walking_param.nbparamank+walking_param.nbparamBb) zeros(2,walking_param.nbparamABCD) [1 1 0 0;0 0 1 1] zeros(2,walking_param.nbparamank) zeros(2,walking_param.nbparamBb)];
            Aeq=[Aeq;ypA];
            beq=[beq;0;0];

            obj.Aeq=Aeq;
            obj.beq=beq;
        end
        
        function [cost gradient]=compute_cost(obj,x)
            cost=1/2*x'*obj.H*x+obj.G*x+obj.C;
            gradient=x'*obj.H+obj.G;
        end
    end
    methods (Static) 
        %% intern methods of compute_ZMP
        function [A_Ca B_Ca]=compute_zmp_Ca(tpassage,psa_zmpinit,psa_zmpfin,nbphases,nbparamABCD)
        % Hypo : We consider ZMP trajectory as 5th order polynomials between via-points.
        % This function compute A and B as: C=A*x+B
        % Where : C are give the polynomials coefficient of ZMP trajectory in one direction.
        % x are the via-point boundary conditions in position, speed and acceleration.
        % Initial and final via-point boundary conditions are known and in psa_zmpinit and psa_zmpfin

            A_Ca=zeros(6*nbphases,nbparamABCD+6);
            for i=1:nbphases
                dt=tpassage(i+1)-tpassage(i);
                m=[1 0 0 0 0 0;
                   0 1 0 0 0 0;
                   0 0 2 0 0 0;
                   1 dt dt^2 dt^3 dt^4 dt^5;
                   0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                   0 0 2 6*dt 12*dt^2 20*dt^3];
                m=inv(m);
                A_Ca(1+(i-1)*6:6+(i-1)*6,:)=[zeros(6,3*(i-1)) m zeros(6,(nbphases-i)*3)];
            end

            B_Ca=[A_Ca(:,1:3) A_Ca(:,end-2:end)]*[psa_zmpinit;psa_zmpfin];

            A_Ca(:,1:3)=[];
            A_Ca(:,end-2:end)=[];
            A_Ca=[A_Ca zeros(size(A_Ca,1),4)];
        end
        %% intern methods of compute_COM
        function [M_cs]=compute_com_M_cs(discretization,w,frequency,nbphases)
            %Compute the matrix M_cs
            M_cs= zeros(1+sum(discretization),(nbphases)*2);
            M_cs(1,1:2)=[cosh(0) sinh(0)];
            nb=0;
            rowinitphase=1;
            for j=1:nbphases
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency;%dt is delta(tj)=t-Tj
                    M_cs(rowinitphase+i,(j-1)*2+1:(j-1)*2+2)=[cosh(w*dt) sinh(w*dt)];
                end
            end
        end
        function [A_A]=compute_com_A_A(nbphases,w)
        %Compute the matrix A_A use to compute the A coefficient from zmp
        %polynomial coefficients
            A_A=[];
            for j=1:nbphases
                Aj=[1 0 2/w^2 0 24/w^4 0;
                    0 1 0 6/w^2 0 120/w^4;
                    0 0 1 0 12/w^2 0;
                    0 0 0 1 0 20/w^2;
                    0 0 0 0 1 0;
                    0 0 0 0 0 1];
                A_A=[A_A;zeros(6,(j-1)*6) Aj zeros(6,(nbphases-j)*6)];
            end
        end
        function [A_y B_y]=compute_com_y(A_Ca,B_Ca,A_A,tpassage,w)
            %put a system G*y=N*x0+H*l as y=G^-1*N*x+G^-1*H*l
            %x0 are the com initial and final positions and optimization
            %variables
            %l=A_l*X+B_l
            %A_l=[A_A.A_Ca 0]
            %B_l=[A_A.B_Ca]
            
            G=wpg_qp_problem_Adrien.compute_com_G(tpassage,w);
            
            H=wpg_qp_problem_Adrien.compute_com_H(tpassage);
            
            N=[1 0 0 0; zeros(size(H,1)-2,4); 0 1 0 0]; %pcominit pconfin scominit scomfin
            
            A_l=A_A*A_Ca;
            B_l=A_A*B_Ca;

            A_global=G\(H*A_A);
            
            A_y=[zeros(size(A_global,1),size(A_Ca,2)-4) G\(N)]+G\(H*A_l);
            B_y=G\(H*B_l);
        end
        function [G]=compute_com_G(tpassage,w)
        %Compute the G matrix
            nbphases=length(tpassage)-1;

            G=[1 0 zeros(1,2*nbphases-2)];

            for i=1:nbphases-1
                Dt=tpassage(i+1)-tpassage(i);
                G=[G;zeros(1,2*(i-1)) cosh(w*Dt) sinh(w*Dt) -1 0 zeros(1,2*(nbphases-i)-2);zeros(1,2*(i-1)) w*sinh(w*Dt) w*cosh(w*Dt) 0 -w zeros(1,2*(nbphases-i)-2)];
            end

            Dt=tpassage(nbphases+1)-tpassage(nbphases);
            G=[G;zeros(1,2*nbphases-2) cosh(w*Dt) sinh(w*Dt)];
        end
        function [H]=compute_com_H(tpassage)
        %compute the H matrix
            nbphases=length(tpassage)-1;

            H=[-1 0 zeros(1,6*nbphases-2)];

            for i=1:nbphases-1
                Dt=tpassage(i+1)-tpassage(i);
                H=[H;zeros(1,6*(i-1)) -Dt^0 -Dt^1 -Dt^2 -Dt^3 -Dt^4 -Dt^5 +1 0 zeros(1,6*(nbphases-i)-2);zeros(1,6*(i-1)) 0 -1*Dt^0 -2*Dt^1 -3*Dt^2 -4*Dt^3 -5*Dt^4 0 +1 zeros(1,6*(nbphases-i)-2)];
            end

            Dt=tpassage(nbphases+1)-tpassage(nbphases);
            H=[H;zeros(1,6*(nbphases-1)) -Dt^0 -Dt^1 -Dt^2 -Dt^3 -Dt^4 -Dt^5];
        end
        function [A_com B_com]=compute_com(A_Ca,B_Ca,A_y,B_y,A_A,M_t,M_cs)
        %%%generate the COM trajectory gradient%%%
        %Use Morisawa algorithm with 5th zmp polynomial coeff known
            A_com=M_cs*A_y+M_t*A_A*A_Ca;
            B_com=M_cs*B_y+M_t*A_A*B_Ca;
        end
        function [A_cons_scom B_cons_scom] = cons_scom_initfin(A_Ca,B_Ca,A_y,B_y,A_A,tpassage,nbphases,w)
            dti=0;%initial time of the first phase
            dtf=tpassage(nbphases+1)-tpassage(nbphases);%fial time step of the last phase
            M_cs=[w*sinh(w*dti) w*cosh(w*dti) zeros(1,(nbphases-1)*2);
                zeros(1,(nbphases-1)*2) w*sinh(w*dtf) w*cosh(w*dtf)];
            M_t=[0 1 2*dti 3*dti^2 4*dti^3 5*dti^4 zeros(1,(nbphases-1)*6);
                zeros(1,(nbphases-1)*6) 0 1 2*dtf 3*dtf^2 4*dtf^3 5*dtf^4];

            A_cons_scom=M_cs*A_y+M_t*A_A*A_Ca-[zeros(2,size(A_Ca,2)-2) [1 0;0 1]];
            B_cons_scom=M_cs*B_y+M_t*A_A*B_Ca;
        end
        %% intern methods of compute_cost_force
        function [A_fcom B_fcom] = compute_force(A_zmp,B_zmp,A_com,B_com,mg,h)
            A_fcom=(A_zmp-A_com)*(-mg)/h;
            B_fcom=(B_zmp-B_com)*(-mg)/h;
        end
        %% intern methods of compute_cost_torque_ssp
        function [A_t_ssp B_t_ssp] = compute_torque_ssp(Apzmp,Bpzmp,Afcom,Bfcom,Apankle,Bpankle,mg,ha)
        %le couple en x est fonction des y
        %le couple en y est fonction des x

            A1=ha*Afcom+mg*Apzmp;
            A2=-mg*Apankle;
            A_t_ssp=[A1 zeros(size(A2,1),size(A2,2))]+[zeros(size(A1,1),size(A1,2)) A2];
            B_t_ssp=ha*Bfcom+mg*Bpzmp-mg*Bpankle;

            A_t_ssp=A_t_ssp(any(any(Apankle,2)+any(Bpankle,2),2),:);
            B_t_ssp=B_t_ssp(any(any(Apankle,2)+any(Bpankle,2),2),:);
        end
        function [A_an B_an] = compute_ankle_positions_SSP(pankinit,pankfin,discretization)
            nbphases=length(discretization);
            A_an=zeros(1+sum(discretization),ceil(nbphases/3));
            %pankle(1)=0;
            nb=0;
            rowinitphase=1;
            for j=1:nbphases
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                if(mod(j,3)==0 || j==1 || j==nbphases || j==2 || j==nbphases-1)
            %         for i=1:nb
            %             pankle=[pankle;0];
            %         end
                else
                    for i=1:nb
                    A_an(rowinitphase+i,ceil(j/3))=1;
                    end
                end
            end

            A_an=A_an(:,any(A_an,1));

            B_an=zeros(1+sum(discretization),1);
            B_an(1+discretization(1)+1:1+sum(discretization(1:2)))=pankinit*ones(discretization(2),1);
            B_an(1+sum(discretization(1:end-2))+1:1+sum(discretization(1:end-1)))=pankfin*ones(discretization(end-1),1);

        end
        %% intern methods of compute_cons_zmp_com
        function [Acons Bcons] = cons_xy_zmp_stability_SSP(xApzmp,xBpzmp,yApzmp,yBpzmp,xApankle,xBpankle,yApankle,yBpankle,discretization,backtoankle,fronttoankle,exttoankle,inttoankle,sole_margin,theta,type_phase,nbparamBb,firstSS)
            xApzmp_=[xApzmp zeros(size(xApankle)) zeros(size(xApzmp,1),nbparamBb)];blabla.toto
            yApzmp_=[yApzmp zeros(size(yApankle)) zeros(size(yApzmp,1),nbparamBb)];
            xApankle=[zeros(size(xApzmp)) xApankle zeros(size(xApzmp,1),nbparamBb) zeros(size(yApzmp_))];
            yApankle=[zeros(size(xApzmp_)) zeros(size(yApzmp)) yApankle zeros(size(yApzmp,1),nbparamBb)];

            xApzmp_=[xApzmp_ zeros(size(xApzmp_))];
            yApzmp_=[zeros(size(yApzmp_)) yApzmp_];

            Acons1=zeros(size(xApankle));
            Bcons1=zeros(size(xBpankle));
            Acons2=zeros(size(xApankle));
            Bcons2=zeros(size(xBpankle));
            Acons3=zeros(size(xApankle));
            Bcons3=zeros(size(xBpankle));
            Acons4=zeros(size(xApankle));
            Bcons4=zeros(size(xBpankle));

            j=1;
            %constraint direction inverse clock-wise
            t=theta(j);
            dx=[cos(t) -sin(t) -cos(t) sin(t)];
            dy=[sin(t) cos(t) -sin(t) -cos(t)];


            xA=xApzmp_-xApankle;
            yA=yApzmp_-yApankle;
            xB=xBpzmp-xBpankle;
            yB=yBpzmp-yBpankle;

            for i=1:sum(discretization)
                if(i>sum(discretization(1:j)))
                        j=j+1;
                        %constraint direction inverse clock-wise
                        t=theta(j);
                        dx=[cos(t) -sin(t) -cos(t) sin(t)];
                        dy=[sin(t) cos(t) -sin(t) -cos(t)];
                end

                if type_phase(j)~=0
                    Acons1(i,:)=dx(1)*(xA(i,:))+dy(1)*(yA(i,:));
                    Bcons1(i)=-dx(1)*(xB(i,:))-dy(1)*(yB(i,:))+fronttoankle-sole_margin;
                    Acons2(i,:)=dx(2)*(xA(i,:))+dy(2)*(yA(i,:));
                    Bcons2(i)=-dx(2)*(xB(i,:))-dy(2)*(yB(i,:));
                    Acons3(i,:)=dx(3)*(xA(i,:))+dy(3)*(yA(i,:));
                    Bcons3(i)=-dx(3)*(xB(i,:))-dy(3)*(yB(i,:))+backtoankle-sole_margin;
                    Acons4(i,:)=dx(4)*(xA(i,:))+dy(4)*(yA(i,:));
                    Bcons4(i)=-dx(4)*(xB(i,:))-dy(4)*(yB(i,:));
                    switch (j)
                        case 1
                            Bcons2(i)=Bcons2(i)+inttoankle-sole_margin;
                            Bcons4(i)=Bcons4(i)+exttoankle-sole_margin;
                        case 3
                            Bcons2(i)=Bcons2(i)+exttoankle-sole_margin;
                            Bcons4(i)=Bcons4(i)+inttoankle-sole_margin;
                        case 4
                            Bcons2(i)=Bcons2(i)+exttoankle-sole_margin;
                            Bcons4(i)=Bcons4(i)+inttoankle-sole_margin;
                        case 6
                            Bcons2(i)=Bcons2(i)+inttoankle-sole_margin;
                            Bcons4(i)=Bcons4(i)+exttoankle-sole_margin;
                    end
                end
            end

            Acons=[Acons1(any(any(Acons1,2)+any(Bcons1,2),2),:);
                   Acons2(any(any(Acons2,2)+any(Bcons2,2),2),:);
                   Acons3(any(any(Acons3,2)+any(Bcons3,2),2),:);
                   Acons4(any(any(Acons4,2)+any(Bcons4,2),2),:)];
            Bcons=[Bcons1(any(any(Acons1,2)+any(Bcons1,2),2),:);
                   Bcons2(any(any(Acons2,2)+any(Bcons2,2),2),:);
                   Bcons3(any(any(Acons3,2)+any(Bcons3,2),2),:);
                   Bcons4(any(any(Acons4,2)+any(Bcons4,2),2),:)];
        end
        %% intern methods of compute_ZMP1
        function [A_Ca1 B_Ca1]=compute_zmp1_Ca1(tpassage,psa_zmp1init,psa_zmp1fin,cutting)
            nbpoly=length(tpassage)-1;
            nbpoly1=(length(tpassage)-1-2)/3+2;
            g=[];

            m22=zeros(6,nbpoly1-3);

            %zmp1 init
            dt=tpassage(2)-tpassage(1);
            dt=dt*cutting;
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            m=inv(m);
            m0=zeros(6,0);
            m1=m(1:6,1:3);
            m2=zeros(6,(nbpoly)*3);

            m3=zeros(6,0);
            m4=m(1:6,4:6);
            m5=zeros(6,3*(nbpoly1-1)+3+3);

            g=[g;m0 m1 m2 m22 m3 m4 m5];

            %decoupage de zmp1 init en deux parties
            dt=tpassage(2)-tpassage(1);
            dt=dt*(1-cutting);
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            m=inv(m);
            m0=zeros(6,0);
            m1=zeros(6,0);
            m2=zeros(6,(nbpoly)*3+3);

            m3=zeros(6,0);
            m4=m;
            m5=zeros(6,3*(nbpoly1-1)+3);

            g=[g;m0 m1 m2 m22 m3 m4 m5];

            for i=2:nbpoly-1
                dt=tpassage(i+1)-tpassage(i);
                    m=[1 0 0 0 0 0;
                       0 1 0 0 0 0;
                       0 0 2 0 0 0;
                       1 dt dt^2 dt^3 dt^4 dt^5;
                       0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                       0 0 2 6*dt 12*dt^2 20*dt^3];
                m=inv(m);

                m0=zeros(6,3*(i-1));
                m1=m(1:6,1:3);
                m2=zeros(6,(nbpoly)*3-(i-1)*3);

                m3=zeros(6,(i+1-mod(i+1,3))+3);
                m4=m(1:6,4:6);
                m5=zeros(6,(nbpoly1-1)*3-(i+1-mod(i+1,3))+3);

                g=[g;m0 m1 m2 m22 m3 m4 m5];
            end

            %zmp1 fin
            dt=tpassage(end)-tpassage(end-1);
            dt=dt*cutting;
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            m=inv(m);

            m0=zeros(6,(nbpoly)*3-3);
            m1=m(1:6,1:3);
            m2=zeros(6,3);

            m3=zeros(6,3*(nbpoly1-1)+3);
            m4=m(1:6,4:6);
            m5=zeros(6,3);

            g=[g;m0 m1 m2 m22 m3 m4 m5];

            %decoupage de zmp1 fin en deux parties
            dt=tpassage(end)-tpassage(end-1);
            dt=dt*(1-cutting);
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            m=inv(m);

            m0=zeros(6,(nbpoly)*3-3);
            m1=zeros(6,0);
            m2=zeros(6,3+3);

            m3=zeros(6,3*(nbpoly1-1)+3);
            m4=m;
            m5=zeros(6,0);

            g=[g;m0 m1 m2 m22 m3 m4 m5];

            A_Ca1=[g(:,4:(nbpoly)*3+(nbpoly1-3)) zeros(size(g,1),4) g(:,(nbpoly)*3+(nbpoly1-3)+1+3:end-3)];
            B_Ca1=[g(:,1:3)]*psa_zmp1init+[g(:,end-2:end)]*psa_zmp1fin;
        end
        %% intern methods of compute_cost_torque_zmp1
        function [Ac]=compute_force_repartition(tpassage,cutting)
            nbpoly=length(tpassage)-1;
            g=[];

            %%%cutting fcom1 init in two%%%
            dt=tpassage(2)-tpassage(1);
            dt=dt*cutting;
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            %m=inv(m);
            m=m\[0.5;0;0;0.75;0;0];
            g=[g;m];

            %second part
            dt=tpassage(2)-tpassage(1);
            dt=dt*(1-cutting);
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            % m=inv(m);
            m=m\[0.75;0;0;0;0;0];
            g=[g;m];


            for i=2:nbpoly-1
                dt=tpassage(i+1)-tpassage(i);
                m=[1 0 0 0 0 0;
                   0 1 0 0 0 0;
                   0 0 2 0 0 0;
                   1 dt dt^2 dt^3 dt^4 dt^5;
                   0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                   0 0 2 6*dt 12*dt^2 20*dt^3];
                %m=inv(m);
                m=m\[1;0;0;0;0;0];

                g=[g;m];
            end

            %%%cutting fcom1 fin in two%%%
            dt=tpassage(end)-tpassage(end-1);
            dt=dt*(1-cutting);
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            % m=inv(m);
            m=m\[1;0;0;0.25;0;0];
            g=[g;m];

            %second part
            dt=tpassage(end)-tpassage(end-1);
            dt=dt*cutting;
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            % m=inv(m);
            m=m\[0.25;0;0;0.5;0;0];
            g=[g;m];

            Ac=g;
            % A=g(:,4:size(g,2)-3);
            % B=[g(:,1:3) g(:,size(g,2)-2:size(g,2))]*[psa_zmpinit;psa_zmpfin];

        end
        function [A_t_zmp1 B_t_zmp1] = compute_torque_zmp1(A_zmp1,B_zmp1,A_fcom,B_fcom,A_an,B_an,mg,ha,Afb)
        %le couple en x est fonction des y
        %le couple en y est fonction des x
            A1=ha*Afb*([A_fcom zeros(size(A_fcom,1),size(A_zmp1,2)-size(A_fcom,2))]);
            A2=+Afb*mg*A_zmp1;
            A3=[zeros(size(A_fcom)) -Afb*mg*A_an zeros(size(A_fcom,1),size(A_zmp1,2)-size(A_fcom,2)-size(A_an,2))];
            A_t_zmp1=A1+A2+A3;
            B_t_zmp1=ha*Afb*(B_fcom)+Afb*mg*B_zmp1-Afb*mg*B_an;

            A_t_zmp1=A_t_zmp1(any(any(A_an,2)+any(B_an,2),2),:);
            B_t_zmp1=B_t_zmp1(any(any(A_an,2)+any(B_an,2),2),:);
        end
        function [A_an B_an] = compute_ankle_positions_zmp1(pankinit1,pankinit2,pankfin1,discretization)
            nbphases=length(discretization);

            A_an=zeros(1+sum(discretization),ceil(nbphases/3));
            nb=0;
            rowinitphase=1;
            for j=1:nbphases-2
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                if(mod(j,3)==0&&j>3)
                    for i=1:nb
                         A_an(rowinitphase+i,ceil(j/3))=1;
                    end
                else
            %         for i=1:nb
            %         pankle(rowinitphase+i)=pstep(ceil(j/3));
            %         end
                end
            end
            A_an=A_an(:,any(A_an,1));

            B_an=zeros(1+sum(discretization),1);
            B_an(1)=pankinit1;
            B_an(1+1:1+discretization(1))=pankinit1*ones(discretization(1),1);
            B_an(1+sum(discretization(1:2))+1:1+sum(discretization(1:3)))=pankinit2*ones(discretization(3),1);
            B_an(1+sum(discretization(1:end-1))+1:1+sum(discretization(1:end)))=pankfin1*ones(discretization(end),1);

        end
        %% intern methods of compute_ZMP2
        function [A_zmp2 B_zmp2]=compute_zmp2(A_zmp,B_zmp,A_zmp1,B_zmp1,k_diag,mg)
            xf2=(eye(length(k_diag))-k_diag)*mg;
            A_zmp2=-mg*(xf2\(A_zmp1-[A_zmp zeros(size(A_zmp,1),size(A_zmp1,2)-size(A_zmp,2))]))+A_zmp1;
            B_zmp2=-mg*(xf2\(B_zmp1-B_zmp))+B_zmp1;
        end
        %% intern methods of compute_cost_torque_zmp2
        function [A_an B_an] = compute_ankle_positions_zmp2(pankinit2,pankfin1,pankfin2,discretization)
            nbphases=length(discretization);
            A_an=zeros(1+sum(discretization),ceil(nbphases/3));
            nb=0;
            rowinitphase=1;
            for j=1:nbphases-3
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                if(mod(j,3)==0&&j>=3)
                    for i=1:nb
                         A_an(rowinitphase+i,ceil(j/3))=1;
                    end
                else
            %         for i=1:nb
            %         pankle(rowinitphase+i)=pstep(ceil(j/3));
            %         end
                end
            end
            A_an=A_an(:,any(A_an,1));

            B_an=zeros(1+sum(discretization),1);
            B_an(1)=pankinit2;
            B_an(1+1:1+discretization(1))=pankinit2*ones(discretization(1),1);
            B_an(1+sum(discretization(1:end-3))+1:1+sum(discretization(1:end-2)))=pankfin1*ones(discretization(end-2),1);
            B_an(1+sum(discretization(1:end-1))+1:1+sum(discretization(1:end)))=pankfin2*ones(discretization(end),1);
        end
        function [A_t_zmp2 B_t_zmp2] = compute_torque_zmp2(A_zmp2,B_zmp2,A_fcom,B_fcom,A_an,B_an,mg,ha,Afb)
        %le couple en x est fonction des y
        %le couple en y est fonction des x
            IAfb=eye(size(Afb,1))-Afb;

            A1=ha*IAfb*([A_fcom zeros(size(A_fcom,1),size(A_zmp2,2)-size(A_fcom,2))]);
            A2=+IAfb*mg*A_zmp2;
            A3=[zeros(size(A_fcom)) -IAfb*mg*A_an zeros(size(A_fcom,1),size(A_zmp2,2)-size(A_fcom,2)-size(A_an,2))];
            A_t_zmp2=A1+A2+A3;

            B_t_zmp2=ha*IAfb*(B_fcom)+IAfb*mg*B_zmp2-IAfb*mg*B_an;
        end
        %% intern methods of compute_cost_acc_zmp2
        function [A_zmp2_acc B_zmp2_acc]=compute_zmp2_acc(A_zmp,B_zmp,Apzmp1,B_zmp1,A_zmp_s,B_zmp_s,A_zmp1_s,B_zmp1_s,A_zmp_a,B_zmp_a,A_zmp1_a,B_zmp1_a,k_diag,k_diag_d,k_diag_dd)
        %Compute the gradient matrix of zmp2 acceleration
            Afb_=(eye(length(k_diag))-k_diag);
            dAfb_=-k_diag_d;
            ddAfb_=-k_diag_dd;

            A1=-(Afb_\(A_zmp1_a-[A_zmp_a zeros(size(A_zmp_a,1),size(A_zmp1_a,2)-size(A_zmp_a,2))]));
            A2=+2*dAfb_*((Afb_)\((Afb_)\(A_zmp1_s-[A_zmp_s zeros(size(A_zmp_s,1),size(A_zmp1_s,2)-size(A_zmp_s,2))])));
            A3=+ddAfb_*((Afb_)\((Afb_)\(Apzmp1-[A_zmp zeros(size(A_zmp,1),size(Apzmp1,2)-size(A_zmp,2))])))-2*dAfb_*dAfb_*((Afb_)\((Afb_)\((Afb_)\(Apzmp1-[A_zmp zeros(size(A_zmp,1),size(Apzmp1,2)-size(A_zmp,2))]))));

            B1=-(Afb_\(B_zmp1_a-B_zmp_a));
            B2=+2*dAfb_*((Afb_)\((Afb_)\(B_zmp1_s-B_zmp_s)));
            B3=+ddAfb_*((Afb_)\((Afb_)\(B_zmp1-B_zmp)))-2*dAfb_*dAfb_*((Afb_)\((Afb_)\((Afb_)\(B_zmp1-B_zmp))));

            A_zmp2_acc=A1+A2+A3+A_zmp1_a;
            B_zmp2_acc=B1+B2+B3+B_zmp1_a;
        end
        %% intern methods of compute_cons_zmp12
        function [Acons Bcons] = cons_xy_zmp1_stability(xApzmp,xBpzmp,yApzmp,yBpzmp,xApankle,xBpankle,yApankle,yBpankle,discretization,backtoankle,fronttoankle,exttoankle,inttoankle,sole_margin,theta,type_phase,nbparamABCD,nbparamBb,firstSS)
            xApzmp_=[xApzmp];
            yApzmp_=[yApzmp];
            xApankle_=[zeros(size(xApzmp_,1),nbparamABCD+4) xApankle zeros(size(xApzmp,1),nbparamBb) zeros(size(yApzmp_))];
            yApankle_=[zeros(size(xApzmp_)) zeros(size(yApzmp_,1),nbparamABCD+4) yApankle zeros(size(yApzmp,1),nbparamBb)];

            xApzmp_=[xApzmp_ zeros(size(xApzmp_))];
            yApzmp_=[zeros(size(yApzmp_)) yApzmp_];

            Acons1=zeros(size(xApankle_));
            Bcons1=zeros(size(xBpankle));
            Acons2=zeros(size(xApankle_));
            Bcons2=zeros(size(xBpankle));
            Acons3=zeros(size(xApankle_));
            Bcons3=zeros(size(xBpankle));
            Acons4=zeros(size(xApankle_));
            Bcons4=zeros(size(xBpankle));

            j=1;
            %constraint direction inverse clock-wise
            theta_=theta(1:end-1);
            theta_=theta_(any(type_phase==0,1));

            t=theta_(j);
            dx=[cos(t) -sin(t) -cos(t) sin(t)];
            dy=[sin(t) cos(t) -sin(t) -cos(t)];

            xA=xApzmp_-xApankle_;
            yA=yApzmp_-yApankle_;
            xB=xBpzmp-xBpankle;
            yB=yBpzmp-yBpankle;

            discretization_=discretization(any(type_phase==0,1));

            for i=1:sum(discretization_)
                if(i>sum(discretization_(1:j)))
                        j=j+1;
                        %constraint direction inverse clock-wise
                        t=theta_(j);
                        dx=[cos(t) -sin(t) -cos(t) sin(t)];
                        dy=[sin(t) cos(t) -sin(t) -cos(t)];
                end

                    Acons1(i,:)=dx(1)*(xA(i,:))+dy(1)*(yA(i,:));
                    Bcons1(i)=-dx(1)*(xB(i,:))-dy(1)*(yB(i,:))+fronttoankle-sole_margin;
                    Acons2(i,:)=dx(2)*(xA(i,:))+dy(2)*(yA(i,:));
                    Bcons2(i)=-dx(2)*(xB(i,:))-dy(2)*(yB(i,:));
                    Acons3(i,:)=dx(3)*(xA(i,:))+dy(3)*(yA(i,:));
                    Bcons3(i)=-dx(3)*(xB(i,:))-dy(3)*(yB(i,:))+backtoankle-sole_margin;
                    Acons4(i,:)=dx(4)*(xA(i,:))+dy(4)*(yA(i,:));
                    Bcons4(i)=-dx(4)*(xB(i,:))-dy(4)*(yB(i,:));

                    if(firstSS==0)
                        Bcons2(i)=Bcons2(i)+(inttoankle*(mod(j,2)==1)+exttoankle*(mod(j,2)==0)-sole_margin);
                        Bcons4(i)=Bcons4(i)+(inttoankle*(mod(j,2)==0)+exttoankle*(mod(j,2)==1)-sole_margin);
                    elseif(firstSS==1)
                        Bcons2(i)=Bcons2(i)+(inttoankle*(mod(j,2)==0)+exttoankle*(mod(j,2)==1)-sole_margin);
                        Bcons4(i)=Bcons4(i)+(inttoankle*(mod(j,2)==1)+exttoankle*(mod(j,2)==0)-sole_margin);
                    else
                        'Choose a first SS foot'
                    end
            end

             Acons=[Acons1;
                 Acons2;
                 Acons3;
                 Acons4];
             Bcons=[Bcons1;
                 Bcons2;
                 Bcons3;
                 Bcons4];
        end
        function [Acons Bcons] = cons_xy_zmp2_stability(xApzmp,xBpzmp,yApzmp,yBpzmp,xApankle,xBpankle,yApankle,yBpankle,discretization,backtoankle,fronttoankle,exttoankle,inttoankle,sole_margin,theta,type_phase,nbparamABCD,nbparamBb,firstSS)
            xApzmp_=[xApzmp];
            yApzmp_=[yApzmp];
            xApankle=[zeros(size(xApzmp_,1),nbparamABCD+4) xApankle zeros(size(xApzmp,1),nbparamBb) zeros(size(yApzmp_))];
            yApankle=[zeros(size(xApzmp_)) zeros(size(yApzmp_,1),nbparamABCD+4) yApankle zeros(size(yApzmp,1),nbparamBb)];

            xApzmp_=[xApzmp_ zeros(size(xApzmp_))];
            yApzmp_=[zeros(size(yApzmp_)) yApzmp_];

            Acons1=zeros(size(xApankle));
            Bcons1=zeros(size(xBpankle));
            Acons2=zeros(size(xApankle));
            Bcons2=zeros(size(xBpankle));
            Acons3=zeros(size(xApankle));
            Bcons3=zeros(size(xBpankle));
            Acons4=zeros(size(xApankle));
            Bcons4=zeros(size(xBpankle));

            j=1;
            %constraint direction inverse clock-wise
            theta_=theta(2:end);
            theta_=theta_(any(type_phase==0,1));

            t=theta_(j);
            dx=[cos(t) -sin(t) -cos(t) sin(t)];
            dy=[sin(t) cos(t) -sin(t) -cos(t)];

            xA=xApzmp_-xApankle;
            yA=yApzmp_-yApankle;
            xB=xBpzmp-xBpankle;
            yB=yBpzmp-yBpankle;

            discretization_=discretization(any(type_phase==0,1));

            for i=1:sum(discretization_)
                if(i>sum(discretization_(1:j)))
                        j=j+1;
                        %constraint direction inverse clock-wise
                        t=theta_(j);
                        dx=[cos(t) -sin(t) -cos(t) sin(t)];
                        dy=[sin(t) cos(t) -sin(t) -cos(t)];
                end

                    Acons1(i,:)=dx(1)*(xA(i,:))+dy(1)*(yA(i,:));
                    Bcons1(i)=-dx(1)*(xB(i,:))-dy(1)*(yB(i,:))+fronttoankle-sole_margin;
                    Acons2(i,:)=dx(2)*(xA(i,:))+dy(2)*(yA(i,:));
                    Bcons2(i)=-dx(2)*(xB(i,:))-dy(2)*(yB(i,:));
                    Acons3(i,:)=dx(3)*(xA(i,:))+dy(3)*(yA(i,:));
                    Bcons3(i)=-dx(3)*(xB(i,:))-dy(3)*(yB(i,:))+backtoankle-sole_margin;
                    Acons4(i,:)=dx(4)*(xA(i,:))+dy(4)*(yA(i,:));
                    Bcons4(i)=-dx(4)*(xB(i,:))-dy(4)*(yB(i,:));

                    if(firstSS==0)
                        Bcons2(i)=Bcons2(i)+(inttoankle*(mod(j,2)==0)+exttoankle*(mod(j,2)==1)-sole_margin);
                        Bcons4(i)=Bcons4(i)+(inttoankle*(mod(j,2)==1)+exttoankle*(mod(j,2)==0)-sole_margin);
                    elseif(firstSS==1)
                        Bcons2(i)=Bcons2(i)+(inttoankle*(mod(j,2)==1)+exttoankle*(mod(j,2)==0)-sole_margin);
                        Bcons4(i)=Bcons4(i)+(inttoankle*(mod(j,2)==0)+exttoankle*(mod(j,2)==1)-sole_margin);
                    else
                        'Choose a first SS foot'
                    end
            end

             Acons=[Acons1;
                 Acons2;
                 Acons3;
                 Acons4];
             Bcons=[Bcons1;
                 Bcons2;
                 Bcons3;
                 Bcons4];
        end
        function [Acons Bcons] = cons_stretching(pankinit,pankfin,backtoankle,fronttoankle,exttoankle,inttoankle,xankmax,xankmin,yankmax,yankmin,theta,nbparamABCD,nbparamank,nbparamBb)
            nbparamtotal=nbparamABCD+4+nbparamank+nbparamBb;
            xApankle=[zeros(1,nbparamtotal*2);
                zeros(nbparamank,nbparamABCD+4) eye(nbparamank) zeros(nbparamank,nbparamBb+nbparamtotal);
                zeros(1,nbparamtotal*2)];
            xBpankle=[pankinit(1);zeros(nbparamank,1);pankfin(1)];
            yApankle=[zeros(1,nbparamtotal*2);
                zeros(nbparamank,nbparamtotal+nbparamABCD+4) eye(nbparamank) zeros(nbparamank,nbparamBb);
                zeros(1,nbparamtotal*2)];
            yBpankle=[pankinit(2);zeros(nbparamank,1);pankfin(2)];

            %clock-wise turn from upper right if foot in x direction
            ABCDleft=[fronttoankle -backtoankle -backtoankle fronttoankle;
                       exttoankle exttoankle -inttoankle -inttoankle];

            ABCDright=[fronttoankle -backtoankle -backtoankle fronttoankle;
                        inttoankle inttoankle -exttoankle -exttoankle];

            Acons=zeros((nbparamank+1)*16,nbparamtotal*2);
            Bcons=zeros((nbparamank+1)*16,1);

            for i=1:nbparamank+1
                t1=theta(i*3-1);
                t2=theta(i*3+1);

                rot=[cos(t2) -sin(t2);sin(t2) cos(t2)];
                dx=[cos(t1) -sin(t1) -cos(t1) sin(t1)];
                dy=[sin(t1) cos(t1) -sin(t1) -cos(t1)];

                if(mod(i,2))
                    ABCD_=ABCDright;
                else
                    ABCD_=ABCDleft;
                end
                [Acons((i-1)*16+1:(i-1)*16+4,:)       Bcons((i-1)*16+1:(i-1)*16+4)]         = wpg_qp_problem_Adrien.cons_builder (rot,dx,dy,ABCD_(:,1),xApankle(i,:),xBpankle(i),yApankle(i,:),yBpankle(i),xApankle(i+1,:),xBpankle(i+1),yApankle(i+1,:),yBpankle(i+1),xankmax,xankmin,yankmax,yankmin,i);
                [Acons((i-1)*16+1+4:(i-1)*16+4+4,:)   Bcons((i-1)*16+1+4:(i-1)*16+4+4)]     = wpg_qp_problem_Adrien.cons_builder (rot,dx,dy,ABCD_(:,2),xApankle(i,:),xBpankle(i),yApankle(i,:),yBpankle(i),xApankle(i+1,:),xBpankle(i+1),yApankle(i+1,:),yBpankle(i+1),xankmax,xankmin,yankmax,yankmin,i);
                [Acons((i-1)*16+1+8:(i-1)*16+4+8,:)   Bcons((i-1)*16+1+8:(i-1)*16+4+8)]     = wpg_qp_problem_Adrien.cons_builder (rot,dx,dy,ABCD_(:,3),xApankle(i,:),xBpankle(i),yApankle(i,:),yBpankle(i),xApankle(i+1,:),xBpankle(i+1),yApankle(i+1,:),yBpankle(i+1),xankmax,xankmin,yankmax,yankmin,i);
                [Acons((i-1)*16+1+12:(i-1)*16+4+12,:) Bcons((i-1)*16+1+12:(i-1)*16+4+12)]   = wpg_qp_problem_Adrien.cons_builder (rot,dx,dy,ABCD_(:,4),xApankle(i,:),xBpankle(i),yApankle(i,:),yBpankle(i),xApankle(i+1,:),xBpankle(i+1),yApankle(i+1,:),yBpankle(i+1),xankmax,xankmin,yankmax,yankmin,i);
            end
        end
        function [Acons Bcons] = cons_builder (rot,dx,dy,ABCD,xApankle1,xBpankle1,yApankle1,yBpankle1,xApankle2,xBpankle2,yApankle2,yBpankle2,xankmax,xankmin,yankmax,yankmin,j)
            ABCD_=rot*ABCD;
            xApankle2_=xApankle2;
            xBpankle2_=xBpankle2+ABCD_(1);
            yApankle2_=yApankle2;
            yBpankle2_=yBpankle2+ABCD_(2);

            Acons1=dx(1)*(xApankle2_-xApankle1)+dy(1)*(yApankle2_-yApankle1);
            Bcons1=-dx(1)*(xBpankle2_-xBpankle1)-dy(1)*(yBpankle2_-yBpankle1)+xankmax;
            Acons2=dx(2)*(xApankle2_-xApankle1)+dy(2)*(yApankle2_-yApankle1);
            Bcons2=-dx(2)*(xBpankle2_-xBpankle1)-dy(2)*(yBpankle2_-yBpankle1);
            Acons3=dx(3)*(xApankle2_-xApankle1)+dy(3)*(yApankle2_-yApankle1);
            Bcons3=-dx(3)*(xBpankle2_-xBpankle1)-dy(3)*(yBpankle2_-yBpankle1)-xankmin;
            Acons4=dx(4)*(xApankle2_-xApankle1)+dy(4)*(yApankle2_-yApankle1);
            Bcons4=-dx(4)*(xBpankle2_-xBpankle1)-dy(4)*(yBpankle2_-yBpankle1);

            Bcons2=Bcons2+yankmax*mod(j-1,2)-yankmin*mod(j,2);
            Bcons4=Bcons4+yankmax*mod(j,2)-yankmin*mod(j-1,2);

            Acons=[Acons1;Acons2;Acons3;Acons4];
            Bcons=[Bcons1;Bcons2;Bcons3;Bcons4];
        end
        %% intern methods of compute_Aeq_beq
        function [A B] = cons_ankle_fixed_path(nbparamABCD,nbparamank,nbparamBb,step_number_pankle_fixed)
            Ascomeq_path=[];
            Bscomeq_path=[];
            if size(step_number_pankle_fixed,1)~=0
                for i=1:size(step_number_pankle_fixed)
                    [Ascomeq_path Bscomeq_path]=wpg_qp_problem_Adrien.cons_pankle_fixed(Ascomeq_path,Bscomeq_path,nbparamABCD,nbparamank,nbparamBb,step_number_pankle_fixed(i,1),step_number_pankle_fixed(i,2:3));
                end
            end
            A=Ascomeq_path;
            B=Bscomeq_path;
        end
        function [A B] = cons_pankle_fixed(AscomeqDSP,Bscomeq,nbparamABCD,nbparamank,nbparamBb,step_number,pankle_fixed)
            nbparamtotal=nbparamABCD+nbparamank+nbparamBb;

            Ascomeq_path=[AscomeqDSP;zeros(1,nbparamABCD+4+step_number-1) 1 zeros(1,nbparamtotal+4+nbparamBb+(nbparamank-step_number))];
            Bscomeq_path=[Bscomeq;-pankle_fixed(1)];

            Ascomeq_path=[Ascomeq_path;zeros(1,nbparamtotal+4+nbparamABCD+4+step_number-1) 1 zeros(1,nbparamBb+(nbparamank-step_number))];
            Bscomeq_path=[Bscomeq_path;-pankle_fixed(2)];

            A=Ascomeq_path;
            B=Bscomeq_path;
        end
        %% %intern methods of qp_generating_problem_matrix
        function [A B] = cost_viapoint_gradient(Afcom2,Bfcom2,Atorque2ssp,Btorque2ssp,lambda,mu)
            A=lambda*[Afcom2 zeros(size(Afcom2,1),size(Atorque2ssp,2)-size(Afcom2,2));zeros(size(Atorque2ssp,1)-size(Afcom2,1),size(Afcom2,2)) zeros(size(Atorque2ssp,1)-size(Afcom2,1),size(Atorque2ssp,2)-size(Afcom2,2))]+mu*Atorque2ssp;
            B=lambda*[Bfcom2 zeros(1,size(Btorque2ssp,2)-size(Bfcom2,2))] +mu*Btorque2ssp;
        end
        function [A B] = cost_viapoint_gradient_nl(Afcom2,Bfcom2,Atorque2ssp,Btorque2ssp,lambda,mu)
            A=lambda*[Afcom2 zeros(size(Afcom2,1),size(Atorque2ssp,2)-size(Afcom2,2));zeros(size(Atorque2ssp,1)-size(Afcom2,1),size(Afcom2,2)) zeros(size(Atorque2ssp,1)-size(Afcom2,1),size(Atorque2ssp,2)-size(Afcom2,2))];
            B=lambda*[Bfcom2 zeros(1,size(Btorque2ssp,2)-size(Bfcom2,2))];
        end
        function [A B] = cost_viapoint_gradient_withDSP_azmp1(Aviapoint,Bviapoint,Atorque2_1,Btorque2_1,Atorque2_2,Btorque2_2,Aa2zmp1,Ba2zmp1,lambda,mu)
            A=[Aviapoint zeros(size(Aviapoint,1),size(Atorque2_1,2)-size(Aviapoint,2));zeros(size(Atorque2_1,1)-size(Aviapoint,1),size(Atorque2_1,2))]+mu*Atorque2_1+mu*Atorque2_2+lambda*Aa2zmp1;
            B=[Bviapoint zeros(1,size(Btorque2_1,2)-size(Bviapoint,2))]+mu*Btorque2_1+mu*Btorque2_2+lambda*Ba2zmp1;
            %A=[Aviapoint zeros(size(Aviapoint,1),size(Atorque2_1,2)-size(Aviapoint,2));zeros(size(Atorque2_1,1)-size(Aviapoint,1),size(Atorque2_1,2))]+mu*Atorque2_1+mu*Atorque2_2;
            %B=[Bviapoint zeros(1,size(Btorque2_1,2)-size(Bviapoint,2))]+mu*Btorque2_1+mu*Btorque2_2;            
        end
        function [A B] = cost_viapoint_gradient_withDSP_nl(Aviapoint,Bviapoint,Atorque2_1,Btorque2_1,Atorque2_2,Btorque2_2,Aa2zmp1,Ba2zmp1,lambda,mu)
            A=[Aviapoint zeros(size(Aviapoint,1),size(Atorque2_1,2)-size(Aviapoint,2));zeros(size(Atorque2_1,1)-size(Aviapoint,1),size(Atorque2_1,2))];
            B=[Bviapoint zeros(1,size(Btorque2_1,2)-size(Bviapoint,2))];
            %A=[Aviapoint zeros(size(Aviapoint,1),size(Atorque2_1,2)-size(Aviapoint,2));zeros(size(Atorque2_1,1)-size(Aviapoint,1),size(Atorque2_1,2))]+mu*Atorque2_1+mu*Atorque2_2;
            %B=[Bviapoint zeros(1,size(Btorque2_1,2)-size(Bviapoint,2))]+mu*Btorque2_1+mu*Btorque2_2;            
        end        
        function [A B] = pankle_fixed_path(AscomeqDSP,Bscomeq,nbparamABCD,nbparamank,nbparamBb,step_number_pankle_fixed)
            AscomeqDSP_path=AscomeqDSP;
            Bscomeq_path=Bscomeq;
            if size(step_number_pankle_fixed,1)~=0
                for i=1:size(step_number_pankle_fixed)
                    [AscomeqDSP_path Bscomeq_path]=wpg_qp_problem_Adrien.pankle_fixed(AscomeqDSP_path,Bscomeq_path,nbparamABCD,nbparamank,nbparamBb,step_number_pankle_fixed(i,1),step_number_pankle_fixed(i,2:3));
                end
            end
            A=AscomeqDSP_path;
            B=Bscomeq_path;
        end
        function [A B] = pankle_fixed(AscomeqDSP,Bscomeq,nbparamABCD,nbparamank,nbparamBb,step_number,pankle_fixed)
            nbparamtotal=nbparamABCD+nbparamank+nbparamBb;

            Ascomeq_path=[AscomeqDSP;zeros(1,nbparamABCD+4+step_number-1) 1 zeros(1,nbparamtotal+4+nbparamBb+(nbparamank-step_number))];
            Bscomeq_path=[Bscomeq;-pankle_fixed(1)];

            Ascomeq_path=[Ascomeq_path;zeros(1,nbparamtotal+4+nbparamABCD+4+step_number-1) 1 zeros(1,nbparamBb+(nbparamank-step_number))];
            Bscomeq_path=[Bscomeq_path;-pankle_fixed(2)];

            A=Ascomeq_path;
            B=Bscomeq_path;
        end
        %% %intern methods of qp_generating_DSP
        function [mdt] = compute_coeff_dt1_matrix(discretization,frequency,cutting)
            nbphases=length(discretization);
            nbpoints=sum(discretization)+1;
            mdt=zeros(nbpoints,nbphases*6+6+6);
            % mdt=[1 0 0 0 0 0 zeros(1,(nbphases-1)*6)];
            mdt(1,1:6)=[1 0 0 0 0 0];

            %%%zmp1 init cut in two%%%
            nb=discretization(1);
            rowinitphase=1;
            for i=1:ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[dt^0 dt^1 dt^2 dt^3 dt^4 dt^5];
                mdt(rowinitphase+i,(1-1)*6+1:(1-1)*6+6)=Dt;
            end
            rowinitphase=1+ceil(nb*cutting);
            for i=1:nb-ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[dt^0 dt^1 dt^2 dt^3 dt^4 dt^5];
                mdt(rowinitphase+i,(2-1)*6+1:(2-1)*6+6)=Dt;
            end

            rowinitphase=1;
            %%%generation of dt%%%
            for j=2:nbphases-1
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency;%dt is delta(tj)=t-Tj
                    Dt=[dt^0 dt^1 dt^2 dt^3 dt^4 dt^5];
                        mdt(rowinitphase+i,(j+1-1)*6+1:(j+1-1)*6+6)=Dt;
                end
            end

            %%%zmp1 fin cut in two%%%
            rowinitphase=rowinitphase+nb;
            nb=discretization(end);
            for i=1:ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[dt^0 dt^1 dt^2 dt^3 dt^4 dt^5];
                mdt(rowinitphase+i,(nbphases+1-1)*6+1:(nbphases+1-1)*6+6)=Dt;
            end
            rowinitphase=rowinitphase+ceil(nb*cutting);
            for i=1:nb-ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[dt^0 dt^1 dt^2 dt^3 dt^4 dt^5];
                mdt(rowinitphase+i,(nbphases+2-1)*6+1:(nbphases+2-1)*6+6)=Dt;
            end

        end
        function [A B]=compute_coeff_zmp1_gradient(tpassage,psa_zmp1init,psa_zmp1fin,cutting)
            nbpoly=length(tpassage)-1;
            nbpoly1=(length(tpassage)-1-2)/3+2;
            g=[];

            m22=zeros(6,nbpoly1-3);

            %zmp1 init
            dt=tpassage(2)-tpassage(1);
            dt=dt*cutting;
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            m=inv(m);
            m0=zeros(6,0);
            m1=m(1:6,1:3);
            m2=zeros(6,(nbpoly)*3);

            m3=zeros(6,0);
            m4=m(1:6,4:6);
            m5=zeros(6,3*(nbpoly1-1)+3+3);

            g=[g;m0 m1 m2 m22 m3 m4 m5];

            %decoupage de zmp1 init en deux parties
            dt=tpassage(2)-tpassage(1);
            dt=dt*(1-cutting);
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            m=inv(m);
            m0=zeros(6,0);
            m1=zeros(6,0);
            m2=zeros(6,(nbpoly)*3+3);

            m3=zeros(6,0);
            m4=m;
            m5=zeros(6,3*(nbpoly1-1)+3);

            g=[g;m0 m1 m2 m22 m3 m4 m5];

            for i=2:nbpoly-1
                dt=tpassage(i+1)-tpassage(i);
                    m=[1 0 0 0 0 0;
                       0 1 0 0 0 0;
                       0 0 2 0 0 0;
                       1 dt dt^2 dt^3 dt^4 dt^5;
                       0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                       0 0 2 6*dt 12*dt^2 20*dt^3];
                m=inv(m);

                m0=zeros(6,3*(i-1));
                m1=m(1:6,1:3);
                m2=zeros(6,(nbpoly)*3-(i-1)*3);

                m3=zeros(6,(i+1-mod(i+1,3))+3);
                m4=m(1:6,4:6);
                m5=zeros(6,(nbpoly1-1)*3-(i+1-mod(i+1,3))+3);

                g=[g;m0 m1 m2 m22 m3 m4 m5];
            end

            %zmp1 fin
            dt=tpassage(end)-tpassage(end-1);
            dt=dt*cutting;
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            m=inv(m);

            m0=zeros(6,(nbpoly)*3-3);
            m1=m(1:6,1:3);
            m2=zeros(6,3);

            m3=zeros(6,3*(nbpoly1-1)+3);
            m4=m(1:6,4:6);
            m5=zeros(6,3);

            g=[g;m0 m1 m2 m22 m3 m4 m5];

            %decoupage de zmp1 fin en deux parties
            dt=tpassage(end)-tpassage(end-1);
            dt=dt*(1-cutting);
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            m=inv(m);

            m0=zeros(6,(nbpoly)*3-3);
            m1=zeros(6,0);
            m2=zeros(6,3+3);

            m3=zeros(6,3*(nbpoly1-1)+3);
            m4=m;
            m5=zeros(6,0);

            g=[g;m0 m1 m2 m22 m3 m4 m5];

            A=[g(:,4:(nbpoly)*3+(nbpoly1-3)) zeros(size(g,1),4) g(:,(nbpoly)*3+(nbpoly1-3)+1+3:end-3)];
            B=[g(:,1:3)]*psa_zmp1init+[g(:,end-2:end)]*psa_zmp1fin;
        end
        function [A]=compute_repartition_matrix(tpassage,cutting)
            nbpoly=length(tpassage)-1;
            g=[];

            %%%cutting fcom1 init in two%%%
            dt=tpassage(2)-tpassage(1);
            dt=dt*cutting;
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            %m=inv(m);
            m=m\[0.5;0;0;0.75;0;0];
            g=[g;m];

            %second part
            dt=tpassage(2)-tpassage(1);
            dt=dt*(1-cutting);
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            % m=inv(m);
            m=m\[0.75;0;0;0;0;0];
            g=[g;m];


            for i=2:nbpoly-1
                dt=tpassage(i+1)-tpassage(i);
                m=[1 0 0 0 0 0;
                   0 1 0 0 0 0;
                   0 0 2 0 0 0;
                   1 dt dt^2 dt^3 dt^4 dt^5;
                   0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                   0 0 2 6*dt 12*dt^2 20*dt^3];
                %m=inv(m);
                m=m\[1;0;0;0;0;0];

                g=[g;m];
            end

            %%%cutting fcom1 fin in two%%%
            dt=tpassage(end)-tpassage(end-1);
            dt=dt*(1-cutting);
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            % m=inv(m);
            m=m\[1;0;0;0.25;0;0];
            g=[g;m];

            %second part
            dt=tpassage(end)-tpassage(end-1);
            dt=dt*cutting;
            m=[1 0 0 0 0 0;
                0 1 0 0 0 0;
                0 0 2 0 0 0;
                1 dt dt^2 dt^3 dt^4 dt^5;
                0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                0 0 2 6*dt 12*dt^2 20*dt^3];
            % m=inv(m);
            m=m\[0.25;0;0;0.5;0;0];
            g=[g;m];

            A=g;
            % A=g(:,4:size(g,2)-3);
            % B=[g(:,1:3) g(:,size(g,2)-2:size(g,2))]*[psa_zmpinit;psa_zmpfin];

        end
        function [Apankle Bpankle] = torque_ankle_positions_DSP_1(pankinit1,pankinit2,pankfin1,discretization)
            nbphases=length(discretization);

            Apankle=zeros(1+sum(discretization),ceil(nbphases/3));
            nb=0;
            rowinitphase=1;
            for j=1:nbphases-2
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                if(mod(j,3)==0&&j>3)
                    for i=1:nb
                         Apankle(rowinitphase+i,ceil(j/3))=1;
                    end
                else
            %         for i=1:nb
            %         pankle(rowinitphase+i)=pstep(ceil(j/3));
            %         end
                end
            end
            Apankle=Apankle(:,any(Apankle,1));

            Bpankle=zeros(1+sum(discretization),1);
            Bpankle(1)=pankinit1;
            Bpankle(1+1:1+discretization(1))=pankinit1*ones(discretization(1),1);
            Bpankle(1+sum(discretization(1:2))+1:1+sum(discretization(1:3)))=pankinit2*ones(discretization(3),1);
            Bpankle(1+sum(discretization(1:end-1))+1:1+sum(discretization(1:end)))=pankfin1*ones(discretization(end),1);

        end
        function [A B] = torque_ankle_DSP_gradient_1(Apzmp1,Bpzmp1,Afcom,Bfcom,Apankle,Bpankle,mg,ha,Afb)
        %le couple en x est fonction des y
        %le couple en y est fonction des x
            A1=ha*Afb*([Afcom zeros(size(Afcom,1),size(Apzmp1,2)-size(Afcom,2))]);
            A2=+Afb*mg*Apzmp1;
            A3=[zeros(size(Afcom)) -Afb*mg*Apankle zeros(size(Afcom,1),size(Apzmp1,2)-size(Afcom,2)-size(Apankle,2))];
            A=A1+A2+A3;
            B=ha*Afb*(Bfcom)+Afb*mg*Bpzmp1-Afb*mg*Bpankle;

            A=A(any(any(Apankle,2)+any(Bpankle,2),2),:);
            B=B(any(any(Apankle,2)+any(Bpankle,2),2),:);
        end
        function [A B]=compute_pzmp2(Apzmp,Bpzmp,Apzmp1,Bpzmp1,Afb,mg)
            xf2=(eye(length(Afb))-Afb)*mg;
            A=-mg*(xf2\(Apzmp1-[Apzmp zeros(size(Apzmp,1),size(Apzmp1,2)-size(Apzmp,2))]))+Apzmp1;
            B=-mg*(xf2\(Bpzmp1-Bpzmp))+Bpzmp1;
        end
        function [Apankle Bpankle] = torque_ankle_positions_DSP_2(pankinit2,pankfin1,pankfin2,discretization)
            nbphases=length(discretization);
            Apankle=zeros(1+sum(discretization),ceil(nbphases/3));
            nb=0;
            rowinitphase=1;
            for j=1:nbphases-3
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                if(mod(j,3)==0&&j>=3)
                    for i=1:nb
                         Apankle(rowinitphase+i,ceil(j/3))=1;
                    end
                else
            %         for i=1:nb
            %         pankle(rowinitphase+i)=pstep(ceil(j/3));
            %         end
                end
            end
            Apankle=Apankle(:,any(Apankle,1));

            Bpankle=zeros(1+sum(discretization),1);
            Bpankle(1)=pankinit2;
            Bpankle(1+1:1+discretization(1))=pankinit2*ones(discretization(1),1);
            Bpankle(1+sum(discretization(1:end-3))+1:1+sum(discretization(1:end-2)))=pankfin1*ones(discretization(end-2),1);
            Bpankle(1+sum(discretization(1:end-1))+1:1+sum(discretization(1:end)))=pankfin2*ones(discretization(end),1);
        end
        function [A B] = torque_ankle_DSP_gradient_2(Apzmp2,Bpzmp2,Afcom,Bfcom,Apankle,Bpankle,mg,ha,Afb)
        %le couple en x est fonction des y
        %le couple en y est fonction des x
            IAfb=eye(size(Afb,1))-Afb;

            A1=ha*IAfb*([Afcom zeros(size(Afcom,1),size(Apzmp2,2)-size(Afcom,2))]);
            A2=+IAfb*mg*Apzmp2;
            A3=[zeros(size(Afcom)) -IAfb*mg*Apankle zeros(size(Afcom,1),size(Apzmp2,2)-size(Afcom,2)-size(Apankle,2))];
            A=A1+A2+A3;

            B=ha*IAfb*(Bfcom)+IAfb*mg*Bpzmp2-IAfb*mg*Bpankle;
        end
        function [mdt] = compute_coeff_dddt1_matrix(discretization,frequency,cutting)
            nbphases=length(discretization);
            nbpoints=sum(discretization)+1;
            mdt=zeros(nbpoints,nbphases*6+6+6);

            mdt(1,1:6)=[0 0 2 0 0 0];

            %%%zmp1 init cut in two%%%
            nb=discretization(1);
            rowinitphase=1;
            for i=1:ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 0 2 6*dt^1 12*dt^2 20*dt^3];
                mdt(rowinitphase+i,(1-1)*6+1:(1-1)*6+6)=Dt;
            end
            rowinitphase=1+ceil(nb*cutting);
            for i=1:nb-ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 0 2 6*dt^1 12*dt^2 20*dt^3];
                mdt(rowinitphase+i,(2-1)*6+1:(2-1)*6+6)=Dt;
            end

            rowinitphase=1;
            %%%generation of dt%%%
            for j=2:nbphases-1
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency;%dt is delta(tj)=t-Tj
                    Dt=[0 0 2 6*dt^1 12*dt^2 20*dt^3];
                        mdt(rowinitphase+i,(j+1-1)*6+1:(j+1-1)*6+6)=Dt;
                end
            end

            %%%zmp1 fin cut in two%%%
            rowinitphase=rowinitphase+nb;
            nb=discretization(end);
            for i=1:ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 0 2 6*dt^1 12*dt^2 20*dt^3];
                mdt(rowinitphase+i,(nbphases+1-1)*6+1:(nbphases+1-1)*6+6)=Dt;
            end
            rowinitphase=rowinitphase+ceil(nb*cutting);
            for i=1:nb-ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 0 2 6*dt^1 12*dt^2 20*dt^3];
                mdt(rowinitphase+i,(nbphases+2-1)*6+1:(nbphases+2-1)*6+6)=Dt;
            end
        end
        function [mdt] = compute_coeff_ddt1_matrix(discretization,frequency,cutting)
            nbphases=length(discretization);
            nbpoints=sum(discretization)+1;
            mdt=zeros(nbpoints,nbphases*6+6+6);

            mdt(1,1:6)=[0 1 0 0 0 0];

            %%%zmp1 init cut in two%%%
            nb=discretization(1);
            rowinitphase=1;
            for i=1:ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 1 2*dt^1 3*dt^2 4*dt^3 5*dt^4];
                mdt(rowinitphase+i,(1-1)*6+1:(1-1)*6+6)=Dt;
            end
            rowinitphase=1+ceil(nb*cutting);
            for i=1:nb-ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 1 2*dt^1 3*dt^2 4*dt^3 5*dt^4];
                mdt(rowinitphase+i,(2-1)*6+1:(2-1)*6+6)=Dt;
            end

            rowinitphase=1;
            %%%generation of dt%%%
            for j=2:nbphases-1
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency;%dt is delta(tj)=t-Tj
                    Dt=[0 1 2*dt^1 3*dt^2 4*dt^3 5*dt^4];
                        mdt(rowinitphase+i,(j+1-1)*6+1:(j+1-1)*6+6)=Dt;
                end
            end

            %%%zmp1 fin cut in two%%%
            rowinitphase=rowinitphase+nb;
            nb=discretization(end);
            for i=1:ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 1 2*dt^1 3*dt^2 4*dt^3 5*dt^4];
                mdt(rowinitphase+i,(nbphases+1-1)*6+1:(nbphases+1-1)*6+6)=Dt;
            end
            rowinitphase=rowinitphase+ceil(nb*cutting);
            for i=1:nb-ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 1 2*dt^1 3*dt^2 4*dt^3 5*dt^4];
                mdt(rowinitphase+i,(nbphases+2-1)*6+1:(nbphases+2-1)*6+6)=Dt;
            end
        end
        function [mdt] = compute_coeff_ddt_matrix(discretization,frequency)
            nbphases=length(discretization);
            nbpoints=sum(discretization)+1;
            mdt=zeros(nbpoints,nbphases*6);
            % mdt=[1 0 0 0 0 0 zeros(1,(nbphases-1)*6)];
            mdt(1,1:6)=[0 1 0 0 0 0];
            nb=0;
            rowinitphase=1;
            for j=1:nbphases
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency;%dt is delta(tj)=t-Tj
                    Dt=[0 dt^0 2*dt^1 3*dt^2 4*dt^3 5*dt^4];
            %         mdt=[mdt;zeros(1,(j-1)*6) Dt zeros(1,(nbphases-j)*6)];
                     mdt(rowinitphase+i,(j-1)*6+1:(j-1)*6+6)=Dt;
                end
            end
        end
        function [mdt] = compute_coeff_dddt_matrix(discretization,frequency)
            nbphases=length(discretization);
            nbpoints=sum(discretization)+1;
            mdt=zeros(nbpoints,nbphases*6);
            % mdt=[1 0 0 0 0 0 zeros(1,(nbphases-1)*6)];
            mdt(1,1:6)=[0 0 2 0 0 0];
            nb=0;
            rowinitphase=1;
            for j=1:nbphases
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency;%dt is delta(tj)=t-Tj
                    Dt=[0 0 2*dt^0 6*dt^1 12*dt^2 20*dt^3];
            %         mdt=[mdt;zeros(1,(j-1)*6) Dt zeros(1,(nbphases-j)*6)];
                        mdt(rowinitphase+i,(j-1)*6+1:(j-1)*6+6)=Dt;
                end
            end
        end
        function [A B]=compute_azmp2_gradient(Apzmp,Bpzmp,Apzmp1,Bpzmp1,Aszmp,Bszmp,Aszmp1,Bszmp1,Aazmp,Bazmp,Aazmp1,Bazmp1,Afb,dAfb,ddAfb)
        %Compute the gradient matrix of zmp2 acceleration
            Afb_=(eye(length(Afb))-Afb);
            dAfb_=-dAfb;
            ddAfb_=-ddAfb;

            A1=-(Afb_\(Aazmp1-[Aazmp zeros(size(Aazmp,1),size(Aazmp1,2)-size(Aazmp,2))]));
            A2=+2*dAfb_*((Afb_)\((Afb_)\(Aszmp1-[Aszmp zeros(size(Aszmp,1),size(Aszmp1,2)-size(Aszmp,2))])));
            A3=+ddAfb_*((Afb_)\((Afb_)\(Apzmp1-[Apzmp zeros(size(Apzmp,1),size(Apzmp1,2)-size(Apzmp,2))])))-2*dAfb_*dAfb_*((Afb_)\((Afb_)\((Afb_)\(Apzmp1-[Apzmp zeros(size(Apzmp,1),size(Apzmp1,2)-size(Apzmp,2))]))));

            B1=-(Afb_\(Bazmp1-Bazmp));
            B2=+2*dAfb_*((Afb_)\((Afb_)\(Bszmp1-Bszmp)));
            B3=+ddAfb_*((Afb_)\((Afb_)\(Bpzmp1-Bpzmp)))-2*dAfb_*dAfb_*((Afb_)\((Afb_)\((Afb_)\(Bpzmp1-Bpzmp))));

            A=A1+A2+A3+Aazmp1;
            B=B1+B2+B3+Bazmp1;
        end
        function [Acons Bcons] = zmp_constraint_stability_DSP_xy_1(xApzmp,xBpzmp,yApzmp,yBpzmp,xApankle,xBpankle,yApankle,yBpankle,discretization,backtoankle,fronttoankle,exttoankle,inttoankle,sole_margin,theta,type_phase,nbparamABCD,nbparamBb,firstSS)
            xApzmp_=[xApzmp];
            yApzmp_=[yApzmp];
            xApankle_=[zeros(size(xApzmp_,1),nbparamABCD+4) xApankle zeros(size(xApzmp,1),nbparamBb) zeros(size(yApzmp_))];
            yApankle_=[zeros(size(xApzmp_)) zeros(size(yApzmp_,1),nbparamABCD+4) yApankle zeros(size(yApzmp,1),nbparamBb)];

            xApzmp_=[xApzmp_ zeros(size(xApzmp_))];
            yApzmp_=[zeros(size(yApzmp_)) yApzmp_];

            Acons1=zeros(size(xApankle_));
            Bcons1=zeros(size(xBpankle));
            Acons2=zeros(size(xApankle_));
            Bcons2=zeros(size(xBpankle));
            Acons3=zeros(size(xApankle_));
            Bcons3=zeros(size(xBpankle));
            Acons4=zeros(size(xApankle_));
            Bcons4=zeros(size(xBpankle));

            j=1;
            %constraint direction inverse clock-wise
            theta_=theta(1:end-1);
            theta_=theta_(any(type_phase==0,1));

            t=theta_(j);
            dx=[cos(t) -sin(t) -cos(t) sin(t)];
            dy=[sin(t) cos(t) -sin(t) -cos(t)];

            xA=xApzmp_-xApankle_;
            yA=yApzmp_-yApankle_;
            xB=xBpzmp-xBpankle;
            yB=yBpzmp-yBpankle;

            discretization_=discretization(any(type_phase==0,1));

            for i=1:sum(discretization_)
                if(i>sum(discretization_(1:j)))
                        j=j+1;
                        %constraint direction inverse clock-wise
                        t=theta_(j);
                        dx=[cos(t) -sin(t) -cos(t) sin(t)];
                        dy=[sin(t) cos(t) -sin(t) -cos(t)];
                end

                    Acons1(i,:)=dx(1)*(xA(i,:))+dy(1)*(yA(i,:));
                    Bcons1(i)=-dx(1)*(xB(i,:))-dy(1)*(yB(i,:))+fronttoankle-sole_margin;
                    Acons2(i,:)=dx(2)*(xA(i,:))+dy(2)*(yA(i,:));
                    Bcons2(i)=-dx(2)*(xB(i,:))-dy(2)*(yB(i,:));
                    Acons3(i,:)=dx(3)*(xA(i,:))+dy(3)*(yA(i,:));
                    Bcons3(i)=-dx(3)*(xB(i,:))-dy(3)*(yB(i,:))+backtoankle-sole_margin;
                    Acons4(i,:)=dx(4)*(xA(i,:))+dy(4)*(yA(i,:));
                    Bcons4(i)=-dx(4)*(xB(i,:))-dy(4)*(yB(i,:));

                    if(firstSS==0)
                        Bcons2(i)=Bcons2(i)+(inttoankle*(mod(j,2)==1)+exttoankle*(mod(j,2)==0)-sole_margin);
                        Bcons4(i)=Bcons4(i)+(inttoankle*(mod(j,2)==0)+exttoankle*(mod(j,2)==1)-sole_margin);
                    elseif(firstSS==1)
                        Bcons2(i)=Bcons2(i)+(inttoankle*(mod(j,2)==0)+exttoankle*(mod(j,2)==1)-sole_margin);
                        Bcons4(i)=Bcons4(i)+(inttoankle*(mod(j,2)==1)+exttoankle*(mod(j,2)==0)-sole_margin);
                    else
                        'Choose a first SS foot'
                    end
            end

             Acons=[Acons1;
                 Acons2;
                 Acons3;
                 Acons4];
             Bcons=[Bcons1;
                 Bcons2;
                 Bcons3;
                 Bcons4];
        end
        function [Acons Bcons] = zmp_constraint_stability_DSP_xy_2(xApzmp,xBpzmp,yApzmp,yBpzmp,xApankle,xBpankle,yApankle,yBpankle,discretization,backtoankle,fronttoankle,exttoankle,inttoankle,sole_margin,theta,type_phase,nbparamABCD,nbparamBb,firstSS)
            xApzmp_=[xApzmp];
            yApzmp_=[yApzmp];
            xApankle=[zeros(size(xApzmp_,1),nbparamABCD+4) xApankle zeros(size(xApzmp,1),nbparamBb) zeros(size(yApzmp_))];
            yApankle=[zeros(size(xApzmp_)) zeros(size(yApzmp_,1),nbparamABCD+4) yApankle zeros(size(yApzmp,1),nbparamBb)];

            xApzmp_=[xApzmp_ zeros(size(xApzmp_))];
            yApzmp_=[zeros(size(yApzmp_)) yApzmp_];

            Acons1=zeros(size(xApankle));
            Bcons1=zeros(size(xBpankle));
            Acons2=zeros(size(xApankle));
            Bcons2=zeros(size(xBpankle));
            Acons3=zeros(size(xApankle));
            Bcons3=zeros(size(xBpankle));
            Acons4=zeros(size(xApankle));
            Bcons4=zeros(size(xBpankle));

            j=1;
            %constraint direction inverse clock-wise
            theta_=theta(2:end);
            theta_=theta_(any(type_phase==0,1));

            t=theta_(j);
            dx=[cos(t) -sin(t) -cos(t) sin(t)];
            dy=[sin(t) cos(t) -sin(t) -cos(t)];

            xA=xApzmp_-xApankle;
            yA=yApzmp_-yApankle;
            xB=xBpzmp-xBpankle;
            yB=yBpzmp-yBpankle;

            discretization_=discretization(any(type_phase==0,1));

            for i=1:sum(discretization_)
                if(i>sum(discretization_(1:j)))
                        j=j+1;
                        %constraint direction inverse clock-wise
                        t=theta_(j);
                        dx=[cos(t) -sin(t) -cos(t) sin(t)];
                        dy=[sin(t) cos(t) -sin(t) -cos(t)];
                end

                    Acons1(i,:)=dx(1)*(xA(i,:))+dy(1)*(yA(i,:));
                    Bcons1(i)=-dx(1)*(xB(i,:))-dy(1)*(yB(i,:))+fronttoankle-sole_margin;
                    Acons2(i,:)=dx(2)*(xA(i,:))+dy(2)*(yA(i,:));
                    Bcons2(i)=-dx(2)*(xB(i,:))-dy(2)*(yB(i,:));
                    Acons3(i,:)=dx(3)*(xA(i,:))+dy(3)*(yA(i,:));
                    Bcons3(i)=-dx(3)*(xB(i,:))-dy(3)*(yB(i,:))+backtoankle-sole_margin;
                    Acons4(i,:)=dx(4)*(xA(i,:))+dy(4)*(yA(i,:));
                    Bcons4(i)=-dx(4)*(xB(i,:))-dy(4)*(yB(i,:));

                    if(firstSS==0)
                        Bcons2(i)=Bcons2(i)+(inttoankle*(mod(j,2)==0)+exttoankle*(mod(j,2)==1)-sole_margin);
                        Bcons4(i)=Bcons4(i)+(inttoankle*(mod(j,2)==1)+exttoankle*(mod(j,2)==0)-sole_margin);
                    elseif(firstSS==1)
                        Bcons2(i)=Bcons2(i)+(inttoankle*(mod(j,2)==1)+exttoankle*(mod(j,2)==0)-sole_margin);
                        Bcons4(i)=Bcons4(i)+(inttoankle*(mod(j,2)==0)+exttoankle*(mod(j,2)==1)-sole_margin);
                    else
                        'Choose a first SS foot'
                    end
            end

             Acons=[Acons1;
                 Acons2;
                 Acons3;
                 Acons4];
             Bcons=[Bcons1;
                 Bcons2;
                 Bcons3;
                 Bcons4];
        end
        function [Acons Bcons] = zmp_constraint_ankle_pos(pankinit,pankfin,backtoankle,fronttoankle,exttoankle,inttoankle,xankmax,xankmin,yankmax,yankmin,theta,nbparamABCD,nbparamank,nbparamBb)
            nbparamtotal=nbparamABCD+4+nbparamank+nbparamBb;
            xApankle=[zeros(1,nbparamtotal*2);
                zeros(nbparamank,nbparamABCD+4) eye(nbparamank) zeros(nbparamank,nbparamBb+nbparamtotal);
                zeros(1,nbparamtotal*2)];
            xBpankle=[pankinit(1);zeros(nbparamank,1);pankfin(1)];
            yApankle=[zeros(1,nbparamtotal*2);
                zeros(nbparamank,nbparamtotal+nbparamABCD+4) eye(nbparamank) zeros(nbparamank,nbparamBb);
                zeros(1,nbparamtotal*2)];
            yBpankle=[pankinit(2);zeros(nbparamank,1);pankfin(2)];

            %clock-wise turn from upper right if foot in x direction
            ABCDleft=[fronttoankle -backtoankle -backtoankle fronttoankle;
                       exttoankle exttoankle -inttoankle -inttoankle];

            ABCDright=[fronttoankle -backtoankle -backtoankle fronttoankle;
                        inttoankle inttoankle -exttoankle -exttoankle];

            Acons=zeros((nbparamank+1)*16*4,nbparamtotal*2);
            Bcons=zeros((nbparamank+1)*16*4,1);

            for i=1:nbparamank+1
                t1=theta(i*3-1);
                t2=theta(i*3+1);

                rot=[cos(t2) -sin(t2);sin(t2) cos(t2)];
                dx=[cos(t1) -sin(t1) -cos(t1) sin(t1)];
                dy=[sin(t1) cos(t1) -sin(t1) -cos(t1)];

                if(mod(i,2))
                    ABCD_=ABCDright;
                else
                    ABCD_=ABCDleft;
                end
                [Acons((i-1)*16+1:(i-1)*16+4,:)       Bcons((i-1)*16+1:(i-1)*16+4)]         = wpg_qp_problem_Adrien.constraint_builder (rot,dx,dy,ABCD_(:,1),xApankle(i,:),xBpankle(i),yApankle(i,:),yBpankle(i),xApankle(i+1,:),xBpankle(i+1),yApankle(i+1,:),yBpankle(i+1),xankmax,xankmin,yankmax,yankmin,i);
                [Acons((i-1)*16+1+4:(i-1)*16+4+4,:)   Bcons((i-1)*16+1+4:(i-1)*16+4+4)]     = wpg_qp_problem_Adrien.constraint_builder (rot,dx,dy,ABCD_(:,2),xApankle(i,:),xBpankle(i),yApankle(i,:),yBpankle(i),xApankle(i+1,:),xBpankle(i+1),yApankle(i+1,:),yBpankle(i+1),xankmax,xankmin,yankmax,yankmin,i);
                [Acons((i-1)*16+1+8:(i-1)*16+4+8,:)   Bcons((i-1)*16+1+8:(i-1)*16+4+8)]     = wpg_qp_problem_Adrien.constraint_builder (rot,dx,dy,ABCD_(:,3),xApankle(i,:),xBpankle(i),yApankle(i,:),yBpankle(i),xApankle(i+1,:),xBpankle(i+1),yApankle(i+1,:),yBpankle(i+1),xankmax,xankmin,yankmax,yankmin,i);
                [Acons((i-1)*16+1+12:(i-1)*16+4+12,:) Bcons((i-1)*16+1+12:(i-1)*16+4+12)]   = wpg_qp_problem_Adrien.constraint_builder (rot,dx,dy,ABCD_(:,4),xApankle(i,:),xBpankle(i),yApankle(i,:),yBpankle(i),xApankle(i+1,:),xBpankle(i+1),yApankle(i+1,:),yBpankle(i+1),xankmax,xankmin,yankmax,yankmin,i);
            end
        end
        function [Acons Bcons]=constraint_builder (rot,dx,dy,ABCD,xApankle1,xBpankle1,yApankle1,yBpankle1,xApankle2,xBpankle2,yApankle2,yBpankle2,xankmax,xankmin,yankmax,yankmin,j)
            ABCD_=rot*ABCD;
            xApankle2_=xApankle2;
            xBpankle2_=xBpankle2+ABCD_(1);
            yApankle2_=yApankle2;
            yBpankle2_=yBpankle2+ABCD_(2);

            Acons1=dx(1)*(xApankle2_-xApankle1)+dy(1)*(yApankle2_-yApankle1);
            Bcons1=-dx(1)*(xBpankle2_-xBpankle1)-dy(1)*(yBpankle2_-yBpankle1)+xankmax;
            Acons2=dx(2)*(xApankle2_-xApankle1)+dy(2)*(yApankle2_-yApankle1);
            Bcons2=-dx(2)*(xBpankle2_-xBpankle1)-dy(2)*(yBpankle2_-yBpankle1);
            Acons3=dx(3)*(xApankle2_-xApankle1)+dy(3)*(yApankle2_-yApankle1);
            Bcons3=-dx(3)*(xBpankle2_-xBpankle1)-dy(3)*(yBpankle2_-yBpankle1)-xankmin;
            Acons4=dx(4)*(xApankle2_-xApankle1)+dy(4)*(yApankle2_-yApankle1);
            Bcons4=-dx(4)*(xBpankle2_-xBpankle1)-dy(4)*(yBpankle2_-yBpankle1);

            Bcons2=Bcons2+yankmax*mod(j-1,2)-yankmin*mod(j,2);
            Bcons4=Bcons4+yankmax*mod(j,2)-yankmin*mod(j-1,2);

            Acons=[Acons1;Acons2;Acons3;Acons4];
            Bcons=[Bcons1;Bcons2;Bcons3;Bcons4];
        end
        %% %intern methods of qp_generating_SSP
        function [mdt] = compute_coeff_dt_matrix(discretization,frequency,nbphases,nbpointdiscret)
            mdt=zeros(nbpointdiscret,nbphases*6);
            mdt(1,1:6)=[1 0 0 0 0 0];
            nb=0;
            rowinitphase=1;
            for j=1:nbphases
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency; %dt is delta(tj)=t-Tj
                    Dt=[dt^0 dt^1 dt^2 dt^3 dt^4 dt^5];
            %         mdt=[mdt;zeros(1,(j-1)*6) Dt zeros(1,(nbphases-j)*6)];
                    mdt(rowinitphase+i,(j-1)*6+1:(j-1)*6+6)=Dt;
                end
            end
        end
        function [mscwt] = compute_coeff_swt_cwt_matrix(discretization,w,frequency,nbphases)
            mscwt= zeros(1+sum(discretization),(nbphases)*2);
            mscwt(1,1:2)=[cosh(0) sinh(0)];
            nb=0;
            rowinitphase=1;
            for j=1:nbphases
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency;%dt is delta(tj)=t-Tj
                    mscwt(rowinitphase+i,(j-1)*2+1:(j-1)*2+2)=[cosh(w*dt) sinh(w*dt)];
                end
            end
        end
        function [A B]=compute_coeff_zmp_gradient(tpassage,psa_zmpinit,psa_zmpfin,nbphases,nbparamABCD)
        % Hypo : We consider ZMP trajectory as 5th order polynomials between via-points.
        % This function compute A and B as: C=A*x+B
        % Where : C are give the polynomials coefficient of ZMP trajectory in one direction.
        % x are the via-point boundary conditions in position, speed and acceleration.
        % Initial and final via-point boundary conditions are known and in psa_zmpinit and psa_zmpfin

            A=zeros(6*nbphases,nbparamABCD+6);
            for i=1:nbphases
                dt=tpassage(i+1)-tpassage(i);
                m=[1 0 0 0 0 0;
                   0 1 0 0 0 0;
                   0 0 2 0 0 0;
                   1 dt dt^2 dt^3 dt^4 dt^5;
                   0 1 2*dt 3*dt^2 4*dt^3 5*dt^4;
                   0 0 2 6*dt 12*dt^2 20*dt^3];
                m=inv(m);
                A(1+(i-1)*6:6+(i-1)*6,:)=[zeros(6,3*(i-1)) m zeros(6,(nbphases-i)*3)];
            end

            B=[A(:,1:3) A(:,end-2:end)]*[psa_zmpinit;psa_zmpfin];

            A(:,1:3)=[];
            A(:,end-2:end)=[];
            A=[A zeros(size(A,1),4)];
        end
        function [A]=com_morisawa_A_gradient(nbphases,w)
        %Compute the matrix gradient of A vector
            A=[];
            for j=1:nbphases
                Aj=[1 0 2/w^2 0 24/w^4 0;
                    0 1 0 6/w^2 0 120/w^4;
                    0 0 1 0 12/w^2 0;
                    0 0 0 1 0 20/w^2;
                    0 0 0 0 1 0;
                    0 0 0 0 0 1];
                A=[A;zeros(6,(j-1)*6) Aj zeros(6,(nbphases-j)*6)];
            end
        end
        function [A B]=com_morisawa_VW_gradient(A_zmp_gradient,B_zmp_gradient,A_gradient,tpassage,w)
            [A_global B_global]=wpg_qp_problem_Adrien.com_morisawa_VW_gradient_global(A_gradient,tpassage,w);

            A=[zeros(size(A_global,1),size(A_zmp_gradient,2)-4) B_global]+A_global*A_zmp_gradient;
            B=A_global*B_zmp_gradient;
        end
        function [A B]=com_morisawa_VW_gradient_global(A_gradient,tpassage,w)
        %put a system Z1*y=X*x+Z2*L as y=A*a+B*x
        %y is a vector of V and W scalar coeff of COM trajectory
        %x is a vector with the initial and final known COM position
        %a is a vector build of ZMP polynomial coeff
            Z1=wpg_qp_problem_Adrien.compute_Z1(tpassage,w);
            Z2=wpg_qp_problem_Adrien.compute_Z2(tpassage);

            X=[1 0 0 0; zeros(size(Z2,1)-2,4); 0 1 0 0]; %pcominit pconfin scominit scomfin

            A=Z1\(Z2*A_gradient); %%%%% Why A_gradient
            B=Z1\(X);
        end
        function [Z1]=compute_Z1(tpassage,w)
        %Compute the Z1 matrix
            nbphases=length(tpassage)-1;

            Z1=[1 0 zeros(1,2*nbphases-2)];

            for i=1:nbphases-1
                Dt=tpassage(i+1)-tpassage(i);
                Z1=[Z1;zeros(1,2*(i-1)) cosh(w*Dt) sinh(w*Dt) -1 0 zeros(1,2*(nbphases-i)-2);zeros(1,2*(i-1)) w*sinh(w*Dt) w*cosh(w*Dt) 0 -w zeros(1,2*(nbphases-i)-2)];
            end

            Dt=tpassage(nbphases+1)-tpassage(nbphases);
            Z1=[Z1;zeros(1,2*nbphases-2) cosh(w*Dt) sinh(w*Dt)];
        end
        function [Z2]=compute_Z2(tpassage)
        %compute the Z2 matrix
            nbphases=length(tpassage)-1;

            Z2=[-1 0 zeros(1,6*nbphases-2)];

            for i=1:nbphases-1
                Dt=tpassage(i+1)-tpassage(i);
                Z2=[Z2;zeros(1,6*(i-1)) -Dt^0 -Dt^1 -Dt^2 -Dt^3 -Dt^4 -Dt^5 +1 0 zeros(1,6*(nbphases-i)-2);zeros(1,6*(i-1)) 0 -1*Dt^0 -2*Dt^1 -3*Dt^2 -4*Dt^3 -5*Dt^4 0 +1 zeros(1,6*(nbphases-i)-2)];
            end

            Dt=tpassage(nbphases+1)-tpassage(nbphases);
            Z2=[Z2;zeros(1,6*(nbphases-1)) -Dt^0 -Dt^1 -Dt^2 -Dt^3 -Dt^4 -Dt^5];
        end
        function [A B] = pcom_generator_morisawa_gradient(A_zmp_gradient,B_zmp_gradient,AVW,BVW,A_gradient,mdt,mscwt)
        %%%generate the COM trajectory gradient%%%
        %Use Morisawa algorithm with 5th zmp polynomial coeff known
            A=mscwt*AVW+mdt*A_gradient*A_zmp_gradient;
            B=mscwt*BVW+mdt*A_gradient*B_zmp_gradient;
        end
        function [A B] = fcom_gradient(Apzmp,Bpzmp,Apcom,Bpcom,mg,h)
            A=(Apzmp-Apcom)*(-mg)/h;
            B=(Bpzmp-Bpcom)*(-mg)/h;
        end
        function [A B] = torque_ankle_SSP_gradient(Apzmp,Bpzmp,Afcom,Bfcom,pankinit,pankfin,discretization,mg,ha,T)
        %le couple en x est fonction des y
        %le couple en y est fonction des x
            [Apankle Bpankle]=wpg_qp_problem_Adrien.torque_ankle_positions_SSP(pankinit,pankfin,discretization);
            Apankle=Apankle(T,:);
            Bpankle=Bpankle(T,:);

            A1=ha*Afcom+mg*Apzmp;
            A2=-mg*Apankle;
            A=[A1 zeros(size(A2,1),size(A2,2))]+[zeros(size(A1,1),size(A1,2)) A2];
            B=ha*Bfcom+mg*Bpzmp-mg*Bpankle;

            A=A(any(any(Apankle,2)+any(Bpankle,2),2),:);
            B=B(any(any(Apankle,2)+any(Bpankle,2),2),:);

        end
        function [Apankle Bpankle] = torque_ankle_positions_SSP(pankinit,pankfin,discretization)
            nbphases=length(discretization);
            Apankle=zeros(1+sum(discretization),ceil(nbphases/3));
            %pankle(1)=0;
            nb=0;
            rowinitphase=1;
            for j=1:nbphases
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                if(mod(j,3)==0 || j==1 || j==nbphases || j==2 || j==nbphases-1)
            %         for i=1:nb
            %             pankle=[pankle;0];
            %         end
                else
                    for i=1:nb
                    Apankle(rowinitphase+i,ceil(j/3))=1;
                    end
                end
            end

            Apankle=Apankle(:,any(Apankle,1));

            Bpankle=zeros(1+sum(discretization),1);
            Bpankle(1+discretization(1)+1:1+sum(discretization(1:2)))=pankinit*ones(discretization(2),1);
            Bpankle(1+sum(discretization(1:end-2))+1:1+sum(discretization(1:end-1)))=pankfin*ones(discretization(end-1),1);

        end
        function [Acons Bcons] = zmp_constraint_stability_SSP_xy(xApzmp,xBpzmp,yApzmp,yBpzmp,xApankle,xBpankle,yApankle,yBpankle,discretization,backtoankle,fronttoankle,exttoankle,inttoankle,sole_margin,theta,type_phase,nbparamBb,firstSS)
            xApzmp_=[xApzmp zeros(size(xApankle)) zeros(size(xApzmp,1),nbparamBb)];
            yApzmp_=[yApzmp zeros(size(yApankle)) zeros(size(yApzmp,1),nbparamBb)];
            xApankle=[zeros(size(xApzmp)) xApankle zeros(size(xApzmp,1),nbparamBb) zeros(size(yApzmp_))];
            yApankle=[zeros(size(xApzmp_)) zeros(size(yApzmp)) yApankle zeros(size(yApzmp,1),nbparamBb)];

            xApzmp_=[xApzmp_ zeros(size(xApzmp_))];
            yApzmp_=[zeros(size(yApzmp_)) yApzmp_];

            Acons1=zeros(size(xApankle));
            Bcons1=zeros(size(xBpankle));
            Acons2=zeros(size(xApankle));
            Bcons2=zeros(size(xBpankle));
            Acons3=zeros(size(xApankle));
            Bcons3=zeros(size(xBpankle));
            Acons4=zeros(size(xApankle));
            Bcons4=zeros(size(xBpankle));

            j=1;
            %constraint direction inverse clock-wise
            t=theta(j);
            dx=[cos(t) -sin(t) -cos(t) sin(t)];
            dy=[sin(t) cos(t) -sin(t) -cos(t)];


            xA=xApzmp_-xApankle;
            yA=yApzmp_-yApankle;
            xB=xBpzmp-xBpankle;
            yB=yBpzmp-yBpankle;

            for i=1:sum(discretization)
                if(i>sum(discretization(1:j)))
                        j=j+1;
                        %constraint direction inverse clock-wise
                        t=theta(j);
                        dx=[cos(t) -sin(t) -cos(t) sin(t)];
                        dy=[sin(t) cos(t) -sin(t) -cos(t)];
                end

                if type_phase(j)~=0
                    Acons1(i,:)=dx(1)*(xA(i,:))+dy(1)*(yA(i,:));
                    Bcons1(i)=-dx(1)*(xB(i,:))-dy(1)*(yB(i,:))+fronttoankle-sole_margin;
                    Acons2(i,:)=dx(2)*(xA(i,:))+dy(2)*(yA(i,:));
                    Bcons2(i)=-dx(2)*(xB(i,:))-dy(2)*(yB(i,:));
                    Acons3(i,:)=dx(3)*(xA(i,:))+dy(3)*(yA(i,:));
                    Bcons3(i)=-dx(3)*(xB(i,:))-dy(3)*(yB(i,:))+backtoankle-sole_margin;
                    Acons4(i,:)=dx(4)*(xA(i,:))+dy(4)*(yA(i,:));
                    Bcons4(i)=-dx(4)*(xB(i,:))-dy(4)*(yB(i,:));
                    switch (j)
                        case 1 %pied droit
                            Bcons2(i)=Bcons2(i)+inttoankle-sole_margin;
                            Bcons4(i)=Bcons4(i)+exttoankle-sole_margin;
                        case 3 %pied gauche
                            Bcons2(i)=Bcons2(i)+exttoankle-sole_margin;
                            Bcons4(i)=Bcons4(i)+inttoankle-sole_margin;
%                         case 4
%                             Bcons2(i)=Bcons2(i)+exttoankle-sole_margin;
%                             Bcons4(i)=Bcons4(i)+inttoankle-sole_margin;
%                         case 6
%                             Bcons2(i)=Bcons2(i)+inttoankle-sole_margin;
%                             Bcons4(i)=Bcons4(i)+exttoankle-sole_margin;
                    end
                end
            end

            Acons=[Acons1(any(any(Acons1,2)+any(Bcons1,2),2),:);
                   Acons2(any(any(Acons2,2)+any(Bcons2,2),2),:);
                   Acons3(any(any(Acons3,2)+any(Bcons3,2),2),:);
                   Acons4(any(any(Acons4,2)+any(Bcons4,2),2),:)];
            Bcons=[Bcons1(any(any(Acons1,2)+any(Bcons1,2),2),:);
                   Bcons2(any(any(Acons2,2)+any(Bcons2,2),2),:);
                   Bcons3(any(any(Acons3,2)+any(Bcons3,2),2),:);
                   Bcons4(any(any(Acons4,2)+any(Bcons4,2),2),:)];
        end
        function [A B] = scom_eqcontraint(A_zmp_gradient,B_zmp_gradient,AVW,BVW,A_gradient,tpassage,nbphases,w)
            dti=0;
            dtf=tpassage(nbphases+1)-tpassage(nbphases);
            mscwt=[w*sinh(w*dti) w*cosh(w*dti) zeros(1,(nbphases-1)*2);
                zeros(1,(nbphases-1)*2) w*sinh(w*dtf) w*cosh(w*dtf)];
            mdt=[0 1 2*dti 3*dti^2 4*dti^3 5*dti^4 zeros(1,(nbphases-1)*6);
                zeros(1,(nbphases-1)*6) 0 1 2*dtf 3*dtf^2 4*dtf^3 5*dtf^4];

            A=mscwt*AVW+mdt*A_gradient*A_zmp_gradient-[zeros(2,size(A_zmp_gradient,2)-2) [1 0;0 1]];
            B=mscwt*BVW+mdt*A_gradient*B_zmp_gradient;
        end
        %% %intern methods global
        function [A B] = compute_traj_discrete(A_gradient,B_gradient,dt)
            A=dt*A_gradient;
            B=dt*B_gradient;
        end
        function [A B C] = compute_quad_matrix(A_gradient,B_gradient)
            A=transpose(A_gradient)*A_gradient;
            B=transpose(B_gradient)*A_gradient;
            C=transpose(B_gradient)*B_gradient;
        end
        function [M_t] = compute_M_t(discretization,frequency,nbphases,nbpointdiscret)
            M_t=zeros(nbpointdiscret,nbphases*6);
            M_t(1,1:6)=[1 0 0 0 0 0];
            nb=0;
            rowinitphase=1;
            for j=1:nbphases
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency; %dt is delta(tj)=t-Tj
                    Dt=[dt^0 dt^1 dt^2 dt^3 dt^4 dt^5];
            %         mdt=[mdt;zeros(1,(j-1)*6) Dt zeros(1,(nbphases-j)*6)];
                    M_t(rowinitphase+i,(j-1)*6+1:(j-1)*6+6)=Dt;
                end
            end
        end
        function [M_t_d] = compute_M_t_d(discretization,frequency)
            %compute the derivative of M_t
            nbphases=length(discretization);
            nbpoints=sum(discretization)+1;
            M_t_d=zeros(nbpoints,nbphases*6);
            % mdt=[1 0 0 0 0 0 zeros(1,(nbphases-1)*6)];
            M_t_d(1,1:6)=[0 1 0 0 0 0];
            nb=0;
            rowinitphase=1;
            for j=1:nbphases
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency;%dt is delta(tj)=t-Tj
                    Dt=[0 dt^0 2*dt^1 3*dt^2 4*dt^3 5*dt^4];
            %         mdt=[mdt;zeros(1,(j-1)*6) Dt zeros(1,(nbphases-j)*6)];
                     M_t_d(rowinitphase+i,(j-1)*6+1:(j-1)*6+6)=Dt;
                end
            end
        end
        function [M_t_dd] = compute_M_t_dd(discretization,frequency)
            nbphases=length(discretization);
            nbpoints=sum(discretization)+1;
            M_t_dd=zeros(nbpoints,nbphases*6);
            % mdt=[1 0 0 0 0 0 zeros(1,(nbphases-1)*6)];
            M_t_dd(1,1:6)=[0 0 2 0 0 0];
            nb=0;
            rowinitphase=1;
            for j=1:nbphases
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency;%dt is delta(tj)=t-Tj
                    Dt=[0 0 2*dt^0 6*dt^1 12*dt^2 20*dt^3];
            %         mdt=[mdt;zeros(1,(j-1)*6) Dt zeros(1,(nbphases-j)*6)];
                        M_t_dd(rowinitphase+i,(j-1)*6+1:(j-1)*6+6)=Dt;
                end
            end
        end
        function [M_t1] = compute_M_t1(discretization,frequency,cutting)
            %compute the double derivative of M_t
            nbphases=length(discretization);
            nbpoints=sum(discretization)+1;
            M_t1=zeros(nbpoints,nbphases*6+6+6);
            % mdt=[1 0 0 0 0 0 zeros(1,(nbphases-1)*6)];
            M_t1(1,1:6)=[1 0 0 0 0 0];

            %%%zmp1 init cut in two%%%
            nb=discretization(1);
            rowinitphase=1;
            for i=1:ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[dt^0 dt^1 dt^2 dt^3 dt^4 dt^5];
                M_t1(rowinitphase+i,(1-1)*6+1:(1-1)*6+6)=Dt;
            end
            rowinitphase=1+ceil(nb*cutting);
            for i=1:nb-ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[dt^0 dt^1 dt^2 dt^3 dt^4 dt^5];
                M_t1(rowinitphase+i,(2-1)*6+1:(2-1)*6+6)=Dt;
            end

            rowinitphase=1;
            %%%generation of dt%%%
            for j=2:nbphases-1
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency;%dt is delta(tj)=t-Tj
                    Dt=[dt^0 dt^1 dt^2 dt^3 dt^4 dt^5];
                        M_t1(rowinitphase+i,(j+1-1)*6+1:(j+1-1)*6+6)=Dt;
                end
            end

            %%%zmp1 fin cut in two%%%
            rowinitphase=rowinitphase+nb;
            nb=discretization(end);
            for i=1:ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[dt^0 dt^1 dt^2 dt^3 dt^4 dt^5];
                M_t1(rowinitphase+i,(nbphases+1-1)*6+1:(nbphases+1-1)*6+6)=Dt;
            end
            rowinitphase=rowinitphase+ceil(nb*cutting);
            for i=1:nb-ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[dt^0 dt^1 dt^2 dt^3 dt^4 dt^5];
                M_t1(rowinitphase+i,(nbphases+2-1)*6+1:(nbphases+2-1)*6+6)=Dt;
            end

        end
        function [M_t1_d] = compute_M_t1_d(discretization,frequency,cutting)
            %compute the derivative of M_t1
            nbphases=length(discretization);
            nbpoints=sum(discretization)+1;
            M_t1_d=zeros(nbpoints,nbphases*6+6+6);

            M_t1_d(1,1:6)=[0 1 0 0 0 0];

            %%%zmp1 init cut in two%%%
            nb=discretization(1);
            rowinitphase=1;
            for i=1:ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 1 2*dt^1 3*dt^2 4*dt^3 5*dt^4];
                M_t1_d(rowinitphase+i,(1-1)*6+1:(1-1)*6+6)=Dt;
            end
            rowinitphase=1+ceil(nb*cutting);
            for i=1:nb-ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 1 2*dt^1 3*dt^2 4*dt^3 5*dt^4];
                M_t1_d(rowinitphase+i,(2-1)*6+1:(2-1)*6+6)=Dt;
            end

            rowinitphase=1;
            %%%generation of dt%%%
            for j=2:nbphases-1
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency;%dt is delta(tj)=t-Tj
                    Dt=[0 1 2*dt^1 3*dt^2 4*dt^3 5*dt^4];
                        M_t1_d(rowinitphase+i,(j+1-1)*6+1:(j+1-1)*6+6)=Dt;
                end
            end

            %%%zmp1 fin cut in two%%%
            rowinitphase=rowinitphase+nb;
            nb=discretization(end);
            for i=1:ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 1 2*dt^1 3*dt^2 4*dt^3 5*dt^4];
                M_t1_d(rowinitphase+i,(nbphases+1-1)*6+1:(nbphases+1-1)*6+6)=Dt;
            end
            rowinitphase=rowinitphase+ceil(nb*cutting);
            for i=1:nb-ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 1 2*dt^1 3*dt^2 4*dt^3 5*dt^4];
                M_t1_d(rowinitphase+i,(nbphases+2-1)*6+1:(nbphases+2-1)*6+6)=Dt;
            end
        end
        function [M_t1_dd] = compute_M_t1_dd(discretization,frequency,cutting)
            %compute the double derivative of M_t1
            nbphases=length(discretization);
            nbpoints=sum(discretization)+1;
            M_t1_dd=zeros(nbpoints,nbphases*6+6+6);

            M_t1_dd(1,1:6)=[0 0 2 0 0 0];

            %%%zmp1 init cut in two%%%
            nb=discretization(1);
            rowinitphase=1;
            for i=1:ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 0 2 6*dt^1 12*dt^2 20*dt^3];
                M_t1_dd(rowinitphase+i,(1-1)*6+1:(1-1)*6+6)=Dt;
            end
            rowinitphase=1+ceil(nb*cutting);
            for i=1:nb-ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 0 2 6*dt^1 12*dt^2 20*dt^3];
                M_t1_dd(rowinitphase+i,(2-1)*6+1:(2-1)*6+6)=Dt;
            end

            rowinitphase=1;
            %%%generation of dt%%%
            for j=2:nbphases-1
                rowinitphase=rowinitphase+nb;
                nb=discretization(j);
                for i=1:nb
                    dt=i/frequency;%dt is delta(tj)=t-Tj
                    Dt=[0 0 2 6*dt^1 12*dt^2 20*dt^3];
                        M_t1_dd(rowinitphase+i,(j+1-1)*6+1:(j+1-1)*6+6)=Dt;
                end
            end

            %%%zmp1 fin cut in two%%%
            rowinitphase=rowinitphase+nb;
            nb=discretization(end);
            for i=1:ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 0 2 6*dt^1 12*dt^2 20*dt^3];
                M_t1_dd(rowinitphase+i,(nbphases+1-1)*6+1:(nbphases+1-1)*6+6)=Dt;
            end
            rowinitphase=rowinitphase+ceil(nb*cutting);
            for i=1:nb-ceil(nb*cutting)
                dt=i/frequency;%dt is delta(tj)=t-Tj
                Dt=[0 0 2 6*dt^1 12*dt^2 20*dt^3];
                M_t1_dd(rowinitphase+i,(nbphases+2-1)*6+1:(nbphases+2-1)*6+6)=Dt;
            end
        end


    end
end