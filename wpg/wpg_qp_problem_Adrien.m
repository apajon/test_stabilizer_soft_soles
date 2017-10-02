classdef wpg_qp_problem_Adrien<handle

    properties
        %quadratic problem properties for a quadprog optimization
        H
        f
        G
        A
        b
        Aeq
        beq
        C
        Hcom
        Gcom
        Ccom
        %storage of gradients of the problem
        A_xCa
        B_xCa
        A_yCa
        B_yCa
        
        A_xzmp
        B_xzmp
        A_yzmp
        B_yzmp
        
        A_xzmp_spd
        B_xzmp_spd
        A_yzmp_spd
        B_yzmp_spd
        
        A_xzmp_acc
        B_xzmp_acc
        A_yzmp_acc
        B_yzmp_acc
        
        A_xzmp_acc2
        B_xzmp_acc2
        C_xzmp_acc2
        A_yzmp_acc2
        B_yzmp_acc2
        C_yzmp_acc2
        
        A_xcom
        B_xcom
        A_ycom
        B_ycom
        
        A_xcom_spd
        B_xcom_spd
        A_ycom_spd
        B_ycom_spd
        
        A_xfcom
        B_xfcom
        A_yfcom
        B_yfcom
        
        A_xfcom2
        B_xfcom2
        C_xfcom2
        A_yfcom2
        B_yfcom2
        C_yfcom2
        
        A_yt_ssp
        B_yt_ssp
        A_xt_ssp
        B_xt_ssp
        
        A_yt2_ssp
        B_yt2_ssp
        C_yt2_ssp
        A_xt2_ssp
        B_xt2_ssp
        C_xt2_ssp
        
        A_cons_zmp_stability
        B_cons_zmp_stability
        
        A_xcons_scom
        B_xcons_scom
        A_ycons_scom
        B_ycons_scom
        
        A_cons_scom
        B_cons_scom
        
        A_xCa1
        B_xCa1
        A_yCa1
        B_yCa1
        
        A_xzmp1
        B_xzmp1
        A_yzmp1
        B_yzmp1
        
        A_xzmp1_spd
        B_xzmp1_spd
        A_yzmp1_spd
        B_yzmp1_spd
        
        A_xzmp1_acc
        B_xzmp1_acc
        A_yzmp1_acc
        B_yzmp1_acc
        
        A_yt_zmp1
        B_yt_zmp1
        A_xt_zmp1
        B_xt_zmp1
        
        A_yt2_zmp1
        B_yt2_zmp1
        C_yt2_zmp1
        A_xt2_zmp1
        B_xt2_zmp1
        C_xt2_zmp1
        
        Ac
        k_diag
        k_diag_d
        k_diag_dd
        
        A_xzmp1_acc2
        B_xzmp1_acc2
        C_xzmp1_acc2
        A_yzmp1_acc2
        B_yzmp1_acc2
        C_yzmp1_acc2
        
        A_xzmp2
        B_xzmp2
        A_yzmp2
        B_yzmp2
        
        A_yt_zmp2
        B_yt_zmp2
        A_xt_zmp2
        B_xt_zmp2
        
        A_yt2_zmp2
        B_yt2_zmp2
        C_yt2_zmp2
        A_xt2_zmp2
        B_xt2_zmp2
        C_xt2_zmp2
        
        A_xzmp2_acc
        B_xzmp2_acc
        A_yzmp2_acc
        B_yzmp2_acc
        
        A_xzmp2_acc2
        B_xzmp2_acc2
        C_xzmp2_acc2
        A_yzmp2_acc2
        B_yzmp2_acc2
        C_yzmp2_acc2
        
        A_cons_zmp1_stab
        B_cons_zmp1_stab
        A_cons_zmp2_stab
        B_cons_zmp2_stab
        
        A_cons_stretch
        B_cons_stretch
        
        A_cons_ineq
        B_cons_ineq
        
        A_cons_ankle
        B_cons_ankle
    end
    methods
        function obj=wpg_qp_problem_Adrien(wpg_param)
            %%
            %%% compute the matrix of time discretization %%%
            M_t = sparse(obj.compute_M_t(wpg_param));
            M_t_d = sparse(obj.compute_M_t_d(wpg_param));
            M_t_dd = sparse(obj.compute_M_t_dd(wpg_param));
            
            obj.compute_ZMP(wpg_param,M_t,M_t_d,M_t_dd);
            obj.compute_COM(wpg_param,M_t,M_t_d);
            obj.compute_cost_force(wpg_param);        
            obj.compute_cost_torque_ssp(wpg_param);
            obj.compute_cons_zmp(wpg_param);
                                   
            %%
            %where to cut the beginning and ending of zmp1 in percentage
            cutting=1/2;
            %%% compute the matrix of time discretization for ZMP1%%%
            %change from M_t because durig the initial and final phases
            %ZMP1 is cut into 2 polynomials
            M_t1=obj.compute_M_t(wpg_param);
            M_t1=sparse(M_t1(any(wpg_param.dt_type_phase==0,2),:));%suppress SSP discretization steps because ZMP1&2 are only defined in DSP
            
            %matrix of force repartition in DSP
            wpg_param_ = wpg_parameters(wpg_param.robot,wpg_param.type_traj,wpg_param.firstSS,wpg_param.frequency,wpg_param.lambda,wpg_param.epsilon,wpg_param.e,wpg_param.rightorleft,5,2,1,1,1);
            obj.Ac=obj.compute_force_repartition(wpg_param_,cutting);
            if (wpg_param.nbpolypi==0)
                obj.Ac(1:12)=[];
            end
            if (wpg_param.nbpolypf==0)
                obj.Ac(end-11:end)=[];
            end
            M_t1_=obj.compute_M_t1(wpg_param_,cutting);
            M_t1_=M_t1_(any(wpg_param_.dt_type_phase==0,2),:);%suppress SSP discretization steps because ZMP1&2 are only defined in DSP
            if (wpg_param.nbpolypi==0)
                M_t1_(:,1:12)=[];
                M_t1_(1:wpg_param_.discretization(1),:)=[];
            end
            if (wpg_param.nbpolypf==0)
                M_t1_(:,end-11:end)=[];
                M_t1_(end-wpg_param_.discretization(end):end,:)=[];
            end
            obj.k_diag=sparse(diag((M_t1_)*obj.Ac));
            
            M_t1_d=obj.compute_M_t1_d(wpg_param_,cutting);
            M_t1_d=M_t1_d(any(wpg_param_.dt_type_phase==0,2),:);%suppress SSP discretization steps because ZMP1&2 are only defined in DSP
            if (wpg_param.nbpolypi==0)
                M_t1_d(:,1:12)=[];
                M_t1_d(1:wpg_param_.discretization(1),:)=[];
            end
            if (wpg_param.nbpolypf==0)
                M_t1_d(:,end-11:end)=[];
                M_t1_d(end-wpg_param_.discretization(end):end,:)=[];
            end
            obj.k_diag_d=sparse(diag((M_t1_d)*obj.Ac));
            
            M_t1_dd=obj.compute_M_t1_dd(wpg_param_,cutting);
            M_t1_dd=M_t1_dd(any(wpg_param_.dt_type_phase==0,2),:);%suppress SSP discretization steps because ZMP1&2 are only defined in DSP
            if (wpg_param.nbpolypi==0)
                M_t1_dd(:,1:12)=[];
                M_t1_dd(1:wpg_param_.discretization(1),:)=[];
            end
            if (wpg_param.nbpolypf==0)
                M_t1_dd(:,end-11:end)=[];
                M_t1_dd(end-wpg_param_.discretization(end):end,:)=[];
            end
            obj.k_diag_dd=sparse(diag((M_t1_dd)*obj.Ac));
            
            clear('wpg_param_','M_t1_','M_t1_d','M_t1_dd');
            
            obj.compute_ZMP1(wpg_param,M_t1,cutting);
            obj.compute_cost_torque_zmp1(wpg_param);
            obj.compute_cost_acc_zmp1(wpg_param);
            obj.compute_ZMP2(wpg_param);
            obj.compute_cost_torque_zmp2(wpg_param);
            obj.compute_cost_acc_zmp2(wpg_param);
            obj.compute_cons_zmp12(wpg_param);
            
            obj.compute_cost_acc_zmp_SSP(wpg_param);
            
            %%
            obj.compute_H_G_C(wpg_param);
            obj.compute_A_b(wpg_param);
            
            %%
%             M_t = obj.compute_M_t(wpg_param);
%             M_t = M_t([1 end],:);
%             M_t_d = obj.compute_M_t_d(wpg_param);
%             M_t_d = M_t_d([1 end],:);
%             M_t_dd = obj.compute_M_t_dd(wpg_param);
%             M_t_dd = M_t_dd([1 end],:);
            
            obj.compute_Aeq_beq(wpg_param);            
            
            %%
%             wpg_param.nbcontrolpointzmp = wpg_param.nbcontrolpointzmp+4;
        end
        

        function obj=compute_ZMP(obj,wpg_param,M_t,M_t_d,M_t_dd)
            %The ZMP is computed from : zmp=A_zmp*X+B_zmp
            %A_zmp=M_t*A_Ca
            %B_zmp=M_t*B_Ca
            %A_Ca*X+B_Ca gives the zmp polynomial coefficients
            %% %compute matrix Ca
            %Ca is the matrix which multiplied by zmp via-points boundary
            %conditions gives a vector with the zmp polynomial coefficients
            [obj.A_xCa obj.B_xCa]=obj.compute_zmp_Ca(wpg_param);
            [obj.A_yCa obj.B_yCa]=obj.compute_zmp_Ca(wpg_param);
            %% %compute matrix A_zmp and B_zmp
            [obj.A_xzmp obj.B_xzmp] = obj.compute_traj_discrete(obj.A_xCa,obj.B_xCa,M_t);
            [obj.A_yzmp obj.B_yzmp] = obj.compute_traj_discrete(obj.A_yCa,obj.B_yCa,M_t);
            
            [obj.A_xzmp_spd obj.B_xzmp_spd] = obj.compute_traj_discrete(obj.A_xCa,obj.B_xCa,M_t_d);
            [obj.A_yzmp_spd obj.B_yzmp_spd] = obj.compute_traj_discrete(obj.A_yCa,obj.B_yCa,M_t_d);
            
            [obj.A_xzmp_acc obj.B_xzmp_acc] = obj.compute_traj_discrete(obj.A_xCa,obj.B_xCa,M_t_dd);
            [obj.A_yzmp_acc obj.B_yzmp_acc] = obj.compute_traj_discrete(obj.A_yCa,obj.B_yCa,M_t_dd);
        end
        function obj=compute_COM(obj,wpg_param,M_t,M_t_d)
            %the COM is computed from : com=A_com*X+B_com
            %com=M_cs*y+M_t*l
            %com=M_cs*(A_y*X+B_y)+M_t*(A_l*X+B_l)
            %A_y=G^-1*[A_A.A_Ca N 0]
            %B_y=G^-1*[A_A.B_Ca]
            %A_l=[A_A.A_Ca 0]
            %B_l=[A_A.B_Ca]
            %% %compute matrix G composed with sinh(wt) and cosh(wt)
            M_cs = obj.compute_com_M_cs(wpg_param);
            M_cs_d = obj.compute_com_M_cs_d(wpg_param);
           
            %% %A coeff of COM poly
            % A_gradient=com_morisawa_A_gradient(length(wpg_param.tpassage)-1,wpg_param.w);
%             A_A=obj.compute_com_A_A(length(wpg_param.tpassage)-1-8-3,wpg_param.w);
            A_A=obj.compute_com_A_A(wpg_param,wpg_param.nbphases,wpg_param.w);
            %% %VW gradient
            [A_xy B_xy]=obj.compute_com_VM(obj.A_xCa,obj.B_xCa,A_A,wpg_param.tpassage,wpg_param.w,wpg_param.poly_degree);
            [A_yy B_yy]=obj.compute_com_VM(obj.A_yCa,obj.B_yCa,A_A,wpg_param.tpassage,wpg_param.w,wpg_param.poly_degree);
            %% %Compute matrix A_com and B_com
            [obj.A_xcom obj.B_xcom]=obj.compute_com(obj.A_xCa,obj.B_xCa,A_xy,B_xy,A_A,M_t,M_cs);
            [obj.A_ycom obj.B_ycom]=obj.compute_com(obj.A_yCa,obj.B_yCa,A_yy,B_yy,A_A,M_t,M_cs);
            
            %% %Compute matrix A_com and B_com
            [obj.A_xcom_spd obj.B_xcom_spd]=obj.compute_com(obj.A_xCa,obj.B_xCa,A_xy,B_xy,A_A,M_t_d,M_cs_d);
            [obj.A_ycom_spd obj.B_ycom_spd]=obj.compute_com(obj.A_yCa,obj.B_yCa,A_yy,B_yy,A_A,M_t_d,M_cs_d);
            
            %% %Constraint eq for com speed initial and final
            %These constraints matrix are defined here because they are a
            %part of COM definition
            %Warning : A*P+B=0 => put -B in constraint in optimized function
            [obj.A_xcons_scom obj.B_xcons_scom] = obj.cons_scom_initfin(obj.A_xCa,obj.B_xCa,A_xy,B_xy,A_A,wpg_param.tpassage,wpg_param.nbphases,wpg_param.w,wpg_param.poly_degree);
%             obj.A_xcons_scom=obj.A_xcons_scom(1,:)-obj.A_xcons_scom(2,:);
%             obj.B_xcons_scom=obj.B_xcons_scom(1,:)-obj.B_xcons_scom(2,:);
            
            [obj.A_ycons_scom obj.B_ycons_scom] = obj.cons_scom_initfin(obj.A_yCa,obj.B_yCa,A_yy,B_yy,A_A,wpg_param.tpassage,wpg_param.nbphases,wpg_param.w,wpg_param.poly_degree);
%             obj.A_ycons_scom=obj.A_ycons_scom(1,:)+obj.A_ycons_scom(2,:);
%             obj.B_ycons_scom=obj.B_ycons_scom(1,:)+obj.B_ycons_scom(2,:);
            %% %%%concatenation%%%
            obj.B_cons_scom=[obj.B_xcons_scom;obj.B_ycons_scom];
%             obj.A_cons_scom=[obj.A_xcons_scom  zeros(size(obj.A_xcons_scom,1),wpg_param.nbparamank) zeros(size(obj.A_xcons_scom,1),wpg_param.nbcontrolpointzmp1)  zeros(size(obj.A_ycons_scom)) zeros(size(obj.A_ycons_scom,1),wpg_param.nbparamank) zeros(size(obj.A_xcons_scom,1),wpg_param.nbcontrolpointzmp1);
%                       zeros(size(obj.A_ycons_scom)) zeros(size(obj.A_ycons_scom,1),wpg_param.nbparamank) zeros(size(obj.A_ycons_scom,1),wpg_param.nbcontrolpointzmp1) obj.A_ycons_scom  zeros(size(obj.A_ycons_scom,1),wpg_param.nbparamank) zeros(size(obj.A_xcons_scom,1),wpg_param.nbcontrolpointzmp1)];
            obj.A_cons_scom=[obj.A_xcons_scom  zeros(size(obj.A_xcons_scom,1),wpg_param.nbpankle) zeros(size(obj.A_xcons_scom,1),wpg_param.nbcontrolpointzmp1)  zeros(size(obj.A_ycons_scom)) zeros(size(obj.A_ycons_scom,1),wpg_param.nbpankle) zeros(size(obj.A_xcons_scom,1),wpg_param.nbcontrolpointzmp1);
                                  zeros(size(obj.A_ycons_scom)) zeros(size(obj.A_ycons_scom,1),wpg_param.nbpankle) zeros(size(obj.A_ycons_scom,1),wpg_param.nbcontrolpointzmp1) obj.A_ycons_scom  zeros(size(obj.A_ycons_scom,1),wpg_param.nbpankle) zeros(size(obj.A_xcons_scom,1),wpg_param.nbcontrolpointzmp1)];
        end
        function obj=compute_cost_force(obj,wpg_param)
            %% %compute force on com in x and y axis
            %fcom=A_fcom*X+B_fcom=m*com_acceleration=m*com_a
            %com_a=g/z_com*(com-zmp)=g/z_com*[(A_com-A_zmp)*X+B_com-B_zmp]
            [obj.A_xfcom obj.B_xfcom]=obj.compute_force(obj.A_xzmp,obj.B_xzmp,obj.A_xcom,obj.B_xcom,wpg_param.mg,wpg_param.z);
            [obj.A_yfcom obj.B_yfcom]=obj.compute_force(obj.A_yzmp,obj.B_yzmp,obj.A_ycom,obj.B_ycom,wpg_param.mg,wpg_param.z);

            %% %compute the cost function criteria linked to the force on coml
            %compute the square of the force on com
            %fcom²=X^T*A_fcom2*X+2*B_fcom2*X+C_fcom2
            [obj.A_xfcom2 obj.B_xfcom2 obj.C_xfcom2] = obj.compute_quad_matrix(obj.A_xfcom,obj.B_xfcom);
            [obj.A_yfcom2 obj.B_yfcom2 obj.C_yfcom2] = obj.compute_quad_matrix(obj.A_yfcom,obj.B_yfcom);

        end
        function obj=compute_cost_torque_ssp(obj,wpg_param)
            %% %compute torques in ankle during SSP in x and y axis
            
            %% %compute the ankle position for the torque computation
            %an_ssp=A_an_ssp*X+B_an_ssp
            [A_xan_ssp B_xan_ssp]=obj.compute_ankle_positions_SSP(wpg_param,wpg_param.pankinit_firstSS(1),wpg_param.pankfin_lastSS(1));
            [A_yan_ssp B_yan_ssp]=obj.compute_ankle_positions_SSP(wpg_param,wpg_param.pankinit_firstSS(2),wpg_param.pankfin_lastSS(2));
            
            %% %torque in ankle
            %compute the torques in ankle in x and y axis
            %?_SSP=A_t_SSP*X+B_t_SSP
            %Warning : torques in x-axis are function of optimization
            %parameters in y-axis and reverse for the torques in y-axis.
            [obj.A_yt_ssp obj.B_yt_ssp] = obj.compute_torque_ssp(obj.A_xzmp,obj.B_xzmp,obj.A_xfcom,obj.B_xfcom,A_xan_ssp,B_xan_ssp,wpg_param.mg,wpg_param.ha);
%             obj.B_yt_ssp=obj.B_yt_ssp+obj.A_yt_ssp(:,end-2:end)*wpg_param.step_number_pankle_fixed(:,2);
%             obj.A_yt_ssp(:,end-2:end)=[];

            [obj.A_xt_ssp obj.B_xt_ssp] = obj.compute_torque_ssp(obj.A_yzmp,obj.B_yzmp,obj.A_yfcom,obj.B_yfcom,A_yan_ssp,B_yan_ssp,wpg_param.mg,wpg_param.ha);
%             obj.B_xt_ssp=obj.B_xt_ssp+obj.A_xt_ssp(:,end-2:end)*wpg_param.step_number_pankle_fixed(:,3);
%             obj.A_xt_ssp(:,end-2:end)=[];
            
            %% %torque 2 in ankle
            %% %compute the cost function criteria linked to the torques in ankle in SSP
            %compute the square of the force on com
            %T_SSP²=X^T*A_t2_SSP*X+2*B_t2_SSP*X+C_t2_SSP
            %Warning : torques in x-axis are function of optimization
            %parameters in y-axis and reverse for the torques in y-axis.
            [obj.A_yt2_ssp obj.B_yt2_ssp obj.C_yt2_ssp] = obj.compute_quad_matrix(obj.A_yt_ssp,obj.B_yt_ssp);
            [obj.A_xt2_ssp obj.B_xt2_ssp obj.C_xt2_ssp] = obj.compute_quad_matrix(obj.A_xt_ssp,obj.B_xt_ssp);
        end 
        function obj=compute_cons_zmp(obj,wpg_param)
            %% %compute the inequality constraint equations in SSP
            %% %compute the ankle position for the torque computation used in the stability constraint
            %an_ssp=A_an_ssp*X+B_an_ssp
            [A_xan_ssp B_xan_ssp]=obj.compute_ankle_positions_SSP(wpg_param,wpg_param.pankinit_firstSS(1),wpg_param.pankfin_lastSS(1));
            [A_yan_ssp B_yan_ssp]=obj.compute_ankle_positions_SSP(wpg_param,wpg_param.pankinit_firstSS(2),wpg_param.pankfin_lastSS(2));
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
            %sin(-a)*OZ>=-inttoankle
            %if the suppport foot us the right foot:
            %sin(a)*OZ<=inttoankle
            %sin(-a)*OZ>=-exttoankle
            %
            %we must verify :
            %A_cons_zmp_stability*X <= B_cons_zmp_stability
            [obj.A_cons_zmp_stability obj.B_cons_zmp_stability]=obj.cons_xy_zmp_stability_SSP(obj.A_xzmp,obj.B_xzmp,obj.A_yzmp,obj.B_yzmp,A_xan_ssp,B_xan_ssp,A_yan_ssp,B_yan_ssp,wpg_param.discretization,wpg_param.backtoankle,wpg_param.fronttoankle,wpg_param.exttoankle,wpg_param.inttoankle,wpg_param.sole_margin,wpg_param.psi_zmp,wpg_param.type_phase,wpg_param.nbcontrolpointzmp1,wpg_param.firstSS);
        end
        
        function obj=compute_ZMP1(obj,wpg_param,M_t1,cutting)
            %coeff gradient computation of zmp1_poly
            [obj.A_xCa1 obj.B_xCa1]=obj.compute_zmp1_Ca1(wpg_param,wpg_param.tpassage,wpg_param.xpsa_zmp1init,wpg_param.xpsa_zmp1fin,cutting);
            [obj.A_yCa1 obj.B_yCa1]=obj.compute_zmp1_Ca1(wpg_param,wpg_param.tpassage,wpg_param.ypsa_zmp1init,wpg_param.ypsa_zmp1fin,cutting);
            %position of ZMP1
            [obj.A_xzmp1 obj.B_xzmp1]=obj.compute_traj_discrete(obj.A_xCa1,obj.B_xCa1,M_t1);
            [obj.A_yzmp1 obj.B_yzmp1]=obj.compute_traj_discrete(obj.A_yCa1,obj.B_yCa1,M_t1);
            %compute the derivative of M_t1
            M_t1_d=obj.compute_M_t_d(wpg_param);
            M_t1_d=M_t1_d(any(wpg_param.dt_type_phase==0,2),:);
            %acceleration of zmp1
            [obj.A_xzmp1_spd obj.B_xzmp1_spd]=obj.compute_traj_discrete(obj.A_xCa1,obj.B_xCa1,M_t1_d);
            [obj.A_yzmp1_spd obj.B_yzmp1_spd]=obj.compute_traj_discrete(obj.A_yCa1,obj.B_yCa1,M_t1_d);
            %compute the double derivative of M_t1
            M_t1_dd=obj.compute_M_t_dd(wpg_param);
            M_t1_dd=M_t1_dd(any(wpg_param.dt_type_phase==0,2),:);
            %acceleration of zmp1
            [obj.A_xzmp1_acc obj.B_xzmp1_acc]=obj.compute_traj_discrete(obj.A_xCa1,obj.B_xCa1,M_t1_dd);
            [obj.A_yzmp1_acc obj.B_yzmp1_acc]=obj.compute_traj_discrete(obj.A_yCa1,obj.B_yCa1,M_t1_dd);
        end
        function obj=compute_cost_torque_zmp1(obj,wpg_param)
            %% %ankle position uses for torques by zmp1
            [A_xan_dsp B_xan_dsp]=obj.compute_ankle_positions_zmp1(wpg_param,wpg_param.pankinit_firstinair(1),wpg_param.pankinit_firstSS(1),wpg_param.pankfin_lastSS(1));
            [A_yan_dsp B_yan_dsp]=obj.compute_ankle_positions_zmp1(wpg_param,wpg_param.pankinit_firstinair(2),wpg_param.pankinit_firstSS(2),wpg_param.pankfin_lastSS(2));
            %yBpankle1(1:wpg_param.discretization(1)+1)=-yBpankle1(1:wpg_param.discretization(1)+1);%les pieds sont parallèles au debut, La création des pankle commence par le pied gauche or zmp1 est sous le pieds droit
            
            %suppress SSP discretization steps because ZMP1 is only defined in DSP
            A_xan_dsp=A_xan_dsp(any(wpg_param.dt_type_phase==0,2),:);
            B_xan_dsp=B_xan_dsp(any(wpg_param.dt_type_phase==0,2),:);
            A_yan_dsp=A_yan_dsp(any(wpg_param.dt_type_phase==0,2),:);
            B_yan_dsp=B_yan_dsp(any(wpg_param.dt_type_phase==0,2),:);

            %% %torques by zmp1
            %suppress SSP discretization steps because ZMP1 is only defined in DSP
            A_xfcom_=obj.A_xfcom(any(wpg_param.dt_type_phase==0,2),:);
            B_xfcom_=obj.B_xfcom(any(wpg_param.dt_type_phase==0,2),:);
            A_yfcom_=obj.A_yfcom(any(wpg_param.dt_type_phase==0,2),:);
            B_yfcom_=obj.B_yfcom(any(wpg_param.dt_type_phase==0,2),:);
            %Warning : torques in x-axis are function of optimization
            %parameters in y-axis and reverse for the torques in y-axis.
            [obj.A_yt_zmp1 obj.B_yt_zmp1] = obj.compute_torque_zmp1(obj.A_xzmp1,obj.B_xzmp1,A_xfcom_,B_xfcom_,A_xan_dsp,B_xan_dsp,wpg_param.mg,wpg_param.ha,obj.k_diag);
%             obj.B_yt_zmp1=obj.B_yt_zmp1+obj.A_yt_zmp1(:,end-2-3:end-3)*[wpg_param.step_number_pankle_fixed(:,2)];
%             obj.A_yt_zmp1(:,end-2-3:end-3)=[];
            
            [obj.A_xt_zmp1 obj.B_xt_zmp1] = obj.compute_torque_zmp1(obj.A_yzmp1,obj.B_yzmp1,A_yfcom_,B_yfcom_,A_yan_dsp,B_yan_dsp,wpg_param.mg,wpg_param.ha,obj.k_diag);
%             obj.B_xt_zmp1=obj.B_xt_zmp1+obj.A_xt_zmp1(:,end-2-3:end-3)*[wpg_param.step_number_pankle_fixed(:,3)];
%             obj.A_xt_zmp1(:,end-2-3:end-3)=[];
            %Warning : torques in x-axis are function of optimization
            %parameters in y-axis and reverse for the torques in y-axis.
            [obj.A_yt2_zmp1 obj.B_yt2_zmp1 obj.C_yt2_zmp1] = obj.compute_quad_matrix(obj.A_yt_zmp1,obj.B_yt_zmp1);
            [obj.A_xt2_zmp1 obj.B_xt2_zmp1 obj.C_xt2_zmp1] = obj.compute_quad_matrix(obj.A_xt_zmp1,obj.B_xt_zmp1);
        end
        function obj=compute_cost_acc_zmp1(obj,wpg_param)
            %% %%%cost function acceleration zmp1 & 2%%%
%             %compute the double derivative of M_t1
%             M_t1_dd=obj.compute_M_t_dd(wpg_param);
%             M_t1_dd=M_t1_dd(any(wpg_param.dt_type_phase==0,2),:);
%             %% %acceleration of zmp1
%             [obj.A_xzmp1_acc obj.B_xzmp1_acc]=obj.compute_traj_discrete(obj.A_xCa1,obj.B_xCa1,M_t1_dd);
%             [obj.A_yzmp1_acc obj.B_yzmp1_acc]=obj.compute_traj_discrete(obj.A_yCa1,obj.B_yCa1,M_t1_dd);
            %% %acceleration2 of zmp1
            [obj.A_xzmp1_acc2 obj.B_xzmp1_acc2 obj.C_xzmp1_acc2]=obj.compute_quad_matrix(obj.A_xzmp1_acc,obj.B_xzmp1_acc);
            [obj.A_yzmp1_acc2 obj.B_yzmp1_acc2 obj.C_yzmp1_acc2]=obj.compute_quad_matrix(obj.A_yzmp1_acc,obj.B_yzmp1_acc);
        end
        function obj=compute_ZMP2(obj,wpg_param)
            %% %position of ZMP2
            A_xzmp_=obj.A_xzmp(any(wpg_param.dt_type_phase==0,2),:);
            B_xzmp_=obj.B_xzmp(any(wpg_param.dt_type_phase==0,2),:);
            A_yzmp_=obj.A_yzmp(any(wpg_param.dt_type_phase==0,2),:);
            B_yzmp_=obj.B_yzmp(any(wpg_param.dt_type_phase==0,2),:);
%             reduce=any(diag(obj.k_diag)~=1,2);
            
            indices=[];
            for i=0:2
                indices=[indices;find(any(diag(obj.k_diag)==1,2))+i];
            end
            reduce=1:length(diag(obj.k_diag));
            reduce(indices)=[];
            
            [obj.A_xzmp2 obj.B_xzmp2]=obj.compute_zmp2(A_xzmp_(reduce,:),B_xzmp_(reduce,:),obj.A_xzmp1(reduce,:),obj.B_xzmp1(reduce,:),obj.k_diag(reduce,reduce),wpg_param.mg);
            [obj.A_yzmp2 obj.B_yzmp2]=obj.compute_zmp2(A_yzmp_(reduce,:),B_yzmp_(reduce,:),obj.A_yzmp1(reduce,:),obj.B_yzmp1(reduce,:),obj.k_diag(reduce,reduce),wpg_param.mg);
        end
        function obj=compute_cost_torque_zmp2(obj,wpg_param)
            %% %ankle position uses for torques by zmp2
            [A_xan_dsp B_xan_dsp]=obj.compute_ankle_positions_zmp2(wpg_param,wpg_param.pankinit_firstSS(1),wpg_param.pankfin_lastSS(1),wpg_param.pankfin_lastinair(1));
            [A_yan_dsp B_yan_dsp]=obj.compute_ankle_positions_zmp2(wpg_param,wpg_param.pankinit_firstSS(2),wpg_param.pankfin_lastSS(2),wpg_param.pankfin_lastinair(2));
    %         yBpankle2(end-discretization(end):end)=-yBpankle2(end-discretization(end):end);
    %         yBpankle2(end-wpg_param.discretization(end)+1:end)=yBpankle2(end-wpg_param.discretization(end)+1:end)-(-1)^(wpg_param.nbstep)*0.095*2;
            A_xan_dsp=A_xan_dsp(any(wpg_param.dt_type_phase==0,2),:);
            B_xan_dsp=B_xan_dsp(any(wpg_param.dt_type_phase==0,2),:);
            A_yan_dsp=A_yan_dsp(any(wpg_param.dt_type_phase==0,2),:);
            B_yan_dsp=B_yan_dsp(any(wpg_param.dt_type_phase==0,2),:);
            %% %torques by zmp1
            %suppress SSP discretization steps because ZMP1 is only defined in DSP
            A_xfcom_=obj.A_xfcom(any(wpg_param.dt_type_phase==0,2),:);
            B_xfcom_=obj.B_xfcom(any(wpg_param.dt_type_phase==0,2),:);
            A_yfcom_=obj.A_yfcom(any(wpg_param.dt_type_phase==0,2),:);
            B_yfcom_=obj.B_yfcom(any(wpg_param.dt_type_phase==0,2),:);
            %% %torques by zmp2
            %Warning : torques in x-axis are function of optimization
            %parameters in y-axis and reverse for the torques in y-axis.
%             reduce=any(diag(obj.k_diag)~=1,2);
            
            indices=[];
            for i=0:2
                indices=[indices;find(any(diag(obj.k_diag)==1,2))+i];
            end
            reduce=1:length(diag(obj.k_diag));
            reduce(indices)=[];
            
            [obj.A_yt_zmp2 obj.B_yt_zmp2] = obj.compute_torque_zmp2(obj.A_xzmp2,obj.B_xzmp2,A_xfcom_(reduce,:),B_xfcom_(reduce,:),A_xan_dsp(reduce,:),B_xan_dsp(reduce,:),wpg_param.mg,wpg_param.ha,obj.k_diag(reduce,reduce));
            
            [obj.A_xt_zmp2 obj.B_xt_zmp2] = obj.compute_torque_zmp2(obj.A_yzmp2,obj.B_yzmp2,A_yfcom_(reduce,:),B_yfcom_(reduce,:),A_yan_dsp(reduce,:),B_yan_dsp(reduce,:),wpg_param.mg,wpg_param.ha,obj.k_diag(reduce,reduce));
            %Warning : torques in x-axis are function of optimization
            %parameters in y-axis and reverse for the torques in y-axis.
            [obj.A_yt2_zmp2 obj.B_yt2_zmp2 obj.C_yt2_zmp2] = obj.compute_quad_matrix(obj.A_yt_zmp2,obj.B_yt_zmp2);
            [obj.A_xt2_zmp2 obj.B_xt2_zmp2 obj.C_xt2_zmp2] = obj.compute_quad_matrix(obj.A_xt_zmp2,obj.B_xt_zmp2);
        end
        function obj=compute_cost_acc_zmp2(obj,wpg_param)
            %compute the double derivative of M_t1
%             M_t1_dd=obj.compute_M_t_dd(wpg_param,cutting);
%             M_t1_dd=M_t1_dd(any(wpg_param.dt_type_phase==0,2),:);
%             M_t1_dd=M_t1_dd(L,:);
            %% %speed of zmp1
%             M_t1_d=obj.compute_M_t_d(wpg_param);
%             M_t1_d=M_t1_d(any(wpg_param.dt_type_phase==0,2),:);
%             [xAszmp1 xBszmp1]=obj.compute_traj_discrete(obj.A_xCa1,obj.B_xCa1,M_t1_d);
%             [yAszmp1 yBszmp1]=obj.compute_traj_discrete(obj.A_yCa1,obj.B_yCa1,M_t1_d);
%             xAszmp1=obj.A_xzmp1_spd;
%             xBszmp1=obj.B_xzmp1_spd;
%             yAszmp1=obj.A_yzmp1_spd;
%             yBszmp1=obj.B_yzmp1_spd;
            
            %% %position of zmp
            A_xzmp_=obj.A_xzmp(any(wpg_param.dt_type_phase==0,2),:);
            B_xzmp_=obj.B_xzmp(any(wpg_param.dt_type_phase==0,2),:);
            A_yzmp_=obj.A_yzmp(any(wpg_param.dt_type_phase==0,2),:);
            B_yzmp_=obj.B_yzmp(any(wpg_param.dt_type_phase==0,2),:);
            %% %speed of zmp
%             M_t_d=obj.compute_M_t_d(wpg_param);
%             M_t_d=M_t_d(any(wpg_param.dt_type_phase==0,2),:);
%             [xAszmp xBszmp] = obj.compute_traj_discrete(obj.A_xCa,obj.B_xCa,M_t_d);
%             [yAszmp yBszmp] = obj.compute_traj_discrete(obj.A_yCa,obj.B_yCa,M_t_d);
            
            xAszmp=obj.A_xzmp_spd;
            xBszmp=obj.B_xzmp_spd;
            yAszmp=obj.A_yzmp_spd;
            yBszmp=obj.B_yzmp_spd;

            xAszmp_=xAszmp(any(wpg_param.dt_type_phase==0,2),:);
            xBszmp_=xBszmp(any(wpg_param.dt_type_phase==0,2),:);
            yAszmp_=yAszmp(any(wpg_param.dt_type_phase==0,2),:);
            yBszmp_=yBszmp(any(wpg_param.dt_type_phase==0,2),:);
            %% %acceleration of zmp
%             M_t_dd=obj.compute_M_t_dd(wpg_param);
%             M_t_dd=M_t_dd(any(wpg_param.dt_type_phase==0,2),:);
%             [xAazmp xBazmp] = obj.compute_traj_discrete(obj.A_xCa,obj.B_xCa,M_t_dd);
%             [yAazmp yBazmp] = obj.compute_traj_discrete(obj.A_yCa,obj.B_yCa,M_t_dd);
            
            xAazmp=obj.A_xzmp_acc;
            xBazmp=obj.B_xzmp_acc;
            yAazmp=obj.A_yzmp_acc;
            yBazmp=obj.B_yzmp_acc;

            xAazmp_=xAazmp(any(wpg_param.dt_type_phase==0,2),:);
            xBazmp_=xBazmp(any(wpg_param.dt_type_phase==0,2),:);
            yAazmp_=yAazmp(any(wpg_param.dt_type_phase==0,2),:);
            yBazmp_=yBazmp(any(wpg_param.dt_type_phase==0,2),:);
            %% %acceleration of zmp2
            indices=[];
            for i=0:2
                indices=[indices;find(any(diag(obj.k_diag)==1,2))+i];
            end
            reduce=1:length(diag(obj.k_diag));
            reduce(indices)=[];
            
            
%             toto=min([ceil(1/100/(1/wpg_param.frequency)) size(A_xzmp_,1)])+1;
%             toto=2;
            
            [obj.A_xzmp2_acc obj.B_xzmp2_acc]=obj.compute_zmp2_acc(A_xzmp_(reduce,:),B_xzmp_(reduce,:),obj.A_xzmp1(reduce,:),obj.B_xzmp1(reduce,:),xAszmp_(reduce,:),xBszmp_(reduce,:),obj.A_xzmp1_spd(reduce,:),obj.B_xzmp1_spd(reduce,:),xAazmp_(reduce,:),xBazmp_(reduce,:),obj.A_xzmp1_acc(reduce,:),obj.B_xzmp1_acc(reduce,:),obj.k_diag(reduce,reduce),obj.k_diag_d(reduce,reduce),obj.k_diag_dd(reduce,reduce));
%             obj.A_xzmp2_acc(:,end-2-3:end-3)=[];
            [obj.A_yzmp2_acc obj.B_yzmp2_acc]=obj.compute_zmp2_acc(A_yzmp_(reduce,:),B_yzmp_(reduce,:),obj.A_yzmp1(reduce,:),obj.B_yzmp1(reduce,:),yAszmp_(reduce,:),yBszmp_(reduce,:),obj.A_yzmp1_spd(reduce,:),obj.B_yzmp1_spd(reduce,:),yAazmp_(reduce,:),yBazmp_(reduce,:),obj.A_yzmp1_acc(reduce,:),obj.B_yzmp1_acc(reduce,:),obj.k_diag(reduce,reduce),obj.k_diag_d(reduce,reduce),obj.k_diag_dd(reduce,reduce));
%             obj.A_yzmp2_acc(:,end-2-3:end-3)=[];
            %acceleration2 of zmp1
            [obj.A_xzmp2_acc2 obj.B_xzmp2_acc2 obj.C_xzmp2_acc2]=obj.compute_quad_matrix(obj.A_xzmp2_acc,obj.B_xzmp2_acc);
            [obj.A_yzmp2_acc2 obj.B_yzmp2_acc2 obj.C_yzmp2_acc2]=obj.compute_quad_matrix(obj.A_yzmp2_acc,obj.B_yzmp2_acc);
        end
        function obj=compute_cons_zmp12(obj,wpg_param)
            %% %ankle position uses for torques by zmp1
            [A_xan_dsp1 B_xan_dsp1]=obj.compute_ankle_positions_zmp1(wpg_param,wpg_param.pankinit_firstinair(1),wpg_param.pankinit_firstSS(1),wpg_param.pankfin_lastSS(1));
            [A_yan_dsp1 B_yan_dsp1]=obj.compute_ankle_positions_zmp1(wpg_param,wpg_param.pankinit_firstinair(2),wpg_param.pankinit_firstSS(2),wpg_param.pankfin_lastSS(2));
            
            %suppress SSP discretization steps because ZMP1 is only defined in DSP
            A_xan_dsp1=A_xan_dsp1(any(wpg_param.dt_type_phase==0,2),:);
            B_xan_dsp1=B_xan_dsp1(any(wpg_param.dt_type_phase==0,2),:);
            A_yan_dsp1=A_yan_dsp1(any(wpg_param.dt_type_phase==0,2),:);
            B_yan_dsp1=B_yan_dsp1(any(wpg_param.dt_type_phase==0,2),:);

            %% %ankle position uses for torques by zmp2
            [A_xan_dsp2 B_xan_dsp2]=obj.compute_ankle_positions_zmp2(wpg_param,wpg_param.pankinit_firstSS(1),wpg_param.pankfin_lastSS(1),wpg_param.pankfin_lastinair(1));
            [A_yan_dsp2 B_yan_dsp2]=obj.compute_ankle_positions_zmp2(wpg_param,wpg_param.pankinit_firstSS(2),wpg_param.pankfin_lastSS(2),wpg_param.pankfin_lastinair(2));

            %suppress SSP discretization steps because ZMP1 is only defined in DSP
            A_xan_dsp2=A_xan_dsp2(any(wpg_param.dt_type_phase==0,2),:);
            B_xan_dsp2=B_xan_dsp2(any(wpg_param.dt_type_phase==0,2),:);
            A_yan_dsp2=A_yan_dsp2(any(wpg_param.dt_type_phase==0,2),:);
            B_yan_dsp2=B_yan_dsp2(any(wpg_param.dt_type_phase==0,2),:);
            
%             reduce_zmp2=any(diag(obj.k_diag)~=1,2);
            
            indices=[];
            for i=0:2
                indices=[indices;find(any(diag(obj.k_diag)==1,2))+i];
            end
            reduce_zmp2=1:length(diag(obj.k_diag));
            reduce_zmp2(indices)=[];
            
            A_xan_dsp2=A_xan_dsp2(reduce_zmp2,:);
            B_xan_dsp2=B_xan_dsp2(reduce_zmp2,:);
            A_yan_dsp2=A_yan_dsp2(reduce_zmp2,:);
            B_yan_dsp2=B_yan_dsp2(reduce_zmp2,:);
            %% %%%constraint stability zmp1&2
            [obj.A_cons_zmp1_stab obj.B_cons_zmp1_stab]=obj.cons_xy_zmp1_stability(obj.A_xzmp1,obj.B_xzmp1,obj.A_yzmp1,obj.B_yzmp1,A_xan_dsp1,B_xan_dsp1,A_yan_dsp1,B_yan_dsp1,wpg_param.discretization,wpg_param.backtoankle,wpg_param.fronttoankle,wpg_param.exttoankle,wpg_param.inttoankle,wpg_param.sole_margin,wpg_param.psi_zmp1,wpg_param.type_phase,wpg_param.nbcontrolpointzmp,wpg_param.nbcontrolpointzmp1,wpg_param.firstSS);
            
            
            [obj.A_cons_zmp2_stab obj.B_cons_zmp2_stab]=obj.cons_xy_zmp2_stability(obj.A_xzmp2,obj.B_xzmp2,obj.A_yzmp2,obj.B_yzmp2,A_xan_dsp2,B_xan_dsp2,A_yan_dsp2,B_yan_dsp2,wpg_param.discretization,wpg_param.backtoankle,wpg_param.fronttoankle,wpg_param.exttoankle,wpg_param.inttoankle,wpg_param.sole_margin,wpg_param.psi_zmp2,wpg_param.type_phase,wpg_param.nbcontrolpointzmp,wpg_param.nbcontrolpointzmp1,wpg_param.firstSS);      
        end
        
        function obj=compute_cost_acc_zmp_SSP(obj,wpg_param)
            A_xzmp_acc=obj.A_xzmp_acc(any(wpg_param.dt_type_phase~=0,2),:);
            B_xzmp_acc=obj.B_xzmp_acc(any(wpg_param.dt_type_phase~=0,2),:);
            A_yzmp_acc=obj.A_yzmp_acc(any(wpg_param.dt_type_phase~=0,2),:);
            B_yzmp_acc=obj.B_yzmp_acc(any(wpg_param.dt_type_phase~=0,2),:);
            
            [obj.A_xzmp_acc2 obj.B_xzmp_acc2 obj.C_xzmp_acc2]=obj.compute_quad_matrix(A_xzmp_acc,B_xzmp_acc);
            [obj.A_yzmp_acc2 obj.B_yzmp_acc2 obj.C_yzmp_acc2]=obj.compute_quad_matrix(A_yzmp_acc,B_yzmp_acc);
        end
        
        function obj=compute_H_G_C(obj,wpg_param)
            %compute the coefficient of the cost function
            %cost=1/2*X^T*H*X+G*X+1/2*C
            %cost_quad=1/2*X^T*H*X+f'*X=1/2*X^T*H*X+G'*X
            n=1;
            q=1;
            m=1;
            o=1;
            p=1;
            r=1;
            s=1;
            %% Compute H
            xH_ssp=p*wpg_param.lambda*[obj.A_xfcom2 zeros(size(obj.A_xfcom2,1),size(obj.A_yt2_ssp,2)-size(obj.A_xfcom2,2)); ...
                zeros(size(obj.A_yt2_ssp,1)-size(obj.A_xfcom2,1),size(obj.A_xfcom2,2)) ... 
                zeros(size(obj.A_yt2_ssp,1)-size(obj.A_xfcom2,1),size(obj.A_yt2_ssp,2)-size(obj.A_xfcom2,2))] ...
                +p*wpg_param.mu*obj.A_yt2_ssp ...
                +o*wpg_param.epsilon*[obj.A_xzmp_acc2 zeros(size(obj.A_xzmp_acc2,1),size(obj.A_yt2_ssp,2)-size(obj.A_xzmp_acc2,2)); ...
                zeros(size(obj.A_yt2_ssp,1)-size(obj.A_xzmp_acc2,1),size(obj.A_xzmp_acc2,2)) ... 
                zeros(size(obj.A_yt2_ssp,1)-size(obj.A_xzmp_acc2,1),size(obj.A_yt2_ssp,2)-size(obj.A_xzmp_acc2,2))];
            xH=[xH_ssp zeros(size(xH_ssp,1),size(obj.A_yt2_zmp1,2)-size(xH_ssp,2)); ...
                zeros(size(obj.A_yt2_zmp1,1)-size(xH_ssp,1),size(obj.A_yt2_zmp1,2))] ...
                +m*0.5*wpg_param.mu*(s*obj.A_yt2_zmp1+r*obj.A_yt2_zmp2) ...
                +wpg_param.epsilon*(q*obj.A_xzmp1_acc2+n*obj.A_xzmp2_acc2);
            
            yH_ssp=p*wpg_param.lambda*[obj.A_yfcom2 zeros(size(obj.A_yfcom2,1),size(obj.A_xt2_ssp,2)-size(obj.A_yfcom2,2)); ...
                zeros(size(obj.A_xt2_ssp,1)-size(obj.A_yfcom2,1),size(obj.A_yfcom2,2)) ...
                zeros(size(obj.A_xt2_ssp,1)-size(obj.A_yfcom2,1),size(obj.A_xt2_ssp,2)-size(obj.A_yfcom2,2))] ...
                +p*wpg_param.mu*obj.A_xt2_ssp ...
                +o*wpg_param.epsilon*[obj.A_yzmp_acc2 zeros(size(obj.A_yzmp_acc2,1),size(obj.A_xt2_ssp,2)-size(obj.A_yzmp_acc2,2)); ...
                zeros(size(obj.A_xt2_ssp,1)-size(obj.A_yzmp_acc2,1),size(obj.A_yzmp_acc2,2)) ...
                zeros(size(obj.A_xt2_ssp,1)-size(obj.A_yzmp_acc2,1),size(obj.A_xt2_ssp,2)-size(obj.A_yzmp_acc2,2))];
            yH=[yH_ssp zeros(size(yH_ssp,1),size(obj.A_xt2_zmp1,2)-size(yH_ssp,2)); ...
                zeros(size(obj.A_xt2_zmp1,1)-size(yH_ssp,1),size(obj.A_xt2_zmp1,2))] ...
                +m*0.5*wpg_param.mu*(s*obj.A_xt2_zmp1+r*obj.A_xt2_zmp2) ...
                +wpg_param.epsilon*(q*obj.A_yzmp1_acc2+n*obj.A_yzmp2_acc2);
            
            obj.H=[xH zeros(size(yH));zeros(size(xH)) yH];
            
            %% Compute G
            xG_ssp=p*wpg_param.lambda*[obj.B_xfcom2 zeros(1,size(obj.B_yt2_ssp,2)-size(obj.B_xfcom2,2))] ...
                +p*wpg_param.mu*obj.B_yt2_ssp ...
                +o*wpg_param.epsilon*[obj.B_xzmp_acc2 zeros(1,size(obj.B_yt2_ssp,2)-size(obj.B_xzmp_acc2,2))];
            xG=[xG_ssp zeros(1,size(obj.B_yt2_zmp1,2)-size(xG_ssp,2))] ...
                +m*0.5*wpg_param.mu*(s*obj.B_yt2_zmp1+r*obj.B_yt2_zmp2) ...
                +wpg_param.epsilon*(q*obj.B_xzmp1_acc2+n*obj.B_xzmp2_acc2);
            
            yG_ssp=p*wpg_param.lambda*[obj.B_yfcom2 zeros(1,size(obj.B_xt2_ssp,2)-size(obj.B_yfcom2,2))] ...
                +p*wpg_param.mu*obj.B_xt2_ssp ...
                +o*wpg_param.epsilon*[obj.B_yzmp_acc2 zeros(1,size(obj.B_xt2_ssp,2)-size(obj.B_yzmp_acc2,2))];
            yG=[yG_ssp zeros(1,size(obj.B_xt2_zmp1,2)-size(yG_ssp,2))] ...
                +m*0.5*wpg_param.mu*(s*obj.B_xt2_zmp1+r*obj.B_xt2_zmp2) ...
                +wpg_param.epsilon*(q*obj.B_yzmp1_acc2+n*obj.B_yzmp2_acc2);
            
            obj.f=[xG yG]';
            obj.G=[xG yG];
            
            %% Compute C
            xC=p*wpg_param.lambda*(obj.C_xfcom2) ...
                +p*(1-wpg_param.lambda)*(obj.C_yt2_ssp) ...
                +m*0.5*(1-wpg_param.lambda)*(s*obj.C_yt2_zmp1+r*obj.C_yt2_zmp2) ...    
                +wpg_param.epsilon*(q*obj.C_xzmp1_acc2+n*obj.C_xzmp2_acc2) ...
                +o*wpg_param.epsilon*(obj.C_xzmp_acc2);
            
            yC=p*wpg_param.lambda*(obj.C_yfcom2) ...
                +p*(1-wpg_param.lambda)*(obj.C_xt2_ssp) ...
                +m*0.5*(1-wpg_param.lambda)*(s*obj.C_xt2_zmp1+r*obj.C_xt2_zmp2) ...    
                +wpg_param.epsilon*(q*obj.C_yzmp1_acc2+n*obj.C_yzmp2_acc2) ...
                +o*wpg_param.epsilon*(obj.C_yzmp_acc2);
            
            obj.C=xC+yC;
        end
        function obj=compute_A_b(obj,wpg_param)
            %% %%%contrainte de non chevauchement des pieds%%%
            if wpg_param.type_traj==1
                obj.A_cons_stretch=[];
                obj.B_cons_stretch=[];
            else
                [obj.A_cons_stretch obj.B_cons_stretch]=obj.cons_stretching(wpg_param);
%                 obj.A_cons_stretch=[];
%                 obj.B_cons_stretch=[];
            end
            
%             obj.A_cons_zmp1_stab=[];
%             obj.A_cons_zmp2_stab=[];
            
%             obj.B_cons_zmp1_stab=[];
%             obj.B_cons_zmp2_stab=[];
            
%             obj.A_cons_zmp_stability=[];
%             obj.B_cons_zmp_stability=[];
            
            obj.A=[obj.A_cons_zmp_stability;    obj.A_cons_zmp1_stab;   obj.A_cons_zmp2_stab;   obj.A_cons_stretch];
            obj.b=[obj.B_cons_zmp_stability;    obj.B_cons_zmp1_stab;   obj.B_cons_zmp2_stab;   obj.B_cons_stretch];
        end
        function obj=compute_Aeq_beq(obj,wpg_param)
            %% %%%Contrainte sur des ankle positions%%%
            [obj.A_cons_ankle obj.B_cons_ankle]=obj.cons_ankle_fixed_path(wpg_param.nbcontrolpointzmp,wpg_param.nbpankle,wpg_param.nbcontrolpointzmp1,wpg_param.step_number_pankle_fixed);
%             obj.A_cons_ankle=[];
%             obj.B_cons_ankle=[];
            Aeq_=[obj.A_cons_scom;obj.A_cons_ankle];
            beq_=[obj.B_cons_scom;obj.B_cons_ankle];
            
            %% %%%Constraints on ZMP init%%%
            xpA=obj.A_xzmp(1,:);
            Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
            beq_=[beq_;obj.B_xzmp(1,:)-wpg_param.xpsa_zmpinit(1)];
            
            ypA=obj.A_yzmp(1,:);
            Aeq_=[Aeq_;zeros(size(ypA,1),size(Aeq_,2)/2) ypA zeros(size(ypA,1),size(Aeq_,2)/2-size(ypA,2))];
            beq_=[beq_;obj.B_yzmp(1,:)-wpg_param.ypsa_zmpinit(1)];
            
            xpA=obj.A_xzmp_spd(1,:);
            Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
            beq_=[beq_;obj.B_xzmp_spd(1,:)-wpg_param.xpsa_zmpinit(2)];
            
            ypA=obj.A_yzmp_spd(1,:);
            Aeq_=[Aeq_;zeros(size(ypA,1),size(Aeq_,2)/2) ypA zeros(size(ypA,1),size(Aeq_,2)/2-size(ypA,2))];
            beq_=[beq_;obj.B_yzmp_spd(1,:)-wpg_param.ypsa_zmpinit(2)];
            
            xpA=obj.A_xzmp_acc(1,:);
            Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
            beq_=[beq_;obj.B_xzmp_acc(1,:)-wpg_param.xpsa_zmpinit(3)];
            
            ypA=obj.A_yzmp_acc(1,:);
            Aeq_=[Aeq_;zeros(size(ypA,1),size(Aeq_,2)/2) ypA zeros(size(ypA,1),size(Aeq_,2)/2-size(ypA,2))];
            beq_=[beq_;obj.B_xzmp_acc(1,:)-wpg_param.ypsa_zmpinit(3)];
            
            %% %%%Constraints on ZMP final%%%
            n=size(obj.A_xzmp,1);
            xpA=obj.A_xzmp(n,:);
            Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
            beq_=[beq_;obj.B_xzmp(n,:)-wpg_param.xpsa_zmpfin(1)];
            
            ypA=obj.A_yzmp(n,:);
            Aeq_=[Aeq_;zeros(size(ypA,1),size(Aeq_,2)/2) ypA zeros(size(ypA,1),size(Aeq_,2)/2-size(ypA,2))];
            beq_=[beq_;obj.B_yzmp(n,:)-wpg_param.ypsa_zmpfin(1)];
            
            xpA=obj.A_xzmp_spd(n,:);
            Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
            beq_=[beq_;obj.B_xzmp_spd(n,:)-wpg_param.xpsa_zmpfin(2)];
            
            ypA=obj.A_yzmp_spd(n,:);
            Aeq_=[Aeq_;zeros(size(ypA,1),size(Aeq_,2)/2) ypA zeros(size(ypA,1),size(Aeq_,2)/2-size(ypA,2))];
            beq_=[beq_;obj.B_yzmp_spd(n,:)-wpg_param.ypsa_zmpfin(2)];

            xpA=obj.A_xzmp_acc(n,:);
            Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
            beq_=[beq_;obj.B_xzmp_acc(n,:)-wpg_param.xpsa_zmpfin(3)];
            
            ypA=obj.A_yzmp_acc(n,:);
            Aeq_=[Aeq_;zeros(size(ypA,1),size(Aeq_,2)/2) ypA zeros(size(ypA,1),size(Aeq_,2)/2-size(ypA,2))];
            beq_=[beq_;obj.B_yzmp_acc(n,:)-wpg_param.ypsa_zmpfin(3)];
            
            %% %%%Constraints on ZMP1%%%
            if wpg_param.type_phase(1)==0
                %% %%%Constraints on ZMP1 init%%%
                xpA=obj.A_xzmp1(1,:);
                Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
                beq_=[beq_;obj.B_xzmp1(1,:)-wpg_param.xpsa_zmp1init(1)];

                ypA=obj.A_yzmp1(1,:);
                Aeq_=[Aeq_;zeros(size(ypA,1),size(Aeq_,2)/2) ypA zeros(size(ypA,1),size(Aeq_,2)/2-size(ypA,2))];
                beq_=[beq_;obj.B_yzmp1(1,:)-wpg_param.ypsa_zmp1init(1)];

                xpA=obj.A_xzmp1_spd(1,:);
                Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
                beq_=[beq_;obj.B_xzmp1_spd(1,:)-wpg_param.xpsa_zmp1init(2)];

                ypA=obj.A_yzmp1_spd(1,:);
                Aeq_=[Aeq_;zeros(size(ypA,1),size(Aeq_,2)/2) ypA zeros(size(ypA,1),size(Aeq_,2)/2-size(ypA,2))];
                beq_=[beq_;obj.B_yzmp1_spd(1,:)-wpg_param.ypsa_zmp1init(2)];

                xpA=obj.A_xzmp1_acc(1,:);
                Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
                beq_=[beq_;obj.B_xzmp1_acc(1,:)-wpg_param.xpsa_zmp1init(3)];

                ypA=obj.A_yzmp1_acc(1,:);
                Aeq_=[Aeq_;zeros(size(ypA,1),size(Aeq_,2)/2) ypA zeros(size(ypA,1),size(Aeq_,2)/2-size(ypA,2))];
                beq_=[beq_;obj.B_xzmp1_acc(1,:)-wpg_param.ypsa_zmp1init(3)];

                %% %%%Constraints on ZMP1 final%%%
                n=size(obj.A_xzmp1,1);            
                xpA=obj.A_xzmp1(n,:);
                Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
                beq_=[beq_;obj.B_xzmp1(n,:)-wpg_param.xpsa_zmp1fin(1)];

                ypA=obj.A_yzmp1(n,:);
                Aeq_=[Aeq_;zeros(size(ypA,1),size(Aeq_,2)/2) ypA zeros(size(ypA,1),size(Aeq_,2)/2-size(ypA,2))];
                beq_=[beq_;obj.B_yzmp1(n,:)-wpg_param.ypsa_zmp1fin(1)];

                xpA=obj.A_xzmp1_spd(n,:);
                Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
                beq_=[beq_;obj.B_xzmp1_spd(n,:)-wpg_param.xpsa_zmp1fin(2)];

                ypA=obj.A_yzmp1_spd(n,:);
                Aeq_=[Aeq_;zeros(size(ypA,1),size(Aeq_,2)/2) ypA zeros(size(ypA,1),size(Aeq_,2)/2-size(ypA,2))];
                beq_=[beq_;obj.B_yzmp1_spd(n,:)-wpg_param.ypsa_zmp1fin(2)];

                xpA=obj.A_xzmp1_acc(n,:);
                Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
                beq_=[beq_;obj.B_xzmp1_acc(n,:)-wpg_param.xpsa_zmp1fin(3)];

                ypA=obj.A_yzmp1_acc(n,:);
                Aeq_=[Aeq_;zeros(size(ypA,1),size(Aeq_,2)/2) ypA zeros(size(ypA,1),size(Aeq_,2)/2-size(ypA,2))];
                beq_=[beq_;obj.B_yzmp1_acc(n,:)-wpg_param.ypsa_zmp1fin(3)];
            end
            %% %%%Constraints on continuity ZMP-ZMP1%%%
            continuity_zmp_phases=find(wpg_param.type_phase==6);
            continuity_zmp=zeros(1,length(continuity_zmp_phases));
            for i=1:length(continuity_zmp_phases)
                continuity_zmp(i)=sum(wpg_param.discretization(1:continuity_zmp_phases(i)))+1;
            end
            if wpg_param.type_phase(end)~=0
                continuity_zmp(end)=[];
            end
            
            discretization_=wpg_param.discretization(any(wpg_param.type_phase==0,1));
            if wpg_param.type_phase(1)==0 
                type_phase_=[0 wpg_param.type_phase(find(any(wpg_param.type_phase(1:end-1)==0,1))+1)];
            else
                type_phase_=[1 wpg_param.type_phase(find(any(wpg_param.type_phase(1:end-1)==0,1))+1)];
            end
            continuity_zmp1_phases=find(type_phase_)-1;
            continuity_zmp1=zeros(1,length(continuity_zmp_phases));
            for i=1:length(continuity_zmp_phases)
                continuity_zmp1(i)=sum(discretization_(1:continuity_zmp1_phases(i)))+1;
            end
            if wpg_param.type_phase(end)~=0
                continuity_zmp1(end)=[];
            end
            
            xpA=[obj.A_xzmp(continuity_zmp,:) zeros(length(continuity_zmp),size(obj.A_xzmp1,2)-size(obj.A_xzmp,2))]-obj.A_xzmp1(continuity_zmp1,:);
            Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
            beq_=[beq_;obj.B_xzmp(continuity_zmp,:)-obj.B_xzmp1(continuity_zmp1,:)];
            
            ypA=[obj.A_yzmp(continuity_zmp,:) zeros(length(continuity_zmp),size(obj.A_yzmp1,2)-size(obj.A_yzmp,2))]-obj.A_yzmp1(continuity_zmp1,:);
            Aeq_=[Aeq_;ypA zeros(size(ypA,1),size(Aeq_,2)-size(ypA,2))];
            beq_=[beq_;obj.B_yzmp(continuity_zmp,:)-obj.B_yzmp1(continuity_zmp1,:)];
            
            xpA=[obj.A_xzmp_spd(continuity_zmp,:) zeros(length(continuity_zmp),size(obj.A_xzmp1_spd,2)-size(obj.A_xzmp_spd,2))]-obj.A_xzmp1_spd(continuity_zmp1,:);
            Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
            beq_=[beq_;obj.B_xzmp_spd(continuity_zmp,:)-obj.B_xzmp1_spd(continuity_zmp1,:)];
            
            ypA=[obj.A_yzmp_spd(continuity_zmp,:) zeros(length(continuity_zmp),size(obj.A_yzmp1_spd,2)-size(obj.A_yzmp_spd,2))]-obj.A_yzmp1_spd(continuity_zmp1,:);
            Aeq_=[Aeq_;ypA zeros(size(ypA,1),size(Aeq_,2)-size(ypA,2))];
            beq_=[beq_;obj.B_yzmp_spd(continuity_zmp,:)-obj.B_yzmp1_spd(continuity_zmp1,:)];
            
            xpA=[obj.A_xzmp_acc(continuity_zmp,:) zeros(length(continuity_zmp),size(obj.A_xzmp1_acc,2)-size(obj.A_xzmp_acc,2))]-obj.A_xzmp1_acc(continuity_zmp1,:);
            Aeq_=[Aeq_;xpA zeros(size(xpA,1),size(Aeq_,2)-size(xpA,2))];
            beq_=[beq_;obj.B_xzmp_acc(continuity_zmp,:)-obj.B_xzmp1_acc(continuity_zmp1,:)];
            
            ypA=[obj.A_yzmp_acc(continuity_zmp,:) zeros(length(continuity_zmp),size(obj.A_yzmp1_acc,2)-size(obj.A_yzmp_acc,2))]-obj.A_yzmp1_acc(continuity_zmp1,:);
            Aeq_=[Aeq_;ypA zeros(size(ypA,1),size(Aeq_,2)-size(ypA,2))];
            beq_=[beq_;obj.B_yzmp_acc(continuity_zmp,:)-obj.B_yzmp1_acc(continuity_zmp1,:)];
            
            %%
            %force a distance in x of 0.125m between the initial and final COM
            %force an equal speed in x between the initial and final COM
            %force a symetry in position and speed in y between the initial
            %and final COM
            xpA=[zeros(4,wpg_param.nbcontrolpointzmp) [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1] zeros(4,wpg_param.nbpankle) zeros(4,wpg_param.nbcontrolpointzmp1) zeros(4,wpg_param.nbcontrolpointzmp+4+wpg_param.nbpankle+wpg_param.nbcontrolpointzmp1)];
            Aeq_=[Aeq_;xpA];
            beq_=[beq_;-wpg_param.xpcominit;-wpg_param.xpcomfin;-wpg_param.xscominit;-wpg_param.xscomfin];
            ypA=[zeros(4,wpg_param.nbcontrolpointzmp+4+wpg_param.nbpankle+wpg_param.nbcontrolpointzmp1) zeros(4,wpg_param.nbcontrolpointzmp) [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1] zeros(4,wpg_param.nbpankle) zeros(4,wpg_param.nbcontrolpointzmp1)];
            Aeq_=[Aeq_;ypA];
            beq_=[beq_;-wpg_param.ypcominit;-wpg_param.ypcomfin;-wpg_param.yscominit;-wpg_param.yscomfin];
            
            %%
            obj.Aeq=Aeq_;
            obj.beq=beq_;
            
        end
        
        function [cost gradient]=compute_cost(obj,x)
            cost=1/2*x'*obj.H*x+obj.G*x+1/2*obj.C;
            gradient=1/2*x'*obj.H+obj.G;
        end
        
        function []=compute_cost_fcom2(obj,wpg_param)
            %compute the coefficient of the cost function
            %cost=1/2*X^T*H*X+G*X+1/2*C
            %cost_quad=1/2*X^T*H*X+f'*X=1/2*X^T*H*X+G'*X
            %% Compute H
            xH_ssp=wpg_param.lambda*[obj.A_xfcom2 zeros(size(obj.A_xfcom2,1),size(obj.A_yt2_ssp,2)-size(obj.A_xfcom2,2)); ...
                zeros(size(obj.A_yt2_ssp,1)-size(obj.A_xfcom2,1),size(obj.A_xfcom2,2)) ... 
                zeros(size(obj.A_yt2_ssp,1)-size(obj.A_xfcom2,1),size(obj.A_yt2_ssp,2)-size(obj.A_xfcom2,2))];
            xH=[xH_ssp zeros(size(xH_ssp,1),size(obj.A_yt2_zmp1,2)-size(xH_ssp,2)); ...
                zeros(size(obj.A_yt2_zmp1,1)-size(xH_ssp,1),size(obj.A_yt2_zmp1,2))];
            
            yH_ssp=wpg_param.lambda*[obj.A_yfcom2 zeros(size(obj.A_yfcom2,1),size(obj.A_xt2_ssp,2)-size(obj.A_yfcom2,2)); ...
                zeros(size(obj.A_xt2_ssp,1)-size(obj.A_yfcom2,1),size(obj.A_yfcom2,2)) ...
                zeros(size(obj.A_xt2_ssp,1)-size(obj.A_yfcom2,1),size(obj.A_xt2_ssp,2)-size(obj.A_yfcom2,2))];
            yH=[yH_ssp zeros(size(yH_ssp,1),size(obj.A_xt2_zmp1,2)-size(yH_ssp,2)); ...
                zeros(size(obj.A_xt2_zmp1,1)-size(yH_ssp,1),size(obj.A_xt2_zmp1,2))];
            
            obj.Hcom=[xH zeros(size(yH));zeros(size(xH)) yH];
            
            %% Compute G
            xG_ssp=wpg_param.lambda*[obj.B_xfcom2 zeros(1,size(obj.B_yt2_ssp,2)-size(obj.B_xfcom2,2))];
            xG=[xG_ssp zeros(1,size(obj.B_yt2_zmp1,2)-size(xG_ssp,2))];
            
            yG_ssp=wpg_param.lambda*[obj.B_yfcom2 zeros(1,size(obj.B_xt2_ssp,2)-size(obj.B_yfcom2,2))];
            yG=[yG_ssp zeros(1,size(obj.B_xt2_zmp1,2)-size(yG_ssp,2))];
            
            obj.Gcom=[xG yG];
            
            %% Compute C
            xC=wpg_param.lambda*(obj.C_xfcom2);
            
            yC=wpg_param.lambda*(obj.C_yfcom2);
            
            obj.Ccom=xC+yC;
            
            scaling = 10000;
            %cost=(1/2*x'*Hcom*x+Gcom*x+Ccom)/scaling; 
            %gradient=(x'*Hcom+Gcom)/scaling;
        end
        
        function []=reduce_constraint(obj,n)
            m=size(obj.B_cons_zmp_stability,1);
            A_cons_zmp_stability=[];
            B_cons_zmp_stability=[];
            for i=1:n:m
                A_cons_zmp_stability=[A_cons_zmp_stability;obj.A_cons_zmp_stability(i,:)];
                B_cons_zmp_stability=[B_cons_zmp_stability;obj.B_cons_zmp_stability(i,:)];
            end
            obj.A_cons_zmp_stability=A_cons_zmp_stability;
            obj.B_cons_zmp_stability=B_cons_zmp_stability;
            
            m=size(obj.B_cons_zmp1_stab,1);
            A_cons_zmp1_stab=[];
            B_cons_zmp1_stab=[];
            for i=1:n:m
                A_cons_zmp1_stab=[A_cons_zmp1_stab;obj.A_cons_zmp1_stab(i,:)];
                B_cons_zmp1_stab=[B_cons_zmp1_stab;obj.B_cons_zmp1_stab(i,:)];
            end
            obj.A_cons_zmp1_stab=A_cons_zmp1_stab;
            obj.B_cons_zmp1_stab=B_cons_zmp1_stab;
            
            m=size(obj.B_cons_zmp2_stab,1);
            A_cons_zmp2_stab=[];
            B_cons_zmp2_stab=[];
            for i=1:n:m
                A_cons_zmp2_stab=[A_cons_zmp2_stab;obj.A_cons_zmp2_stab(i,:)];
                B_cons_zmp2_stab=[B_cons_zmp2_stab;obj.B_cons_zmp2_stab(i,:)];
            end
            obj.A_cons_zmp2_stab=A_cons_zmp2_stab;
            obj.B_cons_zmp2_stab=B_cons_zmp2_stab;
            
            obj.A=[obj.A_cons_zmp_stability;    obj.A_cons_zmp1_stab;   obj.A_cons_zmp2_stab;   obj.A_cons_stretch];
            obj.b=[obj.B_cons_zmp_stability;    obj.B_cons_zmp1_stab;   obj.B_cons_zmp2_stab;   obj.B_cons_stretch];
        end
        
        function []=extend_constraint_symetry(obj,wpg_param)
            Aeq=obj.Aeq;
            beq=obj.beq;
            
            xpA=[zeros(1,3) [-1 0 0] [1 0 0] zeros(1,3) zeros(1,2+0) zeros(1,wpg_param.nbcontrolpointzmp1) zeros(1,wpg_param.nbcontrolpointzmp+0+wpg_param.nbcontrolpointzmp1)];
            Aeq=[Aeq;xpA];
            beq=[beq;-0.125/2];
            ypA=[zeros(1,wpg_param.nbcontrolpointzmp+0+wpg_param.nbcontrolpointzmp1) zeros(1,3) [1 0 0] [1 0 0] zeros(1,3) zeros(1,2+0) zeros(1,wpg_param.nbcontrolpointzmp1)];
            Aeq=[Aeq;ypA];
            beq=[beq;0];

            xsA=[zeros(1,3) [0 -1 0] [0 1 0] zeros(1,3) zeros(1,2+0) zeros(1,wpg_param.nbcontrolpointzmp1) zeros(1,wpg_param.nbcontrolpointzmp+0+wpg_param.nbcontrolpointzmp1)];
            Aeq=[Aeq;xsA];
            beq=[beq;0];
            ysA=[zeros(1,wpg_param.nbcontrolpointzmp+0+wpg_param.nbcontrolpointzmp1) zeros(1,3) [0 1 0] [0 1 0] zeros(1,3) zeros(1,2+0) zeros(1,wpg_param.nbcontrolpointzmp1) ];
            Aeq=[Aeq;ysA];
            beq=[beq;0];

            xaA=[zeros(1,3) [0 0 -1] [0 0 1] zeros(1,3) zeros(1,2+0) zeros(1,wpg_param.nbcontrolpointzmp1) zeros(1,wpg_param.nbcontrolpointzmp+0+wpg_param.nbcontrolpointzmp1)];
            Aeq=[Aeq;xaA];
            beq=[beq;0];
            yaA=[zeros(1,wpg_param.nbcontrolpointzmp+0+wpg_param.nbcontrolpointzmp1) zeros(1,3) [0 0 1] [0 0 1] zeros(1,3) zeros(1,2+0) zeros(1,wpg_param.nbcontrolpointzmp1) ];
            Aeq=[Aeq;yaA];
            beq=[beq;0];
            
            obj.Aeq=Aeq;
            obj.beq=beq;
        end
    end
    methods (Static) 
        %% intern methods of compute_ZMP
        function [A_Ca B_Ca]=compute_zmp_Ca(wpg_param)
        % Hypo : We consider ZMP trajectory as nth order bspline using DeBoor algorithme.
        % This function compute A and B as: C=A*x+B
        % Where : C are give the polynomials coefficient of ZMP trajectory in one direction.
        % x are the bspline control points.
          
            A_Ca=zeros((wpg_param.poly_degree+1)*wpg_param.nbphases,wpg_param.nbcontrolpointzmp);

%             for i=(wpg_param.poly_degree+1):wpg_param.nbcontrolpointzmp
            for i=1+wpg_param.poly_degree:wpg_param.nbphases+wpg_param.poly_degree
                ti=wpg_param.tpassage_ghost-wpg_param.tpassage_ghost(i);
                switch (wpg_param.poly_degree)
                    case 3
                        v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
                        v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
                            (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
                            -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
                            1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
                        v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
                         -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
                            (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
                            -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
                        v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
                        v=[v1' v2' v3' v4'];
                    case 5
                            v=zeros(wpg_param.poly_degree+1,wpg_param.poly_degree+1);
                            for j=1:wpg_param.poly_degree+1
                                fh=str2func(['v' num2str(j) '_n' num2str(wpg_param.poly_degree) '_f']);
                                v(:,j)=fh(ti,i)';
                            end
                    case 6
                        for j=1:wpg_param.poly_degree+1
                            run(['n6_v' num2str(j)]);
                        end
                        v=[v1' v2' v3' v4' v5' v6' v7'];
                    otherwise
                        error('need wpg_param.poly_degree as 5 or 3 or 6, other are not define')
                        
                end
                A_Ca((i-wpg_param.poly_degree-1)*(wpg_param.poly_degree+1)+1:(i-wpg_param.poly_degree-1)*(wpg_param.poly_degree+1)+1+wpg_param.poly_degree,i-wpg_param.poly_degree:i)=v;
            end
            
            A_Ca=[A_Ca zeros(size(A_Ca,1),4)];
            B_Ca=zeros((wpg_param.poly_degree+1)*wpg_param.nbphases,1);
        end
        %% intern methods of compute_COM
        function M_cs=compute_com_M_cs(wpg_param)
            %Compute the matrix M_cs
            M_cs=zeros(wpg_param.nbpointdiscret,wpg_param.nbphases*2);
            nb=0;
            for j=1:wpg_param.nbphases
                tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-wpg_param.tpassage(j);
                M_cs(nb+1:nb+size(tstep,2),(j-1)*2+1:(j-1)*2+2)=[cosh(wpg_param.w*tstep') sinh(wpg_param.w*tstep')];
                nb=nb+size(tstep,2);
            end

            M_cs(wpg_param.nbpointdiscret,(j-1)*2+1:(j-1)*2+2)=[cosh(wpg_param.w*(wpg_param.tpassage(j+1)-wpg_param.tpassage(j))) sinh(wpg_param.w*(wpg_param.tpassage(j+1)-wpg_param.tpassage(j)))];
        end
        function M_cs_d=compute_com_M_cs_d(wpg_param)
            %Compute the matrix M_cs
            M_cs_d=zeros(wpg_param.nbpointdiscret,wpg_param.nbphases*2);
            nb=0;
            for j=1:wpg_param.nbphases
                tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-wpg_param.tpassage(j);
                M_cs_d(nb+1:nb+size(tstep,2),(j-1)*2+1:(j-1)*2+2)=[wpg_param.w*sinh(wpg_param.w*tstep') wpg_param.w*cosh(wpg_param.w*tstep')];
                nb=nb+size(tstep,2);
            end

            M_cs_d(wpg_param.nbpointdiscret,(j-1)*2+1:(j-1)*2+2)=[wpg_param.w*sinh(wpg_param.w*(wpg_param.tpassage(j+1)-wpg_param.tpassage(j))) wpg_param.w*cosh(wpg_param.w*(wpg_param.tpassage(j+1)-wpg_param.tpassage(j)))];
        end
        function A_A=compute_com_A_A(wpg_param,nbphases,w)
        %Compute the matrix A_A use to compute the A coefficient from zmp
        %polynomial coefficients
            switch (wpg_param.poly_degree)
                case 5
                    Aj=[1 0 2/w^2 0 24/w^4 0;
                        0 1 0 6/w^2 0 120/w^4;
                        0 0 1 0 12/w^2 0;
                        0 0 0 1 0 20/w^2;
                        0 0 0 0 1 0;
                        0 0 0 0 0 1];
                    tmp = repmat({Aj},nbphases,1);
                    A_A=blkdiag(tmp{:});
                case 3
                    Aj=[1 0 2/w^2 0 ;
                        0 1 0 6/w^2 ;
                        0 0 1 0;
                        0 0 0 1];
                    tmp = repmat({Aj},nbphases,1);
                    A_A=blkdiag(tmp{:});
                case 6
                    Aj=[1 0 2/w^2 0 24/w^4 0 720/w^6;
                        0 1 0 6/w^2 0 120/w^4 0;
                        0 0 1 0 12/w^2 0 360/w^4;
                        0 0 0 1 0 20/w^2 0;
                        0 0 0 0 1 0 30/w^2;
                        0 0 0 0 0 1 0;
                        0 0 0 0 0 0 1];
                    tmp = repmat({Aj},nbphases,1);
                    A_A=blkdiag(tmp{:});
                otherwise
                    error('need wpg_param.poly_degree as 5 or 3, other are not define')
            end
        end
        function [A_y B_y]=compute_com_VM(A_Ca,B_Ca,A_A,tpassage,w,poly_degree)
            %put a system G*y=N*x0+H*l as y=G^-1*N*x+G^-1*H*l
            %x0 are the com initial and final positions and optimization
            %variables
            %l=A_l*X+B_l
            %A_l=[A_A.A_Ca 0]
            %B_l=[A_A.B_Ca]
            
            G_=wpg_qp_problem_Adrien.compute_com_G(tpassage,w);
            
            H_=wpg_qp_problem_Adrien.compute_com_H(tpassage,poly_degree);
            
            N=sparse([1 0 0 0; zeros(size(H_,1)-2,4); 0 1 0 0]); %pcominit pcomfin scominit scomfin
            
            A_l=sparse(A_A)*sparse(A_Ca);
            B_l=A_A*B_Ca;
            
            A_y=[zeros(size(G_,1),size(A_Ca,2)-4) G_\(N)]+G_\(H_*A_l);
            B_y=G_\(H_*B_l);
        end
        function G=compute_com_G(tpassage,w)
        %Compute the G matrix
            nbphases=length(tpassage)-1;

%             G=[1 0 zeros(1,2*nbphases-2)];
% 
%             for i=1:nbphases-1
%                 Dt=tpassage(i+1)-tpassage(i);
%                 G=[G;zeros(1,2*(i-1)) cosh(w*Dt) sinh(w*Dt) -1 0 zeros(1,2*(nbphases-i)-2);zeros(1,2*(i-1)) w*sinh(w*Dt) w*cosh(w*Dt) 0 -w zeros(1,2*(nbphases-i)-2)];
%             end
% 
%             Dt=tpassage(nbphases+1)-tpassage(nbphases);
%             G=[G;zeros(1,2*nbphases-2) cosh(w*Dt) sinh(w*Dt)];
            
            tmp_Dt=[];
            for i=1:nbphases-1
                Dt=tpassage(i+1)-tpassage(i);
                tmp_Dt=[tmp_Dt {[cosh(w*Dt) sinh(w*Dt);w*sinh(w*Dt) w*cosh(w*Dt)]}];
            end
            tmp_cons = repmat({[-1 0;0 -w]},nbphases-1,1);
            M_Dt=blkdiag(tmp_Dt{:});
            M_cons=blkdiag(tmp_cons{:});
            G=[M_Dt zeros(size(M_Dt,1),2)]+[zeros(size(M_cons,1),2) M_cons];
            G=[1 0 zeros(1,2*nbphases-2);G];
            Dt=tpassage(nbphases+1)-tpassage(nbphases);
            G=[G;zeros(1,2*nbphases-2) cosh(w*Dt) sinh(w*Dt)];
        end
        function H=compute_com_H(tpassage,poly_degree)
        %compute the H matrix
            nbphases=length(tpassage)-1;
            switch (poly_degree)
                case 5
                    tmp_Dt=[];
                    for i=1:nbphases-1
                        Dt=tpassage(i+1)-tpassage(i);
                        tmp_Dt=[tmp_Dt {[-Dt^0 -Dt^1 -Dt^2 -Dt^3 -Dt^4 -Dt^5;0 -1*Dt^0 -2*Dt^1 -3*Dt^2 -4*Dt^3 -5*Dt^4]}];
                    end
                    tmp_cons = repmat({[0 0 0 0 +1 0;0 0 0 0 0 +1]},nbphases-1,1);
                    M_Dt=blkdiag(tmp_Dt{:});
                    M_cons=blkdiag(tmp_cons{:});
                    H=[M_Dt zeros(size(M_Dt,1),poly_degree+1)]+[zeros(size(M_cons,1),2) M_cons zeros(size(M_cons,1),poly_degree+1-2)];
                    H=[-1 0 zeros(1,6*nbphases-2);H];
                    Dt=tpassage(nbphases+1)-tpassage(nbphases);
                    H=[H;zeros(1,6*(nbphases-1)) -Dt^0 -Dt^1 -Dt^2 -Dt^3 -Dt^4 -Dt^5];
                case 3
                    H=[-1 0 zeros(1,4*nbphases-2)];
                    for i=1:nbphases-1
                        Dt=tpassage(i+1)-tpassage(i);
                        H=[H;zeros(1,4*(i-1)) -Dt^0 -Dt^1 -Dt^2 -Dt^3 +1 0 zeros(1,4*(nbphases-i)-2);zeros(1,4*(i-1)) 0 -1*Dt^0 -2*Dt^1 -3*Dt^2 0 +1 zeros(1,4*(nbphases-i)-2)];
                    end
                    Dt=tpassage(nbphases+1)-tpassage(nbphases);
                    H=[H;zeros(1,4*(nbphases-1)) -Dt^0 -Dt^1 -Dt^2 -Dt^3];
                case 6
                    H=[-1 0 zeros(1,(poly_degree+1)*nbphases-2)];
                    for i=1:nbphases-1
                        dt=tpassage(i+1)-tpassage(i);
                        Dt=[];
                        Dt_d=[];
                        for n=0:poly_degree
                            Dt=[Dt dt'.^(n)];
                            Dt_d=[Dt_d n*dt'.^max([0 n-1])];
                        end
                        H=[H;zeros(1,(poly_degree+1)*(i-1)) -Dt +1 0 zeros(1,(poly_degree+1)*(nbphases-i)-2);zeros(1,(poly_degree+1)*(i-1)) -Dt_d 0 +1 zeros(1,(poly_degree+1)*(nbphases-i)-2)];
                    end
                    Dt=tpassage(nbphases+1)-tpassage(nbphases);
                    H=[H;zeros(1,(poly_degree+1)*(nbphases-1)) -Dt^0 -Dt^1 -Dt^2 -Dt^3 -Dt^4 -Dt^5 -Dt^6];
                otherwise
                    error('need wpg_param.poly_degree as 5 or 3, other are not define')
            end
        end
        function [A_com B_com]=compute_com(A_Ca,B_Ca,A_y,B_y,A_A,M_t,M_cs)
        %%%generate the COM trajectory gradient%%%
        %Use Morisawa algorithm with 5th zmp polynomial coeff known
        
            A_com=sparse(M_cs)*sparse(A_y)+sparse(M_t)*sparse(A_A)*sparse(A_Ca);
            B_com=M_cs*B_y+sparse(M_t)*sparse(A_A)*sparse(B_Ca);
        end
        function [A_cons_scom B_cons_scom] = cons_scom_initfin(A_Ca,B_Ca,A_y,B_y,A_A,tpassage,nbphases,w,poly_degree)
            dti=0;%initial time of the first phase
            dtf=tpassage(end)-tpassage(end-1);%fial time step of the last phase
            M_cs=[w*sinh(w*dti) w*cosh(w*dti) zeros(1,(nbphases-1)*2);
                zeros(1,(nbphases-1)*2) w*sinh(w*dtf) w*cosh(w*dtf)];
            switch (poly_degree)
                case 5
                    M_t=[0 1 2*dti 3*dti^2 4*dti^3 5*dti^4 zeros(1,(nbphases-1)*(poly_degree+1));
                        zeros(1,(nbphases-1)*(poly_degree+1)) 0 1 2*dtf 3*dtf^2 4*dtf^3 5*dtf^4];
                case 3
                    M_t=[0 1 2*dti 3*dti^2 zeros(1,(nbphases-1)*(poly_degree+1));
                        zeros(1,(nbphases-1)*(poly_degree+1)) 0 1 2*dtf 3*dtf^2];
                case 6
                    M_t=[0 1 2*dti 3*dti^2 4*dti^3 5*dti^4 6*dti^5 zeros(1,(nbphases-1)*(poly_degree+1));
                        zeros(1,(nbphases-1)*(poly_degree+1)) 0 1 2*dtf 3*dtf^2 4*dtf^3 5*dtf^4 6*dtf^5];
                otherwise
                    error('need wpg_param.poly_degree as 5 or 3 or 6, other are not define')
            end

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
        function [A_an B_an] = compute_ankle_positions_SSP(wpg_param,pankinit,pankfin)            
            nbphases=length(wpg_param.discretization);
            A_an=zeros(1+sum(wpg_param.discretization),wpg_param.nbpankle);
            nb=0;
            rowinitphase=sum(wpg_param.discretization(1:wpg_param.nbpolypi));
            n=2;
            for j=wpg_param.nbpolypi+1:nbphases-wpg_param.nbpolypf
                rowinitphase=rowinitphase+nb;
                nb=wpg_param.discretization(j);
                if wpg_param.type_phase(j)==0
                else
                    for i=1:nb
                        A_an(rowinitphase+i,n)=1;
                    end
                end
%                 if wpg_param.type_phase(max([j-1 1]))~=wpg_param.type_phase(j)&&wpg_param.type_phase(j)==0
                if  wpg_param.type_phase(j)==6                  
                    n=n+1;
                end
            end

            A_an(end,:)=A_an(end-1,:);

            B_an=zeros(wpg_param.nbpointdiscret,1);
%             B_an(sum(wpg_param.discretization(1:wpg_param.nbpolypi))+1:sum(wpg_param.discretization(1:wpg_param.nbpolypi+wpg_param.nbpolyssp)))=pankinit*ones(sum(wpg_param.discretization(wpg_param.nbpolypi+1:wpg_param.nbpolypi+wpg_param.nbpolyssp)),1);
%             B_an(sum(wpg_param.discretization(1:end-wpg_param.nbpolypf-wpg_param.nbpolyssp))+1:sum(wpg_param.discretization(1:end-wpg_param.nbpolypf)))=pankfin*ones(sum(wpg_param.discretization(end-wpg_param.nbpolyssp-wpg_param.nbpolypf+1:end-wpg_param.nbpolypf)),1);

        end
        %% intern methods of compute_cons_zmp
        function [Acons Bcons] = cons_xy_zmp_stability_SSP(xApzmp,xBpzmp,yApzmp,yBpzmp,xApankle,xBpankle,yApankle,yBpankle,discretization,backtoankle,fronttoankle,exttoankle,inttoankle,sole_margin,theta,type_phase,nbcontrolpointzmp1,firstSS)
            xApzmp_=[xApzmp zeros(size(xApankle)) zeros(size(xApzmp,1),nbcontrolpointzmp1)];
            yApzmp_=[yApzmp zeros(size(yApankle)) zeros(size(yApzmp,1),nbcontrolpointzmp1)];
            xApankle=[zeros(size(xApzmp)) xApankle zeros(size(xApzmp,1),nbcontrolpointzmp1) zeros(size(yApzmp_))];
            yApankle=[zeros(size(xApzmp_)) zeros(size(yApzmp)) yApankle zeros(size(yApzmp,1),nbcontrolpointzmp1)];

            xApzmp_=[xApzmp_ zeros(size(xApzmp_))];
            yApzmp_=[zeros(size(yApzmp_)) yApzmp_];

            Acons1=(zeros(size(xApankle)));
            Bcons1=(zeros(size(xBpankle)));
            Acons2=(zeros(size(xApankle)));
            Bcons2=(zeros(size(xBpankle)));
            Acons3=(zeros(size(xApankle)));
            Bcons3=(zeros(size(xBpankle)));
            Acons4=(zeros(size(xApankle)));
            Bcons4=(zeros(size(xBpankle)));

            j=1;
            if type_phase(1)==0
                footRorL=-firstSS;
            else
                footRorL=firstSS;
            end
            %constraint direction inverse clock-wise
            t=theta(j);
            dx=[cos(t) -sin(t) -cos(t) sin(t)];
            dy=[sin(t) cos(t) -sin(t) -cos(t)];


            xA=(xApzmp_-xApankle);
            yA=(yApzmp_-yApankle);
            xB=(xBpzmp-xBpankle);
            yB=(yBpzmp-yBpankle);

            for i=1:sum(discretization)
                if(i>sum(discretization(1:j)))
                        j=j+1;
                        %constraint direction inverse clock-wise
                        t=theta(j);
                        dx=[cos(t) -sin(t) -cos(t) sin(t)];
                        dy=[sin(t) cos(t) -sin(t) -cos(t)];
                        if type_phase(j)==1
                            footRorL=-footRorL;
                        end
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
                    
                    if footRorL==-1
                        Bcons2(i)=Bcons2(i)+exttoankle-sole_margin;
                        Bcons4(i)=Bcons4(i)+inttoankle-sole_margin;
                    elseif footRorL==+1
                        Bcons2(i)=Bcons2(i)+inttoankle-sole_margin;
                        Bcons4(i)=Bcons4(i)+exttoankle-sole_margin;
                    end
                end
            end
            
            if type_phase(end)~=0
                i=i+1;
                Acons1(i,:)=dx(1)*(xA(i,:))+dy(1)*(yA(i,:));
                Bcons1(i)=-dx(1)*(xB(i,:))-dy(1)*(yB(i,:))+fronttoankle-sole_margin;
                Acons2(i,:)=dx(2)*(xA(i,:))+dy(2)*(yA(i,:));
                Bcons2(i)=-dx(2)*(xB(i,:))-dy(2)*(yB(i,:));
                Acons3(i,:)=dx(3)*(xA(i,:))+dy(3)*(yA(i,:));
                Bcons3(i)=-dx(3)*(xB(i,:))-dy(3)*(yB(i,:))+backtoankle-sole_margin;
                Acons4(i,:)=dx(4)*(xA(i,:))+dy(4)*(yA(i,:));
                Bcons4(i)=-dx(4)*(xB(i,:))-dy(4)*(yB(i,:));

                if footRorL==-1
                    Bcons2(i)=Bcons2(i)+exttoankle-sole_margin;
                    Bcons4(i)=Bcons4(i)+inttoankle-sole_margin;
                elseif footRorL==+1
                    Bcons2(i)=Bcons2(i)+inttoankle-sole_margin;
                    Bcons4(i)=Bcons4(i)+exttoankle-sole_margin;
                end
            end

            Acons=[sparse(Acons1(any(any(Acons1,2)+any(Bcons1,2),2),:));
                   sparse(Acons2(any(any(Acons2,2)+any(Bcons2,2),2),:));
                   sparse(Acons3(any(any(Acons3,2)+any(Bcons3,2),2),:));
                   sparse(Acons4(any(any(Acons4,2)+any(Bcons4,2),2),:))];
            Bcons=[Bcons1(any(any(Acons1,2)+any(Bcons1,2),2),:);
                   Bcons2(any(any(Acons2,2)+any(Bcons2,2),2),:);
                   Bcons3(any(any(Acons3,2)+any(Bcons3,2),2),:);
                   Bcons4(any(any(Acons4,2)+any(Bcons4,2),2),:)];
        end
        %% intern methods of compute_ZMP1
        function [A_Ca1 B_Ca1]=compute_zmp1_Ca1(wpg_param,tpassage,psa_zmp1init,psa_zmp1fin,cutting)
            nbpoly=length(tpassage)-1;
%             nbpoly1=wpg_param.nbstep*wpg_param.nbpolyzmp1+wpg_param.nbpolypi+wpg_param.nbpolypf;
%             nbcontrolpointzmp1=wpg_param.nbstep*(wpg_param.poly_degree+wpg_param.nbpolydsp)+wpg_param.nbpolypi+wpg_param.nbpolypf+2*wpg_param.poly_degree;
            nbcontrolpointzmp1=wpg_param.nbcontrolpointzmp1;
            g=zeros((wpg_param.poly_degree+1)*nbpoly,nbcontrolpointzmp1);
            
            %% zmp1
            n=0;
            m=0;
            for i=1+wpg_param.poly_degree:nbpoly+wpg_param.poly_degree
                ti=wpg_param.tpassage_ghost-wpg_param.tpassage_ghost(i);
                if wpg_param.type_phase(i-wpg_param.poly_degree)==0
                    switch (wpg_param.poly_degree)
                        case 3
                            v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
                            v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
                                (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
                                -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
                                1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
                            v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
                             -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
                                (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
                                -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
                            v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
                            v=[v1' v2' v3' v4'];
                        case 5
                            v=zeros(wpg_param.poly_degree+1,wpg_param.poly_degree+1);
                            for j=1:wpg_param.poly_degree+1
                                fh=str2func(['v' num2str(j) '_n' num2str(wpg_param.poly_degree) '_f']);
                                v(:,j)=fh(ti,i)';
                            end
                        case 6
                            for j=1:wpg_param.poly_degree+1
                                run(['n6_v' num2str(j)]);
                            end
                            v=[v1' v2' v3' v4' v5' v6' v7'];

                        otherwise
                            error('need wpg_param.poly_degree as 5 or 3, other are not define')
                    end
                
%                     m1=zeros(wpg_param.poly_degree+1,i-wpg_param.poly_degree-1+n);
%                     m2=v(:,1:end-1-n);
%                     m3=zeros(wpg_param.poly_degree+1,nbpoly+wpg_param.poly_degree+wpg_param.nbparamank-size(m1,2)-size(m2,2));
%                     m4=zeros(wpg_param.poly_degree+1,wpg_param.poly_degree+1+sum(any(wpg_param.type_phase(1:i-wpg_param.poly_degree)==0,1))-1-n);
%                     m=v(:,end-n:end);
%                     m5=zeros(wpg_param.poly_degree+1,wpg_param.nbparamtotal-(nbpoly+wpg_param.poly_degree+wpg_param.nbparamank)-size(m4,2)-size(m,2));
%                     g((i-wpg_param.poly_degree-1)*(wpg_param.poly_degree+1)+1:(i-wpg_param.poly_degree-1)*(wpg_param.poly_degree+1)+1+wpg_param.poly_degree,+wpg_param.nbparamtotal-nbpoly1+wpg_param.nbpolypi+wpg_param.poly_degree)=v(:,1:end-1);
%                     g((i-wpg_param.poly_degree-1)*(wpg_param.poly_degree+1)+1:(i-wpg_param.poly_degree-1)*(wpg_param.poly_degree+1)+1+wpg_param.poly_degree,:)=[m1 m2 m3 m4 m m5];
%                         m1=zeros(wpg_param.poly_degree+1,wpg_param.nbparamtotal-nbpoly1*(wpg_param.poly_degree+1));
                        m1=[];
                        m2=[];
                        m3=[];
                        if sum(i==1+wpg_param.poly_degree:wpg_param.poly_degree+wpg_param.nbpolypi)
                            m4=zeros(wpg_param.poly_degree+1,n);
                            m5=zeros(wpg_param.poly_degree+1,nbcontrolpointzmp1-(wpg_param.poly_degree+1)-n);
                        elseif sum(i==nbpoly+wpg_param.poly_degree-wpg_param.nbpolypf:nbpoly+wpg_param.poly_degree)
                            m4=zeros(wpg_param.poly_degree+1,nbcontrolpointzmp1-(wpg_param.poly_degree+wpg_param.nbpolypf)+n);
                            m5=zeros(wpg_param.poly_degree+1,wpg_param.nbpolypf-1-n);
                        else
                            m4=zeros(wpg_param.poly_degree+1,(wpg_param.poly_degree+wpg_param.nbpolypi)*(wpg_param.nbpolypi~=0)+(ceil((i-(wpg_param.poly_degree+wpg_param.nbpolypi))/(wpg_param.nbpolyssp+wpg_param.nbpolydsp))-1)*(wpg_param.poly_degree+wpg_param.nbpolydsp)+n);
                            m5=zeros(wpg_param.poly_degree+1,nbcontrolpointzmp1-(wpg_param.poly_degree+wpg_param.nbpolypi)*(wpg_param.nbpolypi~=0)-(ceil((i-(wpg_param.poly_degree+wpg_param.nbpolypi))/(wpg_param.nbpolyssp+wpg_param.nbpolydsp))-1)*(wpg_param.poly_degree+wpg_param.nbpolydsp)-(wpg_param.poly_degree+1+n));
                        end
                        m=v;
%                         m5=zeros(wpg_param.poly_degree+1,(wpg_param.poly_degree+1)*(nbpoly1-ceil((i-wpg_param.poly_degree)/3)));
%                         i
%                         [size(m4,2) size(m,2) size(m5,2)]
                        g((i-wpg_param.poly_degree-1)*(wpg_param.poly_degree+1)+1:(i-wpg_param.poly_degree-1)*(wpg_param.poly_degree+1)+1+wpg_param.poly_degree,:)=[m1 m2 m3 m4 m m5];
                        
                    n=n+1;
                else
                    n=0;
                end
            end
            
%             A_Ca1=[g(:,1:nbpoly+wpg_param.poly_degree+wpg_param.nbparamank) zeros(size(g,1),2) g(:,nbpoly+wpg_param.poly_degree+wpg_param.nbparamank+1:end)];
            A_Ca1=[zeros(size(g,1),nbpoly+wpg_param.poly_degree+4) zeros(size(g,1),wpg_param.nbpankle) g];
            B_Ca1=zeros(size(g,1),1);
        end
        %% intern methods of compute_cost_torque_zmp1
        function [Ac]=compute_force_repartition(wpg_param,cutting)
            nbpoly=wpg_param.nbphases;
            g=[];

            %%%cutting fcom1 init in two%%%
            dt=wpg_param.tpassage(2)-wpg_param.tpassage(1);
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
            dt=wpg_param.tpassage(2)-wpg_param.tpassage(1);
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
                dt=wpg_param.tpassage(i+1)-wpg_param.tpassage(i);
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
            dt=wpg_param.tpassage(end)-wpg_param.tpassage(end-1);
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
            dt=wpg_param.tpassage(end)-wpg_param.tpassage(end-1);
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
%             A3=[zeros(size(A_fcom)) -Afb*mg*A_an zeros(size(A_fcom,1),size(A_zmp1,2)-size(A_fcom,2)-size(A_an,2))];
            A3=-Afb*mg*[zeros(size(A_fcom)) A_an zeros(size(A_fcom,1),size(A_zmp1,2)-size(A_fcom,2)-size(A_an,2))];
            A_t_zmp1=A1+A2+A3;
            B_t_zmp1=ha*Afb*(B_fcom)+Afb*mg*B_zmp1-Afb*mg*B_an;

%             A_t_zmp1=A_t_zmp1(any(any(A_an,2)+any(B_an,2),2),:);
%             B_t_zmp1=B_t_zmp1(any(any(A_an,2)+any(B_an,2),2),:);
        end
        function [A_an B_an] = compute_ankle_positions_zmp1(wpg_param,pankinit1,pankinit2,pankfin1)            
            nbphases=length(wpg_param.discretization);
            A_an=zeros(wpg_param.nbpointdiscret,wpg_param.nbpankle);
            nb=0;
            rowinitphase=0;
            n=1;
            for j=1:nbphases
                rowinitphase=rowinitphase+nb;
                nb=wpg_param.discretization(j);
                if wpg_param.type_phase(j)==0
                    for i=1:nb
                        A_an(rowinitphase+i,n)=1;
                    end
                end
                if wpg_param.type_phase(j)==1
                    n=n+1;
                end
            end

            A_an(end,:)=A_an(end-1,:);

            B_an=zeros(wpg_param.nbpointdiscret,1);

        end
        %% intern methods of compute_ZMP2
        function [A_zmp2 B_zmp2]=compute_zmp2(A_zmp,B_zmp,A_zmp1,B_zmp1,k_diag,mg)
            xf2=(eye(length(k_diag))-k_diag)*mg;
            A_zmp2=-mg*(xf2\(A_zmp1-[A_zmp zeros(size(A_zmp,1),size(A_zmp1,2)-size(A_zmp,2))]))+A_zmp1;
            B_zmp2=-mg*(xf2\(B_zmp1-B_zmp))+B_zmp1;
        end
        %% intern methods of compute_cost_torque_zmp2
        function [A_an B_an] = compute_ankle_positions_zmp2(wpg_param,pankinit2,pankfin1,pankfin2)
            nbphases=length(wpg_param.discretization);
            A_an=zeros(wpg_param.nbpointdiscret,wpg_param.nbpankle);
            nb=0;
            rowinitphase=0;
            n=2;
            for j=1:nbphases
                rowinitphase=rowinitphase+nb;
                nb=wpg_param.discretization(j);
                if wpg_param.type_phase(j)==0
                    for i=1:nb
                        A_an(rowinitphase+i,n)=1;
                    end
                end
                if wpg_param.type_phase(j)==1
                    n=n+1;
                end
            end

            A_an(end,:)=A_an(end-1,:);

            B_an=zeros(wpg_param.nbpointdiscret,1);
        end
        function [A_t_zmp2 B_t_zmp2] = compute_torque_zmp2(A_zmp2,B_zmp2,A_fcom,B_fcom,A_an,B_an,mg,ha,Afb)
        %le couple en x est fonction des y
        %le couple en y est fonction des x
            IAfb=eye(size(Afb,1))-Afb;

            A1=ha*IAfb*([A_fcom zeros(size(A_fcom,1),size(A_zmp2,2)-size(A_fcom,2))]);
            A2=+IAfb*mg*A_zmp2;
%             A3=[zeros(size(A_fcom)) -Afb*mg*A_an zeros(size(A_fcom,1),size(A_zmp2,2)-size(A_fcom,2)-size(A_an,2))];
            A3=-IAfb*mg*[zeros(size(A_fcom)) A_an zeros(size(A_fcom,1),size(A_zmp2,2)-size(A_fcom,2)-size(A_an,2))];
            A_t_zmp2=A1+A2+A3;

            B_t_zmp2=ha*IAfb*(B_fcom)+IAfb*mg*B_zmp2-IAfb*mg*B_an;
        end
        %% intern methods of compute_cost_acc_zmp2
        function [A_zmp2_acc B_zmp2_acc]=compute_zmp2_acc(A_zmp,B_zmp,Apzmp1,B_zmp1,A_zmp_s,B_zmp_s,A_zmp1_s,B_zmp1_s,A_zmp_a,B_zmp_a,A_zmp1_a,B_zmp1_a,k_diag,k_diag_d,k_diag_dd)
        %Compute the gradient matrix of zmp2 acceleration
            Afb_=(eye(length(k_diag))-k_diag);
            dAfb_=-k_diag_d;
            ddAfb_=-k_diag_dd;

%             A1=-(Afb_\(A_zmp1_a-[A_zmp_a zeros(size(A_zmp_a,1),size(A_zmp1_a,2)-size(A_zmp_a,2))]));
%             A2=+2*dAfb_*((Afb_)\((Afb_)\(A_zmp1_s-[A_zmp_s zeros(size(A_zmp_s,1),size(A_zmp1_s,2)-size(A_zmp_s,2))])));
%             A3=+ddAfb_*((Afb_)\((Afb_)\(Apzmp1-[A_zmp zeros(size(A_zmp,1),size(Apzmp1,2)-size(A_zmp,2))])))-2*dAfb_*dAfb_*((Afb_)\((Afb_)\((Afb_)\(Apzmp1-[A_zmp zeros(size(A_zmp,1),size(Apzmp1,2)-size(A_zmp,2))]))));

            v_Afb_=diag(Afb_);
            v_ddAfb_=diag(ddAfb_);
            v_dAfb_=diag(dAfb_);
%             bsxfun(@times,v_Afb_.^2,)
            A1=-bsxfun(@times,v_Afb_.^-1,(A_zmp1_a-[A_zmp_a zeros(size(A_zmp_a,1),size(A_zmp1_a,2)-size(A_zmp_a,2))]));
            A2=+2*bsxfun(@times,v_dAfb_,bsxfun(@times,v_Afb_.^-2,(A_zmp1_s-[A_zmp_s zeros(size(A_zmp_s,1),size(A_zmp1_s,2)-size(A_zmp_s,2))])));
            A3=bsxfun(@times,v_ddAfb_,(bsxfun(@times,v_Afb_.^-2,(Apzmp1-[A_zmp zeros(size(A_zmp,1),size(Apzmp1,2)-size(A_zmp,2))]))))-2*bsxfun(@times,v_dAfb_.^2,(bsxfun(@times,v_Afb_.^-3,(Apzmp1-[A_zmp zeros(size(A_zmp,1),size(Apzmp1,2)-size(A_zmp,2))]))));
            
            B1=-(Afb_\(B_zmp1_a-B_zmp_a));
            B2=+2*dAfb_*((Afb_)\((Afb_)\(B_zmp1_s-B_zmp_s)));
%             B3=+ddAfb_*((Afb_)\((Afb_)\(B_zmp1-B_zmp)))-2*dAfb_*dAfb_*((Afb_)\((Afb_)\((Afb_)\(B_zmp1-B_zmp))));
            B3=bsxfun(@times,v_ddAfb_,(bsxfun(@times,v_Afb_.^-2,(B_zmp1-B_zmp))))-2*bsxfun(@times,v_dAfb_.^2,(bsxfun(@times,v_Afb_.^-3,(B_zmp1-B_zmp))));


            A_zmp2_acc=A1+A2+A3+A_zmp1_a;
            B_zmp2_acc=B1+B2+B3+B_zmp1_a;
        end
        %% intern methods of compute_cons_zmp12
        function [Acons Bcons] = cons_xy_zmp1_stability(xApzmp,xBpzmp,yApzmp,yBpzmp,xApankle,xBpankle,yApankle,yBpankle,discretization,backtoankle,fronttoankle,exttoankle,inttoankle,sole_margin,theta,type_phase,nbcontrolpointzmp,nbcontrolpointzmp1,firstSS)
            xApzmp_=xApzmp;
            yApzmp_=yApzmp;
            xApankle_=[zeros(size(xApzmp_,1),nbcontrolpointzmp+4) xApankle zeros(size(xApzmp,1),nbcontrolpointzmp1) zeros(size(yApzmp_))];
            yApankle_=[zeros(size(xApzmp_)) zeros(size(yApzmp_,1),nbcontrolpointzmp+4) yApankle zeros(size(yApzmp,1),nbcontrolpointzmp1)];

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
%             theta_=theta(1:end-1);
            theta_=theta(any(type_phase==0,1));

            t=theta_(j);
            dx=[cos(t) -sin(t) -cos(t) sin(t)];
            dy=[sin(t) cos(t) -sin(t) -cos(t)];

            xA=xApzmp_-xApankle_;
            yA=yApzmp_-yApankle_;
            xB=xBpzmp-xBpankle;
            yB=yBpzmp-yBpankle;

            discretization_=discretization(any(type_phase==0,1));
            type_phase_=[0 type_phase(find(any(type_phase(1:end-1)==0,1))+1)];
            
            if type_phase(1)==0
                footRorL=-firstSS;
            else
                footRorL=firstSS;
                type_phase_(end)=[];
            end

            for i=1:sum(discretization_)
                if(i>sum(discretization_(1:j)))
                        j=j+1;
                        %constraint direction inverse clock-wise
                        t=theta_(j);
                        dx=[cos(t) -sin(t) -cos(t) sin(t)];
                        dy=[sin(t) cos(t) -sin(t) -cos(t)];
                        if type_phase_(j)==1
                            footRorL=-footRorL;
                        end
                end

                    Acons1(i,:)=dx(1)*(xA(i,:))+dy(1)*(yA(i,:));
                    Bcons1(i)=-dx(1)*(xB(i,:))-dy(1)*(yB(i,:))+fronttoankle-sole_margin;
                    Acons2(i,:)=dx(2)*(xA(i,:))+dy(2)*(yA(i,:));
                    Bcons2(i)=-dx(2)*(xB(i,:))-dy(2)*(yB(i,:));
                    Acons3(i,:)=dx(3)*(xA(i,:))+dy(3)*(yA(i,:));
                    Bcons3(i)=-dx(3)*(xB(i,:))-dy(3)*(yB(i,:))+backtoankle-sole_margin;
                    Acons4(i,:)=dx(4)*(xA(i,:))+dy(4)*(yA(i,:));
                    Bcons4(i)=-dx(4)*(xB(i,:))-dy(4)*(yB(i,:));

                    if footRorL==-1
                        Bcons2(i)=Bcons2(i)+exttoankle-sole_margin;
                        Bcons4(i)=Bcons4(i)+inttoankle-sole_margin;
                    elseif footRorL==+1
                        Bcons2(i)=Bcons2(i)+inttoankle-sole_margin;
                        Bcons4(i)=Bcons4(i)+exttoankle-sole_margin;
                    end
            end
            
            if type_phase(end)==0 
                i=i+1;
                Acons1(i,:)=dx(1)*(xA(i,:))+dy(1)*(yA(i,:));
                Bcons1(i)=-dx(1)*(xB(i,:))-dy(1)*(yB(i,:))+fronttoankle-sole_margin;
                Acons2(i,:)=dx(2)*(xA(i,:))+dy(2)*(yA(i,:));
                Bcons2(i)=-dx(2)*(xB(i,:))-dy(2)*(yB(i,:));
                Acons3(i,:)=dx(3)*(xA(i,:))+dy(3)*(yA(i,:));
                Bcons3(i)=-dx(3)*(xB(i,:))-dy(3)*(yB(i,:))+backtoankle-sole_margin;
                Acons4(i,:)=dx(4)*(xA(i,:))+dy(4)*(yA(i,:));
                Bcons4(i)=-dx(4)*(xB(i,:))-dy(4)*(yB(i,:));

                if footRorL==-1
                    Bcons2(i)=Bcons2(i)+exttoankle-sole_margin;
                    Bcons4(i)=Bcons4(i)+inttoankle-sole_margin;
                elseif footRorL==+1
                    Bcons2(i)=Bcons2(i)+inttoankle-sole_margin;
                    Bcons4(i)=Bcons4(i)+exttoankle-sole_margin;
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
        function [Acons Bcons] = cons_xy_zmp2_stability(xApzmp,xBpzmp,yApzmp,yBpzmp,xApankle,xBpankle,yApankle,yBpankle,discretization,backtoankle,fronttoankle,exttoankle,inttoankle,sole_margin,theta,type_phase,nbcontrolpointzmp,nbcontrolpointzmp1,firstSS)
            xApzmp_=xApzmp;
            yApzmp_=yApzmp;
            xApankle=[zeros(size(xApzmp_,1),nbcontrolpointzmp+4) xApankle zeros(size(xApzmp,1),nbcontrolpointzmp1) zeros(size(yApzmp_))];
            yApankle=[zeros(size(xApzmp_)) zeros(size(yApzmp_,1),nbcontrolpointzmp+4) yApankle zeros(size(yApzmp,1),nbcontrolpointzmp1)];

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
%             theta_=theta(2:end);
            theta_=theta(any(type_phase==0,1));

            t=theta_(j);
            dx=[cos(t) -sin(t) -cos(t) sin(t)];
            dy=[sin(t) cos(t) -sin(t) -cos(t)];

            xA=xApzmp_-xApankle;
            yA=yApzmp_-yApankle;
            xB=xBpzmp-xBpankle;
            yB=yBpzmp-yBpankle;

            discretization_=discretization(any(type_phase==0,1));
            
            if type_phase(1)==0
                type_phase_=[0 type_phase(find(any(type_phase(1:end-1)==0,1))+1)];
            else
                type_phase_=[1 type_phase(find(any(type_phase(1:end-1)==0,1))+1)];
                type_phase_(end)=[];
            end
            
            if type_phase(end)==0
                discretization_(any(type_phase_,1))=discretization_(any(type_phase_,1))-(sum(discretization_)+1-size(xApzmp_,1))/sum(type_phase_);
            else
                discretization_(any(type_phase_,1))=discretization_(any(type_phase_,1))-(sum(discretization_)-size(xApzmp_,1))/sum(type_phase_);
            end
            
            
            if type_phase(1)==0
                footRorL=+firstSS;
            else
                footRorL=-firstSS;
            end
 
            for i=1:sum(discretization_)
                if(i>sum(discretization_(1:j)))
                        j=j+1;
                        %constraint direction inverse clock-wise
                        t=theta_(j);
                        dx=[cos(t) -sin(t) -cos(t) sin(t)];
                        dy=[sin(t) cos(t) -sin(t) -cos(t)];
                        
                        if type_phase_(j)==1
                            footRorL=-footRorL;
                        end
                end

                    Acons1(i,:)=dx(1)*(xA(i,:))+dy(1)*(yA(i,:));
                    Bcons1(i)=-dx(1)*(xB(i,:))-dy(1)*(yB(i,:))+fronttoankle-sole_margin;
                    Acons2(i,:)=dx(2)*(xA(i,:))+dy(2)*(yA(i,:));
                    Bcons2(i)=-dx(2)*(xB(i,:))-dy(2)*(yB(i,:));
                    Acons3(i,:)=dx(3)*(xA(i,:))+dy(3)*(yA(i,:));
                    Bcons3(i)=-dx(3)*(xB(i,:))-dy(3)*(yB(i,:))+backtoankle-sole_margin;
                    Acons4(i,:)=dx(4)*(xA(i,:))+dy(4)*(yA(i,:));
                    Bcons4(i)=-dx(4)*(xB(i,:))-dy(4)*(yB(i,:));
                    
                    if footRorL==-1
                        Bcons2(i)=Bcons2(i)+exttoankle-sole_margin;
                        Bcons4(i)=Bcons4(i)+inttoankle-sole_margin;
                    elseif footRorL==+1
                        Bcons2(i)=Bcons2(i)+inttoankle-sole_margin;
                        Bcons4(i)=Bcons4(i)+exttoankle-sole_margin;
                    end
            end
            
            if type_phase(end)==0
                i=i+1;
                Acons1(i,:)=dx(1)*(xA(i,:))+dy(1)*(yA(i,:));
                Bcons1(i)=-dx(1)*(xB(i,:))-dy(1)*(yB(i,:))+fronttoankle-sole_margin;
                Acons2(i,:)=dx(2)*(xA(i,:))+dy(2)*(yA(i,:));
                Bcons2(i)=-dx(2)*(xB(i,:))-dy(2)*(yB(i,:));
                Acons3(i,:)=dx(3)*(xA(i,:))+dy(3)*(yA(i,:));
                Bcons3(i)=-dx(3)*(xB(i,:))-dy(3)*(yB(i,:))+backtoankle-sole_margin;
                Acons4(i,:)=dx(4)*(xA(i,:))+dy(4)*(yA(i,:));
                Bcons4(i)=-dx(4)*(xB(i,:))-dy(4)*(yB(i,:));

                if footRorL==-1
                    Bcons2(i)=Bcons2(i)+exttoankle-sole_margin;
                    Bcons4(i)=Bcons4(i)+inttoankle-sole_margin;
                elseif footRorL==+1
                    Bcons2(i)=Bcons2(i)+inttoankle-sole_margin;
                    Bcons4(i)=Bcons4(i)+exttoankle-sole_margin;
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
        function [Acons Bcons] = cons_stretching(wpg_param)
            nbparamtotal=wpg_param.nbcontrolpointzmp+4+wpg_param.nbpankle+wpg_param.nbcontrolpointzmp1;
            xApankle=[
                zeros(wpg_param.nbpankle,wpg_param.nbcontrolpointzmp+4) eye(wpg_param.nbpankle) zeros(wpg_param.nbpankle,wpg_param.nbcontrolpointzmp1+nbparamtotal);
                ];
            xBpankle=[zeros(wpg_param.nbpankle,1)];
            yApankle=[
                zeros(wpg_param.nbpankle,nbparamtotal+wpg_param.nbcontrolpointzmp+4) eye(wpg_param.nbpankle) zeros(wpg_param.nbpankle,wpg_param.nbcontrolpointzmp1);
                ];
            yBpankle=[zeros(wpg_param.nbpankle,1)];

            %clock-wise turn from upper right if foot in x direction
            ABCDleft=[wpg_param.fronttoankle -wpg_param.backtoankle -wpg_param.backtoankle wpg_param.fronttoankle;
                       wpg_param.exttoankle wpg_param.exttoankle -wpg_param.inttoankle -wpg_param.inttoankle];

            ABCDright=[wpg_param.fronttoankle -wpg_param.backtoankle -wpg_param.backtoankle wpg_param.fronttoankle;
                        wpg_param.inttoankle wpg_param.inttoankle -wpg_param.exttoankle -wpg_param.exttoankle];

            Acons=zeros((wpg_param.nbpankle-3)*16,nbparamtotal*2);
            Bcons=zeros((wpg_param.nbpankle-3)*16,1);

            for i=2:wpg_param.nbpankle-2
%                 n=wpg_param.nbpolypi+i*(wpg_param.nbpolyssp+wpg_param.nbpolydsp)+1;
%                 t1=wpg_param.psi_zmp(n-wpg_param.nbpolyssp-wpg_param.nbpolydsp);
%                 t2=wpg_param.psi_zmp(n);
                
                t1=wpg_param.psi_pstep(i);
                t2=wpg_param.psi_pstep(i+1);

                rot=[cos(t2) -sin(t2);sin(t2) cos(t2)];
                dx=[cos(t1) -sin(t1) -cos(t1) sin(t1)];
                dy=[sin(t1) cos(t1) -sin(t1) -cos(t1)];

                if(mod(i-1,2))
                    ABCD_=ABCDright;
                else
                    ABCD_=ABCDleft;
                end
                [Acons((i-2)*16+1:(i-2)*16+4,:)       Bcons((i-2)*16+1:(i-2)*16+4)]         = wpg_qp_problem_Adrien.cons_builder (rot,dx,dy,ABCD_(:,1),xApankle(i,:),xBpankle(i),yApankle(i,:),yBpankle(i),xApankle(i+1,:),xBpankle(i+1),yApankle(i+1,:),yBpankle(i+1),wpg_param.xankmax,wpg_param.xankmin,wpg_param.yankmax,wpg_param.yankmin,i-1);
                [Acons((i-2)*16+1+4:(i-2)*16+4+4,:)   Bcons((i-2)*16+1+4:(i-2)*16+4+4)]     = wpg_qp_problem_Adrien.cons_builder (rot,dx,dy,ABCD_(:,2),xApankle(i,:),xBpankle(i),yApankle(i,:),yBpankle(i),xApankle(i+1,:),xBpankle(i+1),yApankle(i+1,:),yBpankle(i+1),wpg_param.xankmax,wpg_param.xankmin,wpg_param.yankmax,wpg_param.yankmin,i-1);
                [Acons((i-2)*16+1+8:(i-2)*16+4+8,:)   Bcons((i-2)*16+1+8:(i-2)*16+4+8)]     = wpg_qp_problem_Adrien.cons_builder (rot,dx,dy,ABCD_(:,3),xApankle(i,:),xBpankle(i),yApankle(i,:),yBpankle(i),xApankle(i+1,:),xBpankle(i+1),yApankle(i+1,:),yBpankle(i+1),wpg_param.xankmax,wpg_param.xankmin,wpg_param.yankmax,wpg_param.yankmin,i-1);
                [Acons((i-2)*16+1+12:(i-2)*16+4+12,:) Bcons((i-2)*16+1+12:(i-2)*16+4+12)]   = wpg_qp_problem_Adrien.cons_builder (rot,dx,dy,ABCD_(:,4),xApankle(i,:),xBpankle(i),yApankle(i,:),yBpankle(i),xApankle(i+1,:),xBpankle(i+1),yApankle(i+1,:),yBpankle(i+1),wpg_param.xankmax,wpg_param.xankmin,wpg_param.yankmax,wpg_param.yankmin,i-1);
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
        function [A B] = cons_ankle_fixed_path(nbcontrolpointzmp,nbparamank,nbcontrolpointzmp1,step_number_pankle_fixed)
            Ascomeq_path=[];
            Bscomeq_path=[];
            if size(step_number_pankle_fixed,1)~=0
                for i=1:size(step_number_pankle_fixed)
                    [Ascomeq_path Bscomeq_path]=wpg_qp_problem_Adrien.cons_pankle_fixed(Ascomeq_path,Bscomeq_path,nbcontrolpointzmp,nbparamank,nbcontrolpointzmp1,step_number_pankle_fixed(i,1),step_number_pankle_fixed(i,2:3));
                end
            end
            A=Ascomeq_path;
            B=Bscomeq_path;
        end
        function [A B] = cons_pankle_fixed(AscomeqDSP,Bscomeq,nbcontrolpointzmp,nbparamank,nbcontrolpointzmp1,step_number,pankle_fixed)
            nbparamtotal=nbcontrolpointzmp+nbparamank+nbcontrolpointzmp1;

            Ascomeq_path=[AscomeqDSP;zeros(1,nbcontrolpointzmp+4+step_number-1) 1 zeros(1,nbparamtotal+4+nbcontrolpointzmp1+(nbparamank-step_number))];
            Bscomeq_path=[Bscomeq;-pankle_fixed(1)];

            Ascomeq_path=[Ascomeq_path;zeros(1,nbparamtotal+4+nbcontrolpointzmp+4+step_number-1) 1 zeros(1,nbcontrolpointzmp1+(nbparamank-step_number))];
            Bscomeq_path=[Bscomeq_path;-pankle_fixed(2)];

            A=Ascomeq_path;
            B=Bscomeq_path;
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
        function [M_t] = compute_M_t(wpg_param)
            M_t=zeros(wpg_param.nbpointdiscret,wpg_param.nbphases*(wpg_param.poly_degree+1));
            nb=0;
            for j=1:wpg_param.nbphases
                tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-wpg_param.tpassage(j);
                Dt=[];
                for n=0:wpg_param.poly_degree
                    Dt=[Dt tstep'.^max([0 n])];
                end
                M_t(nb+1:nb+size(Dt,1),(j-1)*(wpg_param.poly_degree+1)+1:(j-1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
                nb=nb+size(Dt,1);
%                 [size(Dt,1)-wpg_param.discretization(j)]
            end
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt (wpg_param.tpassage(end)-wpg_param.tpassage(end-1))'.^max([0 n])];
            end
            M_t(wpg_param.nbpointdiscret,(j-1)*(wpg_param.poly_degree+1)+1:(j-1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
        end
        function [M_t_d] = compute_M_t_d(wpg_param)
            %compute the derivative of M_t
            M_t_d=zeros(wpg_param.nbpointdiscret,wpg_param.nbphases*(wpg_param.poly_degree+1));
            nb=0;
            for j=1:wpg_param.nbphases
                tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-wpg_param.tpassage(j);
                Dt=[];
                for n=0:wpg_param.poly_degree
                    Dt=[Dt n*tstep'.^max([0 n-1])];
                end
                M_t_d(nb+1:nb+size(Dt,1),(j-1)*(wpg_param.poly_degree+1)+1:(j-1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
                nb=nb+size(Dt,1);
%                 [size(Dt,1)-wpg_param.discretization(j)]
            end
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt n*(wpg_param.tpassage(end)-wpg_param.tpassage(end-1))'.^max([0 n-1])];
            end
            M_t_d(wpg_param.nbpointdiscret,(j-1)*(wpg_param.poly_degree+1)+1:(j-1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
        end
        function [M_t_dd] = compute_M_t_dd(wpg_param)
            %compute the second derivative of M_t
            M_t_dd=zeros(wpg_param.nbpointdiscret,wpg_param.nbphases*(wpg_param.poly_degree+1));
            nb=0;
            for j=1:wpg_param.nbphases
                tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-wpg_param.tpassage(j);
                Dt=[];
                for n=0:wpg_param.poly_degree
                    Dt=[Dt n*(n-1)*tstep'.^max([0 n-2])];
                end
                M_t_dd(nb+1:nb+size(Dt,1),(j-1)*(wpg_param.poly_degree+1)+1:(j-1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
                nb=nb+size(Dt,1);
%                 [size(Dt,1)-wpg_param.discretization(j)]
            end
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt n*(n-1)*(wpg_param.tpassage(end)-wpg_param.tpassage(end-1))'.^max([0 n-2])];
            end
            M_t_dd(wpg_param.nbpointdiscret,(j-1)*(wpg_param.poly_degree+1)+1:(j-1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
        end
        function [M_t1] = compute_M_t1(wpg_param,cutting)
            M_t1=zeros(wpg_param.nbpointdiscret,(wpg_param.nbphases+2)*(wpg_param.poly_degree+1));

            %%%zmp1 init cut in two%%%
            nb=0;
            j=1;
            tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*cutting-1/wpg_param.frequency]-wpg_param.tpassage(j);
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt tstep'.^max([0 n])];
            end
            M_t1(nb+1:nb+size(Dt,1),(j-1)*(wpg_param.poly_degree+1)+1:(j-1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
            nb=nb+size(Dt,1);
            
            tstep=[(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*cutting:1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*cutting;
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt tstep'.^max([0 n])];
            end
            M_t1(nb+1:nb+size(Dt,1),(j)*(wpg_param.poly_degree+1)+1:(j)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
            nb=nb+size(Dt,1);

            %%%generation of dt%%%
            for j=2:wpg_param.nbphases-1
                tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-wpg_param.tpassage(j);
                Dt=[];
                for n=0:wpg_param.poly_degree
                    Dt=[Dt tstep'.^max([0 n])];
                end
                M_t1(nb+1:nb+size(Dt,1),(j)*(wpg_param.poly_degree+1)+1:(j+1)*(wpg_param.poly_degree+1))=Dt;
                nb=nb+size(Dt,1);
%                 [size(Dt,1)-wpg_param.discretization(j)]
            end

            %%%zmp1 fin cut in two%%%
            j=wpg_param.nbphases;
            tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*(1-cutting)-1/wpg_param.frequency]-wpg_param.tpassage(j);
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt tstep'.^max([0 n])];
            end
            M_t1(nb+1:nb+size(Dt,1),(j)*(wpg_param.poly_degree+1)+1:(j)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
            nb=nb+size(Dt,1);
            
            tstep=[(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*(1-cutting):1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*(1-cutting);
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt tstep'.^max([0 n])];
            end
            M_t1(nb+1:nb+size(Dt,1),(j+1)*(wpg_param.poly_degree+1)+1:(j+1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
            nb=nb+size(Dt,1);
            
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt (wpg_param.tpassage(end)-(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*(1-cutting))'.^max([0 n])];
            end
            M_t1(wpg_param.nbpointdiscret,(j+1)*(wpg_param.poly_degree+1)+1:(j+1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;

        end
        function [M_t1_d] = compute_M_t1_d(wpg_param,cutting)
            M_t1_d=zeros(wpg_param.nbpointdiscret,(wpg_param.nbphases+2)*(wpg_param.poly_degree+1));

            %%%zmp1 init cut in two%%%
            nb=0;
            j=1;
            tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*cutting-1/wpg_param.frequency]-wpg_param.tpassage(j);
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt n*tstep'.^max([0 n-1])];
            end
            M_t1_d(nb+1:nb+size(Dt,1),(j-1)*(wpg_param.poly_degree+1)+1:(j-1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
            nb=nb+size(Dt,1);
            
            tstep=[(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*cutting:1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*cutting;
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt n*tstep'.^max([0 n-1])];
            end
            M_t1_d(nb+1:nb+size(Dt,1),(j)*(wpg_param.poly_degree+1)+1:(j)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
            nb=nb+size(Dt,1);

            %%%generation of dt%%%
            for j=2:wpg_param.nbphases-1
                tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-wpg_param.tpassage(j);
                Dt=[];
                for n=0:wpg_param.poly_degree
                    Dt=[Dt n*tstep'.^max([0 n-1])];
                end
                M_t1_d(nb+1:nb+size(Dt,1),(j)*(wpg_param.poly_degree+1)+1:(j+1)*(wpg_param.poly_degree+1))=Dt;
                nb=nb+size(Dt,1);
%                 [size(Dt,1)-wpg_param.discretization(j)]
            end

            %%%zmp1 fin cut in two%%%
            j=wpg_param.nbphases;
            tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*(1-cutting)-1/wpg_param.frequency]-wpg_param.tpassage(j);
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt n*tstep'.^max([0 n-1])];
            end
            M_t1_d(nb+1:nb+size(Dt,1),(j)*(wpg_param.poly_degree+1)+1:(j)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
            nb=nb+size(Dt,1);
            
            tstep=[(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*(1-cutting):1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*(1-cutting);
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt n*tstep'.^max([0 n-1])];
            end
            M_t1_d(nb+1:nb+size(Dt,1),(j+1)*(wpg_param.poly_degree+1)+1:(j+1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
            nb=nb+size(Dt,1);
            
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt n*(wpg_param.tpassage(end)-(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*(1-cutting))'.^max([0 n-1])];
            end
            M_t1_d(wpg_param.nbpointdiscret,(j+1)*(wpg_param.poly_degree+1)+1:(j+1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
        end
        function [M_t1_dd] = compute_M_t1_dd(wpg_param,cutting)
            M_t1_dd=zeros(wpg_param.nbpointdiscret,(wpg_param.nbphases+2)*(wpg_param.poly_degree+1));

            %%%zmp1 init cut in two%%%
            nb=0;
            j=1;
            tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*cutting-1/wpg_param.frequency]-wpg_param.tpassage(j);
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt n*(n-1)*tstep'.^max([0 n-2])];
            end
            M_t1_dd(nb+1:nb+size(Dt,1),(j-1)*(wpg_param.poly_degree+1)+1:(j-1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
            nb=nb+size(Dt,1);
            
            tstep=[(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*cutting:1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*cutting;
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt n*(n-1)*tstep'.^max([0 n-2])];
            end
            M_t1_dd(nb+1:nb+size(Dt,1),(j)*(wpg_param.poly_degree+1)+1:(j)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
            nb=nb+size(Dt,1);

            %%%generation of dt%%%
            for j=2:wpg_param.nbphases-1
                tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-wpg_param.tpassage(j);
                Dt=[];
                for n=0:wpg_param.poly_degree
                    Dt=[Dt n*(n-1)*tstep'.^max([0 n-2])];
                end
                M_t1_dd(nb+1:nb+size(Dt,1),(j)*(wpg_param.poly_degree+1)+1:(j+1)*(wpg_param.poly_degree+1))=Dt;
                nb=nb+size(Dt,1);
%                 [size(Dt,1)-wpg_param.discretization(j)]
            end

            %%%zmp1 fin cut in two%%%
            j=wpg_param.nbphases;
            tstep=[wpg_param.tpassage(j):1/wpg_param.frequency:(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*(1-cutting)-1/wpg_param.frequency]-wpg_param.tpassage(j);
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt n*(n-1)*tstep'.^max([0 n-2])];
            end
            M_t1_dd(nb+1:nb+size(Dt,1),(j)*(wpg_param.poly_degree+1)+1:(j)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
            nb=nb+size(Dt,1);
            
            tstep=[(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*(1-cutting):1/wpg_param.frequency:wpg_param.tpassage(j+1)-1/wpg_param.frequency]-(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*(1-cutting);
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt n*(n-1)*tstep'.^max([0 n-2])];
            end
            M_t1_dd(nb+1:nb+size(Dt,1),(j+1)*(wpg_param.poly_degree+1)+1:(j+1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
            nb=nb+size(Dt,1);
            
            Dt=[];
            for n=0:wpg_param.poly_degree
                Dt=[Dt n*(wpg_param.tpassage(end)-(wpg_param.tpassage(j)+wpg_param.tpassage(j+1))*(1-cutting))'.^max([0 n-2])];
            end
            M_t1_dd(wpg_param.nbpointdiscret,(j+1)*(wpg_param.poly_degree+1)+1:(j+1)*(wpg_param.poly_degree+1)+(wpg_param.poly_degree+1))=Dt;
        end
    end
end