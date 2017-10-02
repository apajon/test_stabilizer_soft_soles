function [time_ pzmp pcom fzmp] = zmp_under_foot(wpg_param,zmp,trajectories_zmp,foot_step_wanted,time,zfzmp1,discretization_)
% (foot_step_wanted,wpg_param.nbpankle,time,trajectories_zmp.xpzmp,trajectories_zmp.ypzmp,trajectories_zmp.xpzmp1,trajectories_zmp.ypzmp1,trajectories_zmp.xpzmp2,trajectories_zmp.ypzmp2,trajectories_zmp.xpcom,trajectories_zmp.ypcom,zpcom,zfzmp1,zmp.A_xfcom,zmp.B_xfcom,zmp.A_yfcom,zmp.B_yfcom,wpg_param.psa_abcd,wpg_param.discretization,discretization_,wpg_param.mg)
zpcom=wpg_param.z*ones(size(trajectories_zmp.xpcom));
reduce_zmp2=3;
%% % adapt zmp and ZMP1&2 to have the CoP under the wanted foot step
if foot_step_wanted == 1
    indice_zmp1=1:sum(discretization_(1:wpg_param.nbpolypi));
    
    xpzmp1_=trajectories_zmp.xpzmp1(indice_zmp1);
    ypzmp1_=trajectories_zmp.ypzmp1(indice_zmp1);
    
    xpzmp_=[xpzmp1_];
    ypzmp_=[ypzmp1_];

    time_=time(indice_zmp1);

    xpcom_=trajectories_zmp.xpcom(indice_zmp1);
    ypcom_=trajectories_zmp.ypcom(indice_zmp1);
    zpcom_=zpcom(indice_zmp1);

    zfzmp1_=zfzmp1(indice_zmp1)*wpg_param.mg;
    zfzmp_=[zfzmp1_;];

    xfzmp=-(zmp.A_xfcom*wpg_param.psa_abcd(1:length(wpg_param.psa_abcd)/2)+zmp.B_xfcom);
    yfzmp=-(zmp.A_yfcom*wpg_param.psa_abcd(length(wpg_param.psa_abcd)/2+1:end)+zmp.B_yfcom);
    xfzmp_=xfzmp(indice_zmp1).*(zfzmp_/wpg_param.mg);
    yfzmp_=yfzmp(indice_zmp1).*(zfzmp_/wpg_param.mg);
elseif foot_step_wanted == 2
    indice_zmp=sum(wpg_param.discretization(1:wpg_param.nbpolypi))+1:sum(wpg_param.discretization(1:wpg_param.nbpolypi+wpg_param.nbpolyssp));
    indice_zmp1=sum(discretization_(1:wpg_param.nbpolypi))+1:sum(discretization_(1:wpg_param.nbpolypi+wpg_param.nbpolydsp));
    indice_zmp2=1:sum(discretization_(1:wpg_param.nbpolypi));
    indice_tot=1:sum(wpg_param.discretization(1:wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp));

    xpzmp1_=trajectories_zmp.xpzmp1(indice_zmp1);
    ypzmp1_=trajectories_zmp.ypzmp1(indice_zmp1);
    xpzmp2_=trajectories_zmp.xpzmp2(indice_zmp2);
    ypzmp2_=trajectories_zmp.ypzmp2(indice_zmp2);
    xpzmp_=trajectories_zmp.xpzmp(indice_zmp);
    ypzmp_=trajectories_zmp.ypzmp(indice_zmp);

    xpzmp_=[xpzmp2_;xpzmp_;xpzmp1_];
    ypzmp_=[ypzmp2_;ypzmp_;ypzmp1_];

    time_=time(indice_tot);

    xpcom_=trajectories_zmp.xpcom(indice_tot);
    ypcom_=trajectories_zmp.ypcom(indice_tot);
    zpcom_=zpcom(indice_tot);

    zfzmp1_=zfzmp1(indice_zmp1)*wpg_param.mg;
    zfzmp2_=(1-zfzmp1(indice_zmp2))*wpg_param.mg;
    zfzmp_=ones(length(indice_zmp),1)*wpg_param.mg;
    zfzmp_=[zfzmp2_;zfzmp_;zfzmp1_;];

    xfzmp=-(zmp.A_xfcom*wpg_param.psa_abcd(1:length(wpg_param.psa_abcd)/2)+zmp.B_xfcom);
    yfzmp=-(zmp.A_yfcom*wpg_param.psa_abcd(length(wpg_param.psa_abcd)/2+1:end)+zmp.B_yfcom);
    xfzmp_=xfzmp(indice_tot).*(zfzmp_/wpg_param.mg);
    yfzmp_=yfzmp(indice_tot).*(zfzmp_/wpg_param.mg);
elseif foot_step_wanted==wpg_param.nbpankle-1
    indice_zmp=size(trajectories_zmp.xpzmp,1)-sum(wpg_param.discretization(end-wpg_param.nbpolypf-wpg_param.nbpolyssp+1:end)):size(trajectories_zmp.xpzmp,1)-sum(wpg_param.discretization(end-wpg_param.nbpolypf+1:end))-1;
    indice_zmp1=size(trajectories_zmp.xpzmp1,1)-sum(discretization_(end-wpg_param.nbpolypf+1:end)):size(trajectories_zmp.xpzmp1,1);
    indice_zmp2=size(trajectories_zmp.xpzmp2,1)-sum(discretization_(end-wpg_param.nbpolypf-wpg_param.nbpolydsp+1:end))+2*reduce_zmp2:size(trajectories_zmp.xpzmp2,1)-sum(discretization_(end-wpg_param.nbpolypf+1:end))+reduce_zmp2-1;
    indice_tot=size(trajectories_zmp.xpzmp,1)-sum(wpg_param.discretization(end-wpg_param.nbpolypf-wpg_param.nbpolyssp-wpg_param.nbpolydsp+1:end))+reduce_zmp2:size(trajectories_zmp.xpzmp,1);

    xpzmp1_=trajectories_zmp.xpzmp1(indice_zmp1);
    ypzmp1_=trajectories_zmp.ypzmp1(indice_zmp1);
    xpzmp2_=trajectories_zmp.xpzmp2(indice_zmp2);
    ypzmp2_=trajectories_zmp.ypzmp2(indice_zmp2);
    xpzmp_=trajectories_zmp.xpzmp(indice_zmp);
    ypzmp_=trajectories_zmp.ypzmp(indice_zmp);

    xpzmp_=[xpzmp2_;xpzmp_;xpzmp1_];
    ypzmp_=[ypzmp2_;ypzmp_;ypzmp1_];

    time_=time(indice_tot);

    xpcom_=trajectories_zmp.xpcom(indice_tot);
    ypcom_=trajectories_zmp.ypcom(indice_tot);
    zpcom_=zpcom(indice_tot);

    zfzmp1_=zfzmp1(indice_zmp1)*wpg_param.mg;
%     zfzmp2_=(1-zfzmp1(end-sum(discretization_(end-wpg_param.nbpolypf-wpg_param.nbpolydsp+1:end))+2*reduce_zmp2:end-sum(discretization_(end-wpg_param.nbpolypf+1:end))+reduce_zmp2-1))*wpg_param.mg;
    zfzmp2_=(1-zfzmp1(indice_zmp2+reduce_zmp2*(foot_step_wanted-2)))*wpg_param.mg;
    zfzmp_=ones(length(indice_zmp),1)*wpg_param.mg;
    zfzmp_=[zfzmp2_;zfzmp_;zfzmp1_;];

    xfzmp=-(zmp.A_xfcom*wpg_param.psa_abcd(1:length(wpg_param.psa_abcd)/2)+zmp.B_xfcom);
    yfzmp=-(zmp.A_yfcom*wpg_param.psa_abcd(length(wpg_param.psa_abcd)/2+1:end)+zmp.B_yfcom);
    xfzmp_=xfzmp(indice_tot).*(zfzmp_/wpg_param.mg);
    yfzmp_=yfzmp(indice_tot).*(zfzmp_/wpg_param.mg);
elseif foot_step_wanted==wpg_param.nbpankle
    indice_zmp2=-sum(wpg_param.discretization(end-wpg_param.nbpolypf+1:end))+reduce_zmp2;

    xpzmp2_=trajectories_zmp.xpzmp2(end+indice_zmp2:end);
    ypzmp2_=trajectories_zmp.ypzmp2(end+indice_zmp2:end);

    xpzmp_=[xpzmp2_;];
    ypzmp_=[ypzmp2_;];

    time_=time(end+indice_zmp2:end);

    xpcom_=trajectories_zmp.xpcom(end+indice_zmp2:end);
    ypcom_=trajectories_zmp.ypcom(end+indice_zmp2:end);
    zpcom_=zpcom(end+indice_zmp2:end);

    zfzmp2_=(1-zfzmp1(end+indice_zmp2:end))*wpg_param.mg;
    zfzmp_=[zfzmp2_;];

    xfzmp=-(zmp.A_xfcom*wpg_param.psa_abcd(1:length(wpg_param.psa_abcd)/2)+zmp.B_xfcom);
    yfzmp=-(zmp.A_yfcom*wpg_param.psa_abcd(length(wpg_param.psa_abcd)/2+1:end)+zmp.B_yfcom);
    xfzmp_=xfzmp(end+indice_zmp2:end).*(zfzmp_/wpg_param.mg);
    yfzmp_=yfzmp(end+indice_zmp2:end).*(zfzmp_/wpg_param.mg);
else
    indice_zmp=sum(wpg_param.discretization(1:wpg_param.nbpolypi+wpg_param.nbpolyssp*(foot_step_wanted-2)+wpg_param.nbpolydsp*(foot_step_wanted-2)))+1:sum(wpg_param.discretization(1:wpg_param.nbpolypi+wpg_param.nbpolyssp*(foot_step_wanted-1)+wpg_param.nbpolydsp*(foot_step_wanted-2)));
    indice_zmp1=sum(discretization_(1:wpg_param.nbpolypi+wpg_param.nbpolydsp*(foot_step_wanted-2)))+1:sum(discretization_(1:wpg_param.nbpolypi+wpg_param.nbpolydsp*(foot_step_wanted-1)));
    indice_zmp2=sum(discretization_(1:wpg_param.nbpolypi+wpg_param.nbpolydsp*(foot_step_wanted-3)))+1-reduce_zmp2*(foot_step_wanted-3):sum(discretization_(1:wpg_param.nbpolypi+wpg_param.nbpolydsp*(foot_step_wanted-2)))-reduce_zmp2*(foot_step_wanted-2);
    indice_tot=sum(wpg_param.discretization(1:wpg_param.nbpolypi+wpg_param.nbpolyssp*(foot_step_wanted-2)+wpg_param.nbpolydsp*(foot_step_wanted-2)))+1-length(indice_zmp2):sum(wpg_param.discretization(1:wpg_param.nbpolypi+wpg_param.nbpolyssp*(foot_step_wanted-1)+wpg_param.nbpolydsp*(foot_step_wanted-2)))+length(indice_zmp1);
    
    xpzmp1_=trajectories_zmp.xpzmp1(indice_zmp1);
    ypzmp1_=trajectories_zmp.ypzmp1(indice_zmp1);
    xpzmp2_=trajectories_zmp.xpzmp2(indice_zmp2);
    ypzmp2_=trajectories_zmp.ypzmp2(indice_zmp2);
    xpzmp_=trajectories_zmp.xpzmp(indice_zmp);
    ypzmp_=trajectories_zmp.ypzmp(indice_zmp);

    xpzmp_=[xpzmp2_;xpzmp_;xpzmp1_];
    ypzmp_=[ypzmp2_;ypzmp_;ypzmp1_];

    time_=time(indice_tot);

    xpcom_=trajectories_zmp.xpcom(indice_tot);
    ypcom_=trajectories_zmp.ypcom(indice_tot);
    zpcom_=zpcom(indice_tot);

    zfzmp1_=zfzmp1(indice_zmp1)*wpg_param.mg;
    zfzmp2_=(1-zfzmp1(indice_zmp1(reduce_zmp2+1:end)))*wpg_param.mg;
    zfzmp_=ones(length(indice_zmp),1)*wpg_param.mg;
    zfzmp_=[zfzmp2_;zfzmp_;zfzmp1_;];

    xfzmp=-(zmp.A_xfcom*wpg_param.psa_abcd(1:length(wpg_param.psa_abcd)/2)+zmp.B_xfcom);
    yfzmp=-(zmp.A_yfcom*wpg_param.psa_abcd(length(wpg_param.psa_abcd)/2+1:end)+zmp.B_yfcom);
    xfzmp_=xfzmp(indice_tot).*(zfzmp_/wpg_param.mg);
    yfzmp_=yfzmp(indice_tot).*(zfzmp_/wpg_param.mg);
end

pzmp=[xpzmp_ ypzmp_];
pcom=[xpcom_ ypcom_ zpcom_];
fzmp=[xfzmp_ yfzmp_ zfzmp_];

% %%
% trajectories=[time_' xpzmp_ ypzmp_ xpcom_ ypcom_ xfzmp yfzmp zfzmp];
% 
% zmpcom=fopen('exemple_trajectoire.txt','w');
% for i=1:size(trajectories,1)
%     fprintf(zmpcom,'%f %f %f %f %f %f %f %f\n',trajectories(i,:));
% end
% fclose(zmpcom);
% %%%%%%%%%%%
end