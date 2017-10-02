#include <iostream>
#include <math.h>
#include "mex.h"
#include "Gausstot2Angles.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <fstream>

const double GaussSeidel::epsiGauss = 1e-6;
const double GaussSeidel::epsiZMP = 1e-3;
const double GaussSeidel::epsiSig = 1e-6;
const double GaussSeidel::epsiCou = 1e-6;
const int GaussSeidel::D = 3;
using namespace std;
using namespace Eigen;

extern void _main();

typedef Matrix<double,5,1> Vector5d;
typedef Matrix<double,6,1> Vector6d;

/****************************/

GaussSeidel::GaussSeidel(){}

GaussSeidel::~GaussSeidel(){}

void GaussSeidel::setStaticDataSim(const MatrixXd& Ccc_in,const Vector5d& FtotZMPdes_in,double friction,int m1){

    Ccc = Ccc_in;
    FtotZMPdes = FtotZMPdes_in;
    fric = friction;
    m = m1;// number of contact surface nodes
	FtotZMP.Zero();
	//init
    J.resize(2*D-1,2*D);
    J2.resize(2*D-1,2*D);
	displ_m.resize(D,m);
	displ_ini_m.resize(D,m);
    F3.resize(D*m);
    Fold3.resize(D*m);
	Rpfree.resize(D,m);
	RpabsOld.resize(D,m);
	W.resize(m);
	inv_W.resize(m);
    RtF.resize(D,m);
    RtF3.resize(D*m);
	deltaIni.resize(D,m);
	qt.resize(m);
	qn.resize(m);
	ci.resize(m);
	di.resize(m);
	ind_slip.reserve(m);
	ind_c_slip.reserve(m);
	ind_stick.reserve(m);
	ind_c_stick.reserve(m);
	ind_cont.reserve(m);
	ind_c_cont.reserve(m);
	G1.resize(2*D,2*D);
    Kcart.resize(2*D,2*D);
}

void GaussSeidel::setStaticDataStep(const Vector3d& displ_in,const Vector2d& angleact_in,const VectorXd& F3_in,const MatrixXd& Pfree_in,const MatrixXd& PabsOld_in){
	displ = displ_in;
    angleact = angleact_in;
    F3 = F3_in;
    PabsOld = PabsOld_in;
    Pfree = Pfree_in;
}

void GaussSeidel::run(int contAngle){
	FtotZMP = VectorXd::Zero(5);
	theta = angleact(0);
	phi = angleact(1);
	displangle << displ(0),displ(1),displ(2),theta,phi;
    ZMPdes << FtotZMPdes(3),FtotZMPdes(4);
	//mexPrintf("FtotZMPdes = %g %g %g %g %g %g\n", FtotZMPdes(0), FtotZMPdes(1), FtotZMPdes(2), FtotZMPdes(3), FtotZMPdes(4), FtotZMPdes(5));
	int ite = 0;
    int ite_tot = 0;
	double sum_gauss = 1;
	double norm_delta_displ;
	double norm_delta_angle;
	double max_norm_ddispl = 0.005;/// mm
	double max_norm_dangle = 10.*3.14159/180.;/// °
	//// Compute Rotation
	Rtheta << 1, 0, 0,
		    0, cos(theta), -sin(theta),
			0, sin(theta), cos(theta);
	Rphi << cos(phi), 0, sin(phi),
		    0, 1, 0,
			-sin(phi), 0, cos(phi);
	R = Rtheta * Rphi;
	Rpfree = R * Pfree;
	//RpabsOld = R * PabsOld;
	for (int i = 0; i < m; ++i)
	{
		i3 = 3 * i;
		A = R * Ccc.block<3, 3>(i3, i3) * R.transpose();
		A22 = A.block<2,2>(0,0);
		W[i] = A;
		inv_W[i] = A.inverse();
		//Find the eigen values
		ad = A22(0,0)/2 + A22(1,1)/2;
		ad12 = sqrt(A22(0,0)*A22(0,0) - 2*A22(0,0)*A22(1,1) + A22(1,1)*A22(1,1) + 4*A22(0,1)*A22(1,0))/2;
		lambda_min = ad - ad12;
		lambda_max = ad + ad12;
		qt(i) = lambda_min/(lambda_max*lambda_max);
		qn(i) = 1/A(2,2);
	}
	if (contAngle==1){
		Vector3d displ_ini(0, 0, 0);
		displ_ini_m = displ_ini.replicate(1,m);
		displ_m = displ.replicate(1,m);
		//deltaIni.row(0) = Rpfree.row(0) - RpabsOld.row(0) + (displ_ini_m.row(0) + displ_m.row(0));
		//deltaIni.row(1) = Rpfree.row(1) - RpabsOld.row(1) + (displ_ini_m.row(1) + displ_m.row(1));
		deltaIni.row(0) = displ_ini_m.row(0) + displ_m.row(0);
		deltaIni.row(1) = displ_ini_m.row(1) + displ_m.row(1);
		deltaIni.row(2) = Rpfree.row(2) + (displ_ini_m.row(2) + displ_m.row(2));
	}
	else{
		displ_m = displ.replicate(1,m);
		deltaIni.row(0) = Rpfree.row(0) - PabsOld.row(0) + displ_m.row(0);
		deltaIni.row(1) = Rpfree.row(1) - PabsOld.row(1) + displ_m.row(1);
		deltaIni.row(2) = Rpfree.row(2) + displ_m.row(2);
	}
	//while((sum_gauss > epsiGauss) || ((FtotZMPdes.segment<5>(0)-FtotZMP.segment<5>(0)).norm() > epsiGauss))
	critFZMP = (FtotZMPdes-FtotZMP).norm();
    while((sum_gauss > epsiGauss) || (critFZMP > epsiZMP))
	{
		int cont_stick = 0;
		int cont_slip = 0;
		int cont_c = 0;
		//ci = MatrixXd::Zero(D,D*m);
		//di = MatrixXd::Zero(D*m,D*m);
		Fold3 = F3;
		Map<MatrixXd> F = Map<MatrixXd>(F3.data(),D,m); //Reshape

		RtF = R.transpose() * F;
		RtF3 = Map<VectorXd>(RtF.data(),D*m); //Reshape
        /*if (ite==0){
            for (int i = 0; i < D*m; ++i){
                 //for (int j = 0; j < m; ++j){
                    mexPrintf("Value1 = %g\n", RtF3(i));
                //}
            }
        }  */
		ind_slip.clear();
		ind_c_slip.clear();
		ind_stick.clear();
		ind_c_stick.clear();
		ind_cont.clear();
		ind_c_cont.clear();
		ite = ite + 1;
		for (int i = 0; i < m; ++i)
		{

			deltaTest = deltaIni.col(i);
			int i3 = 3 * i;

			deltaTest = deltaTest  + R*(Ccc.middleCols<3>(i3).transpose()*RtF3 - Ccc.block<3, 3>(i3, i3)*RtF3.segment<3>(i3));
			if (deltaTest(2)<epsiSig){
				F3.segment<3>(i3) = -inv_W[i] * deltaTest;
				if (F3.segment<2>(i3).norm() > (fric * F3(i3+2))){
					Fcfric = Vector3d::Zero();
					if (fric>0.2){
						Fcfric = F3.segment<3>(i3);
						deltaNew = Vector3d::Zero();
					}
					else{
						Fcfric(2) = -qn(i)*deltaTest(2); //deltaTest/W(3,3)
						deltaNew = W[i] * Fcfric + deltaTest;
					}

					//Newton Contact for sliding contacts
					diffFric = 1;
					kn = 0;
					G21 = Matrix3d::Zero();
					G22 = Matrix3d::Zero();
					fricTest;
					while(diffFric > epsiCou){
						phi1 = deltaNew-deltaTest-W[i]*Fcfric;
						ds = (Fcfric.head<2>()-qt(i)*deltaNew.head<2>());
						phi2.head<2>() = Fcfric.head<2>()-(fric*Fcfric(2)*ds/ds.norm());
						phi2(2) = qn(i)*deltaNew(2);
						dsMat << pow(ds(1),2), -ds(0)*ds(1),
								-ds(0)*ds(1), pow(ds(0),2);
						Bigpi = pow(1/ds.norm(),3)*dsMat;
						alpha = fric*Fcfric(2)*Bigpi;
						G21.block<2, 2>(0, 0) = qt(i)*alpha;
						G21(2,2) = qn(i);
						G22.block<2, 2>(0, 0) = MatrixXd::Identity(2,2)-alpha;
						G22.block<2, 1>(0, 2) = -fric*ds/ds.norm();
						d2 = (G21*W[i]+G22).inverse() * (phi2-G21*phi1);
						d1 = phi1 + W[i]*d2;
						deltaNew = deltaNew - d1;
						Fcfric = Fcfric - d2;
						fricTest = Fcfric.head<2>().norm()/Fcfric(2);
						diffFric = fabs(fricTest-fric);
						kn = kn+1;
                        if (kn>500){diffFric = 1e-12;}
					}
					// Find the gradient
					ds = (Fcfric.head<2>()-qt(i)*deltaNew.head<2>());
					Bigpi = pow(1/ds.norm(),3)*dsMat;
					alpha = fric*Fcfric(2)*Bigpi;
					Matrix3d Ci0 = Matrix3d::Zero();
					Ci0.block<2, 2>(0, 0) = MatrixXd::Identity(2,2)-alpha;
					Ci0.block<2, 1>(0, 2) = -fric*ds/ds.norm();
					ci[cont_slip] = Ci0;
					Matrix3d Di0 = Matrix3d::Zero();
					Di0.block<2, 2>(0, 0) = qt(i)*alpha;
					Di0(2, 2) = qn(i);
					di[cont_slip] = Di0;

					F3.segment<3>(i3) = Fcfric;
					ind_slip.push_back(i);
					//ind_slip.push_back(cont_slip);
					ind_c_slip.push_back(cont_c);
					cont_slip = cont_slip + 1;
				}
				else{
					ind_stick.push_back(i);
					ind_c_stick.push_back(cont_c);
					cont_stick = cont_stick + 1;
				}
				ind_cont.push_back(i);
				ind_c_cont.push_back(cont_c);
				cont_c = cont_c + 1;
			}
			else{
				F3.segment<3>(i3) = Vector3d::Zero();
			}

			RtF3.segment<3>(i3) = R.transpose() * F3.segment<3>(i3);
		}

		sum_gauss = 0;
		for (int i = 0; i < m; ++i){
			norm_f_diff = (F3.segment<3>(D*i)-Fold3.segment<3>(D*i)).norm();
			norm_f_new = F3.segment<3>(D*i).norm();
			if (norm_f_new!=0){
				sum_gauss = sum_gauss + norm_f_diff/norm_f_new;
			}
		}

		mslip = (int)ind_slip.size();
		mstick = (int)ind_stick.size();
		mc = (int)ind_cont.size();
		Pfree_mc.resize(D,mc);
		Pfree_mc3.resize(D*mc);
		Fc_mc3.resize(D*mc);
		P_mc3.resize(D*mc);
		RtFc_mc3.resize(D*mc);
		Ftot = Vector3d::Zero();
//		W_mc = MatrixXd::Zero(D*mc,D*mc);
		for (int i = 0; i < mc; ++i){
			Pfree_mc3.segment<3>(D*ind_c_cont[i]) = Pfree.col(ind_cont[i]);
			Pfree_mc.col(ind_c_cont[i]) = Pfree.col(ind_cont[i]);
			Fc_mc3.segment<3>(D*ind_c_cont[i]) = F3.segment<3>(D*ind_cont[i]);
			RtFc_mc3.segment<3>(D*ind_c_cont[i]) = RtF3.segment<3>(D*ind_cont[i]);
			Ftot = Ftot + Fc_mc3.segment<3>(D*ind_c_cont[i]);
//			for (int j = 0; j < mc; ++j){
//                /// ************************* Not optimized, W_mc only necessary when computing J, that is when changing pos and rot of sole ***************************
//                /// ************************* Could be replaced by first computation of Fc in sole frame, then deformation in sole frame, then changed in absolute frame ***************************
//                /// ************************* and Could be recomputed only when R is changed ***************************
//				W_mc.block<3, 3>(D*ind_c_cont[i],D*ind_c_cont[j]) = R * Ccc.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose();
//			}
		}
		RPfree_mc = R * Pfree_mc;
		RPfree_mc3 = Map<VectorXd>(RPfree_mc.data(),D*mc);
		for (int i = 0; i < mc; ++i){
			deltaInt = Vector3d::Zero();
			for (int j = 0; j < mc; ++j){
                /// On recalcule tous les déplacements avec toutes les nouvelles forces. Ce calcul est partiellement fait au dessus avec une partie des anciennes forces
                /// Est-il vraiment nécessaire de le refaire ?
                deltaInt += Ccc.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * RtFc_mc3.segment<3>(D*ind_c_cont[j]);
			}
			P_mc3.segment<3>(D*ind_c_cont[i]) = displ + RPfree_mc3.segment<3>(D*ind_c_cont[i]) + R*deltaInt;
		}

		displ_mc3 = displ.replicate(mc, 1);

        /// Checked if ind_c_cont == 1:mc
        /// Use of ind_c_cont in the following is useless
        /*bool is_seq = true;
        for (int i = 0; i < mc; ++i){
            //mexPrintf("%i ", ind_c_cont[i]);
            if (ind_c_cont[i] != i)
                is_seq = false;
        }
        if (!is_seq) {
            mexPrintf("ind_c_cont = ");
            for (int i = 0; i < mc; ++i){
                //mexPrintf("%i ", ind_c_cont[i]);
                if (ind_c_cont[i] != i)
                    is_seq = false;
            }
            mexPrintf("\n");
        }*/

//		P_mc3 = displ_mc3 + RPfree_mc3 + W_mc * Fc_mc3;
//        mexPrintf("P = \n");
//        for (int i = 0; i < mc; ++i){
//            mexPrintf("%g %g \n", P_mc3(i), P_mc3_bis(i));
//        }
		Mo = Vector3d::Zero();
		for (int i = 0; i < mc; ++i){
			Mo = Mo + P_mc3.segment<3>(D*i).cross(Fc_mc3.segment<3>(D*i));
		}
		Z << -Mo(1)/Ftot(2), Mo(0)/Ftot(2);
        //Vector3d Acio;
        //Acio << -Mo(1)/Ftot(2), Mo(0)/Ftot(2),0.;
        Vector3d ZMPdes3d;
        ZMPdes3d(0)=ZMPdes(0);
        ZMPdes3d(1)=ZMPdes(1);
        ZMPdes3d(2)=0.;
        Mzmp = Mo - ZMPdes3d.cross(Ftot);
		/*if (contAngle==1){
			FtotZMP << Ftot(0), Ftot(1), Ftot(2), Z(0), Z(1), 0;
		}
		else{*/
			FtotZMP << Ftot(0), Ftot(1), Ftot(2), Z(0), Z(1);
		//}
		//FtotZMP << Ftot(0), Ftot(1), Ftot(2), Z(0), Z(1), Mzmp(2); //??????? - ou + Mzmp(2)
        //mexPrintf("FtotZMP    = %g %g %g %g %g %g\n", FtotZMP(0), FtotZMP(1), FtotZMP(2), FtotZMP(3), FtotZMP(4), FtotZMP(5));
		//if (sum_gauss<(0.0001*(FtotZMPdes.segment<5>(0)-FtotZMP.segment<5>(0)).norm())){
        /*if (ite==1){
            gradFtotZ(contAngle);
        	mexPrintf("J = \n");
            for (int i = 0; i < 6; ++i){
                for (int j = 0; j < 6; ++j){
                    mexPrintf("%g ", J(i,j));
                }
                mexPrintf("\n");
            }
        }*/
        FtotZMPerr = FtotZMPdes-FtotZMP;
        if ((sum_gauss<((0.0001*FtotZMPerr).norm())) || (sum_gauss<epsiGauss)){
            ite_tot = ite_tot + 1;
			gradFtotZ(contAngle);

			//displangle = displangle + J.colPivHouseholderQr().solve(FtotZMPdes-FtotZMP);
			//delta_displangle = J.householderQr().solve(FtotZMPerr);
			delta_displangle = J.colPivHouseholderQr().solve(FtotZMPerr);
            //displangle = displangle + J.inverse() *(FtotZMPerr);
			norm_delta_displ = delta_displangle.segment<3>(0).norm();
			norm_delta_angle = delta_displangle.segment<2>(3).norm();
			if (norm_delta_displ>max_norm_ddispl) {
                if (norm_delta_displ*max_norm_dangle>max_norm_ddispl*norm_delta_angle) {
                    delta_displangle *= max_norm_ddispl/norm_delta_displ;
                    //mexPrintf("Reduction of pos-rot step by %g\n",norm_delta_displ/max_norm_ddispl);
                } else {
                    delta_displangle *= max_norm_dangle/norm_delta_angle;
                    //mexPrintf("Reduction of pos-rot step by %g\n",norm_delta_angle/max_norm_dangle);
                }
			} else if (norm_delta_angle>max_norm_dangle) {
                    delta_displangle *= max_norm_dangle/norm_delta_angle;
                    //mexPrintf("Reduction of pos-rot step by %g\n",norm_delta_angle/max_norm_dangle);
			}
			displangle = displangle + delta_displangle;
			displ = displangle.segment<3>(0);
			angleact = displangle.segment<2>(3);
			theta = angleact(0);
			phi = angleact(1);
			//// Compute Rotation
			Rtheta << 1, 0, 0,
					0, cos(theta), -sin(theta),
					0, sin(theta), cos(theta);
			Rphi << cos(phi), 0, sin(phi),
					0, 1, 0,
					-sin(phi), 0, cos(phi);
			R = Rtheta * Rphi;
			Rpfree = R * Pfree;
			//RpabsOld = R * PabsOld;
			for (int i = 0; i < m; ++i)
			{
				i3 = 3 * i;
				A = R * Ccc.block<3, 3>(i3, i3) * R.transpose();
				A22 = A.block<2,2>(0,0);
				W[i] = A;
				inv_W[i] = A.inverse();
				//Find the eigen values
				ad = A22(0,0)/2 + A22(1,1)/2;
				ad12 = sqrt(A22(0,0)*A22(0,0) - 2*A22(0,0)*A22(1,1) + A22(1,1)*A22(1,1) + 4*A22(0,1)*A22(1,0))/2;
				lambda_min = ad - ad12;
				lambda_max = ad + ad12;
				qt(i) = lambda_min/(lambda_max*lambda_max);
				qn(i) = 1/A(2,2);
			}
			if (contAngle==1){
				Vector3d displ_ini(0, 0, 0);
				displ_ini_m = displ_ini.replicate(1,m);
				displ_m = displ.replicate(1,m);
				//deltaIni.row(0) = Rpfree.row(0) - RpabsOld.row(0) + (displ_ini_m.row(0) + displ_m.row(0));
				//deltaIni.row(1) = Rpfree.row(1) - RpabsOld.row(1) + (displ_ini_m.row(1) + displ_m.row(1));
				deltaIni.row(0) = displ_ini_m.row(0) + displ_m.row(0);
				deltaIni.row(1) = displ_ini_m.row(1) + displ_m.row(1);
				deltaIni.row(2) = Rpfree.row(2) + (displ_ini_m.row(2) + displ_m.row(2));
			}
			else{
				displ_m = displ.replicate(1,m);
				deltaIni.row(0) = Rpfree.row(0) - PabsOld.row(0) + displ_m.row(0);
				deltaIni.row(1) = Rpfree.row(1) - PabsOld.row(1) + displ_m.row(1);
				deltaIni.row(2) = Rpfree.row(2) + displ_m.row(2);
			}
            //mexPrintf("FtotZMP    = %g %g %g %g %g %g %g %g %g %g %g %g\n", FtotZMP(0), FtotZMP(1), FtotZMP(2), FtotZMP(3), FtotZMP(4), FtotZMP(5), displ(0), displ(1), displ(2), angleact(0), angleact(1), angleact(2));
		} else {
            //mexPrintf("FtotZMP    = %g %g %g %g %g %g\n", FtotZMP(0), FtotZMP(1), FtotZMP(2), FtotZMP(3), FtotZMP(4), FtotZMP(5));
		}
        critFZMP = (FtotZMPdes-FtotZMP).norm();
        if (ite_tot>250 || ite>25000){
            critFZMP = 1e-12;
            sum_gauss = 1e-12;
            displ << 0,0,0;
			/*if (ite_tot>500){
				mexPrintf("More than 500 iterations \n");
			}
			else{
				mexPrintf("More than 50000 iterations of GaussSeidel \n");
			}*/
        }

	}
	//mexPrintf("%i iterations\n",ite);
    //mexPrintf("Mzmp_z    = %g\n", Mzmp(2));
//	gradFtotZ(contAngle);
//    StiffCart(contAngle);
    StiffCart_Omega(contAngle);
	//Kcart = MatrixXd::Zero(6,6);
}

void GaussSeidel::gradFtotZ(int contAngle)
{

	A_grad = MatrixXd::Zero(D*mc+D*mslip,D*mc+D*mslip);
	B_grad = MatrixXd::Zero(D*mc+D*mslip,2*D-1);
	D_grad = MatrixXd::Zero(2*D-1,D*mc+D*mslip);
	int cont_slip1 = 0;
	/* Compute derivative of rotation */
	dRdtheta << 0, 0, 0,
				 0, -sin(theta), -cos(theta),
				 0, cos(theta), -sin(theta);
	dRdtheta =	dRdtheta*Rphi;
	dRdphi << -sin(phi), 0, cos(phi),
				0, 0, 0,
				-cos(phi), 0, -sin(phi);
	dRdphi = Rtheta*dRdphi;
    if (contAngle==1){
        PfreeFt1 << 0, 0, 0;
        PfreeFt2 << 0, 0, 0;
        PfreeFn << 0, 0, 0;
    }
	// find A_grad
    W_mc = MatrixXd::Zero(D*mc,D*mc);
    for (int i = 0; i < mc; ++i){
        for (int j = 0; j < mc; ++j){
            W_mc.block<3, 3>(D*ind_c_cont[i],D*ind_c_cont[j]) = R * Ccc.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose();
        }
    }
	A_grad.block(0,0,D*mc,D*mc) = -W_mc;

	for (int i = 0; i < mc; ++i){
		B_grad.block<3, 3>(D*i,0) = Matrix3d::Identity();
		D_grad.block<3, 3>(0,D*i) = Matrix3d::Identity();
		D_grad.block<2, 1>(D,D*i+2) = (1/Ftot(2)) * (P_mc3.segment<2>(D*i) - Z.segment<2>(0));
        //D_grad(D+2,D*i) = -P_mc3(D*i+1) + ZMPdes(1);//Z(1);
        //D_grad(D+2,D*i+1) = P_mc3(D*i) - ZMPdes(0);//Z(0);
		if (contAngle==1){
            btheta << 0,0,dRdtheta.row(2)*Pfree_mc3.segment<3>(D*i);
            bphi << 0,0,dRdphi.row(2)*Pfree_mc3.segment<3>(D*i);
            PfreeFt1 = PfreeFt1 + Pfree_mc3.segment<3>(D*i)*Fc_mc3(D*i);
            PfreeFt2 = PfreeFt2 + Pfree_mc3.segment<3>(D*i)*Fc_mc3(D*i+1);
            PfreeFn = PfreeFn + Pfree_mc3.segment<3>(D*i)*Fc_mc3(D*i+2);
		} else {
            btheta = dRdtheta*Pfree_mc3.segment<3>(D*i);
            bphi = dRdphi*Pfree_mc3.segment<3>(D*i);
		}
		for (int j = 0; j < mc; ++j){
			btheta = btheta + (dRdtheta * Ccc.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose() + R * Ccc.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * dRdtheta.transpose()) * Fc_mc3.segment<3>(D*j);
			bphi = bphi + (dRdphi * Ccc.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose() + R * Ccc.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * dRdphi.transpose()) * Fc_mc3.segment<3>(D*j);
		}
		B_grad.block<3, 1>(D*i,D) = btheta;
		B_grad.block<3, 1>(D*i,D+1) = bphi;
		if (cont_slip1 < ind_c_slip.size()){
			if (i == ind_c_slip[cont_slip1]){
				A_grad.block<3, 3>(D*mc+D*cont_slip1,D*i) = ci[cont_slip1];
				A_grad.block<3, 3>(D*mc+D*cont_slip1,D*mc+D*cont_slip1) = di[cont_slip1];
				A_grad.block<3, 3>(D*i,D*mc+D*cont_slip1) = MatrixXd::Identity(D,D);
				D_grad(D,D*mc+D*cont_slip1) = (1/Ftot(2)) * Fc_mc3(D*i+2);
				D_grad(D+1,D*mc+D*cont_slip1+1) = (1/Ftot(2)) * Fc_mc3(D*i+2);
				//D_grad(D+2,D*mc+D*cont_slip1) = Fc_mc3(D*i+1);
				//D_grad(D+2,D*mc+D*cont_slip1+1) = -Fc_mc3(D*i);
				cont_slip1 = cont_slip1 + 1;
			}
		}
	}
//    A_1B = A_grad.inverse() * B_grad;
    A_1B = A_grad.colPivHouseholderQr().solve(B_grad);
	J = D_grad * A_1B;
    if (contAngle==1){
        J2 = MatrixXd::Zero(2*D-1,2*D-1);
        J2(D,D) = dRdtheta.row(0).dot(PfreeFn)/Ftot(2);
        J2(D,D+1) = dRdphi.row(0).dot(PfreeFn)/Ftot(2);
        J2(D+1,D) = dRdtheta.row(1).dot(PfreeFn)/Ftot(2);
        J2(D+1,D+1) = dRdphi.row(1).dot(PfreeFn)/Ftot(2);
        //J2(D+2,D) = dRdtheta.row(0).dot(PfreeFt2) - dRdtheta.row(1).dot(PfreeFt1);
        //J2(D+2,D+1) = dRdphi.row(0).dot(PfreeFt2) - dRdphi.row(1).dot(PfreeFt1);
        J = J + J2;
//        mexPrintf("J2 = \n");
//        for (int i = 3; i < 6; ++i){
//            for (int j = 0; j < 6; ++j){
//                mexPrintf("%g ", J2(i,j));
//            }
//            mexPrintf("\n");
//        }
    }
//	mexPrintf("J = \n");
//    for (int i = 0; i < 6; ++i){
//        for (int j = 0; j < 6; ++j){
//            mexPrintf("%g ", J(i,j));
//        }
//        mexPrintf("\n");
//    }
}

/*void GaussSeidel::StiffCart(int contAngle)
{
	int cont_slip1 = 0;
	EZ = MatrixXd::Zero(D+3,D*mc+D*mslip);
	for (int i = 0; i < mc; ++i){
		EZ.block<3, 3>(0,D*i) = Matrix3d::Identity();
		EZ(D,D*i+2) = (P_mc3(D*i+1)-Z(1));
		EZ(D+1,D*i+2) = -(P_mc3(D*i)-Z(0));
		EZ(D+2,D*i) = (-P_mc3(D*i+1)+Z(1));
		EZ(D+2,D*i+1) = (P_mc3(D*i)-Z(0));
		if (cont_slip1 < ind_c_slip.size()){
			if (i == ind_c_slip[cont_slip1]){
				EZ(D,D*mc+D*cont_slip1+1) = Fc_mc3(D*i+2);
				EZ(D+1,D*mc+D*cont_slip1) = Fc_mc3(D*i+2);
				EZ(D+2,D*mc+D*cont_slip1) = Fc_mc3(D*i+1);
				EZ(D+2,D*mc+D*cont_slip1+1) = -Fc_mc3(D*i);
				cont_slip1 = cont_slip1 + 1;
			}
		}

	}
	//MatrixXd G1;
	G1 << 1, 0, 0,      0,      0,  -Z(1),
		0, 1, 0,        0,      0,  Z(0),
		0, 0, 1,        Z(1), -Z(0), 0,
		0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 0, 1;
	Kcart = EZ * A_1B;// * G1;
    if (contAngle==1){
        J2 = MatrixXd::Zero(2*D,2*D);
        J2(D,D) = dRdtheta.row(1).dot(PfreeFn);
        J2(D,D+1) = dRdphi.row(1).dot(PfreeFn);
        J2(D,D+2) = dRdpsi.row(1).dot(PfreeFn);
        J2(D+1,D) = -dRdtheta.row(0).dot(PfreeFn);
        J2(D+1,D+1) = -dRdphi.row(0).dot(PfreeFn);
        J2(D+1,D+2) = -dRdpsi.row(0).dot(PfreeFn);
        J2(D+2,D) = dRdtheta.row(0).dot(PfreeFt2) - dRdtheta.row(1).dot(PfreeFt1);
        J2(D+2,D+1) = dRdphi.row(0).dot(PfreeFt2) - dRdphi.row(1).dot(PfreeFt1);
        J2(D+2,D+2) = dRdpsi.row(0).dot(PfreeFt2) - dRdpsi.row(1).dot(PfreeFt1);
        Kcart += J2;
    }
}*/

void GaussSeidel::StiffCart_Omega(int contAngle)
{
	EZ = MatrixXd::Zero(D+3,D*mc+D*mslip);
	A_grad = MatrixXd::Zero(D*mc+D*mslip,D*mc+D*mslip);
	B_grad = MatrixXd::Zero(D*mc+D*mslip,2*D);
	Fc_mc_hat.resize(D*mc,D);
	delta_mc_hat.resize(D*mc,D);
	delta_mc.resize(D*mc);
	int cont_slip1 = 0;
	// find A_grad
    W_mc = MatrixXd::Zero(D*mc,D*mc);
    for (int i = 0; i < mc; ++i){
        for (int j = 0; j < mc; ++j){
            W_mc.block<3, 3>(D*ind_c_cont[i],D*ind_c_cont[j]) = R * Ccc.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose();
        }
    }
	A_grad.block(0,0,D*mc,D*mc) = -W_mc;

    delta_mc = W_mc * Fc_mc3;
    for (int i = 0; i < mc; ++i){
//        cross(Fc_mc3.segment<3>(D*i), Fc_mc_hat.block<3,3>(D*i,0));
//        cross(delta_mc.segment<3>(D*i), delta_mc_hat.block<3,3>(D*i,0));
        Fc_mc_hat.block<3,3>(D*i,0) << 0, -Fc_mc3(D*i+2), Fc_mc3(D*i+1),
                    Fc_mc3(D*i+2), 0, -Fc_mc3(D*i),
                    -Fc_mc3(D*i+1), Fc_mc3(D*i), 0;
        delta_mc_hat.block<3,3>(D*i,0) << 0, -delta_mc(D*i+2), delta_mc(D*i+1),
                    delta_mc(D*i+2), 0, -delta_mc(D*i),
                    -delta_mc(D*i+1), delta_mc(D*i), 0;
    }

    B_grad.block(0,D,D*mc,D) = W_mc * Fc_mc_hat;
    B_grad.block(0,D,D*mc,D) -= delta_mc_hat;

	for (int i = 0; i < mc; ++i){
		B_grad.block<3, 3>(D*i,0) = Matrix3d::Identity();
		EZ.block<3, 3>(0,D*i) = Matrix3d::Identity();
		EZ(D,D*i+2) = (P_mc3(D*i+1)-Z(1));
		EZ(D+1,D*i+2) = -(P_mc3(D*i)-Z(0));
		EZ(D+2,D*i) = (-P_mc3(D*i+1)+Z(1));
		EZ(D+2,D*i+1) = (P_mc3(D*i)-Z(0));
        /// Finally do not consider special problem of step 1 for cartesian stiffness
        /// Once contact state obtained, compute considering nodes stuck as other steps
//		if (contAngle==1){
//            B_grad(D*i+2,D) += RPfree_mc3(D*i+1);
//            B_grad(D*i+2,D+1) -= RPfree_mc3(D*i);
//		} else {
            B_grad(D*i,D+1) += RPfree_mc3(D*i+2);
            B_grad(D*i,D+2) -= RPfree_mc3(D*i+1);
            B_grad(D*i+1,D) -= RPfree_mc3(D*i+2);
            B_grad(D*i+1,D+2) += RPfree_mc3(D*i);
            B_grad(D*i+2,D) += RPfree_mc3(D*i+1);
            B_grad(D*i+2,D+1) -= RPfree_mc3(D*i);

//		}
		if (cont_slip1 < ind_c_slip.size()){
			if (i == ind_c_slip[cont_slip1]){
				A_grad.block<3, 3>(D*mc+D*cont_slip1,D*i) = ci[cont_slip1];
				A_grad.block<3, 3>(D*mc+D*cont_slip1,D*mc+D*cont_slip1) = di[cont_slip1];
				A_grad.block<3, 3>(D*i,D*mc+D*cont_slip1) = MatrixXd::Identity(D,D);
				EZ(D,D*mc+D*cont_slip1+1) = Fc_mc3(D*i+2);
				EZ(D+1,D*mc+D*cont_slip1) = Fc_mc3(D*i+2);
				EZ(D+2,D*mc+D*cont_slip1) = Fc_mc3(D*i+1);
				EZ(D+2,D*mc+D*cont_slip1+1) = -Fc_mc3(D*i);
				cont_slip1 = cont_slip1 + 1;
			}
		}

	}
	//MatrixXd G1;
	//mexPrintf("Compute G1 \n");
	G1 << 1, 0, 0,      0,      displ(2),  -displ(1)+ZMPdes(1),
		0, 1, 0,        -displ(2),      0,   displ(0)-ZMPdes(0),
		0, 0, 1,        displ(1)-ZMPdes(1), -displ(0)+ZMPdes(0), 0,
		0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 0, 1;
//    A_1B = A_grad.inverse() * B_grad;
	//mexPrintf("Compute A^-1*B \n");
    A_1B = A_grad.householderQr().solve(B_grad);
	//mexPrintf("Compute Kcart \n");
    /// Finally do not consider special problem of step 1 for cartesian stiffness
    /// Once contact state obtained, compute considering nodes stuck as other steps
//    J2 = MatrixXd::Zero(2*D,2*D);
//    if (contAngle==1){
//		RPfreet_hat_mc = Matrix3d::Zero();
//        for (int i = 0; i < mc; ++i){
//            RPfreet_hat_mc(0,1) = -RPfree_mc3(D*i+2);
//            RPfreet_hat_mc(0,2) = RPfree_mc3(D*i+1);
//            RPfreet_hat_mc(1,0) = RPfree_mc3(D*i+2);
//            RPfreet_hat_mc(1,2) = -RPfree_mc3(D*i);
//            J2.block<3,3>(D,D) += Fc_mc_hat.block<3,3>(D*i,0) * RPfreet_hat_mc;
//		}
//    }
//	Kcart = (EZ * A_1B + J2) * G1;
	Kcart = (EZ * A_1B) * G1;
//        mexPrintf("J2 = \n");
//        for (int i = 3; i < 6; ++i){
//            for (int j = 0; j < 6; ++j){
//                mexPrintf("%g ", J2(i,j));
//            }
//            mexPrintf("\n");
//        }
}

//template <typename DerivedA,typename DerivedB>
//void cross(const EigenBase<DerivedA>& vec, EigenBase<DerivedB>& vec_hat)
//void cross(const Ref<MatrixXd>& vec, Ref<MatrixXd> vec_hat)
//{
////    vec_hat << 0, -vec(2), vec(1),
////                vec(2), 0, -vec(0),
////                -vec(1), vec(0), 0;
//    vec_hat(0,0) = 0;
//    vec_hat(0,1) = -vec(2);
//    vec_hat(0,2) = vec(1);
//    vec_hat(1,0) = vec(2);
//    vec_hat(1,1) = 0;
//    vec_hat(1,2) = -vec(0);
//    vec_hat(2,0) = -vec(1);
//    vec_hat(2,1) = vec(0);
//    vec_hat(2,2) = 0;
//}

void GaussSeidel::set_data_output(double *displ_out,double *angleact_out,double *FtotZMP_out,double *Fc_mc3_out,double *ind_cont_out,double *Kcart_out)
{
    Map<Vector5d>(FtotZMP_out,FtotZMP.rows()) = FtotZMP;
    Map<Vector3d>(displ_out,displ.rows()) = displ;
    Map<Vector2d>(angleact_out,angleact.rows()) = angleact;
    Map<VectorXd>(Fc_mc3_out,Fc_mc3.rows()) = Fc_mc3;
    for (mwSize i=0; i<mc; ++i){
        ind_cont_out[i] = ind_cont[i];
    }
    Map<MatrixXd>(Kcart_out,Kcart.rows(),Kcart.cols()) = Kcart;
}

/*********************/
