#include <iostream>
#include <math.h>
#include "mex.h"
#include "GaussSeidel.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <fstream>

const double GaussSeidel::epsiGauss = 1e-12;
const double GaussSeidel::epsiSig = 1e-12;
const double GaussSeidel::epsiCou = 1e-12;
const int GaussSeidel::D = 3;
using namespace std;
using namespace Eigen;

extern void _main();

typedef Matrix<double,5,1> Vector5d;
typedef Matrix<double,6,1> Vector6d;

/****************************/

GaussSeidel::GaussSeidel(){}

GaussSeidel::~GaussSeidel(){}

void GaussSeidel::setStaticDataSim(const MatrixXd& Cs_in,double friction,int m1){
    Cs = Cs_in;
    fric = friction;
    m = m1;
	FtotZMP.Zero();
	Rtheta.resize(D,D);
	Rphi.resize(D,D);
	R.resize(D,D);
	RPfree.resize(D,m);
	RpabsOld.resize(D,m);
	displ_ini_m.resize(D,m);
	displ_m.resize(D,m);
	W.resize(m);
	inv_W.resize(m);
	ci.resize(m);
	di.resize(m);
	displangle.resize(6);
	deltaIni.resize(D,m);
	qt.resize(m);
	qn.resize(m);
	RtF.resize(D,m);
	ind_slip.reserve(m);
	ind_c_slip.reserve(m);
	ind_stick.reserve(m);
	ind_c_stick.reserve(m);
	ind_cont.reserve(m);
	ind_c_cont.reserve(m);
    P3.resize(D*m);
}

void GaussSeidel::setStaticDataStep(const Vector3d& displ_in,const Vector3d& angleact_in,const VectorXd& F3_in,const MatrixXd& Pfree_in,const MatrixXd& PabsOld_in){
	displ = displ_in;
    angleact = angleact_in;
    F3 = F3_in;
    PabsOld = PabsOld_in;
    Pfree = Pfree_in;
}

void GaussSeidel::run(int contAngle){
	theta = angleact(0);
	phi = angleact(1);
	psi = angleact(2);
	displangle << displ(0),displ(1),displ(2),theta,phi,psi;
	int ite = 0;
	double sum_gauss = 1;
	//// Compute Rotation
	Rtheta << 1, 0, 0,
		    0, cos(theta), -sin(theta),
			0, sin(theta), cos(theta);
	Rphi << cos(phi), 0, sin(phi),
		    0, 1, 0,
			-sin(phi), 0, cos(phi);
	Rpsi << cos(psi), -sin(psi), 0,
		    sin(psi), cos(psi), 0,
			0, 0, 1;
	R = Rtheta * Rphi * Rpsi;
	RPfree = R * Pfree;
	//RpabsOld = R * PabsOld;
	for (int i = 0; i < m; ++i)
	{
		i3 = 3 * i;
		A = R * Cs.block<3, 3>(i3, i3) * R.transpose();
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
        ///Later define displ_ini and psiini as inputs, when optimizing foot position and orientation
		double psiini = 0.;
        Rpsiini << cos(psiini), -sin(psiini), 0,
                sin(psiini), cos(psiini), 0,
                0, 0, 1;
        Rini = Rtheta * Rphi * Rpsiini;
        RiniPfree = Rini * Pfree;
		displ_ini_m = displ_ini.replicate(1,m);
		displ_m = displ.replicate(1,m);
		deltaIni.row(0) = displ_m.row(0) + RPfree.row(0) - displ_ini_m.row(0) - RiniPfree.row(0);
		deltaIni.row(1) = displ_m.row(1) + RPfree.row(1) - displ_ini_m.row(1) - RiniPfree.row(1);
		deltaIni.row(2) = displ_m.row(2) + RPfree.row(2);
	}
	else{
		displ_m = displ.replicate(1,m);
		deltaIni.row(0) = RPfree.row(0) - PabsOld.row(0) + displ_m.row(0);
		deltaIni.row(1) = RPfree.row(1) - PabsOld.row(1) + displ_m.row(1);
		deltaIni.row(2) = RPfree.row(2) + displ_m.row(2);
	}
    while(sum_gauss > epsiGauss)
	{
		int cont_stick = 0;
		int cont_slip = 0;
		int cont_c = 0;
		Fold3 = F3;
		Map<MatrixXd> F = Map<MatrixXd>(F3.data(),D,m);

		RtF = R.transpose() * F;
		RtF3 = Map<VectorXd>(RtF.data(),D*m);
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

			deltaTest = deltaTest  + R*(Cs.middleCols<3>(i3).transpose()*RtF3 - Cs.block<3, 3>(i3, i3)*RtF3.segment<3>(i3));
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
                    //Vector2d acc;
                    //acc = deltaNew.head<2>().norm()*Fcfric.head<2>() + fric * Fcfric(2) * deltaNew.head<2>();
                    //mexPrintf("FtotZMP = %g\n", acc(1));
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

		mslip = ind_slip.size();
		mstick = ind_stick.size();
		mc = ind_cont.size();
		Pfree_mc.resize(D,mc);
		Pfree_mc3.resize(D*mc);
		Fc_mc3.resize(D*mc);
		Ftot = Vector3d::Zero();
		W_mc = MatrixXd::Zero(D*mc,D*mc);
		for (int i = 0; i < mc; ++i){
			Pfree_mc3.segment<3>(D*ind_c_cont[i]) = Pfree.col(ind_cont[i]);
			Pfree_mc.col(ind_c_cont[i]) = Pfree.col(ind_cont[i]);
			Fc_mc3.segment<3>(D*ind_c_cont[i]) = F3.segment<3>(D*ind_cont[i]);
			Ftot = Ftot + Fc_mc3.segment<3>(D*ind_c_cont[i]);
			for (int j = 0; j < mc; ++j){
				W_mc.block<3, 3>(D*ind_c_cont[i],D*ind_c_cont[j]) = R * Cs.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose();
			}
		}

		RPfree_mc = R * Pfree_mc;
		displ_mc3 = displ.replicate(mc, 1);
		RPfree_mc3 = Map<VectorXd>(RPfree_mc.data(),D*mc);

		P_mc3 = displ_mc3 + RPfree_mc3 + W_mc * Fc_mc3;
		Mo = Vector3d::Zero();
		for (int i = 0; i < mc; ++i){
			Mo = Mo + P_mc3.segment<3>(D*i).cross(Fc_mc3.segment<3>(D*i));
		}
		Z << -Mo(1)/Ftot(2), Mo(0)/Ftot(2);
		FtotZMP << Ftot(0), Ftot(1), Ftot(2), Z(0), Z(1), 0;
        
        //Find position of the surface points
        RPfree = R * Pfree;
        RPfree3 = Map<VectorXd>(RPfree.data(),D*m);
        for (int i = 0; i < m; ++i){
            deltaInt = Vector3d::Zero();
            for (int j = 0; j < m; ++j){
                deltaInt += Cs.block<3, 3>(D*i,D*j) * RtF3.segment<3>(D*j);
            }
            P3.segment<3>(D*i) = displ + RPfree3.segment<3>(D*i) + R*deltaInt;
        }     
	}
}

void GaussSeidel::set_data_output(double *displ_out,double *angleact_out,double *FtotZMP_out,double *Fc_mc3_out,double *ind_cont_out,double *ind_slip_out,double *Mzmp_out,double *PSurf3)
{
    Map<Vector6d>(FtotZMP_out,FtotZMP.rows()) = FtotZMP;
    Map<Vector3d>(displ_out,displ.rows()) = displ;
    Map<Vector3d>(angleact_out,angleact.rows()) = angleact;
    Map<VectorXd>(Fc_mc3_out,Fc_mc3.rows()) = Fc_mc3;
    Map<Vector3d>(Mzmp_out,Mzmp.rows()) = Mzmp;
    for (mwSize i=0; i<mc; ++i){
        ind_cont_out[i] = ind_cont[i];
    }
    for (mwSize i=0; i<mslip; ++i){
        ind_slip_out[i] = ind_slip[i];
    }
    Map<VectorXd>(PSurf3,P3.rows()) = P3;
}