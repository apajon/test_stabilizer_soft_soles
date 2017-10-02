/*==========================================================
 * Gauss.cpp - example in MATLAB External Interfaces
 *
 * Illustrates how to use some C++ language features in a MEX-file.
//  * It makes use of member functions, constructors, destructors, and the
 * iostream.
 *
 * The routine simply defines a class, constructs a simple object,
 * and displays the initial values of the internal variables.  It
 * then sets the data members of the object based on the input given
 * to the MEX-file and displays the changed values.
 *
 * This file uses the extension .cpp.  Other common C++ extensions such
 * as .C, .cc, and .cxx are also supported.
 *
 * The calling syntax is:
 *
 *		Gauss( num1, num2 )
 *
 * Limitations:
 * On Windows, this example uses mexPrintf instead cout.  Iostreams
 * (such as cout) are not supported with MATLAB with all C++ compilers.
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2009 The MathWorks, Inc.
 *
 *========================================================*/
/* $Revision: 1.5.4.4 $ */

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

const double GaussSeidel::epsiGauss = 1e-6;
const double GaussSeidel::epsiSig = 1e-6;
const double GaussSeidel::epsiCou = 1e-6;
const int GaussSeidel::D = 3;
using namespace std;
using namespace Eigen;

extern void _main();

typedef Matrix<double,5,1> Vector5d;
typedef Matrix<double,6,1> Vector6d;

/****************************/

GaussSeidel::GaussSeidel()
{
}

GaussSeidel::~GaussSeidel(){
}

void GaussSeidel::setStaticDataSim(const MatrixXd& Ccc_in,const Vector6d& FtotZMPdes_in,double friction,int m1){

    //MatrixXd Ccc1 = read("Ccc.txt");
    Ccc = Ccc_in;
    FtotZMPdes = FtotZMPdes_in;
	/*FtotZMPdes(0) = FtotZMPdes_in[0];
	FtotZMPdes(1) = FtotZMPdes_in[1];
	FtotZMPdes(2) = FtotZMPdes_in[2];
	FtotZMPdes(3) = FtotZMPdes_in[3];
	FtotZMPdes(4) = FtotZMPdes_in[4]; */
    fric = friction;
    m = m1;// number of contact surface nodes
	FtotZMP.Zero();
	//init
	Rtheta.resize(D,D);
	Rphi.resize(D,D);
	R.resize(D,D);
	Rpfree.resize(D,m);
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
    ZMPdes << FtotZMPdes(3),FtotZMPdes(4),0.;
	//mexPrintf("FtotZMPdes = %g %g %g %g %g %g\n", FtotZMPdes(0), FtotZMPdes(1), FtotZMPdes(2), FtotZMPdes(3), FtotZMPdes(4), FtotZMPdes(5));
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
	while (sum_gauss > epsiGauss)
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
		//std::cout << *ind_slip.begin() << std::endl;
		//std::cout << F << std::endl;
        //for (int i = 0; i < D*m; ++i){
//             for (int j = 0; j < 3*m; ++j){
        //        mexPrintf("Value1 = %g\n", RtF3(i));
        //}
//     }
		ind_slip.clear();
		ind_c_slip.clear();
		ind_stick.clear();
		ind_c_stick.clear();
		ind_cont.clear();
		ind_c_cont.clear();
		//std::cout << deltaIni << std::endl;
		ite = ite + 1;
		for (int i = 0; i < m; ++i)
		{

			deltaTest = deltaIni.col(i);
            /*if (ite==1){
                mexPrintf("Value1 = %g\n", deltaTest(0));
                mexPrintf("Value1 = %g\n", deltaTest(1));
                mexPrintf("Value1 = %g\n", deltaTest(2));
            }*/
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
					while(diffFric > epsiCou && kn < 500){
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
			//std::cout << Pfree_mc3 << std::endl << std::endl;
			Pfree_mc.col(ind_c_cont[i]) = Pfree.col(ind_cont[i]);
			Fc_mc3.segment<3>(D*ind_c_cont[i]) = F3.segment<3>(D*ind_cont[i]);
			Ftot = Ftot + Fc_mc3.segment<3>(D*ind_c_cont[i]);
			for (int j = 0; j < mc; ++j){
				W_mc.block<3, 3>(D*ind_c_cont[i],D*ind_c_cont[j]) = R * Ccc.block<3, 3>(D*ind_cont[i],D*ind_cont[j]) * R.transpose();
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
		//Mzmp = Mo - ZMPdes.cross(Ftot);
		FtotZMP << Ftot(0), Ftot(1), Ftot(2), Z(0), Z(1), 0;
//        mexPrintf("FtotZMP    = %g %g %g %g %g %g\n", FtotZMP(0), FtotZMP(1), FtotZMP(2), FtotZMP(3), FtotZMP(4), FtotZMP(5));
//             for (int i = 0; i < 5; ++i){
//             for (int j = 0; j < 3*m; ++j){
  //              mexPrintf("Value1 = %g\n", FtotZMP(i));
    //         }
//     }

		//mexPrintf("SumGauss %g\n",sum_gauss);
	}
    //gradFtotZ(contAngle);
}

/*void GaussSeidel::display()
{
#ifdef _WIN32
	mexPrintf("Value1 = %g\n", val1);
	mexPrintf("Value2 = %g\n\n", val2);
#else
  cout << "Value1 = " << val1 << "\n";
  cout << "Value2 = " << val2 << "\n\n";
#endif
}*/

void GaussSeidel::set_data_output(double *displ_out,double *angleact_out,double *FtotZMP_out,double *Fc_mc3_out,double *ind_cont_out, double *Mzmp_out, mwSize n,mwSize n1)
{
    //Matrix3d Fio = R;
    //Vector5d FtotZMP1 = FtotZMP;
    for (mwSize i=0; i<6; ++i){
        FtotZMP_out[i] = FtotZMP(i);
    }


    for (mwSize i=0; i<3; ++i){
        displ_out[i] = displ(i);
        angleact_out[i] = angleact(i);
    }

    for (mwSize i=0; i<3; ++i){
        displ_out[i] = displ(i);
        angleact_out[i] = angleact(i);
        Mzmp_out[i] = Mzmp(i);
    }
    for (mwSize i=0; i<3*mc; ++i){
        Fc_mc3_out[i] = Fc_mc3(i);
    }
    for (mwSize i=0; i<mc; ++i){
        ind_cont_out[i] = ind_cont[i];
    }
    /*for (mwSize i=0; i<n; ++i){
        for (mwSize j=0; j<n1; ++j){
            z[i*n1 + j] = Fio(j,i);
        }
    }*/
}

/*********************/

static void Gauss(int contAngle,double fric,int m,Vector3d& displ_in,Vector3d& angleact_in,Vector6d& FtotZMPdes_in, const VectorXd& F3_in,const MatrixXd& Pfree_in, const MatrixXd& PabsOld_in, const MatrixXd& Ccc_in,double *displ_out,double *angleact_out,double *FtotZMPOut,double *Fc_mc3_out,double *ind_cont_out, double *Mzmp_out,mwSize n,mwSize n1)
{
/*#ifdef _WIN32
	mexPrintf("\nThe initialized data in object:\n");
#else
  cout << "\nThe initialized data in object:\n";
#endif*/
  GaussSeidel *GS = new GaussSeidel; // Create a  GaussSeidel object

  GS->setStaticDataSim(Ccc_in,FtotZMPdes_in,fric,m);

  GS->setStaticDataStep(displ_in,angleact_in,F3_in,Pfree_in,PabsOld_in);

  GS->run(contAngle);

  GS->set_data_output(displ_out,angleact_out,FtotZMPOut,Fc_mc3_out,ind_cont_out,Mzmp_out,n,n1); // Set data members to incoming

/*#ifdef _WIN32
  mexPrintf("After setting the object's data to your input:\n");
#else
  cout << "After setting the object's data to your input:\n";
#endif*/
  //GS->display();           // Make sure the set_data() worked
  delete(GS);
  flush(cout);
  return;
}

void mexFunction(int num_output, mxArray *output[], int num_input, const mxArray *input[])
{
  int contAngle;
  double fric;
  int m;                        /* input scalar */
  double *displ_in;               /* 3x1 input matrix */
  double *angleact_in;               /* 3x1 input matrix */
  double *FtotZMPdes_in;           /* 6x1 input matrix */
  double *F3_in;                    /* mx1 input matrix */
  double *Pfree_in;                    /* mx1 input matrix */
  double *PabsOld_in;                    /* mx1 input matrix */
  double *Ccc_in;
  size_t ncols;                   /* size of matrix */
  double *displ_out;              /* output matrix */
  double *angleact_out;              /* output matrix */
  double *FtotZMP_out;              /* output matrix */
  double *Fc_mc3_out;           /* ?x1 output matrix */
  double *Mzmp_out;           /* 3x1 output matrix */
  double *ind_cont_out;
  /* Check for proper number of arguments */

  if (num_input != 10) {
    mexErrMsgIdAndTxt("MATLAB:Gauss:nargin",
            "Gauss requires 10 input arguments.");
  } else if (num_output >= 7) {
    mexErrMsgIdAndTxt("MATLAB:Gauss:nargout",
            "Gauss requires 6 output argument.");
  }

  //mexPrintf("Entering Gauss Seidel in C++\n");

  contAngle = mxGetScalar(input[0]);
  fric = mxGetScalar(input[1]);
  //MatrixXd Ccc = read("Ccc.txt");

  /* get the value of the scalar input  */
  m = mxGetScalar(input[2]);

  displ_in = mxGetPr(input[3]);
  Vector3d displ_in1 = Map<Vector3d>(displ_in);

  angleact_in = mxGetPr(input[4]);
  Vector3d angleact_in1 = Map<Vector3d>(angleact_in);
  /* create a pointer to the real data in the input matrix  */

  FtotZMPdes_in = mxGetPr(input[5]);
  Vector6d FtotZMPdes_in1 = Map<Vector6d>(FtotZMPdes_in);
  //FtotZMPdes_in = mxGetPr(input[5]);

  int rows = mxGetM(input[6]); // # rows of X
  int cols = mxGetN(input[6]); // # cols of X
  F3_in = mxGetPr(input[6]);
  //input[6] = mxCreateDoubleMatrix(3*m,1,mxREAL);  //

  VectorXd F3_in1 = Map<VectorXd>(F3_in, rows);

  rows = mxGetM(input[7]); // # rows of X
  cols = mxGetN(input[7]); // # cols of X
  Pfree_in = mxGetPr(input[7]);
  MatrixXd Pfree_in1 = Map<MatrixXd>(Pfree_in, rows, cols);

  rows = mxGetM(input[8]); // # rows of X
  cols = mxGetN(input[8]); // # cols of X
  PabsOld_in = mxGetPr(input[8]);
  MatrixXd PabsOld_in1 = Map<MatrixXd>(PabsOld_in, rows, cols);
  //MexMat F3_in (mxGetPr(input[6]), N, D);

  rows = mxGetM(input[9]); // # rows of X
  cols = mxGetN(input[9]); // # cols of X
  Ccc_in = mxGetPr(input[9]);
  MatrixXd Ccc_in1 = Map<MatrixXd>(Ccc_in, rows, cols);

  /* get dimensions of the input matrix */
  //ncols = mxGetN(input[3]);
  //int ncol = 3
  /* create the output matrix */
  output[0] = mxCreateDoubleMatrix(3,1,mxREAL);  //
  /* get a pointer to the real data in the output matrix */
  displ_out = mxGetPr(output[0]);

  output[1] = mxCreateDoubleMatrix(3,1,mxREAL);  //
  /* get a pointer to the real data in the output matrix */
  angleact_out = mxGetPr(output[1]);

  output[2] = mxCreateDoubleMatrix(6,1,mxREAL);  //FtotZMP
  FtotZMP_out = mxGetPr(output[2]); //FtotZMP

  output[3] = mxCreateDoubleMatrix(3*m,1,mxREAL);  //FtotZMP
  Fc_mc3_out = mxGetPr(output[3]); //FtotZMP

  output[4] = mxCreateDoubleMatrix(m,1,mxREAL);  //FtotZMP
  ind_cont_out = mxGetPr(output[4]); //FtotZMP

  output[5] = mxCreateDoubleMatrix(3,1,mxREAL);  //Mzmp
  Mzmp_out = mxGetPr(output[5]); //FtotZMP

  Gauss(contAngle,fric,m,displ_in1,angleact_in1,FtotZMPdes_in1,F3_in1,Pfree_in1,PabsOld_in1,Ccc_in1,displ_out,angleact_out,FtotZMP_out,Fc_mc3_out,ind_cont_out,Mzmp_out,3,3);

  //mexPrintf("Exiting Gauss Seidel in C++\n");

  return;
}
