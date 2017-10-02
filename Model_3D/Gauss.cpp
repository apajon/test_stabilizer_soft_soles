/*==========================================================
 * Gauss.cpp - Gauss-Seidel contact for fixed orientation and displacement
 *==========================================================
*/

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

using namespace std;
using namespace Eigen;

typedef Matrix<double,5,1> Vector5d;
typedef Matrix<double,6,1> Vector6d;

static void Gauss(int contAngle,double fric,int m,Vector3d& displ_in,Vector3d& angleact_in,const VectorXd& F3_in,const MatrixXd& Pfree_in, const MatrixXd& PabsOld_in, const MatrixXd& Cs_in,double *displ_out,double *angleact_out,double *FtotZMPOut,double *Fc_mc3_out,double *ind_cont_out,double *ind_slip_out,double *Mzmp_out,double *PSurf3)
{
  GaussSeidel *GS = new GaussSeidel; // Create a  GaussSeidel object

  GS->setStaticDataSim(Cs_in,fric,m);

  GS->setStaticDataStep(displ_in,angleact_in,F3_in,Pfree_in,PabsOld_in);

  GS->run(contAngle);

  GS->set_data_output(displ_out,angleact_out,FtotZMPOut,Fc_mc3_out,ind_cont_out,ind_slip_out,Mzmp_out,PSurf3);

  delete(GS);
  flush(cout);
  return;
}

void mexFunction(int num_output, mxArray *output[], int num_input, const mxArray *input[])
{
  int contAngle;
  double fric;
  int m;
  double *displ_in;
  double *angleact_in;
  double *F3_in;
  double *Pfree_in;
  double *PabsOld_in;
  double *Cs_in;
  size_t ncols;
  double *displ_out;
  double *angleact_out;
  double *FtotZMP_out;
  double *Fc_mc3_out;
  double *ind_cont_out;
  double *ind_slip_out;
  double *Mzmp_out; 
  double *PSurf3;
  /* Check for proper number of arguments */
  if (num_input != 9) {
    mexErrMsgIdAndTxt("MATLAB:Gauss:nargin",
            "Gauss requires 8 input arguments.");
  } else if (num_output >= 9) {
    mexErrMsgIdAndTxt("MATLAB:Gauss:nargout",
            "Gauss requires 8 output argument.");
  }

  contAngle = mxGetScalar(input[0]);
  fric = mxGetScalar(input[1]);

  m = mxGetScalar(input[2]);

  displ_in = mxGetPr(input[3]);
  Vector3d displ_in1 = Map<Vector3d>(displ_in);

  angleact_in = mxGetPr(input[4]);
  Vector3d angleact_in1 = Map<Vector3d>(angleact_in);

  int rows = mxGetM(input[5]);
  int cols = mxGetN(input[5]);
  F3_in = mxGetPr(input[5]);

  VectorXd F3_in1 = Map<VectorXd>(F3_in, rows);

  rows = mxGetM(input[6]);
  cols = mxGetN(input[6]);
  Pfree_in = mxGetPr(input[6]);
  MatrixXd Pfree_in1 = Map<MatrixXd>(Pfree_in, rows, cols);

  rows = mxGetM(input[7]);
  cols = mxGetN(input[7]);
  PabsOld_in = mxGetPr(input[7]);
  MatrixXd PabsOld_in1 = Map<MatrixXd>(PabsOld_in, rows, cols);

  rows = mxGetM(input[8]);
  cols = mxGetN(input[8]);
  Cs_in = mxGetPr(input[8]);
  MatrixXd Cs_in1 = Map<MatrixXd>(Cs_in, rows, cols);

  output[0] = mxCreateDoubleMatrix(3,1,mxREAL);
  displ_out = mxGetPr(output[0]);

  output[1] = mxCreateDoubleMatrix(3,1,mxREAL);
  angleact_out = mxGetPr(output[1]);

  output[2] = mxCreateDoubleMatrix(6,1,mxREAL);
  FtotZMP_out = mxGetPr(output[2]);

  output[3] = mxCreateDoubleMatrix(3*m,1,mxREAL);
  Fc_mc3_out = mxGetPr(output[3]);

  output[4] = mxCreateDoubleMatrix(m,1,mxREAL);
  ind_cont_out = mxGetPr(output[4]);

  output[5] = mxCreateDoubleMatrix(m,1,mxREAL);
  ind_slip_out = mxGetPr(output[5]); 
  
  output[6] = mxCreateDoubleMatrix(3,1,mxREAL);
  Mzmp_out = mxGetPr(output[6]);
  
  output[7] = mxCreateDoubleMatrix(3*m,1,mxREAL);
  PSurf3 = mxGetPr(output[7]);  

  Gauss(contAngle,fric,m,displ_in1,angleact_in1,F3_in1,Pfree_in1,PabsOld_in1,Cs_in1,displ_out,angleact_out,FtotZMP_out,Fc_mc3_out,ind_cont_out,ind_slip_out,Mzmp_out,PSurf3);
  
  return;
}
