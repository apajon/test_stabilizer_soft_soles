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
#include "Gausstot.h"
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

extern void _main();

typedef Matrix<double,5,1> Vector5d;
typedef Matrix<double,6,1> Vector6d;

static void Gauss(int contAngle,double fric,int m,Vector3d& displ_in,Vector3d& angleact_in,Vector6d& FtotZMPdes_in, const VectorXd& F3_in,const MatrixXd& Pfree_in, const MatrixXd& PabsOld_in, const MatrixXd& Ccc_in,double *displ_out,double *angleact_out,double *FtotZMPOut,double *Fc_mc3_out,double *Pc_mc3_out,double *ind_cont_out,double *Kcart_out,double *J_out)
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

  GS->set_data_output(displ_out,angleact_out,FtotZMPOut,Fc_mc3_out,Pc_mc3_out,ind_cont_out,Kcart_out,J_out); // Set data members to incoming

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
  int m;                       
  double *displ_in;               
  double *angleact_in;               
  double *FtotZMPdes_in;          
  double *F3_in;                    
  double *Pfree_in;                    
  double *PabsOld_in;                    
  double *Ccc_in;
  double *displ_out;             
  double *angleact_out;             
  double *FtotZMP_out;          
  double *Fc_mc3_out;           
  double *Pc_mc3_out;           
  double *ind_cont_out;
  double *Kcart_out;
  double *J_out;
  /* Check for proper number of arguments */

  if (num_input != 10) {
    mexErrMsgIdAndTxt("MATLAB:Gauss:nargin",
            "Gauss requires 10 input arguments.");
  } else if (num_output >= 9) {
    mexErrMsgIdAndTxt("MATLAB:Gauss:nargout",
            "Gauss requires 8 output argument.");
  }

  //mexPrintf("Entering Gauss Seidel in C++\n");

  contAngle = (int)mxGetScalar(input[0]);
  fric = mxGetScalar(input[1]);
  //MatrixXd Ccc = read("Ccc.txt");

  /* get the value of the scalar input  */
  m = (int)mxGetScalar(input[2]);

  displ_in = mxGetPr(input[3]);
  Vector3d displ_in1 = Map<Vector3d>(displ_in);

  angleact_in = mxGetPr(input[4]);
  Vector3d angleact_in1 = Map<Vector3d>(angleact_in);
  /* create a pointer to the real data in the input matrix  */

  FtotZMPdes_in = mxGetPr(input[5]);
  Vector6d FtotZMPdes_in1 = Map<Vector6d>(FtotZMPdes_in);
  //FtotZMPdes_in = mxGetPr(input[5]);

  size_t rows = mxGetM(input[6]); // # rows of X
  size_t cols = mxGetN(input[6]); // # cols of X
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
  displ_out = mxGetPr(output[0]);

  output[1] = mxCreateDoubleMatrix(3,1,mxREAL);  //
  angleact_out = mxGetPr(output[1]);

  output[2] = mxCreateDoubleMatrix(6,1,mxREAL);  
  FtotZMP_out = mxGetPr(output[2]); 

  output[3] = mxCreateDoubleMatrix(3*m,1,mxREAL);  
  Fc_mc3_out = mxGetPr(output[3]);

  output[4] = mxCreateDoubleMatrix(3*m,1,mxREAL);
  Pc_mc3_out = mxGetPr(output[4]);
  
  output[5] = mxCreateDoubleMatrix(m,1,mxREAL);
  ind_cont_out = mxGetPr(output[5]);

  output[6] = mxCreateDoubleMatrix(6,6,mxREAL);
  Kcart_out = mxGetPr(output[6]);

  output[7] = mxCreateDoubleMatrix(6,6,mxREAL);
  J_out = mxGetPr(output[7]);  
  
  Gauss(contAngle,fric,m,displ_in1,angleact_in1,FtotZMPdes_in1,F3_in1,Pfree_in1,PabsOld_in1,Ccc_in1,displ_out,angleact_out,FtotZMP_out,Fc_mc3_out,Pc_mc3_out,ind_cont_out,Kcart_out,J_out);

  //mexPrintf("Exiting Gauss Seidel in C++\n");

  return;
}
