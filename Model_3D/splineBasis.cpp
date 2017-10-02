#include <unsupported/Eigen/Splines>
#include <Eigen/StdVector>
#include "mex.h"

const int DEGREE = 3;

typedef Eigen::Matrix<int, DEGREE+1, 1> VectorDi;
typedef Eigen::Matrix<double, DEGREE+1, DEGREE+1> MatrixDd;

struct PreComputedSpline
{
  VectorDi idx_u;
  VectorDi idx_v;
  MatrixDd basisValues;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

typedef std::vector<PreComputedSpline, Eigen::aligned_allocator<PreComputedSpline> > PreComputedSplineVector;

//basis for 2d uniform clamped spline for each (u,v) = UVR(:,i)
void computeBasis(PreComputedSplineVector& out, const Eigen::MatrixXd& UVR, const Eigen::ArrayXd& ku, const Eigen::ArrayXd& kv)
{
  typedef Eigen::Spline<double, 1, DEGREE> SplineD;
  
  Eigen::DenseIndex mu = ku.size() + 2*DEGREE;
  Eigen::ArrayXd knots_u(mu);
  knots_u << Eigen::ArrayXd::Constant(DEGREE, ku(0)), ku, Eigen::ArrayXd::Constant(DEGREE, ku(ku.size()-1));
  SplineD su(knots_u, Eigen::ArrayXd::Zero(mu));

  Eigen::DenseIndex mv = kv.size() + 2*DEGREE;
  Eigen::ArrayXd knots_v(mv);
  knots_v << Eigen::ArrayXd::Constant(DEGREE, kv(0)), kv, Eigen::ArrayXd::Constant(DEGREE, kv(kv.size()-1));
  SplineD sv(knots_v, Eigen::ArrayXd::Zero(mv));

  const Eigen::DenseIndex n = UVR.cols();
  out.resize(n);
  for (Eigen::DenseIndex i = 0; i < n; ++i)
  {
    double u = UVR(0, i);
    double v = UVR(1, i);
	double r = UVR(2, i);

    int iu = static_cast<int>(su.span(u));
    out[i].idx_u = VectorDi::LinSpaced(iu - DEGREE, iu + 1);
    SplineD::BasisVectorType bu = su.basisFunctions(u);

    int iv = static_cast<int>(sv.span(v));
    out[i].idx_v = VectorDi::LinSpaced(iv - DEGREE, iv + 1);
    SplineD::BasisVectorType bv = sv.basisFunctions(v);

    out[i].basisValues = bu.matrix().transpose()*bv.matrix();
  }
}

//Gateway function
void mexFunction(int num_output, mxArray *output[],
                  int num_input, const mxArray *input[])
{
  if (num_input != 3)
    mexErrMsgTxt("Invalid number of inputs. 3 inputs are required.");
  else if (num_output>1)
    mexErrMsgTxt("Too many outut arguments.");
  
  if (!mxIsDouble(input[0]))
    mexErrMsgTxt("First input argument must be a matrix.");
  size_t m0 = mxGetM(input[0]);
  size_t n0 = mxGetN(input[0]);
  if (m0 != 3)
    mexErrMsgTxt("First input argument must have 3 rows.");
  
  if (!mxIsDouble(input[1]))
    mexErrMsgTxt("Second input argument must be a vector.");
  size_t m1 = mxGetM(input[1]);
  size_t n1 = mxGetN(input[1]);
  if (m1 != 1 && n1 != 1)
    mexErrMsgTxt("Second input argument must be a row or column vector.");
  size_t s1 = std::max(m1,n1);
  
  if (!mxIsDouble(input[2]))
    mexErrMsgTxt("Second input argument must be a vector.");
  size_t m2 = mxGetM(input[2]);
  size_t n2 = mxGetN(input[2]);
  if (m2 != 1 && n2 != 1)
    mexErrMsgTxt("Second input argument must be a row or column vector.");
  size_t s2 = std::max(m2,n2);
  
  Eigen::MatrixXd arg0 = Eigen::Map<const Eigen::MatrixXd>(mxGetPr(input[0]),m0,n0);
  Eigen::ArrayXd arg1 = Eigen::Map<const Eigen::ArrayXd>(mxGetPr(input[1]),s1);
  Eigen::ArrayXd arg2 = Eigen::Map<const Eigen::ArrayXd>(mxGetPr(input[2]),s2);
  
  //call routine
  PreComputedSplineVector out;
  computeBasis(out, arg0, arg1, arg2);
  
  //output
  mxArray * out_struct;
  int num_info_fields = 6; 
  const char *info_field_names[] = {"idx_u", "idx_v", "basisValues", "u", "v", "r"}; 
  out_struct = mxCreateStructMatrix(n0, 1, num_info_fields, info_field_names);
  
  for (size_t i=0; i<n0; ++i)
  {
    mxArray* idx_u = mxCreateNumericMatrix(1, DEGREE+1, mxINT32_CLASS, mxREAL);
    for (int j=0; j<static_cast<int>(DEGREE)+1; ++j)
      ((INT32_T *)mxGetData(idx_u))[j] = out[i].idx_u[j] + 1;
    mxSetField(out_struct, i, "idx_u", idx_u);
    
    mxArray* idx_v = mxCreateNumericMatrix(1, DEGREE+1, mxINT32_CLASS, mxREAL);
    for (int j=0; j<static_cast<int>(DEGREE)+1; ++j)
      ((INT32_T *)mxGetData(idx_v))[j] = out[i].idx_v[j] + 1;
    mxSetField(out_struct, i, "idx_v", idx_v);
    
    mxArray* basis = mxCreateDoubleMatrix(DEGREE+1, DEGREE+1, mxREAL);
    Eigen::Map<MatrixDd>(mxGetPr(basis)) = out[i].basisValues;
    mxSetField(out_struct, i, "basisValues", basis);

	mxSetField(out_struct, i, "u", mxCreateDoubleScalar(arg0(0,i)));
	mxSetField(out_struct, i, "v", mxCreateDoubleScalar(arg0(1,i)));
	mxSetField(out_struct, i, "r", mxCreateDoubleScalar(arg0(2,i)));
  }
  
  output[0] = out_struct;
}