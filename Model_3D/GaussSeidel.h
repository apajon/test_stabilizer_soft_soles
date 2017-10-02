/*
 * GaussSeidel.h
 */

#include <Eigen/Dense>
#include <Eigen/StdVector>
// #ifndef _GAUSS_SEIDEL_H_
// #define _GAUSS_SEIDEL_H_
//#include "mex.h"

using namespace std;
using namespace Eigen;

class GaussSeidel {
    public:
	  typedef std::vector<Matrix3d,aligned_allocator<Matrix3d> > MatrixArray;
	  typedef Matrix<double,5,1> Vector5d;
	  typedef Matrix<double,6,1> Vector6d;
	  typedef Matrix<double,5,5> Matrix5d;

      GaussSeidel();
      ~GaussSeidel();

      void set_data_output(double *displ_out,double *angleact_out,double *FtotZMPOut,double *Fc_mc3_out,double *ind_cont_out,double *ind_slip_out,double *Mzmp_out,double *PSurf3);
      void setStaticDataSim(const MatrixXd& Cs_in,double friction,int m);
      void setStaticDataStep(const Vector3d& displ_in,const Vector3d& angleact_in,const VectorXd& F3_in,const MatrixXd& Pfree_in,const MatrixXd& PabsOld_in);
	  void run(int contAngle);

    private:
    	//Variables
		int m;
		int mslip;
		int mstick;
		int mc;
		int contAngle;
		const static double epsiGauss; // Desired precision
        const static double epsiZMP;
		const static double epsiSig; // Signorini's law tolerance
		const static double epsiCou; // Coulomb contact precision
		const static int D; //Dimension
		double fric;
		double theta;
		double phi;
		double psi;
        double critFZMP;
		MatrixXd J;
		MatrixXd J2;
		Matrix3d RPfreet_hat_mc;
		Vector6d FtotZMP;
        Vector6d FtotZMPdes;
        Vector6d FtotZMPerr;
		Vector3d displ;
        Vector3d displ_ini;
        Vector3d angleact;
		MatrixXd displ_m;
        MatrixXd displ_mc3;
		MatrixXd displ_ini_m;
		Vector6d displangle;
        Vector6d delta_displangle;
		MatrixXd Cs;
		VectorXd F3;
		VectorXd Fold3;
        Matrix3d R;
        Matrix3d Rini;
		Matrix3d Rtheta;
		Matrix3d Rphi;
		Matrix3d Rpsi;
		Matrix3d Rpsiini;
		MatrixXd PabsOld;
		MatrixXd RpabsOld;
		MatrixXd Pfree;
		MatrixXd RPfree;
		MatrixXd RiniPfree;
        VectorXd RPfree3;
        VectorXd P3;
		MatrixArray ci;
		MatrixArray di;
		MatrixArray W;
		MatrixArray inv_W;
		MatrixXd RtF;
		VectorXd RtF3;
		MatrixXd deltaIni;
		VectorXd qt;
		VectorXd qn;
		Vector3d deltaTest;
		Matrix3d A;
		Matrix2d A22;
		double ad;
		double ad12;
		double lambda_min;
		int i3;
		double norm_f_diff;
		double norm_f_new;
		double lambda_max;
		Matrix3d B;
		std::vector<int> ind_slip;
		std::vector<int> ind_c_slip;
		std::vector<int> ind_stick;
		std::vector<int> ind_c_stick;
		std::vector<int> ind_cont;
		std::vector<int> ind_c_cont;
		Vector3d Fcfric;
		Vector3d deltaNew;
		Matrix2d Bigpi;
		Matrix2d alpha;
		Matrix2d dsMat;
		Vector3d phi1;
		Vector3d d2;
		Vector3d d1;
		Vector3d phi2;
		Vector2d ds;
		Matrix3d G21;
		Matrix3d G22;
		double diffFric;
		int kn;
		double fricTest;
		Vector3d Ftot;
		Vector3d deltaInt;
		Vector3d Mo;
		Vector2d Z;
		Vector3d ZMPdes;
		Vector3d Mzmp;

		MatrixXd Pfree_mc;
		MatrixXd W_mc;
		VectorXd Pfree_mc3;
		MatrixXd RPfree_mc;
		VectorXd RPfree_mc3;
		VectorXd Fc_mc3;
		VectorXd RtFc_mc3;
		VectorXd P_mc3;
		VectorXd P_mc3_bis;

		MatrixXd A_grad;
		MatrixXd B_grad;
        MatrixXd B_omega;
		MatrixXd D_grad;
		MatrixXd Fc_mc_hat;
		MatrixXd delta_mc_hat;
		VectorXd delta_mc;
        MatrixXd A_1B;
        MatrixXd Xi;
        MatrixXd G;
        MatrixXd Kcart;
		Vector3d btheta;
		Vector3d bphi;
		Vector3d bpsi;
		Vector3d PfreeFt1;
		Vector3d PfreeFt2;
		Vector3d PfreeFn;
		Matrix3d dRdtheta;
		Matrix3d dRdphi;
		Matrix3d dRdpsi;
		Matrix3d dRinidtheta;
		Matrix3d dRinidphi;
};
