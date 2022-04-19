// Junction.h: interface for the CJunction class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_JUNCTION_H__994D7993_F853_4ECC_B268_25EDA9DD99C0__INCLUDED_)
#define AFX_JUNCTION_H__994D7993_F853_4ECC_B268_25EDA9DD99C0__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Boundary.h"
#include "Pipe.h"
#include "Properties.h"

class CJunction : public CBoundary  
{
public:
	CJunction(void);
	virtual ~CJunction();
	CJunction(const CJunction& inJunc);				// Copy constructor
	CJunction& operator=(const CJunction& inJunc);	// Overloaded operator=

	void Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rJPIPES, int** &rJPIPES_ENDS, double* &rENDCORR, 
					int i, bool ex, int npipes, std::string param_dir, std::string res_dir, int assyid, string parent_assy_res_dir, string calling_object_str);
	void ConfigureBranchAreas(CProperties* pPpt);
	void RunBoundary(CProperties* pPpt, double time, bool &rRESTORE, int timestep);
	void ReadInput(CProperties* pPpt, char *InputFile);
	void ListProperties(CProperties* pPpt);
	void PrintToScreen(CProperties* pPpt);
	void PrintToFile(double time, double ca);

	// Junction models
	// ----------------------------------------------------------------------------------------------------
	void HCPJ(CProperties* pPpt);
	void NHCPJ(CProperties* pPpt);
	bool BWPLJ(CProperties* pPpt, int timestep);
	void BWPLJContinuityNR(CProperties* pPpt, int timestep, int counter);
	void BWPLJContinuityFP(CProperties* pPpt, int timestep);
	void BWPLJLossTerms(CProperties* pPpt, int timestep, int counter, bool &rESCAPE);

	// Deprecated models
	// ----------------------------------------------------------------------------------------------------
	//void NHPLTJ(CProperties* pPpt, double time, bool &rRESTORE);
	//void NHPLJ(CProperties* pPpt, double time, bool &rRESTORE, int timestep);
	//void NHPLJ2(CProperties* pPpt, double time, bool &rRESTORE, int timestep);
	//void NHPLPC(CProperties* pPpt, double time, bool &rRESTORE);
	//void NHPLJ_Winterbone(CProperties* pPpt, double time, bool &rRESTORE);

	// Functions attaining loss coefficients
	// ----------------------------------------------------------------------------------------------------
	//void LossCoeffs_Interpolation(CProperties* pPpt, int local_flow_type, double time);
	//void LossCoeffs_Equations(CProperties* pPpt, int local_flow_type);	

	// Functions containing code common to various junction models
	// ----------------------------------------------------------------------------------------------------
	void CommonSetup();
	//void CommonContinuityFP(CProperties* pPpt, bool &rRESTORE);
	//void CommonFlowType();
	//void CommonLowFlowBranch();
	//void CommonSetBranches(CProperties* pPpt, int local_flow_type);
	//void CommonMassFlowMachNo(CProperties* pPpt);
	//void CommonEntropyLevels(CProperties* pPpt, int flow_type);
	//void CommonLambdaInStarCorrection();
	//void CommonLambdaOut(CProperties* pPpt);
	void CommonUpdate(CProperties* pPpt);

	// Inline functions
	// ----------------------------------------------------------------------------------------------------
	inline int Get_MAX_FLOW_BRANCH(){return MAX_FLOW_BRANCH;}
	inline double Get_m_dot(int BRANCH){return m_dot[BRANCH];}
	
private:

	// ====================================================================================================
	// Parameter file for Exhaust Junction [0]
	// ====================================================================================================
	
	int num_branches;				// Number of branches; must be at least 2 and agree with number determined by program

	// Model selection (in non-homentropic flow; only constant pressure model exists in homentropic flow)
	// ----------------------------------------------------------------------------------------------------
	bool CONSTP;					// Constant pressure (1=true) or pressure loss junction (0=false) (default)
	double constp_until_time;       // Run constant pressure junction irrespective of model choice until this time (s)

		// If using constant pressure model, i.e., CONSTP == 1 == true
		// ----------------------------------------------------------------------------------------------------
		double tol_A_star;				// Tolerance for A_star loop convergence

		// Else using pressure loss model, i.e., CONSTP == 0 == false
		// ----------------------------------------------------------------------------------------------------
	
		// Branch angles (used only in pressure loss models)
		// ---------------------------------------------------------------------------------------------------- 
		int* id_branch;					// Branch IDs
		double* ref0_branch;			// Rotation angles (degrees) within x-y plane with respect to first reference pipe (branch [0])
		double* ref1_branch;			// Rotation angle (degrees) out of x-y plane with respect to second reference pipe (branch [1])
	
		// Algorithm control
		// ----------------------------------------------------------------------------------------------------
		double tol_mass_sum;			// Tolerance for continuity loop convergence
		double tol_del_star;			// Tolerance for pressure loss loop convergence

	// Deprecated
	// ----------------------------------------------------------------------------------------------------
	//int type;						// Type of junction in use (HCPJID/NHCPJID/NHPLTJID/NHPLJID/NHPLPCID)
	//double angle_deg;				// Branch angle (degrees)
	//double tol_main;				// Tolerance on main loop convergence
	//double tol_cont;				// Tolerance on continuity loop convergence
	//double loop_limit_cont;			// Loop limit; number of continuity loops before continuing
	double loop_limit_switch;		// Loop limit; number of main loops before continuing if flow type is continuously alternating
	
	// ====================================================================================================
	// End of file
	// ====================================================================================================

	// Common variables (irrespective of model)
	// ----------------------------------------------------------------------------------------------------
	int MAX_FLOW_BRANCH;			// Branch carrying greatest positive (towards junction) mass flow

	// Configuration details
	// ---------------------
	double* Fb;					// Branch pipe cross-sectional areas (m^2)

	// Calculation variables
	// ---------------------
	double* lambda_in_star;		// Working lambda_in_stars
	double* lambda_in_star_n;	// Uncorrected lambda_in_stars
	double* lambda_out;			// Corrected lambda_outs
	double* AA;					// Working entropy levels
	double* AA_n;				// Uncorrected entropy levels
	double* A_star;				// Starred N.D. speed of sound
	double* U_star;				// Starred N.D. velocity
	double* del_A_star;			// Pressure difference terms
	double* del_A_star_old;		// Pressure difference terms (old)

	// BWPLJ variables
	// ---------------
	double** ref_datum_branch;	// Stores calculated angle (radians) between each branch and every other branch
	double* C;					// Loss coefficients
	double* del_star;			// Pressure loss terms
	int DATUM;					// Datum branch - carries greatest positive (towards junction) mass flow
	double h_0_sep;				// Calculated enthalpy of all separating flows
	
	// Flow type variables
	// -------------------
	int pipe1;					// Pipe '1'=[0]
	int pipe2;					// Pipe '2'=[1]
	int pipe3;					// Pipe '3'=[2]
	int bpipe1;					// Benson '1' = my [0]
	int bpipe2;					// Benson '2' = my [2]
	int bpipe3;					// Benson '3' = my [1]

	int** opposite;				// Stores which "opposite" flowtype to pick given the actual flowtype and the low_flow_branch
	int low_flow_branch;		// Branch in which the magnitude of the flow is smallest
	int next_lowest_branch;		// Branch with the second lowest flow magnitude

	int a, b;					// Enumerations for non-com branches
	int com, i2, i3;			// Branch labels: common, other 'i2', other 'i3'
	int* pipe_flow;				// Branch flow directions
	int* pipe_flow_old;			// Branch flow directions (old)
	int* coeff;					// Which loss coefficient to lookup [La] [Lb]; used in interpolation (NHPLPC only)
	int flow_type_winterbone, flow_type_winterbone_old, flow_type_winterbone_orig; // Standard flow type numbering
	int flow_type_winterbone_next;
	int flow_type_actual;		// The flowtype found in fresh iterations by the flowtype test
	int flow_type_benson;		// Benson flow type numbering differs from that of Winterbone
	double *G1, *G2;			// Functions G1, G2 in Benson's 'T' junction pressure difference equation (NHPLTJ only)
	double C1,C2,C3,C4,C5,C6;	// Empirical loss coefficients for functions G1, G2 in Benson's 'T' junction (NHPLTJ only)
	double no_flow;				// The values of U_star below which a flow type of 0 is selected

	// Mass flow & Mach No. variables
	// ------------------------------
	double m_dot_total;			// Resultant mass flow rate (kg/s)
	double* m_dot;				// Branch mass flow rates (kg/s)
	double Mn;					// Common branch Mach No.
	double* W;					// Mass flow ratios, for other pipe 'a' and 'b'
	
	// Loss coefficient variables
	// --------------------------
	double** K_loss;			// Loss coefficients K, with labels, e.g. K12
	int label, value;			// Enumerations for K label (e.g. K12) and K value
	double psiT;				// com area/branch area
	//double angle_rad;			// Branch angle (radians)
	double* L;					// Loss coefficients L, for other pipe 'a' and 'b'
	double* x_star;
	double* x_star_old;

	// Algorithm control
	// -----------------
	double loop_limit_main;		// Loop limit; number of main loops before further action
//	double loop_limit_same;		// Loop limit; number of main loops before continuing if stuck in one flowtype
	
	// NHPLJ only:
	bool ALTERNATING;			// Flag; when true signals alternating flow types
	bool FIX_FLOW_TYPE;			// Flag; when true sets original flow type to prevent alternating flow types
	bool OPPOSITE;				// Flag; when true sets 'opposite' flow type depending on the low flow branch
	bool FIRST_OPPOSITE;		// Flag; signals first running of the code which sets the 'opposite' flow type
	bool OTHER;					// Flag; when true sets other flow types, irrespective of actual flow type
	bool ZERO_FLOW;				// Flag; when true sets a zero flow in the branch with the lowest flow magnitude

	bool TESTED_ACTUAL;			// Flag; when true indicates low has been tested with a fixed flow_type_actual
	bool TESTED_ORIGINAL;		// Flag; when true indicates low has been tested with a fixed flow_type_winterbone_orig
	
	int switch_counter;			// Counts number of consecutive changes in flow type
	int fix_counter;			// Counts number of consecutive attempts under FIX_FLOW_TYPE
	int other_counter;			// Counts number of consecutive attempts under OTHER
	int zero_counter;			// Counts number of consecutive attempts under ZERO_FLOW

	int normal_counter_total;	// Number of simulation iterations completed under normal operation
	int fix_counter_total;		// Number of simulation iterations completed under FIX_FLOW_TYPE
	int opposite_counter_total;	// Number of simulation iterations completed under OPPOSITE
	int other_counter_total;	// Number of simulation iterations completed under OTHER
	int zero_counter_total;		// Number of simulation iterations completed under ZERO_FLOW

	bool nhplj_first_loop;		// Flag
	bool first_loop;			// Flag; first execution of loop for this function call
	bool first_loop_ever;		// Flag; first execution of loop during overall simulation
	
	double error1, error2;
	double error1_old, error2_old;

	// HCPJ only variables
	// -------------------
	bool hcpj_first_loop;		// Flag
	double FT;					// Summation of branch areas	
	double* Kq;					// Summation of values
	
	// Other variables
	// ---------------
	char pause;

public:
	//	bool rollback;				// Signals to main function to apply roll back
	bool rolledback;			// Records whether roll back has been applied
	bool stoprollback;			// Signals when roll back should be stopped
};
#endif // !defined(AFX_JUNCTION_H__994D7993_F853_4ECC_B268_25EDA9DD99C0__INCLUDED_)
