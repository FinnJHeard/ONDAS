// APLDev.h: interface for the CAPLDev class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_APLDEV_H__87820053_A2AE_415A_B790_5D4F8BD891B5__INCLUDED_)
#define AFX_APLDEV_H__87820053_A2AE_415A_B790_5D4F8BD891B5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Boundary.h"
#include "EndEnvironment.h"
#include "Junction.h"
//#include "Node.h"
//#include "PathLine.h"
//#include "Pipe.h"
#include "Properties.h"
#include "Time.h"
#include "Transmissive.h"

class CAPLDev : public CBoundary  
{
public:
	CAPLDev();
	virtual ~CAPLDev();

	void Initialise(CTime* pMyTime, CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, /*int** &rAPLDevASSYS,*/ int** &rAPLDevPIPES, int** &rAPLDevPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, CTransmissive* TransmPtr, std::string param_dir, int assyid, bool transm_ptr_exists, string parent_assy_res_dir, CEndEnvironment* &rEndEnv, CTransmissive* &rTransm, CJunction* &rJunc, string calling_object_str);
	void InitialisePost(CTime* pMyTime, CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, /*int** &rAPLDevASSYS,*/ int** &rAPLDevPIPES, int** &rAPLDevPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, CTransmissive* TransmPtr, std::string param_dir, int assyid, bool transm_ptr_exists, string parent_assy_res_dir, CEndEnvironment* &rEndEnv, CTransmissive* &rTransm, CJunction* &rJunc, string calling_object_str);
  void InitialisePostConfig(CProperties* pPpt);
	void ReadInput(CProperties* pPpt, char *InputFile);
	void ListProperties(CProperties* pPpt);

	void RunBoundary(CProperties* pPpt, double time, double vel_grad, double mfr_under_est, CTime &rMyTime, double DELZe, double DELZi, double ca, double del_ca, double TIMEe, double TIMEi, int timestep, bool STEADY);

	void General(CProperties* pPpt, double time, double vel_grad, double mfr_under_est, int timestep, bool STEADY);
	void Loss(CProperties* pPpt, double time, double &rlambda_in_c2, double &rAA_c2, double &rlambda_out1, double &rlambda_out2, double vel_grad, double mfr_under_est, int timestep, bool STEADY);
  double LossEquations(CProperties* pPpt, double time, int timestep, double Ma);
	void Calibrate(CProperties* pPpt, CTime &rMyTime, double DELZe, double DELZi, double ca, double del_ca, double TIMEe, double TIMEi, int timestep, bool STEADY);

		double MeanlineCalc(CProperties* pPpt, double time, int timestep);
	
	void PrintToFile(CProperties* pPpt, int timestep, double time, double ca);

	//void InterpolateGraph(CProperties* pPpt, int device_graph, double S, double Re);
	//void QuadInterpolate(CProperties* pPpt, int device_graph, double S, double Re);
	void ReadResistanceCoefficients(char *InputFile);

	//void LoadSteadyOld(char* InputFile);
	void LoadSteadyCharacteristics(CProperties* pPpt, char* InputFile);
	void LoadLossCoefficients(CTime* pMyTime, CProperties* pPpt, char* InputFile, int admissionMapNo);
	void QuadraticallyInterpolateSteadyCharacteristics(int speedCurve, double PR, double& rMFP, double& rUCs, double& rEta);
	void LinearlyInterpolateSteadyCharacteristics(int speedCurve, double PR, double& rMFP, double& rUCs, double& rEta);
	void InterpolateLossCoefficients(CProperties* pPpt, double time, double curveValue, double xValue, double &rYValue, double &rZValue, bool &rEXTRAPOLATE, int lblCurve, int lblXtemp, int lblYtemp, int lblZtemp, bool SHOW_WARNINGS, int admissionMapNo);
	void RecordTurbineOperatingPoint(CProperties* pPpt);

	void PrintToScreen(CProperties* pPpt);
	
	char* Get_DEVICE();

	inline bool Get_RESET_INITIAL_CONDITIONS() {return RESET_INITIAL_CONDITIONS;}
	inline void Set_RESET_INITIAL_CONDITIONS(bool reset_initial_conditions) {RESET_INITIAL_CONDITIONS = reset_initial_conditions;}

  inline void Set_pUs(double pUs_temp) {pUs = pUs_temp;}
	
private:
	// ====================================================================================================
	// Parameter file for Exhaust APLD[0]
	// ====================================================================================================						
	
	int DEVICE;							// Type of device: gauze (1), throttle (2), EGR valve (3)
	double PARAM_S;						// Curve parameter value, e.g., solidity or speed parameter
	double PARAM_R;						// X-axis parameter value, e.g., upstream Re no. to downstream Mach no.
	double X;							// Dynamic throttle factor to allow for unsteady flow; set to 1 for gauze and egr valve
	bool NAPLDev;						// Non-adiabatic formulation (1=true) or Adiabatic formulation (0=false)
	
	// Configuration
	// ----------------------------------------------------------------------------------------------------
	int nEntries;						// Number of entry end environments
	bool transmEnt;						// Entries controlled by transmissive (1=true) or end environment boundaries (0=false)
	int* entry;							// Holds ID numbers of entry end environments
	int exitEndEnv;						// ID of exit end environment
	int juncInlet;						// ID of rotor inlet junction
	int innerInlet;						// ID of pipe immediately upstream of rotor inlet junction, inner limb
	int innerFollowing;					// ID of pipe immediately downstream of rotor inlet junction, inner limb
	int outerInlet;						// ID of pipe immediately upstream of rotor inlet junction, outer limb
	int outerFollowing;					// ID of pipe immediately downstream of rotor inlet junction, outer limb
	bool averagePR;						// Calculate turbine stage PR by averaging across all entries (1=true) or specify the entry to use (0=false)
	bool massAverageMFP;				// Calculate turbine stage MFP by mass averaged across entries (1=true) or normal average (0=false) (only if nEntries>1)
	
		// If not averaging to calculate PR, i.e. averagePR == 0 == false
		// ----------------------------------------------------------------------------------------------------
		int entryPR;						// ID of entry used to calulate PR (useful under partial admission conditions)

	// Loss coefficient
	// ----------------------------------------------------------------------------------------------------
	bool CONSTANT_K;					// Constant loss coefficient K (1=true) or variable (0=false)

		// If applying constant loss coefficient K, i.e., CONSTANT_K == 1 == true
		// ----------------------------------------------------------------------------------------------------
		double K_VALUE;						// Constant loss coefficient K value
		double etaConst;					// Associated constant eta value
		double D_VALUE;						// Associated constant D value

		// Else variable loss coefficient K, i.e., CONSTANT_K == 0 == false
		// ----------------------------------------------------------------------------------------------------	
		bool INTERP_M2;						// Interpolate for loss coefficient based on downstream Mach number (1=true) or upstream Re (0=false)
		char** LOSS_FILE;					// Names of the file(s) containing the full and partial admission loss data
		bool USE_LAMBDA;					// Use local admission ratio to interpolate loss coefficient (1=true) or use full admission value only (0=false)
		
			// If interpolating K between full and partial admission values, i.e., USE_LAMBDA == 1 == true
			// ----------------------------------------------------------------------------------------------------
			bool ABS_LAMBDA;					// If USE_LAMBDA, using absolute mass flow to calculate local admission ratio (1=true) or use the actual value (0=false). (default=0)
			bool LAMBDA_THRESHOLD;				// If not ABS_LAMBDA, apply admission ratio threshold so that LAMBDA = {0,1} (1=true) or not (0=false). (default=1)
			bool LINEAR_INTERP;					// Linearly interpolating local admission ratio (1=true) or polynomially (0=false). (default=0)
			double weightPartial;				// If LINEAR_INTERP, factor (0-1) to weight interpolation of K in favour (->1) of the partial admission value (default=0)
			double POLY_DEGREE;					// If not LINEAR_INTERP, the degree of polynomial to use. (default=6)

				
		bool PRINT_LOAD_LOSS;				// Print loading of loss coefficients to screen (1=true) or not (0=false)
		bool PRINT_LOSS_FILES;				// Print loss file to screen (1=true) or not (0=false)
		bool CALIBRATE;						// Calibrate the loss curve during this simulaton (1=true) or not (0=false)
		
			// If calibrating, i.e. CALIBRATE == 1 == true
			// ----------------------------------------------------------------------------------------------------
			char* STEADY_FILE;					// File containing the steady PR-MFP characteristic(s)
			double tolPR;						// Permissible difference in PR between consecutive data points (%)
			double tol_MFP;						// Tolerance when attempting to match the desired MFP (%)
			int waitCount;						// Maximum number of time steps to wait for steady conditions
			double tInterval;					// Time interval to wait before recording operating point (s)
			bool CONST_BASELINE;				// Use constant value for initial K, otherwise load from FULL_LOSS_FILE
			double K_baseline;					// If CONST_BASELINE, initial (constant) value of loss coefficient K
			char* RESTART_LOSS_FILE;			// The name of the file containing the restart loss data
			bool PRINT_STEADY_SCREEN;			// Print processed steady mass flow characteristic data to screen (1=true) or not (0=false)
			bool PRINT_STEADY_FILE;				// Print processed steady mass flow characteristic data to file (1=true) or not (0=false)
			double CALIBRATION_DELAY;			// Apply delay for calibration to start (s); set to 0 to start without delay
			double PR_THRESHOLD;				// Set the PR threshold, any points below this PR will not consider the non-adiabatic term during calibration

	// Measurements
	// ----------------------------------------------------------------------------------------------------
	bool USE_DEF_FREQ;					// Use default sampling rate (1=true) or 'freq' below (0=false)
	int freq;							// Record data every freq timesteps
	double freqPulse;					// Frequency to be used to separate cycle data (Hz)
	int degCycle;						// Number of degrees in cycle, e.g., 360 (pulse generator, two-stroke), 720 (four-stroke)
	//bool PRINT_MOVIE_FILE;			// Record data for movie (1=true) or not (0=false) LOCATED IN Boundary.h

	// Turbine geometry
	// ----------------------------------------------------------------------------------------------------
	double tipDiameter;					// Turbine rotor tip diameter (mm)

	// ====================================================================================================
	// End of file
	// ====================================================================================================
	
  // Loop variables
  double lambda_in_n[2];
	double AA_n[2];
	double lambda_in[2];
	double lambda_out[2];
	double A_star[2];
	double U_star[2];
	double AA_1, lambda_in_d1, AA2_over_AA1, AA_c2;
  double lambda_in_star_n2;
  double lambda_in_star_c2;
  double lambda_in_c2;
	double M2_max, temp_M2;
	double del_M2;
	double phi, a, b, alpha, k_ds, k_us;
  bool FIXING_K;

  int nMa;
  int lblMa;
  int lblMaError;
  int MaLevel;
  double** errorMap;
  double lowMa, highMa;
  int lblNeg, lblPos;
  bool FINISHED;

  // General APLD variables
	// ====================================================================================================
  double velHead;             // The upstream velocity head (bar)
  double delP;                // The static pressure loss across the device (bar)
  double pUs;                 // The upstream static pressure (bar)
  double a_ds;				  // The downstream sonic velocity (m/s)
  double T0_us;				  // The upstream stagnation temperature (K)

	// Operating point variables
	// ====================================================================================================

	double timeInCycle, timeInCyclePrev;

	// Separate values for each limb, and a number of locations per limb
	// ----------------------------------------------------------------------------------------------------
	int nLocs;							// Number of locations
	int LIMB_INLET;						// Limb location at the housing entry
	int LIMB_EXIT;						// Limb location closest to the rotor
	int ROTOR;							// Location just upstream of the rotor
	
	double **p, **T, **m, **u, **p0, **T0;	// Flow parameters measured at LIMB_INLET, LIMB_EXIT, and ROTOR
	double* m_stage;						// Stage mass flow rates summed over entries at LIMB_INLET and LIMB_EXIT
	double* m_abs_sum;
	double* p_sum;							// Static pressures summed over entries at LIMB_INLET and LIMB_EXIT
	double* lambda;							// Admission ratio, lambda, at [LIMB_INLET] and [LIMB_EXIT] though only [LIMB_EXIT] is really used, [LIMB_INLET] is for information
	double* tipVel;
	double* PR;
	double* MFP;
	double* UCs;
	double* Ws;
	double* tempEta;

	double mRotorInletJunction;
	double m_dot_inner, T_inner, T0_inner, p_inner, p0_inner, u_inner;
	double m_dot_outer, T_outer, T0_outer, p_outer, p0_outer, u_outer;
	double m_dot_inner_following, m_dot_outer_following;
	double p2_shaft;
	
	CNode* tempNode; // Used to look at rotor inlet junction nodes

	//double a0;							// First polynomial constant for admission ratio interpolating, if not linearly
	//double a1;							// Second polynomial constant for admission ratio interpolating, if not linearly

	bool SKIP;

	// Values for whole stage
	// ----------------------------------------------------------------------------------------------------
	double p_exit;						// MEasured exit static pressure - only one exit station
	double PR_sum;						// Sum PR
	double UCs_sum;						// Sum UCs
	double p0_sum;						// Sum p0
	double mT0_sum;						// Sum mass averaged T0
	double PR_stage_mean;				// Measured mean PR across all entries
	double MFP_stage;					// Measured stage MFP
	double UCs_stage_mean;				// Measured mean velocity ratio, U/Cis, across all entries
	double Ws_stage_inlet;				// Measured stage isentropic power (turbine inlet), Ws
	double UCs_stage_shaft;				// Measured velocity ratio (close to shaft), UCs_stage_shaft
	double Ws_stage_shaft;				// Measured stage isentropic power (close to shaft), Ws_shaft
	double Wact_stage_shaft;			// Stage actual power = eta*Ws_shaft, Wact
	double dh0_turb;					// Stage specific actual stagnation enthalpy change = Wact_stage_shaft/mRotorInletJunction, dh0_turb
	
	int nAdmissionMaps;					// Number of admission maps (e.g., full admission + two partial admission if twin entry)
	double M[2];						// Up and downstream Mach no.
	double* K;							// Interpolated loss coefficients for each admission map, K
	double* eta;						// Instantaneous efficiency, eta, associated to interpolated K for each admission map
	double* D;							// Instantaneous non-adiabatic term, D, associated to interpolated K for each admission map
	bool* OUTSIDE_RANGE;				// Is operation outside of calibrated range on each K loss coefficient map
	double K_OVERALL;					// Final value of K applied at the N/APLD boundary, following interpolation between different admission map values
	double eta_OVERALL;					// Final value of eta found at the N/APLD boundary, following interpolation between different admission map values
	double D_OVERALL;					// Final value of D found at the N/APLD boundary, following interpolation between different admission map values
	bool OUTSIDE_RANGE_OVERALL;			// Is operation outside of calibrated range on any of the K loss coefficient maps (full and partial, if in use)?

	double ****Coeff;					// Matrix of resistance coefficients for this adiabatic pressure loss device
	int graph, num_graphs, curve, point;
	int *num_curves;
	int **num_points;

	int UPSTREAM, DOWNSTREAM;
	bool LEFT_TO_RIGHT;
	int counter;						// Main loop counter
	int counter_between_steady;			// Number of timesteps since last steady conditions

	// Labels for vector MFP_vs_PR
	// ----------------------------------------------------------------------------------------------------
	int lblMFPvsPR_p0_applied, lblMFPvsPR_P0, lblMFPvsPR_PR, lblMFPvsPR_MFP, lblMFPvsPR_M2, lblMFPvsPR_UCs, lblMFPvsPR_eta;

	// Steady mass flow characteristics from file
	// ----------------------------------------------------------------------------------------------------
	double*** steadyPRMFP;				// Array containing read-in steady data [lblSpeed=speed parameter][lblPR=PR][lblMFP=MFP][lblUCs=U/Cis][lblEta=eta]
	int lblSpeed, lblPR, lblMFP, lblUCs, lblEta;	// Labels for steadyPRMFP
	int numSteadySpeeds;				// Number of speeds in steady PR-MFP data
	int* numSteadyPoints;				// Number of (non-zero MFP) data points on each steady mass flow characteristic
	int steadyCurve;					// Current steady mass flow characteristic [0 -> numSteadySpeeds - 1]

	// Loss coefficients from restart file
	// ----------------------------------------------------------------------------------------------------
	double**** K_loss_file;				// Array containing read-in loss coefficient data [lblSpeed=speed parameter][lblX=M2][lblY=K][lblZ=U/Cis]
	int lblSpeedFile, lblM2File, lblKFile, lblKPrevFile, lblKPrevPrevFile, lblErrorFile, lblErrorPrevFile, lblErrorPrevPrevFile, lblFactorFile, lblKUCsFile, lblKEtaFile, lblKDFile; // Labels for K_loss_file
	int* numLossSpeeds;					// Number of speeds in loss coefficient data on each admission map
	int** numLossPoints;				// Number of (non-zero K) values on each loss curve on each admission map
	int lossCurve;						// Current loss curve [0 -> numLossSpeeds[admissionMapNo] - 1]

	// Loss coefficients at runtime
	// ----------------------------------------------------------------------------------------------------
	double*** K_loss_coeffs;			// Vector containing adjustable loss coefficients [speed][point][values]
	int lblM2, lblK, lblKPrev, lblKPrevPrev, lblError, lblErrorPrev, lblErrorPrevPrev, lblFactor, lblKUCs, lblKEta, lblKD; // Labels for K_loss_coeffs
		
	vector<double*> MFP_vs_PR;	// Vector containing recorded steady mass flow characteristic
	CTransmissive* pTransm;	// Pointer to controlling transmissive boundary
	bool TRANSM_PTR_EXISTS; // Whether or not there is a transmissive boundary to point at
	bool REFINING;			// Refine K to give desired MFP
	bool LAST_RUN_LESS;		// Switching flag
	double delK;			// Loss coefficient step size
	//FILE* pNEW_LOSS_FILE;	// Pointer to new loss output file
	FILE* pCalibrationProgressK;	// Pointer to file recording calibration progress in terms of loss coefficient curves
	FILE* pCalibrationFinalK;		// Pointer to file recording calibrated loss coefficient curves
  FILE* pCalibrationFinalKDynasty;		// Pointer to file recording calibrated loss coefficient curves, Dynasty format
	FILE* pCalibrationProgressMFPPR;// Pointer to file recording calibration progress in terms of steady mass flow characteristics
	FILE* pCalibrationFinalMFPPR;	// Pointer to file recording calibrated steady mass flow characteristics
	FILE* pProcessedMap;	// Pointer to file recording processed/trimmed experimental steady map
	FILE* pFinalKRestart;	// Pointer to file recording calibrated loss coefficient curves including error and factor data for use as a restart file
	bool M2_TOO_LOW;		// Flag to remember state of M2
	double tol_M2;			// Tolerance when matching M2
	int K_COUNTER;			// Identifies which element of K_loss_coeffs to calibrate
	double p0_applied;		// Total pressure applied during calibration
	double p0_applied_safe; 
	double del_p0_applied;	// Interval length (bar)
	int row;
	double PR_point, MFP_point, UCs_point;	// Operating point parameter values
	
	// Transmissive, end environment, and junction control
	// ----------------------------------------------------------------------------------------------------
	CEndEnvironment** pEntryEndEnv;	// Pointers to end environments controlling conditions at turbine housing entries
	CTransmissive** pEntryTransm;	// Pointers to transmissive boundaries controlling conditions at turbine housing entries
	CNode** pRotorInletNode;		// Pointers to pipe nodes closest to rotor boundary in each limb
	CEndEnvironment* pExitEndEnv;	// Pointer to exit end environment
	CJunction* pRotorInletJunction;	// Pointer to rotor inlet junction
	CPipe* pInnerInlet;				// Pointer to pipe immediately upstream of rotor inlet junction, inner limb
	CPipe* pInnerFollowing;			// Pointer to pipe immediately downstream of rotor inlet junction, inner limb
	CPipe* pOuterInlet;				// Pointer to pipe immediately upstream of rotor inlet junction, outer limb
	CPipe* pOuterFollowing;			// Pointer to pipe immediately downstream of rotor inlet junction, outer limb
		
	// Labels for pEntryEndEnv
	// ----------------------------------------------------------------------------------------------------
	int SINGLE_ENTRY, FIRST_ENTRY, SECOND_ENTRY;

	int pass_counter;				// Counts number of calibration passes carried out
	bool PASS_INCREASING;			// Flag denoting whether pass increases or decreases the applied stagnation pressure
	bool RESET_INITIAL_CONDITIONS;	// Flag to instruct main function to reset initial gas flow conditions everywhere
	//bool FAKE_STEADY;				// Flag to indicate steady conditions should be assumed

	double tCurrent;				// Current simulation time (s)
	double tRecord;					// Time at which previous operating point was recorded (s)
	
  // Calibration variables
	double currentError;
	double error_max;
	double error_max_percent;
	int max_error_element;

  // Boundary method variables
  double error, minError;

  // Non-adiabatic variables
  double*** D_temp;
  double d_temp, Wact_turb;

  double starttime, endtime, elapsedtime;
  double mstarttime, mendtime, melapsedtime;
  bool COMPUTETIME;

  // MeanlineCalc stuffs
	double Cteta_ds;
	double p_us;
	double T_us;
	double W_us;
	double p0r_us;
	double T0r_us;
	double Wcr_us;
	double p0r_ds;
	double T0r_ds;
	double Wcr_ds;
	double p_ds;
	double T_ds;
	double W_ds;
	double K_meanline;				// Pressure loss coefficient computed by turbomachinary losses consideration
	//double d;							// Additional d_term for Eq. 8.168 
	//double dh0_turb;
	//double a_ds;
	bool DONE;
};

#endif // !defined(AFX_APLDEV_H__87820053_A2AE_415A_B790_5D4F8BD891B5__INCLUDED_)
