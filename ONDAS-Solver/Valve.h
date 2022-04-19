// Valve.h: interface for the CValve class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VALVE_H__51AFE53A_550C_42F9_87DD_3E9F531D1245__INCLUDED_)
#define AFX_VALVE_H__51AFE53A_550C_42F9_87DD_3E9F531D1245__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Boundary.h"
#include "Pipe.h"
#include "Properties.h"

class CValve : public CBoundary  
{
public:

	CValve();
	virtual ~CValve();

	CValve(const CValve& inValve);				// Copy constructor
	CValve& operator=(const CValve& inValve);	// Overloaded operator=

	void Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rVPIPES, int** &rVPIPES_ENDS, double* &rENDCORR,
					int id, bool ex, int npipes, CEngine* EngPtr, std::string param_dir/*, char* vt_dir*/, int assyid, int cyl_id, string parent_assy_res_dir, 
					int nPipesInAssy);

	//void Initialise(CProperties* pPpt, CPipe** pPipes, CPipe*& rPipe, int**& rTRANSPIPES, int**& rTRANSPIPES_ENDS, double*& rENDCORR,
	//	int id, bool ex, int npipes, std::string param_dir, int assyid, string parent_assy_res_dir, int ntransm, int nPipesInAssy, string calling_object_str);

	void InitialiseVolumeValve(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rVPIPES, int** &rVPIPES_ENDS, double* &rENDCORR,
						int id, bool ex, int npipes, std::string param_dir, int assyid, int vol_id, string parent_assy_res_dir);

	void ConfigureExtra(CProperties* pPpt, CPipe** pExhaustSystem, CPipe** pIntakeSystem, int nExPipes, int nInPipes);
	void ListProperties(CProperties* pPpt, char* prefix);

	void PrintToFile(CProperties* pPpt, int timestep, double time, double ca_elapsed, double ca);

	void PrintToScreen(CProperties* pPpt);
	void PrintToScreenVolume(CProperties* pPpt);

	// VVT
	void LoadValveTiming(char* InputFile);
	void LoadFlowArrays(char* InputFile);

	// AVT
	void LoadActiveValveTarget(char* InputFile);
	void RecordAVT(CProperties* pPpt, int timestep, double time_or_CADelapsed, double time_step_length); //void RecordAVT(CProperties* pPpt, int timestep, double time, CPipe* Pipe);
	
	void EffectiveArea(CProperties* pPpt, double ca);
	void EffectiveArea(CProperties* pPpt);

	// Homentropic methods
	void Poppet_H(CProperties* pPpt, int timestep, double time, double &rDMDT, double psi, double AREF, double PREF, double PC, double TC, double Fp);
	
	// Non-homentropic methods
	double Poppet_NH(CProperties* pPpt, double psi, double P0, double T0, double AREF, int cyl_id, bool UPDATE_PIPE, int timestep, double time);

	void Equations(CProperties* pPpt, double &lambda_in_c, double &lambda_out_c, double &AA_c, double &AA_calc, 
				  double lambda_in_n, double AA_n, double Ac, double &U, double &C, double &PpbyPc, 
				  double psi, double rc, bool &rCHOKED, bool &ERROR, double T);	
	
	void NHP(CProperties* pPpt, double &rlambda_in, double &rlambda_out, double lambda_out_old, 
		 double &rAA, CPathLine &rPathLine, double &rDMDT, double psi, double K, double AREF, 
		 double PREF, double PC, double TC, double AC, double Fp, double ca, 
		 int &rpipe_flow, int end, double XPIPE, int timestep, double time);

	inline double	Get_CLIN(){return (*(pCLIN[ONE_SIDE]))[R+1];}
	inline double	Get_CLOUT(){return (*(pCLOUT[ONE_SIDE]))[R+1];}
	inline double	Get_eff_area(){return eff_area;}
	inline void		Set_eff_area(double area_to_be_set){this->eff_area = area_to_be_set;}
	inline double	Get_FP(){return FP;}
	inline void		Set_valve_Cf_or_Cd(double Cf_or_Cd_to_be_set) { this->valve_Cf_or_Cd = Cf_or_Cd_to_be_set; }
	inline double	Get_final_lift(){return final_lift;}
	inline void		Set_final_lift(double lift_to_be_set_mm) { this->final_lift = lift_to_be_set_mm; } // [mm]
	inline double	Get_VC(){return VC;}

	inline bool Get_open(void){return open;}
	inline void Set_open(bool set_open){open = set_open;}

	// Volume valve functions
	// ----------------------
	inline void Set_USE_PIPE_DIA(bool USE_PIPE_DIA_temp){USE_PIPE_DIA = USE_PIPE_DIA_temp;}
	inline void Set_ref_dia(double dia_temp){USE_PIPE_DIA ? ref_dia = this->pBN[ONE_SIDE]->d*1e3 : ref_dia = dia_temp;}
	inline void Set_psi_val(double psi_val_temp){psi_val_temp > 1 ? psi_val = 1 : psi_val = psi_val_temp;}
	inline double Get_psi_val(){return psi_val;}

	inline double Get_dynamic_head(CProperties* pPpt){return 0.5*pBN[ONE_SIDE]->rho*pow(pBN[ONE_SIDE]->U*(EX ? pPpt->AREFe : pPpt->AREFi),2);}

	inline double Get_pBN_p_dash(){return pBN[ONE_SIDE]->p_dash;}
	inline double Get_pBN_T(){return pBN[ONE_SIDE]->T;}
	inline double Get_pBN_U(){return pBN[ONE_SIDE]->U;}
	inline double Get_pBN_Re(){return pBN[ONE_SIDE]->Re;}

private:
	bool open;				// Valve open or closed
	double FP;				// Area of joining pipe (m^2)
	double VO;				// Valve opening angle (degrees)
	double VC;				// Valve closing angle (degrees)

	// Algorithm control
	// -----------------
	int num_loops_current;
	int num_loops_max;

	// Valve timing related variables
	// ------------------------------
	CEngine* pEng;			// Pointer to engine on which valve is located
	int CYL_ID;				// ID of cylinder on which the valve is located

	// Valve parameters
	// ----------------------------------------------------------------------------------------------------
	double ref_dia;			// Valve reference diameter (m)
	//bool L_OR_A;			// Is lift data lift values or area values?
	int L_OR_A;			// Is lift data lift values or area values?

	// If valve lift or area profile is read from file
	// ----------------------------------------------------------------------------------------------------
	double cam_timing_ang;	// Cam-timing angle; angle between TDC firing and angle of first timing point
	int cam_timing_loc;		// Cam timing angle locates start of valve opening (0), maximum opening (1), or end of valve closing (2)

	// If valve data are lift values, i.e., L_OR_A == 1 == true
	// ----------------------------------------------------------------------------------------------------
	double lash;			// Valve lash (m)
	 
//	// Else exhaust valve data are area values, i.e., L_OR_A == 0 == false
//	// ----------------------------------------------------------------------------------------------------

	// Exhaust valve flow or discharge coefficients (can be used with either lift or area values)
	// ----------------------------------------------------------------------------------------------------
	bool const_ref_area;	// Use constant reference area (i.e., values are flow coefficient) (1) or curtain reference area (i.e., values are flow coefficient) (0) to interpret file values
	double onum;			// Multiplier to the flow or discharge coefficient arrays ("def"=1)

	// Variable valve timing and lift
	// ----------------------------------------------------------------------------------------------------
	//double lift_multiplier;	// Lift multiplier; scaling factor for lift values
	//double area_multiplier;	// Area multiplier; scaling factor for area values
	double opening_multiplier; // Valve opening (phi value or lift or area) multiplier
	double ang_multiplier;	// Valve opening duration multiplier; scaling factor for theta values

	// ----------------------------------------------------------------------------------------------------
	// ----------------------------------------------------------------------------------------------------

	int datapoints;			// Number of rows in valve timing array
	double** vt_raw;		// Array containing valve timing raw data prior to adjustment
	double** vt;			// Array containing adjusted valve timing; 2 columns (angle then lift or area)
	int datapoints_flow;	// Number of rows in valve flow arrays
	double** Cd;			// Array containing flow coefficients; 3 columns (reference (lift/diameter), forward Cd, reverse Cd)

	int THETA, LIFT, AREA;	// Column labels for lift or area valve timing arrays
	int REF, FORWARD_CD, REVERSE_CD;	// Column labels for flow arrays

	double raw_lift;		// Lift (mm) before application of multipliers and lash etc.
	double final_lift;		// Actual lift (mm) seen
	double valve_Cf_or_Cd;		// Coefficient of discharge
	double eff_area;		// Effective valve area (m^2)

	double dmdt;			// Mass flow rate (kgs^-1)

	bool SONIC;				// Sonic flow across valve

	// Volume valve variables
	// ----------------------
	int VOL_ID;				// ID of volume on which the valve is located
	double psi_val;			// Area ratio
	bool USE_PIPE_DIA;		// Set valve reference diameter same as adjoining pipe (1=true) (0=false)?

public:

	// Target, applied and recorded values - arrays contain data for one pulse period
	double** target_applied_recorded;

	// Labels
	int TIME;
	int TIME_IN_PERIOD;

	int TARGET_P;
	int APPLIED_P;
	int RECORDED_P;
	int NEXT_P;
	int THIS_P;

	int TARGET_T;
	int TARGET_T0;
	int APPLIED_T;
	int RECORDED_T;
	int NEXT_T;
	int THIS_T;

	int TARGET_V;
	int APPLIED_V;
	int RECORDED_V;
	int NEXT_V;
	int THIS_V;

	int TARGET_MFR;
	int APPLIED_MFR;
	int RECORDED_MFR;
	int NEXT_MFR;
	int THIS_MFR;

	int VEL_GRAD;
	int P_VEL;
	int P_VEL_GRAD;
	int UNDER_EST;  // MFR underestimate
	int UNDER_EST_PREV;  // MFR underestimate from previous cycle
	int UNDER_EST_FACT;  // factor to multiply the underestimate by to get loss value
	int K_LOSS;  // loss value
	int WAVE_TIME;
	int WAVE_PARAMETER;
	int WAVE_MFR;

	// AVT
	// ----------------------------------------------------------------------------------------------------
	bool ACTIVE;			// Valve is used for AVT (1=true) or not (0=false)
	bool TARGET_PIPE_EX;	// Is target location in the exhaust or intake system
	CPipe** pTargetPipeSys; // Pointer to the exhaust or intake system (wherever the target location resides)
	int target_pipe_num;	// Number of the pipe in the exhaust or intake system where the target location resides
	double target_pipe_len;	// Length of the pipe in the exhaust or intake system where the target location resides
	double target_pipe_loc; // Location (as a fraction of pipe length) in the pipe where the target location resides
	CNode* pAVTNodeLeft;	// Pointer to the node on the left hand side of the target location
	CNode* pAVTNodeRight;	// Pointer to the node on the right hand side of the target location
	CNode AVTMeasureNode;	// Interpolated pseudo-node based on the left and right hand side nodes
	bool SINGLE_NODED_PIPE;	// In case pipe containing target location only has a single node
	int numAVTDataPoints;	// Number of rows in active valve array
	double** avt_raw;		// Array containing active valve raw data [CAD, target value] (prior to any adjustment)
	int CAD, TARGET, RECORDED, ERROR, ERROR_GRAD, SIGN_CHANGE, CORR_FACTOR, PSI_UNCORR, PSI_CORR, PSI_CORR_PREV;		// Column labels for active valve arrays in avt_raw

	double target_cyc_av;	// Cycle average of target waveform
	int target_max_row;		// Location (row number) of maximum value in target waveform
	double target_max_ang;	// Location (angle) of maximum value in target waveform
	double target_max_val;	// Maximum value in target waveform

	double recorded_cyc_av;	// Cycle average of recorded waveform at target location
	int recorded_max_row;		// Location (row number) of maximum value in recorded waveform
	double recorded_max_ang;		// Location (angle) of maximum value in recorded waveform
	double recorded_max_val;		// Maximum value in recorded waveform at target location

	double cyc_av_error;	// Error between the target cycle average and that recorded at the target location (= value_cyc_av - target_cyc_av)
	double cyc_av_error_prev;	// Previous value of cyc_av_error
	double phase_lag_between_peaks; // Phase angle between target and recorded peaks (positive values indicate recorded peak location is lagging target peak location)

	double plenum_press;	// Current plenum pressure [bar]
	double plenum_press_prev; //Previous plenum pressure [bar]
	double cyc_av_error_grad; // Rate of change of error
	double cyc_av_tol;		// Tolerance (%-age) within which a match is achieved
	bool CYC_AV_MATCHED;	// Flag to recorded whether a cycle average match has been achieved
	bool CYC_AV_ERR_SIGN_CHANGE; // Flag to check if error in cycle average value has changed sign
	double cyc_av_factor;		// Convergence rate factor (plenum pressure)
	//double AVTFactorSF;		// Convergence rate factor (scale factor)
	double cyc_av_factor_orig;	// Stores the original value of cyc_av_factor
	//double AVTFactorSFOrig; // Stores the original value of AVTFactorSF

	bool PHASE_MATCHING;	// Currently attempting to phase shift to improve match
	bool PHASE_MATCHED;		// Flag to recorded whether a phase match of the peaks has been achieved
	double phase_shift;		// Amount by which to phase shift calculated waveform

	double AVTargetValue;	// Target value for the current point in the cycle (interpolated from read-in array avt_raw) 
	double AVTCurrentValue;	// Value at the target location for the current point in the cycle
	double AVTCumulativeVal;// Cumulative value of AVTCurrentValue over the cycle 
	int AVTInCycleStepCount;// Yimestep counter within a cycle
	double AVTError;		// Instantaneous error between the target value (read from file) and the current value at the target location (= AVTCurrentValue - AVTargetValue)
	double AVTErrorPrev;	// Previous timestep value of AVTError
	double AVTErrorGrad;	// Rate of change of error
	bool FIRST_CYCLE;		// Flag while simulation is within its first cycle

	double time_in_period;	// The time (seconds) elapsed since the start of the current cycle
	//double CAD_in_period;	// Crank angle degree elapsed since the start of the current cycle
	double fraction_cycle;	// Fraction through current cycle (0-1) 
	



	//FILE* FILE_LOC;					// File containing data at specified measuring locations in attached pipe
	//FILE* FILE_p_des_app_seen;
	//bool PRINT_HEADERS;				// Have column headers been printed
	//double* AVTResults;
};

#endif // !defined(AFX_VALVE_H__51AFE53A_550C_42F9_87DD_3E9F531D1245__INCLUDED_)