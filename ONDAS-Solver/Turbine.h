// Turbine.h: interface for the CTurbine class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TURBINE_H__DB25CD3F_14D5_4C53_A761_6CC3A10E13A7__INCLUDED_)
#define AFX_TURBINE_H__DB25CD3F_14D5_4C53_A761_6CC3A10E13A7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Boundary.h"
#include "DataPoint.h"
//#include "Node.h"
#include "PathLine.h"
#include "Pipe.h"
#include "Properties.h"

class CTurbine : public CBoundary  
{
public:
	CTurbine();
	virtual ~CTurbine();

	void Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, /*int** &rTURBASSYS,*/ int** &rTURBPIPES, int** &rTURBPIPES_ENDS, int* &rTURBPIPES_NUMS, double* &rENDCORR, 
		int id, bool ex, std::string param_dir, int assyid, string parent_assy_res_dir, string calling_object_str);

	void InitialiseMap(CProperties* pPpt);
	void ReadInput(CProperties* pPpt, char *InputFile);
	void ListProperties(CProperties* pPpt);

	// Calculation functions
	void RunBoundary(CProperties* pPpt, int timestep, double time, double DELZ);
	void VariablePressure(CProperties* pPpt/*, int& rINLET, int& rOUTLET,*//* bool TURBINE, double* &rC,*//* double time*/); 
	void ConstantPressure(CProperties* pPpt);
	void EquivArea(CProperties* pPpt, int timestep, double time); // For axial turbines and constant outlet pressure only
	void Interpolate(CProperties* pPpt, double speed, double query_temp);
	void InstantaneousWorkAndMassFlow(CProperties* pPpt, double del_t);
	void CycleWorkAndMassFlow(CProperties* pPpt, double Nc);
	void InstantaneousTransmissionLoss(CProperties* pPpt, double del_t);
	void Matching(CProperties* pPpt, double del_t);
	
	// Map processing
	void LoadMap(CProperties* pPpt, char *InputFile, int* &rnum_curves, int** &rnum_points, double**** &rRaw, bool SINGLE_CURVE);
	void LoadMapSAE(CProperties* pPpt, char *InputFile);

	void ProcessMapVar(CProperties* pPpt, int num_extra_reverse, int num_extra_choked, char* proc_map_file_str, bool TURBINE, double T);

	void ProcessMapConst(CProperties* pPpt, double T);
	void ProcessDataPtConst(CProperties* pPpt, double pr, double mfp, double &rlambda_in_star1, double &rlambda_out_star1, double T);	

	// Read/write functions
	void PrintToScreen(CProperties* pPpt);
	void PrintToFile(CProperties* pPpt, int timestep, double time, double ca);
	void PrintToMovieFile(CProperties* pPpt, int timestep, double time, double ca);

private:
	// ====================================================================================================
	// Parameter file for Exhaust Turbine [0]
	// ====================================================================================================
			
/*
	// Type
	// ----------------------------------------------------------------------------------------------------
	bool RADIAL;			// Type of turbine; RADIAL (1=true) or AXIAL (0=false)

		// If this is an axial flow turbine, i.e., RADIAL == 0 == false
		// ----------------------------------------------------------------------------------------------------
		bool EQUIVALENT_AREA;	// If axial turbine, either EQUIVALENT_AREA (1=true) or SIMPLE_UNIQUE_CURVE (0=false) - N.B. equivalent area method uses nozzle so has a constant outlet pressure

		// Else this is a radial turbine, i.e., RADIAL == 1 == true
		// ----------------------------------------------------------------------------------------------------
		bool SINGLE_SPEED;	// If radial, hold turbine at constant speed (only one speed line provided) (1=true) or not (0=false)

	bool FIXED_SPEED;		// Fix rotor rpm (1=true) or not (0=false)

		// If speed is fixed, i.e., FIXED_SPEED == 1 == true
		// ----------------------------------------------------------------------------------------------------
		double fixedSpeed;		// Fixed rotor speed (r/min)
*/

	// Type
	// ----------------------------------------------------------------------------------------------------
	bool EQUIVALENT_AREA;	// Treat turbine as nozzle of equivalent area (implies constant outlet pressure) (1=true) or not (0=false) 
	bool FIXED_SPEED;		// Fix rotor rpm (1=true) or not (0=false)

		// If speed is fixed, i.e., FIXED_SPEED == 1 == true
		// ----------------------------------------------------------------------------------------------------	
		double fixedSpeed;		// Fixed rotor speed (r/min)

	double Pb;				// The value of back pressure (bar) for a constant outlet pressure turbine
	double T0;				// The value of constant temperature (K) for a constant outlet pressure turbine
	double I;				// Turbine inertia (kg.m^2)

	// Map
	// ----------------------------------------------------------------------------------------------------
	char* MAP_FILE;			// The filename of the map for this turbine
	bool SINGLE_SPEED;		// Only one speed line is provided on the map (1=true) or not (0=false) 
	double tolSpeed;		// Percentage variation in speed parameter allowed for points on the same line
	bool ECHO_MAP;			// Print SAE map to screen (1=true) or not (0=false)
	bool PRINT_PROC_MAP;	// Print processed SAE map to screen (1=true) or not (0=false)

	// Measurements
	// ----------------------------------------------------------------------------------------------------
	bool USE_DEF_FREQ;		// Use default sampling rate (1=true) or 'freq' below (0=false)
	int freq;				// Print data to file every freq timesteps

	// Deprecated
	// ----------------------------------------------------------------------------------------------------
	//bool INLET_STAG;		// Map shows stagnation (p01 & T01) (1=true) or static values (p1 & T1) (0=false) at turbine inlet
	//bool OUTLET_STAG;		// Map shows stagnation (p02) (1=true) or static values (p2) (0=false) at turbine outlet
	//int MFP_P;				// MFP (m_dot.sqrt(T(0)1/p(0)1) pressure component has units of Pa (0), kPa (1), or bar (2)

	// ====================================================================================================
	// End of file
	// ====================================================================================================

public:

	// Labels
	// ----------------------------------------------------------------------------------------------------
	int SP, MFP, PR, ETA;
	int INLET, OUTLET;
	
	bool VARIABLE;					// Variable outlet pressure (1) or constant outlet pressure (0) turbine

	// SAE map data
	// ----------------------------------------------------------------------------------------------------
	std::string strNomenclature;	// Nomenclature from first line of SAE file
	double** rawSAE;				// SAE map data prior to separation into speed lines
	int speeds;						// Number of speed lines
	double*** SAE;					// SAE map data split into speed lines
	int* numPointsSAE;				// Number of points on each speed line

	// Map data
	// ----------------------------------------------------------------------------------------------------
	double ****Raw;					// Matrix of raw map data
	int *num_curves;				// Number of curves on each graph
	int **num_points;				// Number of points on each curve on each graph
	CDataPoint** Data;				// Matrix of processed map data point objects
	CDataPoint interp_data_pt;		// The interpolated data point

	char *pr_str, *mfp_str, *sp_str, *eta_str;	// Strings describing the form of the map parameters
	char *mfp_units_str; //, *sp_units_str;			// Strings describing the units of the parameters

	bool ON_MAP;					// Instantaneous speed lies within range of available speed curves (true) or has to be extrapolated (false)
	bool ON_CURVE;					// Instantaneous point lies within range of available curve data (both curves) (true) or has to be extrapolated (false)
	bool REVERSE_FLOW;				// Reverse flow through turbine
	bool CHOKED_FLOW;				// Choked flow through turbine
	
	int curve_above, curve_below;	// Interpolated point lies between these speed curves 

	// Instantaneous variables
	// ----------------------------------------------------------------------------------------------------
	double T01;						// Instantaneous stagnation temperature T01 ahead of the turbine
	double T1;						// Instantaneous static temperature T1 ahead of the turbine

	// Variable outlet pressure turbine parameters
	// ----------------------------------------------------------------------------------------------------
	int num_extra_reverse;			// Number of extra points to include for reverse flow data
	int num_extra_choked;			// Number of extra points to include for choked flow data

		// Zero and reverse flow variables
		// ----------------------------------------------------------------------------------------------------
		double* x_min;
	
		// Choked data variables
		// ----------------------------------------------------------------------------------------------------
		double* Pt;
		double* theta_t;
		double* G1_max;
		double* K1;
		double* K2;
		double* K3;
		double* K4;
		double* x_max;
		
		// Calculation variables
		// ----------------------------------------------------------------------------------------------------
		double* C;

	// Shaft dynamics (common)
	// ----------------------------------------------------------------------------------------------------
	double NT;						// Instantaneous turbine speed (s^-1)
	double dNTdt;					// Instantaneous shaft acceleration (s^-2)

	// Turbine work and mass flow parameters (common)
	// ----------------------------------------------------------------------------------------------------
	double eA_over_m_dot;			// Instantaneous turbine work_per_unit_mass; this is the parameter (eta_TS)*((C_is^2)/2) (J.kg^-1)
	double m_dot;					// Instantaneous turbine mass flow rate (kg^s-1)
	double eT;						// Theoretical turbine instantaneous power (W)
	double eA;						// Actual instantaneous turbine work per unit time, i.e. the actual instantaneous power (W), eA = eta_TS*eT
	double eta_TS;					// Instantaneous turbine efficiency
	
	// Instantaneous matching parameters
	// ----------------------------------------------------------------------------------------------------
	double W_TI;					// Instantaneous turbine power (copy) (W)
	double L_TI;					// Instantaneous turbine torque (N.m)
	
	// Cycle work and mass flow running totals
	// ----------------------------------------------------------------------------------------------------
	double work_per_cycle_Th;
	double work_per_cycle_A;
	double mass_flow_per_cycle;

	// Cycle/total values (calculated once per cycle)
	// ----------------------------------------------------------------------------------------------------
	double W_Th;					// Theoretical turbine cycle/total power (W)
	double W_TA;					// Actual turbine cycle/total power (W)
	double eta_T;					// Average turbine efficiciency  = W_TA/W_Th = (actual turbine total power/total available power)
	double m_dot_T;					// Turbine cycle/total mass flow rate (kg.s^-1)

	// Instantaneous matching parameters
	// ----------------------------------------------------------------------------------------------------
	double W_B;						// Instantaneous bearing power loss (W)
	double L_B;						// Instantaneous bearing torque (N.m)
	double eta_M;					// Instantaneous mechanical efficiency

	double query;					// Instantaneous look-up query value
};

#endif // !defined(AFX_TURBINE_H__DB25CD3F_14D5_4C53_A761_6CC3A10E13A7__INCLUDED_)
