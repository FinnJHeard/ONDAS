// Pipe.h: interface for the CPipe class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PIPE_H__C27D3667_10DC_401D_A116_3140F7A8D453__INCLUDED_)
#define AFX_PIPE_H__C27D3667_10DC_401D_A116_3140F7A8D453__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "Properties.h"
#include "FiniteVolume.h"
#include "Node.h"
#include "PathLine.h"
#include "Engine.h"

class CPipe  
{
public:
	CPipe();
	virtual ~CPipe();

	// Class setup
	// ----------------------------------------------------------------------------------------------------
	void ReadInput(CProperties* pPpt, char *InputFile);
	void Initialise(CProperties* pPpt, int assyid, int id, bool ex, CEngine* EngPtr, std::string param_dir, string parent_assy_res_dir, string calling_object_str);
	void InitialConditions(CProperties* pPpt);

	// Main control
	// ----------------------------------------------------------------------------------------------------
	void RunPropagation(CProperties* pPpt, double DELZ, int timestep, bool &rRESTORE);
	
	// Mesh Method of Characteristics (MMOC)
	// ----------------------------------------------------------------------------------------------------
	double TimeStepMOC(CProperties* pPpt);
	void HomentropicMOC(CProperties* pPpt, double DELZ);
	void NonHomentropicMOC(CProperties* pPpt, double DELZ, int timestep, bool &rRESTORE);
		void PathLines(CProperties* pPpt, double DELZ);
		void PathLinesAtDuctEnds(CProperties* pPpt, int timestep);
		void InterpolateEntropyAtMeshPoints(CProperties* pPpt);
		void MOC_Non_Homentropic(CProperties* pPpt, double DELZ, bool &rRESTORE);
		void MOC_Non_Homentropic_Interior(CProperties* pPpt, double DELZ, bool &rRESTORE);
		void MOC_Non_Homentropic_Boundary(CProperties* pPpt, double DELZ, bool &rRESTORE);
			void MOC_Non_Homentropic_lambdaI(CProperties* pPpt, double DELZ, int S, bool &rRESTORE);
			void MOC_Non_Homentropic_lambdaII(CProperties* pPpt, double DELZ, int S, bool &rRESTORE);
			double MOC_Non_Homentropic_Calcs(CProperties* pPpt, int S, double lambda_A, double lambda_B, 
											 double beta_A, double beta_B, double delXdX, 
											 double AAR, double AAQ, double AAR_dash, double DELX,
											 double XR, double XW_sign, double dDdX_sign, 
											 double d, double Tg, double DELZ, double f, double Re);
		void RemovePathLine(CProperties* pPpt, int SIDE, int timestep);

	// W_alpha_beta schemes
	// ----------------------------------------------------------------------------------------------------
	void W_alpha_beta(CProperties* pPpt, double DELZ);
	void TVD(CProperties* pPpt, double del_t);
	double G_flux_limiter(CProperties* pPpt, double lambda_max, double del_t, double del_x, double r);
	void W_alpha_beta_derive_lambdas(CProperties* pPpt);
	void W_alpha_beta_prep_bcs(CProperties* pPpt);

	// Glimm's method for "exact" solution to Riemann problem
	// ----------------------------------------------------------------------------------------------------
	double phi(CProperties* pPpt, double x, double epsilon);

	// Filling and emptying
	// ----------------------------------------------------------------------------------------------------
	void FillingAndEmptying(CProperties* pPpt, double DELZ);
	double* FillingAndEmptyingEquations(CProperties* pPpt, double T0, double p0, double ps, bool& rREVERSED);
	void FillingAndEmptyingAdjust(CProperties* pPpt, double DELZ, int& rCOUNTER, double TIME);
	void FillingAndEmptyingAdjust2(CProperties* pPpt, double DELZ, int& rCOUNTER, double TIME);
	void FillingAndEmptyingAdjust3(CProperties* pPpt, double DELZ, int& rCOUNTER, double TIME);
	void FillingAndEmptyingSimple(CProperties* pPpt);//, double DELZ, int& rCOUNTER, double TIME);

	// Other functions
	// ----------------------------------------------------------------------------------------------------
	double FrictionFactor(CProperties* pPpt, double d, double Re, char* callingfunc, double cftrans);

	// Colebrook-White (1939) implicit relationship for friction factor
	double Colebrook(double d, double Re);

	// Swamee and Jain (1976) explicit relationship for friction factor; 
	// - for Reynold's No. 5x10^3 <= Re <= 10^8
	// - for relative roughness 10^-6 <= epsilon/d <= 10^-2		
	inline double SwameeJain(double d, double Re){return 0.25/(pow(log10( this->epsilon/(3.7*d) + 5.74/pow(Re, 0.9) ), 2));}

	// Haaland (1983) explicit relationship for friction factor	
	inline double Haaland(double d, double Re){return pow(1 / (-1.8*log10( pow((this->epsilon/d)/3.7, 1.11) + 6.9/Re )), 2);}

	// Ricardo WAVE explicit relationship for friction factor	
	inline double Ricardo(double Re){return 4* /*convert to US ff*/ (0.027*pow(Re*0.1 /*make sure Re is based on 0.1d*/ ,-0.25));}

	double ConvectiveHTCoefficient(CProperties* pPpt, double f, double rho, double u, double T, double Re, double del_x);
	void PrintPathLines(){for(int k=0; k<this->num_pathlines; ++k){PathLine[k].Print(); cout << endl;}}
	void Update(CProperties* pPpt, double DELZ);
	void RecordTappings(CProperties* pPpt, int timestep, double time);
	
	void SetupFiles(CProperties* pPpt, string parent_assy_res_dir);
	//void PrintToFile(CProperties* pPpt, int timestep, double time, double ca);
	void PrintToFileLocations(CProperties* pPpt, int timestep, double time, double ca_elapsed, double ca);
	void PrintToFileFandE(CProperties* pPpt, int timestep, double time, double ca_elapsed, double ca);
	void PrintToFileMovie(CProperties* pPpt, int timestep, double time, double ca);
	void CloseFiles(CProperties* pPpt);
	
	void PrintToScreen(CProperties* pPpt);
	void PrintBoundaries(CProperties* pPpt);
	void ListProperties(CProperties* pPpt);
	
	void Backup(CProperties* pPpt);
	void Restore(CProperties* pPpt);
	bool Steady(CProperties* pPpt);

	char* Identify();
	double d(double x);
	double dddx(CProperties* pPpt, double x);
	double d2ddx2(CProperties* pPpt, double x);
	double f(double x);
	double dfdx(double x);

public:

	bool BUFFER;					// Is this a buffer pipe used for anechoic termination?
	bool DAMPER;					// Is this a damper pipe used for anechoic termination?

	// Pipe propagation method
	// ----------------------------------------------------------------------------------------------------
	int METHOD;						// Propagation method for this pipe: same as default (0), MMOC (1), W_ALPHA_BETA (2), filling and emptying (3)
	double delay;					// Time delay (s) from beginning of simulation to delay propagation

	CNode * Node;					// Vector of nodes
	CPathLine * PathLine;			// Vector of path line characteristics

	CNode * Node_Backup;			// Vector of nodes (backup)
	CPathLine * PathLine_Backup;	// Vector of path line characteristics (backup)

	double* Measurements;
	CNode* MeasureNode;

	// Filling and emptying
	// ----------------------------------------------------------------------------------------------------
	CFiniteVolume * FV;				// Vector of finite volumes

	CEngine * pEng;					// Engine pointer - points to the engine in the simulation if it exists

	int AssyID;						// ID of assembly to which pipe belongs
	char* RES_DIR;					// Results directory for this pipe
	FILE *FILE_LOC;					// Contains data at specified measurement locations
	FILE **FILE_OVERALL_MOV;		// For overall pipe measurements of the required properties
	FILE *FILE_FandE;				// Contains filling and emptying pipe volume data

	double *** Results;
	int rowcounter;
	int sample_factor;
	
	// Read from file
	// ----------------------------------------------------------------------------------------------------
	char** labels;					// Parameter list of labels in no order
	double* values;					// Parameter list of values in same order as labels
	char** strings;					// Parameter list of strings in same order as labels

	// ====================================================================================================
	// Parameter file for Exhaust Pipe [0]
	// ====================================================================================================

	char* strDesc; 					// Optional object description (max. 500 characters - use underscores for spaces)

	// Geometry
	// ----------------------------------------------------------------------------------------------------
	double length;					// The physical length of the system (m)
	double d_odd;					// Diameter at the LH end of the pipe (m)
	double d_even;					// Diameter at the RH end of the pipe (m)
	double bend_angle;				// Turning angle for a curved duct (degrees); range 0-180, default = 0
	int n_int_points;				// Number of internal points (i.e., excluding ends) at which to specify diameter (max. 4)
	double* xi;						// Location of internal diameters measured from odd end (mm)
	double* di;						// Internal pipe diameters (mm)
	bool LINEAR;					// Linear pipe diameter variation (1=true) or quadratic (0=false)
	bool LINEAR_F;					// Linear pipe area variation (1=true) or linear pipe diameter (0=false)
	double x_int;					// Distance from odd end of interior diameter (m)
	double d_int;					// Pipe diameter at some interior point (m)
	
	// Meshing
	// ----------------------------------------------------------------------------------------------------
	bool GLOBAL_MP;					// Use default meshing parameters (1=true) or those specified below (0=false)
	double discret;					// Mesh target discretization length (m)
	int min_meshes;					// Minimum number of meshes (overrides discret)
	
	// Coefficients
	// ----------------------------------------------------------------------------------------------------
	bool GLOBAL_EF;					// Use default roughness height and enhancement factors (1=true) or those specified below (0=false)
	double epsilon;					// Roughness height (m)
	double CFTRANS;					// Friction enhancement factor (common unsteady value = 1.875)				
	double HGTRANS;					// Heat transfer enhancement factor (common unsteady value = 2.625)
	double CPTRANS;					// Pressure loss enhancement factor (set to 0 to switch off)

	// Derived pipe geometry
	double C;						// Gradient of pipe diameter linear variation dD/dX (this is not c, C=dD/dx)
	//double a, b, c;					// Quadratric pipe diameter variation with x, i.e. dia = ax^2 + bx + c
	//double A, B, C;					// Quadratric pipe diameter variation with X, i.e. dia = AX^2 + BX + C
	double vol;						// Pipe volume (m^3)
	double C_p;						// Pressure loss coefficient due to bend loss

	// Pipe initialisation
	// ----------------------------------------------------------------------------------------------------
	double Tw;						// Wall temperature (K)
	int nsections;					// Number of different initialisation sections within pipe (>=1)
	double* xri;					// Location divisions as a fraction of pipe length measured from odd end 
	double* pri;					// Initial pipe pressure (bar)
	double* Tri;					// Initial pipe temperature (K)
	double* vri;					// Initial pipe velocity (K)
						
	// Measurements
	// ----------------------------------------------------------------------------------------------------
	int ntappings;					// Number of measuring locations
	bool USE_DEF_FREQ;				// Use default sampling rate (1=true) or 'freq' below (0=false)
	int freq;						// Print data to files every freq timesteps
//	bool PRINT_TO_FILE_MOV;			// Record data at all pipe nodes for movie (1=true) or not (0=false)
	double print_from_time;			// Only start printing data after this simulation time onwards (s)
	double* loc_measure;			// List of normalized positions along pipe of measuring location
	int max_pts;					// Number of data points per location to record over whole simulation
	int num_props_measured;			// Number of properties to measure, not counting time
	bool DIAMETER;					// Record diameter (1=true) or not (0=false)
	bool AREA;						// Record area (1=true) or not (0=false)
	bool STATIC_PRESSURE;			// Record static pressure (1=true) or not (0=false)
	bool TEMPERATURE;				// Record temperature (1=true) or not (0=false)
	bool DENSITY;					// Record density (1=true) or not (0=false)
	bool VELOCITY; 					// Record velocity (1=true) or not (0=false)
	bool MACH_NUMBER;				// Record Mach number (1=true) or not (0=false) 
	bool MASS_FLOW_RATE;			// Record mass flow rate (1=true) or not (0=false)
	bool REYNOLDS_NO;				// Record nodal Reynold's No. (1=true) or not (0=false)

	// Measuring locations
	// ----------------------------------------------------------------------------------------------------
	double* Pressure;				// Static pressure at each measuring location (bar)
	double* Temperature;			// Static temperature at each measuring location (K)
	double* Velocity;				// Velocity at each measuring location (m/s)
	double* MassFlowRate;			// Mass flow rate at each location (kg/s)
	double* Re;						// Reynolds number at each location
	double* PR;						// PR at each location
	double* MFP;					// MFP at each location

	// Screen output
	// ----------------------------------------------------------------------------------------------------
	bool SHOW_DATA;					// Show pipe data on screen (1=true) or not (0=false)

	// Grid properties
	// ===============
	bool EX;						// Exhaust or intake
	int ID;							// Id number
	int meshes;						// Number of meshes (per pipe)
	int N;							// Number of grid nodes (per pipe)
	int num_pathlines;				// The number of pathlines to track (per pipe)
	int MMOC, W_ALPHA_BETA, FandE, JOINER;	// Labels to represent MMOC (=1), W_ALPHA_BETA (=2), FandE (=3)

	double eff_length;				// The effective length of the system (i.e. including any end corrections) (m)
	double XPIPE;					// The N.D. physical length of the system ()
	double XPIPE_EFF;				// The N.D. effective length of the system ()

	double AREF;					// Reference speed of sound for pipe (ms^-1); will be either AREFe or AREFi
//	double Tref;					// Reference temperature for the pipe (K)
	double xmesh;					// Length of a mesh (m)
	double XMESH;					// N.D. length of a mesh = xmesh/xref ()

	int odd_end_flow, even_end_flow;				// Pipe end flow directions
	int odd_end_flow_backup, even_end_flow_backup;	// Pipe end flow directions (backup)

	double end_corr_odd_p, end_corr_even_p;	// End correction parameters (Lc/d)
	double end_corr_odd, end_corr_even;		// End corrections (Lc/d)*d (m)


	// Filling and emptying
	//double CL1_old, CL2_old;
	//double betaCL1, betaCL2;
//	char* direction_str_odd;
//	char* direction_str_even;

	// Outflow
	//double ps_throat_odd, ps_throat_even;

	// Inflow
	//double p0_throat_odd, p0_throat_even, T0_throat_odd, T0_throat_even;

	//double lowest_error_odd;

	//double tol_odd, tol_even;
	
	
	// New function
	double* A_throat;
	double* U_throat;
	double* A_throat_old;
	double* U_throat_old;
	double* A_error;
	double* U_error;
	double* A_error_old;
	double* U_error_old;
	bool* CONVERGED;
	double* beta;
	double* tol;
	char** direction_str;

	double* ps_throat;
	double* ps_throat_old;
};

#endif // !defined(AFX_PIPE_H__C27D3667_10DC_401D_A116_3140F7A8D453__INCLUDED_)
