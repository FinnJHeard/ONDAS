// Properties.h: interface for the CProperties class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PROPERTIES_H__ED758614_BA0B_4333_A438_18F179BCD211__INCLUDED_)
#define AFX_PROPERTIES_H__ED758614_BA0B_4333_A438_18F179BCD211__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//#include <string.h>
//#include <vector>

#include <string>
#include <vector>

using std::string;
using std::vector;

class CProperties  
{
public:
	CProperties();
	virtual ~CProperties();

	void ReadCaseName(char* InputFile, int num_parameters);
	void ReadInput(char* InputFile, int num_parameters);
	
	void Out(std::string text);
	void Out(const char* text);
	void Out(char* text);
	void Out(char text);
	void Out(int text);
	void Out(double text);
	void Out(double text, int precision);

	void SetParams(char* InputFile);
	void ListProperties(int CODE_VERSION);

	void Configure(char *InputFile, char** &rCONFP, double* &rENDCORR, int** &rASSYS, int** &rPIPES, int** &rPIPES_ENDS, int* &rNUM_PIPES, int number_of, int strlen, char* bc_type, int npipes, vector<string> &rConfStrs, int assyid);
	void ReadConfigFile(char *InputFile, char** &rCONFP, double* &rENDCORR, int number_of, int strlen, char* bc_type, int &rno_identified, int npipes);
	void ConvertStrings(char *InputFile, char** &rCONFP, int** &rPIPES, int** &rASSYS, int** &rPIPES_ENDS, int* &rNUM_PIPES, int number_of, char* bc_type, int no_identified, int npipes, vector<string> &rConfStrs, int assyid);
	void SetupFiles(void);
	void CreateFileNames(char** &rFileName, int num_objects, char* object_type);
	void CreateFileNames2(char*** &rFileName, int num_objects, int num_locs, char* object_type);
	void ReadLossCoefficients(char *InputFile);

	double ViscosityAir(double T);
//	double ConvectiveHeatTransferCoefficient(double f, double rho, double u, double Re, double del_x);
	double cpAir(double T);
	double cvAir(double T);
	double gammaAir(double T);

	// Haemodynamics
	double ViscosityBloodFunc(double T);
	double BLThickness(double diameter);

public:

	// ====================================================================================================
	// Main input parameter file
	// ====================================================================================================

	// Case description
	// ----------------------------------------------------------------------------------------------------
	char* strDesc; 				// Case description (max. 500 characters - use underscores for spaces)

	// Code version
	// ----------------------------------------------------------------------------------------------------
	int	FILE_VERSION;			// File version (yyyymmdd) - allows mismatch between input files and code to be flagged

	// Parameter file control
	// ----------------------------------------------------------------------------------------------------
	int NUM_PARAMS;				// Max. number of parameters that can be read from a parameter file
	bool SHOW_globals;			// Show the global parameters in preamble (1=true) or not (0=false)

	// Simulation control
	// ----------------------------------------------------------------------------------------------------
	double ZMAX;				// The maximum simulation time (s)
	double courant;				// Courant number
	double delt_min;			// Minimum and first time step (s)
	double delt_max;			// Maximum time step (s)
	bool CONSTANT_DELT;			// Constant delt_max timestep (1=true) or program calculated accordinging to COURANT (0=false)
	bool CONTINUOUS;			// Continuous (1) or periodic (0) simulation. If zero cylinders, will be continuous.
	bool SUPERSONIC;			// Permit supersonic flow without restore (1=true) or not (0=false)
	bool BEEP;					// Sound beep for user input/at end of simulation (1=true) or not (0=false)
	bool HAEMODYNAMICS;			// 0 use gas dynamic models, 1 use haemodynamic models.
	
	// Propagation method
	// ----------------------------------------------------------------------------------------------------
	int DEF_METHOD;				// Default propagation method: MMOC (1), W_ALPHA_BETA (2), filling and emptying (3)
	
		// If mesh method of characteristics (MMOC), i.e., DEF_METHOD = 1 
		// ----------------------------------------------------------------------------------------------------
		bool HOMENTROPIC;			// Homentropic simulation (true == 1, false == 0)

			// If non-homentropic MMOC, i.e., HOMENTROPIC = 0
			// ----------------------------------------------------------------------------------------------------
			int	NUM_PATH_MULT;			// Number of pathlines multiplier; number of pathlines = NUM_PATH_MULT*(2*Nodes - 1); min = 2*Nodes - 1

		// If W(alpha,beta) scheme, i.e., DEF_METHOD = 2: MacCormack (1,0), Two-step L-W (0.5,0.5), Lerat & Peyret optimum (1+sqrt(5)/2,0.5)
		// ----------------------------------------------------------------------------------------------------
		bool COMBINED_WAB_MOC;		// Continue to propagate MMOC pathlines (1=true) or extrapolate solution vector to boundary nodes (0=false)
		bool SOURCES;				// Include source terms (should normally be on) (1=true) or not (0=false)
		bool WORK;					// Include work in source terms (1=true) or not (0=false)
		double alpha;				// alpha value
		double beta;				// beta value
		bool ALTERNATE_MAC;			// MacCormack only - alternate forward pred./backward corrector, backward pred./forward corr. (1=true) or not (0=false)
		bool VISC;					// Apply artificial viscosity (1=true) or not (0=false)
		double Cx;					// Artificial viscosity coefficient (pressure based)
		double Cx_alpha;			// Simple artificial viscosity coefficient
		bool TVD;					// Run scheme with (1) or without (0) TVD criterion

		// (3) Filling and emptying
		// ----------------------------------------------------------------------------------------------------

	// Gas properties
	// ----------------------------------------------------------------------------------------------------
	bool COMPRESSIBLE;			// Use compressible (isentropic) or incompressible thermodynamic relations
	double R_universal;			// Universal gas constant (J.mol^-1.K^-1)
	double M_air;				// Molecular mass of air (kg.mol^-1)
	bool constProps;			// Use constant thermodynamic properties below (1=true) or Zucrow & Hoffman function of temperature (0=false)
	double Cp_air;				// Constant value of specific heat at constant pressure for air (J.kg^-1.K^-1)

	// Source term models
	// ----------------------------------------------------------------------------------------------------
	int FRICTION_MODEL;			// Method of calculating friction factor; constant (0), Swamee and Jain (1), Haaland (2), Colebrook-White implicit (3)
	double ff_const;			// Value of friction factor when constant friction factor is specified above; normal range 0.0035-0.008
	int VISCOSITY_AIR;			// Method of calculating absolute viscosity of air; constant (0), Sutherland (1), Blair (2)
	double mu_air_const;		// Value of absolute viscosity to use when constant viscosity is specified (kg.m^-1.s^-1)
	int HEAT_TRANSFER;			// Method of calculating convective heat transfer; Reynolds' analogy (0), Nusselt Relation (1)
	double tol_steady;			// Percentage tolerance on nodal velocities when checking for steady flow (%)

	// Haemodynamic properties (k values from Olufsen (1999), "Structured tree outflow condition for blood flow in larger systemic arteries")
	// ----------------------------------------------------------------------------------------------------
	double rho_blood;			// Density of blood [kg.m-3]
	double ViscosityBlood;		// Dynamic viscosity of blood [Pa.s]
	double k1;					// Fitting coefficient k1 for large vessels [kg.s-1.m-1]
	double k2; 					// Fitting coefficient k2 for large vessels [m-1]
	double k3; 					// Fitting coefficient k3 for large vessels [kg.s-1.m-1]
	bool BL_RELATIVE;			// Boundary layer thickness takes relative (1=true) or fixed value (0=false)

	// Boundary layer thickness 
	// ----------------------------------------------------------------------------------------------------
	double BL_rel;			// Relative boundary layer thickness as a fraction of pipe diameter (when BL_RELATIVE=1)
	double BL_abs;			// Fixed absolute boundary layer thickness [m] (when BL_RELATIVE=0)

	// Default values
	// ----------------------------------------------------------------------------------------------------
	double discret;				// Global/default mesh target discretization length (mm)
	int min_meshes;				// Global/default minimum number of meshes (overrides discret)
	double epsilon;				// Global/default pipe surface roughness height (mm)
	double CFTRANS;				// Global/default friction enhancement factor (common unsteady value = 1.875)				
	double HGTRANS;				// Global/default heat transfer enhancement factor (common unsteady value = 2.625)
	double CPTRANS;				// Global/default pressure loss enhancement factor (set to 0 to switch off)
	int freq;					// Default data sampling rate; once per freq timesteps

	// Reference values
	// ----------------------------------------------------------------------------------------------------
	double PREF; 				// Reference (outside) pressure (bar)
	double xref;				// Reference length (m). Only one xref!
	double TREFe;				// Exhaust manifold reference temperature (K)
	double TREFi;				// Intake manifold reference temperature (K)
	double fref;				// Reference pipe area (m^2)

	// Boundary conditions
	// ----------------------------------------------------------------------------------------------------
	int PARENT;					// Special integer identifier that when used implies boundary is located on parent assembly
	int STRLEN;					// Common length in characters of configuration sub-string (per pipe) (e.g. 4 for 123e124e125o)
	double ZERO_TOL;			// Expressions within this tolerance can be considered as 0 (used when testing equality)

	// Screen output
	// ----------------------------------------------------------------------------------------------------
	double outputInterval;		// Time interval between printing simulation progress output to screen (s)
	bool SHOW_conn;				// Show the odd and even end boundary connections for each pipe (1=true) or not (0=false)
	bool SHOW_config;			// Show the interpretation of the boundary configuration files in preamble (1=true) or not (0=false)
	bool SHOW_params;			// Show the simulation object parameters in preamble (1=true) or not (0=false)
	bool SHOW_calls;			// Show function calls in screen output (1=true) or not (0=false)
	bool SHOW_x;				// Show pipe x (m) values in screen output (1=true) or not (0=false)
	bool SHOW_X;				// Show pipe X () values in screen output (1=true) or not (0=false)
	bool SHOW_DELX;				// Show pipe DELX_L and DELX_R values in screen output (1=true) or not (0=false)
	bool SHOW_d;				// Show pipe d (m) values in screen output (1=true) or not (0=false)
	bool SHOW_dddx;				// Show pipe dddx (m) values in screen output (1=true) or not (0=false)
	bool SHOW_d2ddx2;			// Show pipe d2ddx2 (m) values in screen output (1=true) or not (0=false)
	bool SHOW_f;				// Show pipe f (m^2) values in screen output (1=true) or not (0=false)
	bool SHOW_cfa;				// Show pipe cfa (m^2) values in screen output (1=true) or not (0=false)
	bool SHOW_cfa_delx;			// Show pipe cfa_delx (m) values in screen output (1=true) or not (0=false)
	bool SHOW_dfdx;				// Show domain dfdx (m) values in screen output (1=true) or not (0=false)
	bool SHOW_vol;				// Show domain volume (m) values in screen output (1=true) or not (0=false)
	bool SHOW_W;				// Show domain solution vector W[] values in screen output (1=true) or not (0=false)
	bool SHOW_F;				// Show domain flux vector F[] values in screen output (1=true) or not (0=false)
	bool SHOW_C;				// Show domain source vector C[] values in screen output (1=true) or not (0=false)
	bool SHOW_W_pred;			// Show domain predicted solution vector W_pred[] values in screen output (1=true) or not (0=false)
	bool SHOW_F_pred;			// Show domain predicted flux vector F_pred[] values in screen output (1=true) or not (0=false)
	bool SHOW_C_pred;			// Show domain predicted source vector C_pred[] values in screen output (1=true) or not (0=false)
	bool SHOW_pathlines;		// Show pathlines in screen output (1=true) or not (0=false)

	// File output
	// ----------------------------------------------------------------------------------------------------
	string strFileExt; 			// File extension for results files, e.g., .orf	or .txt

	// Engines, cylinders & valves
	// ----------------------------------------------------------------------------------------------------
	char* CYL_DIR;				// Folder containing cylinder P & T data files - must include the "\" following the directory name
	double POPPET_H_TOL1;		// Homentropic poppet valve code: loop tolerance 1
	double POPPET_H_TOL2;		// Homentropic poppet valve code: loop tolerance 2
	char* VT_DIR;				// Folder containing valve timing files - must include the "\" following the directory name

	// End environments	
	// ----------------------------------------------------------------------------------------------------
	char* ENDENV_DIR;			// Folder containing stagnation value files - must include the "\" following the directory name
	double NOZZ_TOL;			// Homentropic and non-homentropic nozzle code: loop tolerance
	double NHI_TOL;				// Non-homentropic inflow code: loop tolerance
	bool USE_PHI;				// Take account of area ratio (1=true) or treat as fully open (0=false) under INFLOW

	// Junctions
	// ----------------------------------------------------------------------------------------------------
	char* JUNC_DIR;				// Folder containing junction files, e.g. loss coefficients - must include the "\" following the directory name
	int max_branches;			// Maximum number of branches per junction (required in order to calculate max string length)
	double ****Loss;			// Matrix of loss coefficients for pulse converter model

	// Adiabatic pressure loss devices
	// ----------------------------------------------------------------------------------------------------
	char* APLD_DIR;				// Folder containing loss files - must include the "\" following the directory name
	double LOSS_TOL;			// Loss loop tolerance

	// Sudden area changes
	// ----------------------------------------------------------------------------------------------------
	double tolSuddenEnlrgN;		// Convergence tolerance on variable N in sudden enlargement function
	double tolSuddenEnlrgLin;	// Convergence tolerance on variable lambda_in in sudden enlargement function

	// Turbines
	// ----------------------------------------------------------------------------------------------------
	char* TURB_DIR;				// Folder containing turbines maps - must include the "\" following the directory name

	// Transmissive boundaries
	// ----------------------------------------------------------------------------------------------------
	char* TRANSM_DIR;			// Folder containing input wave files - must include the "\" following the directory name

	// Assemblies
	// ----------------------------------------------------------------------------------------------------
	int NASSEMBLY;				// Number of assemblies in simulation

	// ====================================================================================================
	// End of file
	// ====================================================================================================

	// File names and directories
	char* case_name;			// Name of directory containing main input and assembly directory files for the current case
	char* case_name_slash;		// Name of directory containing main input and assembly directory files for the current case plus slash
	std::string	assemblies_dir;	// Path to assemblies directory
	std::string cases_dir;		// Path to cases directory
	std::string param_dir;		// Path to parameter files directory
	std::string case_dir;		// Path to current case directory
	std::string case_res_dir;	// Path to current case directory
	
	FILE* OUTFILE;		// Simulation case .out file

	// Read in arrays for input files
	char** labels;				// Parameter list of labels in no order
	double* values;				// Parameter list of values in same order as labels
	char** strings;				// Parameter list of strings in same order as labels


	int CONSTANT_FF, SWAMEEJAIN_FF, HAALAND_FF, COLEBROOK_FF, RICARDO_FF; // Labels for the above
	int CONSTANT_VISCOSITY, SUTHERLAND_VISCOSITY, BLAIR_VISCOSITY, RICARDO_VISCOSITY; // Labels for the above
	int REYNOLDS_ANALOGY, NUSSELT_RELATION, RICARDO_HT; // Labels for the above

	double DELZe;				// Dimensionless time step - exhaust manifold
	double DELZi;				// Dimensionless time step - intake manifold

	double AREFe;				// Exhaust manifold reference speed of sound (m/s)
	double AREFi;				// Intake manifold reference speed of sound (m/s)	
	double tol_steady_multiplier; // Multiplier that widens steady tolerance (normally 1)

	// Derived gas properties
	double R_air;				// Gas constant for air (J.kg^-1.K^-1)
	double Cv_air;				// Constant value of specific heat at constant volume for air (J.kg^-1.K^-1)
	double gamma_air;			// Constant ratio of specific heats


	// Constants
	// ======
	double AA;
	double BB;
	double Q;		
	double QI;		
	//double Eta;
	double sigma;		// Stefan-Boltzman constant (for cylinder radiative heat transfer)
	
	bool STOP;			// Signals simulation to offer to stop (e.g., in case of a transmissive pulse match)

	// Output files
	FILE **EX_J_FILE;

	double K_entropy;
};

#endif // !defined(AFX_PROPERTIES_H__ED758614_BA0B_4333_A438_18F179BCD211__INCLUDED_)
