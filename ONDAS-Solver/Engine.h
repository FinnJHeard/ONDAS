// Engine.h: interface for the CEngine class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ENGINE_H__D7805F5A_BCEC_4BB7_90AD_32160B1592E0__INCLUDED_)
#define AFX_ENGINE_H__D7805F5A_BCEC_4BB7_90AD_32160B1592E0__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Properties.h"
//#include "Cylinder.h"
#include "Time.h"

class CEngine  
{
public:
	CEngine();
	virtual ~CEngine();
	CEngine(const CEngine& inEng);				// Copy constructor
	CEngine& operator=(const CEngine& inEng);	// Overloaded operator=

//	void ReadProperties(CProperties* pPpt, int id);
	void Initialise(CProperties* pPpt, int id, std::string param_dir, char* eng_file, int ncyls, int ninpipes, int assyid, string parent_assy_res_dir, string calling_object_str);
	void ReadInput(CProperties* pPpt, char *InputFile, int ninpipes);
	void ListProperties(CProperties* pPpt);
	void PrintToScreen();
	void RunBoundary(CProperties* pPpt, CTime &rMyTime, double DELZe, double DELZi, double time);
	void Update(CTime &rMyTime, double DELZe, double DELZi, CProperties* pPpt, double time);
	char* Identify();

	// Inlines
//	inline int Get_NEXVALVES(){return NEXVALVES;}
//	inline int Get_NINVALVES(){return NINVALVES;}

public:

	// Identification
	// ----------------------------------------------------------------------------------------------------
	char* strDesc; 			// Optional object description (max. 500 characters - use underscores for spaces)
	int ID;					// Engine ID
	char* ENG_FILE;			// Engine parameter file
	int NCYLS;				// Number of cylinders
	int AssyID;				// ID of assembly on which boundary belongs

	// Read from file
	// ==============
	char** labels;			// Parameter list of labels in no order
	double* values;			// Parameter list of values in same order as labels
	char** strings;			// Parameter list of strings in same order as labels

	// ====================================================================================================
	// Parameter file for Engine []
	// ====================================================================================================

	// Operation
	// ----------------------------------------------------------------------------------------------------
	double ca_start;		// Starting crank angle relative to TDCF during cycle (0-720 for 4-stroke cycle) (degrees)
	double cycle;			// 2 or 4-stroke cycle
	double reveng;			// Engine speed (rev/min)

	// Cylinder geometry
	// ----------------------------------------------------------------------------------------------------
	double dcyl;			// Bore (mm)
	double stroke;			// Stroke (mm)	
	double conrod;			// Connecting rod length (mm)
	double cr;				// Compression ratio
	double EVO;				// Exhaust valve opening (degrees CA) - taken from first exhaust valve VO
	double IVO;				// Intake valve opening (degrees CA) - taken from first intake valve VO
	double IVC;				// Intake valve closing (degrees CA) - taken from first intake valve VC

	// Special cases
	// ----------------------------------------------------------------------------------------------------
	bool IVOL;				// Variable volume (1=true) or constant volume (0=false) cylinder model
	bool IAIR;				// Intake valve present (1=true) or no intake to cylinder (0=false)
	bool IPIPE;				// Intake valve connects to pipe (1=true) or air receiver (0=false)
							
		// Air receiver (IPIPE=false)
		// ----------------------------------------------------------------------------------------------------
		double PA;				// Receiver air pressure (bar)
		double TA;				// Receiver air temperature (K)

	// Cylinder model
	// ----------------------------------------------------------------------------------------------------
	int MODEL;				// Select (0) constant conditions cylinder, (1) constant volume cylinder, (2) P & T read from file, (3) P & T reset at EVO, (4) Watson single-zone

		// (0) Constant conditions (cylinder acts as a reservoir)
		// ----------------------------------------------------------------------------------------------------
		double pres;			// Constant cylinder pressure (bar)
		double Tres;			// Constant cylinder temperature (K)
		double Vres;			// Constant cylinder volume (litres)

		// (1) Constant volume, variable conditions (cylinder acts as a volume)
		// ----------------------------------------------------------------------------------------------------	
		double vconst;			// Constant cylinder volume (litres)	
		double pinit;			// Initial cylinder pressure (bar)
		double Tinit;			// Initial cylinder temperature (K)

		// (2) Read cylinder pressure & temperature from file
		// ----------------------------------------------------------------------------------------------------
		char* cyl_file;			// The file name from which to read the cylinder P & T trace

		// (3) Cylinder pressure & temperature reset at every EVO or TDC event
		// ----------------------------------------------------------------------------------------------------
		bool RESET_EVO;			// Reset cylinder conditions at EVO (1=true) or at TDC (0=false)
		double pcr;				// Cylinder pressure at reset (e.g. EVO or TDC) (bar)
		double Tcr;				// Cylinder temperature at reset (e.g. EVO or TDC) (K)

		// (4) Watson single-zone heat release model
		// ----------------------------------------------------------------------------------------------------
		int x;					// For a fuel composition CxHy
		int y;					// For a fuel composition CxHy
		double phiFAstoich;		// Stoichiometric fuel/air ratio
		double lhv;				// Fuel LHV (J/kg)
		double ign_ca;			// Ignition (degrees CA)
		double comb_stop;		// Point at which combustion is assumed to be complete (degrees CA)
		double inj_start;		// Start of fuel injection (degrees CA)
		double inj_stop;		// End of fuel injection (degrees CA)
		double m_f0;			// Total fuel mass injected per cycle per cylinder (kg)
		double comb_a;			// 0.8 to 0.95
		double comb_b;			// 0.25 to 0.45
		double comb_c;			// 0.25 to 0.5
		bool ADD_MFB;			// Add mass of fuel burned to cylinder mass (1=true) (0=false)
	
	// Cylinder heat transfer model
	// ----------------------------------------------------------------------------------------------------
	int HT_MODEL;			// Select (0) none, (1) Annand convective heat transfer, (2) Woschni correlation

		// (1) Annand convective heat transfer
		// ----------------------------------------------------------------------------------------------------
		double a;				// 0.2 - 0.8 depending on engine type
		double b;				// Usual to set b = 0.7
		double Tw;				// Cylinder wall temperature (K)

		// (2) Woschni correlation
		// ----------------------------------------------------------------------------------------------------
		double T_c;				// Cylinder coolant temperature (K)
		double h_cc;			// Coolant heat transfer coefficient (W/(m^2.K))
		double k;				// Liner wall thermal conductivity (W/(m.K))
		double tw;				// Cylinder liner wall thickness (m)
		double epsilon;			// Emissivity; approx. empirical value at peak cylinder conditions
		double coeffs;			// A combination of heat transfer coefficients = k*h_cc/(k + h_cc*tw)
		double WOS_TOL;			// Woschni heat transfer radiative loop tolerance

	// Valve and general AVT parameters
	// ----------------------------------------------------------------------------------------------------
	bool PRE_vt;			// Show valve timing in preamble (1=true) or not (0=false)
	int nSamples;			// The number of data points to record to represent a pulse
	double tol;				// The tolerance (%) within which a pulse match will occur
	bool SUSPEND;			// Suspend simulation when match is achieved (1==true) or not (0==false)
	
	int	nExPipesInAssy;		// Number of exhaust pipes in corresponding assembly
	int	nInPipesInAssy;		// Number of intake pipes in corresponding assembly
	double period;			// Time period of an engine cycle

	// Exhaust valve parameters
	// ----------------------------------------------------------------------------------------------------
	int NEXVALVES;			// Number of exhaust valves per cylinder
	double REF_DIA_EXH;		// Exhaust valve reference diameter (m)
	int EV_WAVEFORM;		// Waveform: switch (open/close) = 0, square = 1, triangular = 2, sinusoidal  = 3, Fourier series = 4, file = 5
	int EV_L_OR_A;			// Valve opening profile (whether calculated or from file) is: lift (0), geometric area (1), psi (valve to pipe area ratio) (2)
	double EV_MIN_OPENING;		// Exhaust valve minimum opening: lift (mm) or area (m^2) or psi value (valve to pipe area ratio)
	double EV_MAX_OPENING;		// Exhaust valve maximum opening: lift (mm) or area (m^2) or psi value (valve to pipe area ratio)
	
	// If exhaust valve acts as a switch (open/close), i.e., EV_WAVEFORM == 0
	// ----------------------------------------------------------------------------------------------------
	bool EV_OPEN;			// Exhaust valve is fully open (1=true) or fully closed (0=false)
	
	// If exhaust valve opening profile is square or triangular or sinusoidal or a Fourier series, i.e., EV_WAVEFORM == 1 or 2 or 3 or 4
	// ----------------------------------------------------------------------------------------------------
	double EV_BEGIN_OPENING;	// Exhaust valve begin opening (degrees)
	int EV_END_or_phi;			// Use EV_END angle (0) or phi (1) to determine closing point 
	double EV_END_CLOSING;		// Exhaust valve end closing (degrees)
	double EV_phi;					// Valve opening duration as fraction of wavelength (0-1)

	// If exhaust valve opening profile is read from file, i.e., EV_WAVEFORM == 5
	// ----------------------------------------------------------------------------------------------------
	char* EV_FILE;			// Exhaust valve data file
	double CAM_TNG_ANG_EXH;	// Cam timing angle (degrees)
	int CAM_TNG_LOC_EXH;	// Cam timing angle locates start of valve opening (0), maximum opening (1), or end of valve closing (2)

	// If exhaust valve data are lift values, i.e., EV_L_OR_A == 0
	// ----------------------------------------------------------------------------------------------------
	double LASH_EXH;		// Exhaust valve lash (m)
	char* EVFA_FILE;		// Exhaust valve flow array file
	bool CONST_REF_AREA_EXH;// Use constant reference area (i.e., values are flow coefficient) (1) or curtain reference area (i.e., values are flow coefficient) (0) to interpret file values
	double ONUM_EXH;		// Multiplier to the flow or discharge coefficient arrays ("def"=1)

	// Else exhaust valve data are area values, i.e., EV_L_OR_A == 1
	// ----------------------------------------------------------------------------------------------------
	
//	// Exhaust valve flow or discharge coefficients (can only be used with lift data since l/d is the lookup)
//	// ----------------------------------------------------------------------------------------------------
	
	// Variable valve timing and lift
	// ----------------------------------------------------------------------------------------------------
	double EV_OPENING_MULT;	// Valve opening (phi value or lift or area) multiplier
	double ANG_MULT_EXH;	// Valve opening duration multiplier; scaling factor for theta values
	double scale_factor_ev;	// Valve opening duration compression factor

	// Active Valve Train
	// ----------------------------------------------------------------------------------------------------
	bool EV_ACTIVE;			// Active exhaust valve (1=true) or not (0=false)
	bool EV_TARGET_PIPE_EX; // Target location is in exhaust system (1=true) or intake system (0=false)
	int EV_target_pipe_num;		// Number of pipe to target
	double EV_target_pipe_loc;	// Location of target tapping along target pipe as a fraction of its pipe length
	char* EV_ACTIVE_FILE_P; // Target pressure file
	
	// Intake valve parameters
	// ----------------------------------------------------------------------------------------------------
	int NINVALVES;			// Number of intake valves per cylinder
	double REF_DIA_INT;		// Intake valve reference diameter (m)
	int IV_WAVEFORM;		// Waveform: switch (open/close) = 0, square = 1, triangular = 2, sinusoidal  = 3, Fourier series = 4, file = 5
	int IV_L_OR_A;			// Valve opening profile (whether calculated or from file) is: lift (0), effective flow area (1), psi (valve to pipe area ratio) (2)
	double IV_MIN_OPENING;		// Intake valve minimum opening: lift (mm) or effective flow area (m^2) or psi value (valve to pipe area ratio)
	double IV_MAX_OPENING;		// Intake valve maximum opening: lift (mm) or effective flow area (m^2) or psi value (valve to pipe area ratio)

	// If intake valve acts as a switch (open/close), i.e., EV_WAVEFORM == 0
	// ----------------------------------------------------------------------------------------------------
	bool IV_OPEN;			// Intake valve is fully open (1=true) or fully closed (0=false)
	
	// If intake valve opening profile is square or triangular or sinusoidal or a Fourier series, i.e., IV_WAVEFORM == 1 or 2 or 3 or 4
	// ----------------------------------------------------------------------------------------------------
	double IV_BEGIN_OPENING;	// Intake valve begin opening (degrees)
	int IV_END_or_phi;			// Use IV_END angle (0) or phi (1) to determine closing point 
	double IV_END_CLOSING;		// Intake valve end closing (degrees)
	double IV_phi;				// Intake valve opening duration as fraction of wavelength (0-1)

	// If intake valve opening profile is read from file, i.e., EV_WAVEFORM == 5
	// ----------------------------------------------------------------------------------------------------
	char* IV_FILE;				// Intake valve data file
	double CAM_TNG_ANG_INT;		// Cam timing angle (degrees)
	int CAM_TNG_LOC_INT;		// Cam timing angle locates start of valve opening (0), maximum opening (1), or end of valve closing (2)

	// If intake valve data are lift values, i.e., IV_L_OR_A == 1
	// ----------------------------------------------------------------------------------------------------
	double LASH_INT;		// Intake valve lash (m)
	char* IVFA_FILE;		// Intake valve flow array file
	bool CONST_REF_AREA_INT;// Use constant reference area (i.e., values are flow coefficient) (1) or curtain reference area (i.e., values are flow coefficient) (0) to interpret file values
	double ONUM_INT;		// Multiplier to the flow or discharge coefficient arrays ("def"=1)

	// Else intake valve data are area values, i.e., EV_L_OR_A == 0
	// ----------------------------------------------------------------------------------------------------

//	// Intake valve flow or discharge coefficients (can only be used with lift data since l/d is the lookup)
//	// ----------------------------------------------------------------------------------------------------

	// Variable valve timing and lift
	// ----------------------------------------------------------------------------------------------------
	double IV_OPENING_MULT;	// Valve opening (phi value or lift or area) multiplier
	double ANG_MULT_INT;	// Valve opening duration multiplier; scaling factor for theta values
	double scale_factor_iv;	// Valve opening duration compression factor

	// Active Valve Train
	// ----------------------------------------------------------------------------------------------------
	bool IV_ACTIVE;			// Active intake valve (1=true) or not (0=false)
	bool IV_TARGET_PIPE_EX; // Target location is in exhaust system (1=true) or intake system (0=false)
	int IV_target_pipe_num;		// Number of pipe to target
	double IV_target_pipe_loc;	// Location of target tapping along target pipe as a fraction of its pipe length
	char* IV_ACTIVE_FILE_P; // Target pressure file
	
	// Simulation variables
	// ----------------------------------------------------------------------------------------------------
	double del_theta;		// Increment in degrees CA
	double rotation;		// Angle during revolution (will vary 0-360) relative to TDC (top-dead-centre) (degrees)
	double ca;				// Crank angle during cycle (will vary 0-720 for 4-stroke cycle) relative to TDCF (top-dead-centre-firing) (degrees)
	double ca_old;			// Stores previous value of ca (degrees)
	double ca_elapsed;		// Number of degrees elapsed relative to TDCF at start of simulation, e.g. 14400 degs = 20 full four-stroke cycles (degrees)
	double degrees_elapsed; // Number of degrees elapsed since start of simulation = ca_elapsed - start_ca (degrees)
	int rev;				// Ordinal of current rev. 1st revolution => rev = 1
	bool NEW_CYCLE;			// Flag set to true during the first iteration of a new cycle

	// Overall engine work and power
	// ----------------------------------------------------------------------------------------------------
	double del_WK_total;	// Work done during the present iteration (J)
	double P_inst;			// Instantaneous power (W)
	double CYCLE_START_TIME;// Current cycle start time (s)
	double CYCLE_WK;		// Cumulative cycle work (J)
	double PREV_CYCLE_WK;	// Previous cycle work (J)
	double PREV_CYCLE_POWER;// Power calulated over most recently complete cycle (W)
	double CYCLE_FUEL;		// Cumulative cycle fuel (injected) (kg)
	double PREV_CYCLE_FUEL;	// Previous cycle fuel (injected) (kg)
	double eta;				// Thermal efficiency
	double PREV_CYCLE_SFC;	// (Indicated) specific fuel consumption (g/kWh)
};

#endif // !defined(AFX_ENGINE_H__D7805F5A_BCEC_4BB7_90AD_32160B1592E0__INCLUDED_)

// REMEMBER TO ADD NEW VARIABLES TO THE OVERLOADED OPERATORS