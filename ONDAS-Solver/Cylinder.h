// Cylinder.h: interface for the CCylinder class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CYLINDER_H__8C8B53E5_114C_4A94_BAB5_3B763BA669B8__INCLUDED_)
#define AFX_CYLINDER_H__8C8B53E5_114C_4A94_BAB5_3B763BA669B8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Engine.h"
#include "PathLine.h"
#include "Pipe.h"
#include "Properties.h"
#include "Node.h"
#include "Valve.h"

class CCylinder : public CBoundary
{
public:
	CCylinder();
	virtual ~CCylinder();
	CCylinder(const CCylinder& inCyl);				// Copy constructor
	CCylinder& operator=(const CCylinder& inCyl);	// Overloaded operator=

	void Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rExhaustPipe, /*CPipe* &rIntakePipe,*/ 
		/*int** &rCYLASSYS,*/ int** &rCYLPIPES, int** &rCYLPIPES_ENDS, double* &rENDCORR, 
		int id, bool ex, int npipes, CEngine* EngPtr, std::string param_dir/*, char* cyl_dir*//*, char* vt_dir*/, 
		int assyid, string parent_assy_res_dir, string calling_object_str, int nExPipesInAssy, int nInPipesInAssy);
	
	void Configure(CProperties* pPpt, CPipe** pExhaustSystem, CPipe** pIntakeSystem, int nExPipes, int nInPipes);
	
	void ReadInput(CProperties* pPpt, char *InputFile);
	void ListProperties(CProperties* pPpt);

	void RunBoundary(CProperties* pPpt, CTime &rMyTime, double DELZe, double DELZi, double ca, double del_ca, int timestep, double time);
	void Update(CTime &rMyTime, double DELZe, double DELZi, CProperties* pPpt, /*double ca, double del_ca,*/ int timestep, double time);
	
	double CombustionDQDT(CProperties* pPpt, double DT);
	double SingleZoneWatson(CProperties* pPpt);
	bool AllClosed(void);

	void LoadCylinderData(char* InputFile);
	void Interpolate(CProperties* pPpt, double ca);

	void AirReceiver(CProperties* pPpt, double& rDMIDT, double PAIR, double TAIR, double PCYL, double TCYL, double AAIR, double ACYL, double F2);

	void PrintToScreen(CProperties* pPpt);
	void PrintToFile(CProperties* pPpt, int timestep, double time, double ca);
	void PrintToFileMovie(CProperties* pPpt, int timestep, double time, double ca);
	void SetupFiles(CProperties* pPpt, string parent_assy_res_dir);
	void CloseFiles();
	char* Identify();

//	void CylinderUpdate(double dt,int cylnum,CPipe* InPipe,CPipe* ExPipe, CProperties* pPpt, bool stopped);

public:
	CEngine* EnginePtr;	// Pointer to engine on which cylinder is located

private:
	int ID;
	
	int AssyID;				// ID of assembly on which boundary belongs
	char* RES_DIR;			// Results directory for this boundary
	
	// Read from file
	// ==============
	char** labels;		// Parameter list of labels in no order
	double* values;		// Parameter list of values in same order as labels
	char** strings;		// Parameter list of strings in same order as labels
	
	double offset;		// Cylinder offset from TDCF (degrees CA)
	int datapoints_cyl;	// Number of rows in cylinder data array
	double** cyl_data;	// Array containing cylinder pressure and temperature data

	// Working variables
	// =================
	double THETA;		// Cylinder angle after TDCF (0-720 for 4-stroke) (degrees CA)
	double PC;			// Cylinder pressure (bar)
	double TC;			// Cylinder temperature (K)
	double VC;			// Cylinder volume (m^3)
	double MC;			// Mass in cylinder (kg)
	double AC;			// Speed of sound in cylinder (m/s)
	double CONRAT;		// Ratio of (con-rod length)/(crank radius). CONRAT = 2.0*EnginePtr->CONROD/EnginePtr->STROKE;
	double FNN;			// FNN = sqrt(pow(CONRAT, 2) - pow(sin(THETAR), 2));
	double X;			// Piston displacement (m). X = 0.5*EnginePtr->STROKE*(1.0 + CONRAT - FNN - cos(THETAR));
	double FCYL;		// Cylinder cross-sectional area (m^2). FCYL = 0.25*PI*pow(EnginePtr->DCYL, 2);
	double DXDT;		// Change in X with time
	double DTHDT;		// Change in THETA with time (degreesCA.s^-1)
	double DXDTH;		// Change in X with THETA
	double DVCDT;		// Cylinder volume time rate of change
	double DPCDT;		// Cylinder pressure time rate of change

	// Special cases
	// =============
	bool EVO_SET;			// Has PCR and TCR values been reset at EVO this cycle?

	// Combustion
	// ==========
	bool COMB_RESET;		// Have the combustion variables been reset yet this cycle?
	bool INJ_PREV_ITER;		// Injection taking place previous iteration?
	double m_finj;			// Fuel mass injected so far this cycle (kg)
	double m_fub;			// Mass of unburnt fuel in cylinder (kg)
	double m_fb;			// Mass of fuel burned so far this cycle (kg)
	double del_m_fb;		// Mass of fuel burned this iteration (kg)
	double phiFAoverall;	// Overall fuel/air equivalence ratio
	double QLHVrate;		// Combustion heat release rate (J/s)

	// Cycle variables
	// ===============	
	bool WAIT;				// Do not run the cylinder if starting while a valve is open!
	int NCYCLES;			// Cycle counter
	double W_cig;			// Gross indicated work per cycle (J)
	double W_cig_prev;		// Gross indicated work, previous cycle (J)
	double W_cin;			// Net indicated work per cycle (J)
	double W_cin_prev;		// Net indicated work, previous cycle (J)

public:	
	// Valves
	// ======
	CValve* ExhaustValve;	// Exhaust valve array
	CValve* IntakeValve;	// Intake valve array
	double *DMEDT, *DMIDT;	// Mass flow rate arrays, exhaust and intake valves

private:	
	// Motored variables (Woschni heat transfer only)
	// ==============================================
	double PC_mot;			// Motored cylinder pressure (bar)
	double TC_mot;			// Motored cylinder temperature (K)
	double MC_mot;			// Motored mass in cylinder (kg)
	double AC_mot;			// Motored speed of sound in cylinder (m/s)
	double DPC_motDT;		// Motored cylinder pressure time rate of change
	double *DMEDT_mot, *DMIDT_mot;	// Motored mass flow rate arrays, exhaust and intake valves
	bool MOT_SET;			// Have the motored reference variables been set yet this cycle?
	double PC_mot_ref, TC_mot_ref, VC_mot_ref; // Motored reference values (e.g. recorded at IVC)
	double W_cig_mot;		// Motored gross indicated work per cycle (J)
	double W_cig_mot_prev;	// Motored gross indicated work, previous cycle (J)
	double W_cin_mot;		// Motored net indicated work per cycle (J)
	double W_cin_mot_prev;	// Motored net indicated work, previous cycle (J)

	// Measurements
	// ============
	FILE *OUTPUT_FILE;
	FILE *FILE_CYCLE;
	FILE *FILE_OVERALL_MOV;		// For overall measurements of the required properties - movie

	int mov_freq;					// Record data every x timesteps, for movie
	bool PRINT_TO_FILE_MOV;			// Record data for movie (1=true) or not (0=false)
	int print_from_cycle;			// Only start printing data from the start of this cycle onwards

	int num_props_measured;			// Number of properties to measure, not counting time
	bool CYL_PRESSURE;				// Record cylinder pressure (1=true) or not (0=false)
	bool CYL_TEMPERATURE;			// Record cylinder temperature (1=true) or not (0=false)
};

#endif // !defined(AFX_CYLINDER_H__8C8B53E5_114C_4A94_BAB5_3B763BA669B8__INCLUDED_)
