// Volume.h: interface for the CVolume class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VOLUME_H__F075B625_EE86_43A1_8F56_7F61274E4725__INCLUDED_)
#define AFX_VOLUME_H__F075B625_EE86_43A1_8F56_7F61274E4725__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//#include "PathLine.h"
#include "Pipe.h"
#include "Properties.h"
#include "Node.h"
#include "Valve.h"

class CVolume  
{
public:
	CVolume();
	virtual ~CVolume();

	void Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rVOLPIPES, int** &rVOLPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, std::string param_dir, std::string res_dir, int assyid, string parent_assy_res_dir, string calling_object_str);
	void Configure(CProperties* pPpt, CPipe** pExhaustSystem, CPipe** pIntakeSystem, int nExPipes, int nInPipes);
	void RunBoundary(CProperties* pPpt, double DELZ, double time, int timestep);
	void ReadInput(CProperties* pPpt, char *InputFile);
	void ListProperties(CProperties* pPpt);
	void PrintToScreen(CProperties* pPpt);
	void PrintToFile(CProperties* pPpt, int timestep, double time, double ca);
	void SetupFiles(CProperties* pPpt, string parent_assy_res_dir);
	void CloseFiles();
	char* Identify();

private:
	int ID;	
	int AssyID;				// ID of assembly on which boundary belongs
	bool EX;
	char* RES_DIR;					// Results directory for this boundary

	// Read from file
	// --------------
	char** labels;		// Parameter list of labels in no order
	double* values;		// Parameter list of values in same order as labels
	char** strings;		// Parameter list of strings in same order as labels

	// Volume geometry
	// --------------------------------------------------------------------------------
	double V_INIT;					// Volume (m^3)
	bool USE_V_INIT;				// Use V_INIT directly (1=true) or derive volume from pipe geometry below (0=false)
	double length;					// The physical length of the system (m)
	bool LINEAR;					// Linear pipe diameter variation (1=true) or quadratic (0=false)
	double C;						// Gradient of pipe diameter linear variation dD/dX (this is not c, C=dD/dx)
	double a, b, c;					// Quadratric pipe diameter variation with X, i.e. dia = aX^2 + bX + c
	double D_odd;					// Diameter at the LH end of the pipe (m)
	double D_even;					// Diameter at the RH end of the pipe (m)
	double D_int;					// Pipe diameter at some interior point (m)
	double x_int;					// Distance from odd end of interior diameter (m)
	double vol;						// Pipe volume (m^3)

	// Initial conditions
	// ------------------
	double p0_INIT;		// Initial volume pressure (bar)
	double T0_INIT; 	// Initial volume temeperature (K)

	// Valve parameters
	// ----------------
	int NVOLUMEVALVES;	// Number of valves/pipes joined at the volume
	bool USE_PIPE_DIA;	// Set valve reference diameter same as adjoining pipe (1=true) (0=false)?
	double ref_dia;		// Valve reference diameter (mm)

	// Working variables
	// -----------------
	double p0;			// Volume pressure (bar)
	double T0;			// Volume temperature (K)
	double V;			// Volume volume (m^3)
	double m;			// Mass in volume (kg)
	double a0;			// Speed of sound in volume (m/s)
	double dVdt;		// Volume volume time rate of change
	double dp0dt;		// Volume pressure time rate of change

	// Temp
	double* pX_temp;
	double* TX_temp;
	double*	regime;
	bool* FORWARD;

	// Measurements
	// ------------
	FILE *OUTPUT_FILE;
	bool USE_DEF_FREQ;				// Use default sampling rate (1=true) or 'freq' below (0=false)
	int freq;						// Print data to file every freq timesteps

public:	
	// Valves
	// ------
	CValve* VolumeValve;// Volume valve array
	double *dmdt;		// Volume valves mass flow rate array
	double* av;			// Speed of sound for mass flow dmdt
};

#endif // !defined(AFX_VOLUME_H__F075B625_EE86_43A1_8F56_7F61274E4725__INCLUDED_)
