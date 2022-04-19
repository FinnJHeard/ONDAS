// Assembly.h: interface for the CAssembly class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ASSEMBLY_H__8456099F_6D9D_491C_9AD0_A227DAA186AD__INCLUDED_)
#define AFX_ASSEMBLY_H__8456099F_6D9D_491C_9AD0_A227DAA186AD__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//#include "Properties.h"

#include "StdAfx.h"
#include "Anechoic.h"
#include "APLDev.h"
#include "Boundary.h"
#include "EndCap.h"
#include "Cylinder.h"
#include "EndEnvironment.h"
#include "Engine.h"
#include "Junction.h"
#include "Pipe.h"
#include "Sudden.h"
#include "Transmissive.h"
#include "Turbine.h"
#include "Volume.h"
#include "Time.h"

#include <string>
#include <vector>

using std::string;
using std::vector;

class CAssembly
{
public:
	CAssembly();
	virtual ~CAssembly();

	void Initialise(CProperties* pPpt/*, CAssembly* rAssy*/, int id/*, bool ex, int npipes*//*, std::string param_dir_temp*/);
	void InitialiseBCs(CTime* pMyTime, CProperties* pPpt, CAssembly* &rAssy);
	CPipe** PipePtrs(int NPIPES, CAssembly* &rAssy, int** ASSYS, int** PIPES, int i, bool EX);
	CPipe** PipePtrsCyl(int NPIPES, CAssembly* &rAssy, int** ASSYS, int** PIPES, int i, bool* EX);
	void ConfigureBCs(CProperties* pPpt);
	void ConfigureAll(CProperties* pPpt, vector<string> &rConfStrs);
	void ListPropertiesAll(CProperties* pPpt);
	void ReadAssyDir(CProperties* pPpt, char *InputFile);
	void ReadInput(CProperties* pPpt, char *InputFile);
	void ListProperties(CProperties* pPpt);
	void PrintBoundaryConnections(CProperties* pPpt);
	char* Identify();

	void RunBoundary(CProperties* pPpt, CTime MyTime, double DELZe, double DELZi, double* TIMEe, double* TIMEi, int timestep, bool RESTORE);

	// Declare pipe and boundary objects
	// ----------------------------------------------------------------------------------------------------
	CPipe *Exhaust;
	CPipe *Intake;
	CAnechoic *EX_ANEC;
	CAnechoic *IN_ANEC;
	CEndEnvironment *EX_END;
	CEndEnvironment *IN_END;
	CEndCap *CE;
	CEndCap *CI;
	CSudden *SE;
	CSudden *SI;
	CJunction *JE;
	CJunction *JI;
	CJunction* JE_backup;
	//CJunction* JI_backup;
	CEngine *Eng;
	//CEngine *Eng_backup;
	CCylinder *Cyl;
	//CCylinder *Cyl_backup;
	CAPLDev *APLDevE;
	CAPLDev *APLDevI;
	CTurbine *Turbine;
	CTransmissive *TransmE;
	CTransmissive *TransmI;
//	CAssembly *Assembly;
	CVolume *ExVolume;
	CVolume *InVolume;

	// Identification
	// ----------------------------------------------------------------------------------------------------
	int ID;					// Assembly number

	// Read from file
	// ----------------------------------------------------------------------------------------------------
	char** labels;			// Parameter list of labels in no order
	double* values;			// Parameter list of values in same order as labels
	char** strings;			// Parameter list of strings in same order as labels

	// Directories
	// ----------------------------------------------------------------------------------------------------
	string param_dir;		// Path to parameter files directory
	string res_dir;			// Path to numbered results directory
	string assy_res_dir;	// Path to assembly results directory within numbered results directory
	char* assy_dir;			// Name of directory containing assembly configuration and parameter files
	char* assy_dir_slash;	// Name of directory containing assembly configuration and parameter files plus slash

//	int LEVEL;				// Level of sub-assembly within assembly structure (0 is the top level)

	// Pipe parameters
	// ----------------------------------------------------------------------------------------------------
	int NEXPIPES;				// Number of exhaust pipes
	int NINPIPES;				// Number of intake pipes
	double exVolume;			// Total volume of exhaust system pipes
	double inVolume;			// Total volume of intake system pipes
	
	// Engine, cylinder & valves
	// ----------------------------------------------------------------------------------------------------
	char* ENG_FILE;				// Engine parameter file
	int NENGINES;				// Number of engines
	int NCYLS;					// Number of cylinders
//	char* CYL_DIR;				// Folder containing cylinder P & T data files - must include the "\" following the directory name
	char** CYLCONFP;			// Cylinder configuration strings
	double* CYLENDCORR;			// Cylinder end correction parameters
	int** CYLASSYS;				// Pair of pipe numbers representing the exhaust and intake pipe assemblies
	int** CYLPIPES;				// Pair of pipe numbers representing the exhaust and intake pipe [ID][EXHAUST or INTAKE]
	int** CYLPIPES_ENDS;		// List of ODD/EVEN tags giving the end of each pipe to use
	int* CYLPIPES_NUMS;			// The number of pipes joined at each cylinder

	// End cap parameters
	// ----------------------------------------------------------------------------------------------------
	int NEXENDCAP;				// Number of endcaps in exhaust manifold
	char** EXEndCapCONFP;		// Exhaust manifold end cap configuration strings
	double* EXEndCapENDCORR;	// Exhaust manifold end cap end correction parameters
	int** EXEndCapASSYS;		// Lists the assembly number for the pipe joining each exhaust end cap
	int** EXEndCapPIPES;		// Lists the pipe number joining each exhaust end cap
	int** EXEndCapPIPES_ENDS;	// Lists the ODD/EVEN tag referring to which end of the pipe to use
	int* EXEndCapPIPES_NUMS;	// The number of pipes joined at the end cap
	
	int NINENDCAP;				// Number of endcaps in intake manifold
	char** INEndCapCONFP;		// Intake manifold end cap configuration strings
	double* INEndCapENDCORR;	// Intake manifold end cap end correction parameters
	int** INEndCapASSYS;		// Lists the assembly number for the pipe joining each intake end cap
	int** INEndCapPIPES;		// Lists the pipe number joining each intake end cap
	int** INEndCapPIPES_ENDS;	// Lists the ODD/EVEN tag referring to which end of the pipe to use
	int* INEndCapPIPES_NUMS;	// The number of pipes joined at the end cap

	// End environment parameters
	// ----------------------------------------------------------------------------------------------------
	int NEXEND;					// Number of end environments in exhaust manifold
	char** EXECONFP;			// Exhaust manifold end environment configuration strings
	double* EXEENDCORR;			// Exhaust manifold end environment end end correction parameters
	int** EXEASSYS;				// Lists the assembly number for the pipe joining each exhaust end environment
	int** EXEPIPES;				// Lists the pipe number joining each exhaust end environment
	int** EXEPIPES_ENDS;		// Lists the ODD/EVEN tag referring to which end of the pipe to use
	int* EXEPIPES_NUMS;			// The number of pipes joined at the end environment
	
	int NINEND;					// Number of end environments in intake manifold
	char** INECONFP;			// Intake manifold end environment configuration strings
	double* INEENDCORR;			// Intake manifold end environment end correction parameters
	int** INEASSYS;				// Lists the assembly number for the pipe joining each intake end environment
	int** INEPIPES;				// Lists the pipe number joining each intake end environment
	int** INEPIPES_ENDS;		// Lists the ODD/EVEN tag referring to which end of the pipe to use
	int* INEPIPES_NUMS;			// The number of pipes joined at the end environment

	// Anechoic end parameters
	// ----------------------------------------------------------------------------------------------------
	int NEXANEC;				// Number of anechoic ends in exhaust manifold
	char** EXANECCONFP;			// Exhaust manifold anechoic end configuration strings
	double* EXANECENDCORR;		// Exhaust manifold anechoic end end correction parameters
	int** EXANECASSYS;			// Lists the assembly number for the pipe joining each exhaust anechoic end
	int** EXANECPIPES;			// Lists the pipe number joining each exhaust anechoic end
	int** EXANECPIPES_ENDS;		// Lists the ODD/EVEN tag referring to which end of the pipe to use
	int* EXANECPIPES_NUMS;			// The number of pipes joined at the anechoic end

	int NINANEC;				// Number of anechoic ends in intake manifold
	char** INANECCONFP;			// Intake manifold anechoic end configuration strings
	double* INANECENDCORR;		// Intake manifold anechoic end end correction parameters
	int** INANECASSYS;			// Lists the assembly number for the pipe joining each intake anechoic end
	int** INANECPIPES;			// Lists the pipe number joining each intake anechoic end
	int** INANECPIPES_ENDS;		// Lists the ODD/EVEN tag referring to which end of the pipe to use
	int* INANECPIPES_NUMS;		// The number of pipes joined at the anechoic end

	// Junction parameters
	// ----------------------------------------------------------------------------------------------------
	int NEXJUNCS;				// Number of junctions in exhaust manifold
	char** EXJCONFP;			// Exhaust manifold junction configuration strings
	double* EXJENDCORR;			// Exhaust manifold junction end correction parameters
	int** EXJASSYS;				// List the assembly numbers for the pipes constituting each junction
	int** EXJPIPES;				// List of pipe numbers constituting each junction
	int** EXJPIPES_ENDS;		// List of ODD/EVEN tags giving the end of each pipe to use
	int* EXJPIPES_NUMS;			// The number of pipes joined at the junction

	int NINJUNCS;				// Number of junctions in intake manifold
	char** INJCONFP;			// Intake manifold junction configuration strings
	double* INJENDCORR;			// Intake manifold junction end correction parameters
	int** INJASSYS;				// List the assembly numbers for the pipes constituting each junction
	int** INJPIPES;				// List of pipe numbers constituting each junction
	int** INJPIPES_ENDS;		// List of ODD/EVEN tags giving the end of each pipe to use
	int* INJPIPES_NUMS;			// The number of pipes joined at the junction

	//double ****Loss;			// Matrix of loss coefficients for pulse converter model

	// Adiabatic pressure loss device parameters
	// ----------------------------------------------------------------------------------------------------
	int NEXAPLDev;				// Number of adiabatic pressure loss devices in exhaust manifold
	char** EXAPLDevCONFP;		// Exhaust manifold adiabatic pressure loss device configuration strings
	double* EXAPLDevENDCORR;	// Exhaust manifold adiabatic pressure loss device end correction parameters
	int** EXAPLDevASSYS;		// Exhaust manifold pair of assembly numbers for the pipes being joined at the adiabatic pressure loss device
	int** EXAPLDevPIPES;		// Exhaust manifold pair of pipe numbers being joined at the adiabatic pressure loss device
	int** EXAPLDevPIPES_ENDS;	// List of ODD/EVEN tags giving the end of each pipe to use
	int* EXAPLDevPIPES_NUMS;	// The number of pipes joined at adiabatic pressure loss device

	int NINAPLDev;				// Number of adiabatic pressure loss devices in intake manifold
	char** INAPLDevCONFP;		// Intake manifold adiabatic pressure loss device configuration strings
	double* INAPLDevENDCORR;	// Intake manifold adiabatic pressure loss device end correction parameters
	int** INAPLDevASSYS;		// Intake manifold pair of assembly numbers for the pipes being joined at the adiabatic pressure loss device
	int** INAPLDevPIPES;		// Intake manifold pair of pipe numbers being joined at the adiabatic pressure loss device
	int** INAPLDevPIPES_ENDS;	// List of ODD/EVEN tags giving the end of each pipe to use
	int* INAPLDevPIPES_NUMS;	// The number of pipes joined at adiabatic pressure loss device
//	int * INAPLDevTYPES;		// List of the TYPE of each device; intake manifold

//	char* APLD_DIR;				// Folder containing loss files - must include the "\" following the directory name

	// Sudden area change parameters
	// ----------------------------------------------------------------------------------------------------
	int NEXSUD;					// Number of sudden area changes in exhaust manifold
	char** EXSCONFP;			// Exhaust manifold sudden area change configuration strings
	double* EXSENDCORR;			// Exhaust manifold sudden area change end correction parameters
	int** EXSASSYS;				// Exhaust manifold pair of assembly numbers for the pipes being joined at the sudden area change
	int** EXSPIPES;				// Exhaust manifold pair of pipe numbers being joined at the sudden area change
	int** EXSPIPES_ENDS;		// List of ODD/EVEN tags giving the end of each pipe to use
	int* EXSPIPES_NUMS;			// The number of pipes joined at the sudden area change

	int NINSUD;					// Number of sudden area changes in intake manifold
	char** INSCONFP;			// Intake manifold sudden area change configuration strings
	double* INSENDCORR;			// Intake manifold sudden area change end correction parameters
	int** INSASSYS;				// Intake manifold pair of assembly numbers for the pipes being joined at the sudden area change
	int** INSPIPES;				// Intake manifold pair of pipe numbers being joined at the sudden area change
	int** INSPIPES_ENDS;		// List of ODD/EVEN tags giving the end of each pipe to use
	int* INSPIPES_NUMS;			// The number of pipes joined at the sudden area change

	// Turbine parameters
	// ----------------------------------------------------------------------------------------------------
	int NTURBINE;				// Number of turbines
	char** TURBCONFP;			// Turbine configuration strings
	double* TURBENDCORR;		// Turbine end correction parameters
	int** TURBASSYS;			// Turbine pair of pair of assembly numbers for the pipes being joined at the turbine
	int** TURBPIPES;			// Turbine pair of pair of pipe numbers being joined at the turbine
	int** TURBPIPES_ENDS;		// List of ODD/EVEN tags giving the end of each pipe to use
	int* TURBPIPES_NUMS;		// The number of pipes joined at the turbine

//	char* TURB_DIR;				// Folder containing turbines maps - must include the "\" following the directory name

	// Transmissive boundary parameters
	// ----------------------------------------------------------------------------------------------------
	int NEXTRANSM;				// Number of transmissive boundaries in exhaust manifold
	char** EXTCONFP;			// Exhaust manifold transmissive boundary configuration strings
	double* EXTENDCORR;			// Exhaust manifold transmissive boundary end correction parameters
	int** EXTASSYS;				// Lists the assembly numbers for the pipes joining each exhaust transmissive boundary
	int** EXTPIPES;				// Lists the pipe number joining each exhaust transmissive boundary
	int** EXTPIPES_ENDS;		// Lists the ODD/EVEN tag referring to which end of the pipe to use
	int* EXTPIPES_NUMS;			// The number of pipes joined at the transmissive boundary

	int NINTRANSM;				// Number of transmissive boundaries in intake manifold
	char** INTCONFP;			// Intake manifold transmissive boundary configuration strings
	double* INTENDCORR;			// Intake manifold transmissive boundary end correction parameters
	int** INTASSYS;				// Lists the assembly numbers for the pipes joining each intake transmissive boundary
	int** INTPIPES;				// Lists the pipe number joining each intake transmissive boundary
	int** INTPIPES_ENDS;		// Lists the ODD/EVEN tag referring to which end of the pipe to use
	int* INTPIPES_NUMS;			// The number of pipes joined at the transmissive boundary
	
//	char* TRANSM_DIR;			// Folder containing input wave files - must include the "\" following the directory name

	// Assembly boundary parameters
	// ----------------------------------------------------------------------------------------------------
//	int NASSEMBLY;				// Number of sub-assemblies

	// Volume boundary parameters
	// ----------------------------------------------------------------------------------------------------
	int NEXVOLUMES;				// Number of volumes in exhaust manifold
	char** EXVOLCONFP;			// Exhaust manifold volume configuration strings
	double* EXVOLENDCORR;		// Exhaust manifold volume end correction parameters
	int** EXVOLASSYS;			// List the assembly numbers for the pipes joining each volume
	int** EXVOLPIPES;			// List of pipe numbers joining each volume
	int** EXVOLPIPES_ENDS;		// List of ODD/EVEN tags giving the end of each pipe to use
	int* EXVOLPIPES_NUMS;		// The number of pipes joined to the volume

	int NINVOLUMES;				// Number of volumes in intake manifold
	char** INVOLCONFP;			// Intake manifold volume configuration strings
	double* INVOLENDCORR;		// Intake manifold volume end correction parameters
	int** INVOLASSYS;			// List the assembly numbers for the pipes joining each volume
	int** INVOLPIPES;			// List of pipe numbers joining each volume
	int** INVOLPIPES_ENDS;		// List of ODD/EVEN tags giving the end of each pipe to use
	int* INVOLPIPES_NUMS;		// The number of pipes joined to the volume
};

#endif // !defined(AFX_ASSEMBLY_H__8456099F_6D9D_491C_9AD0_A227DAA186AD__INCLUDED_)
