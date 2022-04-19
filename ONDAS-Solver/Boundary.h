// Boundary.h: interface for the CBoundary class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BOUNDARY_H__61719847_1887_4688_B688_8AD8EAEB9502__INCLUDED_)
#define AFX_BOUNDARY_H__61719847_1887_4688_B688_8AD8EAEB9502__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Pipe.h"
#include "Properties.h"

class CBoundary  
{
public:
	CBoundary();
	virtual ~CBoundary();

	// Common administrative routines
	// ----------------------------------------------------------------------------------------------------
	void InitialiseGen(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rPIPES, int** &rPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, int assyid, string calling_function_str, string parent_assy_res_dir);
	void Configure(CProperties* pPpt);
	void PrintConnections(CProperties* pPpt);
	void SetupFiles(CProperties* pPpt, std::string res_str);
	void CloseFiles(CProperties* pPpt);
	char* Identify();

	// Common homentropic boundary conditions
	// ----------------------------------------------------------------------------------------------------
	void common_HI_code(CProperties* pPpt, double lambda_in, double &rlambda_out, bool &rCHOKED, double A0, double T);
	double common_HN_code(CProperties* pPpt, int timestep, double time, double lambda_in, double phi, double Pb, bool &rCHOKED, double T);

	// Common non-homentropic boundary conditions
	// ----------------------------------------------------------------------------------------------------
	void common_NHI_code(CProperties* pPpt, double lambda_in_n, double lambda_out_n, double AA_n, double &rlambda_in_c, double &rlambda_out_c, double &rAA_c, double psi, double P0, double T0, bool &rCHOKED, bool &rSONIC, double T, int timestep, double time, bool INFLOW_CONST_P, double tolLoop);
	double common_NHN_code(CProperties* pPpt, double lambda_in, double AA, double phi, double Pb, bool &rCHOKED, double T, double time);
	
	// Common pipe update for homentropic and non-homentropic boundary conditions
	// ----------------------------------------------------------------------------------------------------
	void common_UPDATE_H(CProperties* pPpt, double* lambda_in_c, double* lambda_out_c, int* pipe_flow, bool* CHOKED);
	void common_UPDATE_NH(CProperties* pPpt, double* lambda_in_c, double* lambda_out_c, double* AA_c, int* pipe_flow, bool* CHOKED);

	// Inlines
	// ----------------------------------------------------------------------------------------------------
	inline int Get_NPIPES(){return NPIPES;}
	
protected:

public:
	// Identification
	// ----------------------------------------------------------------------------------------------------
	char* strDesc; 					// Optional object description (max. 500 characters - use underscores for spaces)
	int ID;							// Boundary number
	int* NAME;						// Boundary name(s) e.g. INFLOW, or (TURBINLET, TURBOUTLET)
	bool EX;						// Exhaust or intake boundary
	int AssyID;						// ID of assembly on which boundary belongs
	char* RES_DIR;					// Results directory for this boundary // e.g., ..dat/case/case-name/case-name.001.res/assembly0.res/

	// Configuration
	// ----------------------------------------------------------------------------------------------------
	int NPIPES;						// Number of pipes joined to this boundary
	int* end;						// Vector of ODD or EVEN end of the pipes at the boundary
	CPipe** pPipe;					// Vector of pipe pointers
	CNode** pBN;					// Vector of boundary node pointers
	double*** pCLIN;				// Vector of pointers pointing to CL1 or CL2 depending on the end
	double*** pCLOUT;				// Vector of pointers pointing to CL1 or CL2 depending on the end
	double** CLIN_STAR;				// Vector of starred values of CL1 or CL2 depending on the end
	double** CLOUT_STAR;			// Vector of starred values of CL1 or CL2 depending on the end
	int** pend_flow;				// Pointers to the relevant pipe end flows
	CPathLine** pPathLine;			// Pointers to the relevant pipe end pathlines

	// Reading into these arrays from parameter file
	char** labels;					// Parameter list of labels in no order
	double* values;					// Parameter list of values in same order as labels
	char** strings;					// Parameter list of strings in same order as labels

	// Printing to files, screen
	FILE *OUTPUT_FILE;				// General boundary data output file
	FILE *OUTPUT_FILE_DEBUG;		// Boundary data output file for debugging
	FILE *MOVIE_FILE;				// Boundary data ready for animating
	double print_from_time;			// Only start printing data after this simulation time onwards (s)
	bool PRINT_DEBUG_FILE;			// Record data for debugging (1=true) or not (0=false)
	bool PRINT_MOVIE_FILE;			// Record data for movie (1=true) or not (0=false)
	
	// Working variables
	int* pipe_flow;					// Record the current flow direction at the boundary
	int* pipe_flow_old;				// Record the old flow direction at the boundary

public:
	// Printing to files, screen
	bool PRINT_TO_SCREEN;			// Print data to screen (1=true) or not (0=false)
};

#endif // !defined(AFX_BOUNDARY_H__61719847_1887_4688_B688_8AD8EAEB9502__INCLUDED_)
