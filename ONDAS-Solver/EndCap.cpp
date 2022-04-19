// EndCap.cpp: implementation of the CEndCap class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "EndCap.h"
#include "Tools.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CEndCap::CEndCap()
{

}

CEndCap::~CEndCap()
{
	delete [] pPipe;
	delete [] pBN;
	delete [] pCLIN;
	delete [] pCLOUT;
	delete [] pend_flow;
}

// Copy constructor
CEndCap::CEndCap(const CEndCap& inEndCap)
{
	ID = inEndCap.ID;
//	EX = inEndCap.EX;
//	AREF = inEndCap.AREF;
//	end = inEndCap.end;
	
	// Yes, copy pointers, this is what is desired
	pPipe = inEndCap.pPipe;
	pBN = inEndCap.pBN;
	pCLIN = inEndCap.pCLIN;
	pCLOUT = inEndCap.pCLOUT;
	pend_flow = inEndCap.pend_flow;
}

CEndCap& CEndCap::operator=(const CEndCap& inEndCap)
{
	if(this != &inEndCap)
	{
		ID = inEndCap.ID;
//		EX = inEndCap.EX;
//		AREF = inEndCap.AREF;
//		end = inEndCap.end;
	
		// Yes, copy pointers, this is what is desired
		pPipe = inEndCap.pPipe;
		pBN = inEndCap.pBN;
		pCLIN = inEndCap.pCLIN;
		pCLOUT = inEndCap.pCLOUT;
		pend_flow = inEndCap.pend_flow;	
	}
	return *this;
}

void CEndCap::Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, /*int** &rCASSYS,*/ int** &rEndCapPIPES, int** &rEndCapPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, int assyid, string parent_assy_res_dir, string calling_object_str)
{
	if (pPpt->SHOW_calls) { pPpt->Out("CEndCap.Initialise ("); pPpt->Out(calling_object_str); pPpt->Out(")\n"); }

	InitialiseGen(pPpt, pPipes, rPipe, rEndCapPIPES, rEndCapPIPES_ENDS, rENDCORR, id, ex, npipes, assyid, "CEndCap", parent_assy_res_dir);
	NAME = new int [NPIPES];
	NAME[ONE_SIDE] = CLOSED;

	std::string res_str;
	if (EX) res_str = "res_ex_closed"; else res_str = "res_in_closed";
	PRINT_MOVIE_FILE = false;
	SetupFiles(pPpt, res_str); // Open results file
}

void CEndCap::RunBoundary()
{
	(*(pCLOUT[ONE_SIDE]))[R+1] = (*(pCLIN[ONE_SIDE]))[R+1];
	*(pend_flow[ONE_SIDE]) = NOFLOW;
}

void CEndCap::PrintToFile(CProperties* pPpt, int timestep, double time, double ca_elapsed, double ca)
// ============================================================ //
// Prints instantaneous boundary data to file.					//
// This function is called from the main function.				//
// ============================================================ //
{
	if (pPpt->SHOW_calls) { pPpt->Out(Identify()); pPpt->Out(".PrintToFile\n"); }

	if (timestep == 0) {

		std::string temp_str = "Object results file for ";
		temp_str += Identify();
		fprintf(OUTPUT_FILE, "%s\n", Underline(StringToChar(temp_str), "-"));

		fprintf(OUTPUT_FILE, "%s", "Time(s)");
		if (!pPpt->CONTINUOUS) fprintf(OUTPUT_FILE, "\t%s%c%s\t%s%c%s", "Elapsed (", Deg(), "CA)", "Crank angle (", Deg(), "CA)");

		// Print internal boundary conditions to file
		fprintf(OUTPUT_FILE, "\t%s\t%s\t%s", "Static pressure (bar)", "Static temperature (K)", "Density (kg/m^3)");

		// Print lambda_in and lambda_out to file
		fprintf(OUTPUT_FILE, "\t%s\t%s", "lambda_in", "lambda_out");

		// Print lambda_in_star and lambda_out_star to file
		fprintf(OUTPUT_FILE, "\t%s\t%s", "lambda_in_star", "lambda_out_star");

		fprintf(OUTPUT_FILE, "\n");
	}

	if (timestep % pPpt->freq == 0) { // Print data at the specified sampling frequency

		fprintf(OUTPUT_FILE, "%f", time);
		if (!pPpt->CONTINUOUS) fprintf(OUTPUT_FILE, "\t%f\t%f", ca_elapsed, ca);
		
		// Print internal boundary conditions to file
		fprintf(OUTPUT_FILE, "\t%f\t%f\t%f", pBN[ONE_SIDE]->p_dash * pPpt->PREF, pBN[ONE_SIDE]->T, pBN[ONE_SIDE]->rho);

		// Print lambda_in and lambda_out to file
		fprintf(OUTPUT_FILE, "\t%f\t%f", (*(pCLIN[0]))[R + 1], (*(pCLOUT[0]))[R + 1]);

		// Print lambda_in_star and lambda_out_star to file
		fprintf(OUTPUT_FILE, "\t%f\t%f", CLIN_STAR[0][R + 1], CLOUT_STAR[0][R + 1]);

		fprintf(OUTPUT_FILE, "\n");
	}
}