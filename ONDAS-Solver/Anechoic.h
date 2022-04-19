// Anechoic.h: interface for the CAnechoic class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ANECHOIC_H__AFD7243A_7C18_4CCA_B92B_C68BE6542AEC__INCLUDED_)
#define AFX_ANECHOIC_H__AFD7243A_7C18_4CCA_B92B_C68BE6542AEC__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Boundary.h"

class CAnechoic : public CBoundary 
{
public:
	CAnechoic();
	virtual ~CAnechoic();

	void Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rEPIPES, int** &rEPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, std::string param_dir, int assyid, string parent_assy_res_dir, string calling_object_str);
	void ReadInput(CProperties* pPpt, char *InputFile);
	//inline void InitialiseAnechoic(CProperties* pPpt){if(ANECHOIC){p_anechoic = pBN[ONE_SIDE]->p_dash*pPpt->PREF;T_anechoic = pBN[ONE_SIDE]->T;}; 
	//void InitialiseAnechoic(CProperties* pPpt);
	void ListProperties(CProperties* pPpt);
	
	void RunBoundary(CProperties* pPpt, int timestep, double time);
	//void HE(CProperties* pPpt);
	//void NHE(CProperties* pPpt, int timestep, double time);
	void Anechoic(CProperties* pPpt, int timestep, double time);
	void AnechoicOld(CProperties* pPpt, int timestep, double time);
	void InitialiseBufferDamper(CProperties* pPpt);

	void PrintToScreen(CProperties* pPpt);

private:

	// ====================================================================================================
	// Parameter file for Exhaust Anechoic End [0]
	// ====================================================================================================
						
	// Basic operation
	// ----------------------------------------------------------------------------------------------------
	//bool ANECHOIC;			// Anechoic (1=true) or ambient termination (0=false)
											
		// If ambient termination, i.e. ANECHOIC == 0 == false
		// ----------------------------------------------------------------------------------------------------
		double P0;				// Stagnation pressure for constant operation (bar)
		double T0;				// Stagnation temperature for constant operation (K)
		double phi;				// Nozzle area ratio for constant operation
							
		// Variable operation (ignores above phi or P0)
		// ----------------------------------------------------------------------------------------------------
		bool VAR_P0;			// Are reservoir conditions variable (1=true) (0=false)?

			// If ambient conditions variable, i.e. VAR_P0 == 1 == true
			// ----------------------------------------------------------------------------------------------------
			bool P0_COS;			// Stagnation pressure variation is sinusoidal (1=true) or linear/step (0=false)
	
				// If linear/step variation, i.e. P0_COS == 0 == false
				// ----------------------------------------------------------------------------------------------------
				double P0_START;			// Stagnation pressure before step change or at start of variation (bar)
				double P0_END;				// Stagnation pressure for step change or at end of variation (bar)
				double T0_START;			// Stagnation temperature before step change or at start of variation (K)
				double T0_END;				// Stagnation temperature for step change or at end of variation (K)

				// else sinusoidal/file variation (requires VAR_P0==1), i.e. P0_COS == 1 == true
				// ----------------------------------------------------------------------------------------------------
				double FREQ;				// Sinusoidal frequency (Hz)
				double PHI;					// Pulse duty cycle
				double P0_LOW;				// Lowest stagnation pressure (bar)
				double P0_HIGH;				// Highest stagnation pressure (bar)

		bool VAR_PHI;		// Is nozzle area ratio variable (1=true) (0=false)?
	
			// If nozzle area ratio variable, i.e. VAR_PHI == 1 == true
			// ----------------------------------------------------------------------------------------------------
			double PHI_START;			// Nozzle area ratio before step change or at start of variation
			double PHI_END;				// Nozzle area ratio for step change or at end of variation

			// Variable phi or P0 operation (ignores above phi or P0)
			// ----------------------------------------------------------------------------------------------------
			bool SQUARE;				// Is the variation a step change (1=true) or a linear variation (0=false)?
			double START_TIME;			// Time at which to start the step change or variation (s)
			double VAR_TIME;			// Time over which step change exists, or time taken for variation (s)

	// Other
	// ----------------------------------------------------------------------------------------------------
	bool NOTIFY;			// Notify user (true == 1) of nozzle flow direction change or not (false == 0)
	bool PRINT_MOVIE_FILE;	// Record data for movie (1=true) or not (0=false)
						
	// ====================================================================================================
	// End of file
	// ====================================================================================================

public:
	CPipe Buffer;			// Buffer pipe between desired boundary location and stagnation conditions
	CPipe Damper;			// Damper pipe after buffer pipe

private:
	int REAL_PIPE, BUFFER_PIPE; // Labels to the real and buffer pipe boundary nodes

	double p_anechoic;		// Anechoic pressure recorded from boundary node initial pressure
	double T_anechoic;

	double lambda_in_an;
	double lambda_out_an;

	bool SONIC;				// Sonic flow across valve
};

#endif // !defined(AFX_ANECHOIC_H__AFD7243A_7C18_4CCA_B92B_C68BE6542AEC__INCLUDED_)