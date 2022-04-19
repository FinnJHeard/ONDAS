// EndEnvironment.h: interface for the CEndEnvironment class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ENDENVIRONMENT_H__3AABC5A8_A44D_4B6B_A88C_6C0B07A7CD1E__INCLUDED_)
#define AFX_ENDENVIRONMENT_H__3AABC5A8_A44D_4B6B_A88C_6C0B07A7CD1E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Boundary.h"

class CEndEnvironment : public CBoundary 
{
public:
	CEndEnvironment();
	virtual ~CEndEnvironment();

	void Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rEPIPES, int** &rEPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, std::string param_dir, int assyid, string parent_assy_res_dir, string calling_object_str);
	void ReadInput(CProperties* pPpt, char *InputFile);
	//inline void InitialiseAnechoic(CProperties* pPpt){if(ANECHOIC){p_anechoic = pBN[ONE_SIDE]->p_dash*pPpt->PREF;T_anechoic = pBN[ONE_SIDE]->T;}; 
	void InitialiseAnechoic(CProperties* pPpt);
	void ListProperties(CProperties* pPpt);
	void LoadStagnationFile(CProperties* pPpt, char* InputFile, double** &rArray, int &rDatapoints);

	void RunBoundary(CProperties* pPpt, int timestep, double time);
	void HE(CProperties* pPpt, int timestep, double time);
	void NHE(CProperties* pPpt, int timestep, double time);
	void Anechoic(CProperties* pPpt, int timestep, double time);

	void PrintToScreen(CProperties* pPpt);
	void PrintToFile(CProperties* pPpt, int timestep, double time, double ca_elapsed, double ca, CPipe* Pipe);

	//inline void Set_P0(double P0_temp){P0 = P0_temp;}
	void Set_P0_gradual(double P0_temp);
	inline double Get_P0(){return P0;}
	inline void Set_VAR_P0(bool VAR_P0_temp){VAR_P0 = VAR_P0_temp;}
	inline bool Get_VAR_P0(){return VAR_P0;}
	inline void Set_CALIBRATE(bool CALIBRATE_temp){CALIBRATE = CALIBRATE_temp;}
	inline bool Get_CALIBRATE_READY(){return CALIBRATE_READY;}

private:
	// ====================================================================================================
	// Parameter file for Exhaust End Environment [0]
	// ====================================================================================================
						
	// Basic operation
	// ----------------------------------------------------------------------------------------------------
	//bool ANECHOIC;			// Anechoic (1=true) or ambient termination (0=false)
																								
		// If ambient termination, i.e. ANECHOIC == 0 == false
		// ----------------------------------------------------------------------------------------------------
		double P0;				// Stagnation pressure for constant operation (bar)
		double P0_desired;		// Value of stagnation pressure requested by APLDev calibration routine
		double T0;				// Stagnation temperature for constant operation (K)
		bool STATIC;			// Desired pressure and temperature are static values (1=true) or stagnation (0=false)
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
				double phaseAngle;			// Phase angle (0-360 degrees)
				double PHI;					// Pulse duty cycle
				double P0_LOW;				// Lowest stagnation pressure (bar)
				double P0_HIGH;				// Highest stagnation pressure (bar)
				bool VAR_BY_FILE;			// Apply reservoir conditions from file (1=true) or use sinusoial variables above (0=false)?

					// If variation by file, i.e. VAR_BY_FILE == 1 == true
					// ----------------------------------------------------------------------------------------------------
					char* P0_FILE;				// File containing stagnation pressure variation	
					char* T0_FILE;				// File containing stagnation temperature variation

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

	// Measurements
	// ----------------------------------------------------------------------------------------------------
	bool USE_DEF_FREQ;		// Use default sampling rate (1=true) or 'freq' below (0=false)
	int freq;				// Print data to file every freq timesteps
	bool ROTOR_EXIT_TAPPING;// Use rotor exit tapping (1=true) or not (0=false)
	int p2_pipe;			// Pipe number from which to record p2 when calculating PR
	double p2_loc;			// Location of p2 tapping along p2_pipe as a fraction of pipe length
	bool NOTIFY;			// Notify user (true == 1) of nozzle flow direction change or not (false == 0)
						
	// ====================================================================================================
	// End of file
	// ====================================================================================================

	// Stagnation values from files
	// ----------------------------------------------------------------------------------------------------
	double fileP0, fileT0;	// Stores file pressure and temperature values
	double** p0FromFile;	// Array containing read-in stagnation pressure data
	int p0Datapoints;		// Number of datapoints in read-in stagnation pressure data
	double** T0FromFile;	// Array containing read-in stagnation temperature data
	int T0Datapoints;		// Number of datapoints in read-in stagnation temperature data

//	double p_anechoic;		// Anechoic pressure recorded from boundary node initial pressure
//	double T_anechoic;
//	double lambda_in_an;
//	double lambda_out_an;

	bool SONIC;
	bool CALIBRATE;  // Flag set to true if attached APLDev is undergoing calibration, otherwise false
	bool CALIBRATE_READY; // Flag set to true once the desired P0 is achieved, under APLDev calibration
};

#endif // !defined(AFX_ENDENVIRONMENT_H__3AABC5A8_A44D_4B6B_A88C_6C0B07A7CD1E__INCLUDED_)
