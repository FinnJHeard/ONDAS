// EndEnvironment.cpp: implementation of the CEndEnvironment class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "EndEnvironment.h"
#include "Tools.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CEndEnvironment::CEndEnvironment()
{

}

CEndEnvironment::~CEndEnvironment()
{

}

void CEndEnvironment::Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rEPIPES, int** &rEPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, std::string param_dir, int assyid, string parent_assy_res_dir, string calling_object_str)
{
	if (pPpt->SHOW_calls) { pPpt->Out("CEndEnvironment.Initialise ("); pPpt->Out(calling_object_str); pPpt->Out(")\n"); }

	InitialiseGen(pPpt, pPipes, rPipe, rEPIPES, rEPIPES_ENDS, rENDCORR, id, ex, npipes, assyid, calling_object_str, parent_assy_res_dir);

	std::string bcname_str = "END";
	ReadInput(pPpt, ConstructString(pPpt, param_dir, bcname_str, EX, ID));
	
	// Boundary name
	NAME = new int [NPIPES];
	NAME[ONE_SIDE] = END_ENVIRONMENT;
	
	std::string res_str;
	if (EX) res_str = "res_ex_endenv"; else res_str = "res_in_endenv";
	SetupFiles(pPpt, res_str); // Open results file
	
	SONIC = false;

	P0_desired = P0; // This won't change unless object is part of a APLDev calibration routine

	if(VAR_P0 && P0_COS	&& VAR_BY_FILE) // If stagnation conditions vary according to file					
	{
		// Load stagnation pressure values from file
		LoadStagnationFile(pPpt, ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->ENDENV_DIR), P0_FILE), p0FromFile, p0Datapoints);
		LoadStagnationFile(pPpt, ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->ENDENV_DIR), T0_FILE), T0FromFile, T0Datapoints);
	}

	fileP0 = P0; 
	fileT0 = T0;

  CALIBRATE = false;  // Flag set to true if attached APLDev is undergoing calibration, otherwise false
  CALIBRATE_READY = false;  // Flag set to true once the desired P0 is achieved, under APLDev calibration
 }

void CEndEnvironment::InitialiseAnechoic(CProperties* pPpt)
{
//	if(ANECHOIC)
	{
//		p_anechoic = pBN[ONE_SIDE]->p_dash*pPpt->PREF;
//		T_anechoic = pBN[ONE_SIDE]->T;

		// Zero velocity, so lambda_in = lambda_out = A
//		lambda_in_an = sqrt(pPpt->gammaAir()*pPpt->R_air*T_anechoic)/pPipe[ONE_SIDE]->AREF;
//		lambda_out_an = sqrt(pPpt->gammaAir()*pPpt->R_air*T_anechoic)/pPipe[ONE_SIDE]->AREF;
	}
	
} 
void CEndEnvironment::ReadInput(CProperties* pPpt, char *InputFile)
{
	int last_entry;
	CommonReadInput(InputFile, pPpt->NUM_PARAMS, this->labels, this->values, this->strings, last_entry);
	for(int r=0; r<last_entry+1; ++r)
	{
		// ====================================================================================================
		// Parameter file for Exhaust End Environment [0]
		// ====================================================================================================
		
		// Optional object description (max. 500 characters - use underscores for spaces)
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "strDesc") == 0) strDesc = strings[r];

		// Basic operation
		// ----------------------------------------------------------------------------------------------------
		//if(strcmp(labels[r], "ANECHOIC") == 0) ANECHOIC = bool(values[r]);

			// If ambient termination, i.e. ANECHOIC == 0 == false
			// ----------------------------------------------------------------------------------------------------
			if(strcmp(labels[r], "P0") == 0) P0 = values[r];
			if(strcmp(labels[r], "T0") == 0) T0 = values[r];
			if(strcmp(labels[r], "STATIC") == 0) STATIC = DoubleToBool(values[r]);
			if(strcmp(labels[r], "phi") == 0) phi = values[r];

			// Variable operation (ignores above phi or P0)
			// ----------------------------------------------------------------------------------------------------
			if(strcmp(labels[r], "VAR_P0") == 0) VAR_P0 = DoubleToBool(values[r]);

				// If ambient conditions variable, i.e. VAR_P0 == 1 == true
				// ----------------------------------------------------------------------------------------------------
				if(strcmp(labels[r], "P0_COS") == 0) P0_COS = DoubleToBool(values[r]);

				// If linear/step variation, i.e. P0_COS == 0 == false
				// ----------------------------------------------------------------------------------------------------
				if(strcmp(labels[r], "P0_START") == 0) P0_START = values[r];
				if(strcmp(labels[r], "P0_END") == 0) P0_END = values[r];
				if(strcmp(labels[r], "T0_START") == 0) T0_START = values[r];
				if(strcmp(labels[r], "T0_END") == 0) T0_END = values[r];

				// else sinusoidal/file variation (requires VAR_P0==1), i.e. P0_COS == 1 == true
				// ----------------------------------------------------------------------------------------------------
				if(strcmp(labels[r], "FREQ") == 0) FREQ = values[r];
				if(strcmp(labels[r], "phaseAngle") == 0) phaseAngle = values[r];
				if(strcmp(labels[r], "PHI") == 0) PHI = values[r];
				if(strcmp(labels[r], "P0_LOW") == 0) P0_LOW = values[r];
				if(strcmp(labels[r], "P0_HIGH") == 0) P0_HIGH = values[r];
				if(strcmp(labels[r], "VAR_BY_FILE") == 0) VAR_BY_FILE = DoubleToBool(values[r]);

					// If variation by file, i.e. VAR_BY_FILE == 1 == true
					// ----------------------------------------------------------------------------------------------------
					if(strcmp(labels[r], "P0_FILE") == 0) P0_FILE = strings[r];
					if(strcmp(labels[r], "T0_FILE") == 0) T0_FILE = strings[r];

			if(strcmp(labels[r], "VAR_PHI") == 0) VAR_PHI = DoubleToBool(values[r]);

				// If nozzle area ratio variable, i.e. VAR_PHI == 1 == true
				// ----------------------------------------------------------------------------------------------------
				if(strcmp(labels[r], "PHI_START") == 0) PHI_START = values[r];
				if(strcmp(labels[r], "PHI_END") == 0) PHI_END = values[r];

				// Variable phi or P0 operation (ignores above phi or P0)
				// ----------------------------------------------------------------------------------------------------
				if(strcmp(labels[r], "SQUARE") == 0) SQUARE = DoubleToBool(values[r]);
				if(strcmp(labels[r], "START_TIME") == 0) START_TIME = values[r];
				if(strcmp(labels[r], "VAR_TIME") == 0) VAR_TIME = values[r];

		// Measurements
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "USE_DEF_FREQ") == 0) USE_DEF_FREQ = DoubleToBool(values[r]);
		if(strcmp(labels[r], "freq") == 0) freq = int(values[r]);
		if(strcmp(labels[r], "ROTOR_EXIT_TAPPING") == 0) ROTOR_EXIT_TAPPING = DoubleToBool(values[r]);
		if(strcmp(labels[r], "p2_pipe") == 0) p2_pipe = int(values[r]);
		if(strcmp(labels[r], "p2_loc") == 0) p2_loc = values[r];
		if(strcmp(labels[r], "print_from_time") == 0) print_from_time = values[r];
		if(strcmp(labels[r], "PRINT_DEBUG_FILE") == 0) PRINT_DEBUG_FILE = DoubleToBool(values[r]);
		if(strcmp(labels[r], "PRINT_MOVIE_FILE") == 0) PRINT_MOVIE_FILE = DoubleToBool(values[r]);
		if(strcmp(labels[r], "NOTIFY") == 0) NOTIFY = DoubleToBool(values[r]);

		// ====================================================================================================
		// End of file
		// ====================================================================================================
	}
	// Set some derived parameters
	if(USE_DEF_FREQ) freq = pPpt->freq;		// Use the default sampling rate
}
void CEndEnvironment::ListProperties(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ListProperties\n");}

	// ====================================================================================================
	// Parameter file for Exhaust End Environment []
	// ====================================================================================================

	pPpt->Out(Underline(Identify(), "=", "\t", strDesc));
	pPpt->Out("\n");

	// Basic operation
	// ----------------------------------------------------------------------------------------------------
//	if(ANECHOIC)
//	{
//		pPpt->Out("\tAnechoic termination, ANECHOIC\t\t\t=\t"); pPpt->Out(TrueOrFalse(ANECHOIC)); pPpt->Out("\n");
//		pPpt->Out("\tAnechoic ref. pressure, p_anechoic\t\t=\t"); pPpt->Out(p_anechoic); pPpt->Out(" bar\n");
//		pPpt->Out("\tAnechoic ref. temperature, T_anechoic\t\t=\t"); pPpt->Out(T_anechoic); pPpt->Out(" K\n");
//	}
//	else
	{
		if(STATIC){pPpt->Out("\tDesired p & T are static values, STATIC\t\t=\t"); pPpt->Out(TrueOrFalse(STATIC)); pPpt->Out("\n");}
		else{pPpt->Out("\tDesired p & T are stagnation values, STATIC\t=\t"); pPpt->Out(TrueOrFalse(STATIC)); pPpt->Out("\n");}
	
		if(phi!=0)
		{
		if(VAR_P0)
		{
			pPpt->Out("\tVariable stagnation conditions");
			if(P0_COS)
			{
				if(VAR_BY_FILE)
				{
					pPpt->Out(" - variation by file: \n");
					pPpt->Out("\tStagnation pressure file, P0_FILE\t\t=\t"); pPpt->Out(P0_FILE); pPpt->Out("\n");
					pPpt->Out("\tStagnation temperature file, T0_FILE\t\t=\t"); pPpt->Out(T0_FILE); pPpt->Out("\n");
					pPpt->Out("\tFrequency, FREQ\t\t\t\t\t=\t"); pPpt->Out(FREQ); pPpt->Out(" Hz\n");
				}
				else
				{
					pPpt->Out(" - sinusoidal variation:\n");
					pPpt->Out("\tSinusoidal frequency, FREQ\t\t\t=\t"); pPpt->Out(FREQ); pPpt->Out(" Hz\n");
					pPpt->Out("\tPhase angle, phaseAngle\t\t\t\t=\t"); pPpt->Out(phaseAngle); pPpt->Out(" degrees\n");
					pPpt->Out("\tPulse duty cycle, PHI\t\t\t\t=\t"); pPpt->Out(PHI); pPpt->Out("\n");
					pPpt->Out("\tLowest stagnation pressure, P0_LOW\t\t=\t"); pPpt->Out(P0_LOW); pPpt->Out(" bar\n");
					pPpt->Out("\tHighest stagnation pressure, P0_HIGH\t\t=\t"); pPpt->Out(P0_HIGH); pPpt->Out(" bar\n");
				}
			}
			else
			{
				if(!SQUARE)
				{
					pPpt->Out(" - linear ramp:\n");
					pPpt->Out("\tStart stagnation pressure, P0_START\t\t=\t"); pPpt->Out(P0_START); pPpt->Out(" bar\n");
					pPpt->Out("\tEnd stagnation pressure, P0_END\t\t=\t"); pPpt->Out(P0_END); pPpt->Out(" bar\n");
					pPpt->Out("\tStart stagnation temperature, T0_START\t\t=\t"); pPpt->Out(T0_START); pPpt->Out(" K\n");
					pPpt->Out("\tEnd stagnation temperature, T0_END\t\t=\t"); pPpt->Out(T0_END); pPpt->Out(" K\n");
					pPpt->Out("\tVariation starts at, START_TIME\t\t=\t"); pPpt->Out(START_TIME); pPpt->Out(" s\n");
					pPpt->Out("\tTime taken, VAR_TIME\t\t\t\t=\t"); pPpt->Out(VAR_TIME); pPpt->Out(" s\n");
				}
				else
				{
					pPpt->Out(" - step change:\n");
					pPpt->Out("\tStart stagnation pressure, P0_START\t\t=\t"); pPpt->Out(P0_START); pPpt->Out(" bar\n");
					pPpt->Out("\tStep stagnation pressure, P0_END\t\t=\t"); pPpt->Out(P0_END); pPpt->Out(" bar\n");
					pPpt->Out("\tStart stagnation temperature, T0_START\t\t=\t"); pPpt->Out(T0_START); pPpt->Out(" K\n");
					pPpt->Out("\tStep stagnation temperature, T0_END\t\t=\t"); pPpt->Out(T0_END); pPpt->Out(" K\n");
					pPpt->Out("\tStep change starts at, START_TIME\t\t=\t"); pPpt->Out(START_TIME); pPpt->Out(" s\n");
					pPpt->Out("\tTime over which step change exists, VAR_TIME\t=\t"); pPpt->Out(VAR_TIME); pPpt->Out(" s\n");
				}
			}
		}
		else
		{
			pPpt->Out("\tConstant stagnation/back pressure, P0\t\t=\t"); pPpt->Out(P0); pPpt->Out(" bar\n");
			pPpt->Out("\tConstant stagnation/back temperature, T0\t=\t"); pPpt->Out(T0); pPpt->Out(" K\n");
		}
		}
		if(VAR_PHI)
		{
			pPpt->Out("\n");
			pPpt->Out("\tNozzle variable operation");
			if(!SQUARE)
			{
				pPpt->Out(" - linear ramp:\n");
				pPpt->Out("\tStart area ratio, PHI_START\t\t\t=\t"); pPpt->Out(PHI_START); pPpt->Out("\n");
				pPpt->Out("\tEnd area ratio, PHI_END\t\t\t\t=\t"); pPpt->Out(PHI_END); pPpt->Out("\n");
				pPpt->Out("\tVariation starts at, START_TIME\t\t\t=\t"); pPpt->Out(START_TIME); pPpt->Out(" s\n");
				pPpt->Out("\tTime taken, VAR_TIME\t\t\t\t=\t"); pPpt->Out(VAR_TIME); pPpt->Out(" s\n");
			}
			else
			{
				pPpt->Out(" - step change:\n");
				pPpt->Out("\tStart area ratio, PHI_START\t\t\t=\t"); pPpt->Out(PHI_START); pPpt->Out("\n");
				pPpt->Out("\tStep area ratio, PHI_END\t\t\t\t=\t"); pPpt->Out(PHI_END); pPpt->Out("\n");
				pPpt->Out("\tStep change starts at, START_TIME\t\t\t=\t"); pPpt->Out(START_TIME); pPpt->Out(" s\n");
				pPpt->Out("\tTime over which step change exists, VAR_TIME\t=\t"); pPpt->Out(VAR_TIME); pPpt->Out(" s\n");
			}
		}
		else
		{
			if(phi==1)
			{
				pPpt->Out("\tOpen pipe end, phi\t\t\t\t=\t"); pPpt->Out(phi); pPpt->Out("\n");
			}
			else
			{
				pPpt->Out("\tConstant nozzle area ratio, phi\t\t\t=\t"); pPpt->Out(phi); pPpt->Out("\n");
			}
		}
	}

	// Measurements
	// ----------------------------------------------------------------------------------------------------
	if(USE_DEF_FREQ)
	{
		if(freq==1) pPpt->Out("\tUsing default sampling rate, pPpt->freq\t\t=\tonce per timestep\n");
		else{pPpt->Out("\tUsing default sampling rate, pPpt->freq)\t\t=\tonce per "); pPpt->Out(freq); pPpt->Out(" timesteps\n");}
	}
	else
	{
		if(freq==1) pPpt->Out("\tUsing local sampling rate, freq\t\t\t=\tonce per timestep\n");
		else{pPpt->Out("\tUsing local sampling rate, freq\t\t\t=\tonce per "); pPpt->Out(freq); pPpt->Out(" timesteps\n");}
	}

	if(ROTOR_EXIT_TAPPING)
	{
		pPpt->Out("\tRotor exit tapping located in pipe, p2_pipe\t=\t");
		if(EX) pPpt->Out("Exhaust["); else pPpt->Out("Inatke[");
		pPpt->Out(p2_pipe); pPpt->Out("]\n");
		pPpt->Out("\t- at a length fraction of, p2_loc\t\t=\t"); pPpt->Out(p2_loc); pPpt->Out("\n");
	}
	pPpt->Out("\tRecord onwards from, print_from_time\t\t=\t"); pPpt->Out(print_from_time); pPpt->Out(" s\n");	
	if(PRINT_DEBUG_FILE) pPpt->Out("\tPrinting debug file\n"); else pPpt->Out("\tNot printing debug file\n");
	if(PRINT_MOVIE_FILE) pPpt->Out("\tPrinting movie file\n"); else pPpt->Out("\tNot printing movie file\n");
	//if(NOTIFY) pPpt->Out("\tNozzle flow direction changes will be notified\n");
	//else pPpt->Out("\tNozzle flow direction changes will not be notified\n");
	pPpt->Out("\n");
	pPpt->Out("\n");
						
// ====================================================================================================
// End of file
// ====================================================================================================
}
void CEndEnvironment::LoadStagnationFile(CProperties* pPpt, char* InputFile, double** &rArray, int &rDatapoints)
{
	double tempInterval, tempValue;
	int row, col;
	int numColumns = 2;
	
	FILE* stream;
	stream = fopen(InputFile, "r");

	if(stream == NULL)
	{
		pPpt->Out(Identify());
		pPpt->Out(":CEndEnvironment::LoadStagnationFile: Error while trying to open the stagnation values file "); pPpt->Out(InputFile); pPpt->Out("\n");
		exit(1);
	}
	else
	{
		rDatapoints = 0;
		fseek(stream, 0L, SEEK_SET);				// Set pointer to beginning of file
		do
		{	
			fscanf(stream, "%lf", &tempInterval);	// Runs over interval
			//cout << tempInterval << "\t";
			fscanf(stream, "%lf", &tempValue);		// Runs over value
			//cout << tempValue << "\n";
			++rDatapoints;
		}while(fscanf(stream, "%l")!=EOF);
		
		//cout << "dataPoints = " << dataPoints << "\n";
		rArray = new double* [rDatapoints];
		for(row=0; row<rDatapoints; ++row) rArray[row] = new double [numColumns];
		
		fseek(stream, 0L, SEEK_SET);				// Reset pointer to beginning of file
			
		for(row=0; row<rDatapoints; ++row)
		{
			for(col=0; col<numColumns; ++col)
			{
				fscanf(stream, "%lf", &rArray[row][col]);	// Record value into array
//cout << rArray[row][col] << "\t";
			}
//cout << "\n";
		}
	}
	fclose(stream);
}

void CEndEnvironment::RunBoundary(CProperties* pPpt, int timestep, double time)
// ====================================================================================================
// Controls execution of the end environment boundary method
// ====================================================================================================
{
	// Update transient conditions
	// ----------------------------------------------------------------------------------------------------
	if(CALIBRATE) // End environment is being controlled by APLDev calibration procedure
	{
		//if(P0!=P0_desired) // Adjust end environment stagnation pressure to the desired (passed) setting gradually, rather than a sudden step change
		if(fabs(P0 - P0_desired) > pPpt->ZERO_TOL) // Adjust end environment stagnation pressure to the desired (passed) setting gradually, rather than a sudden step change
		{
			CALIBRATE_READY = false; // Don't attempt APLDev adjustment until the desired P0 is reached
			//cout << "P0 before = " << P0 << endl;
			//cout << "P0_desired = " << P0_desired << endl;
			//P0 += (P0_desired - P0)/100;
			P0 += (P0_desired - P0)*pPpt->tol_steady*pPpt->tol_steady_multiplier;
			//cout << "P0 after = " << P0 << endl << endl;
		}
		else
		{
		CALIBRATE_READY = true;
		//cout << "CALIBRATE_READY = true" << endl;
		}
	}
	else // Normal operation
	{
		if(VAR_P0)
		{ // Calculate variable stagnation conditions according to input parameters
			if(P0_COS) // Sinusoidal/file variation
			{
				double time_in_period = fmod(time, (1/FREQ));
				//cout << "time_in_period = " << time_in_period << endl;
				//if(time_in_period_prev>time_in_period) // e.g. at the start of a new pulse
				//	time_in_period_prev = time_in_period_prev - period; // Will give a -ve but that's ok
				double integer_part;
				modf(time/(1/FREQ), &integer_part);
  				
				if(VAR_BY_FILE) {
					int i;
					double xValue;
					double fraction_through_pulse = time_in_period/(1/FREQ);

					  //double exact_no_of_pulses =(time_in_period/(1/FREQ)) + (fmod(phaseAngle, 360)/360); // Adjust by phase angle first
					  //double integer_part_pulse;
					  //modf(exact_no_of_pulses, &integer_part_pulse);
					  //double fraction_through_pulse = exact_no_of_pulses - integer_part_pulse; // Decimal part, i.e., fraction through pulse following phase adjustment
					  //pPpt->Out("\n");
					  //pPpt->Out("\n");
					  //pPpt->Out("FREQ = "); pPpt->Out(FREQ); pPpt->Out("\n");
					  //pPpt->Out("period = "); pPpt->Out(1/FREQ); pPpt->Out("\n");
					  //pPpt->Out("time_in_period = "); pPpt->Out(time_in_period); pPpt->Out("\n");
					  //pPpt->Out("phaseAngle = "); pPpt->Out(phaseAngle); pPpt->Out("\n");
					  //pPpt->Out("exact_no_of_pulses = "); pPpt->Out(exact_no_of_pulses); pPpt->Out("\n");
					  //pPpt->Out("integer_part_pulse = "); pPpt->Out(integer_part_pulse); pPpt->Out("\n");
					  //pPpt->Out("fraction_through_pulse = "); pPpt->Out(fraction_through_pulse); pPpt->Out("\n");
					  //pPpt->Out("\n");
					  //pPpt->Out("\n");

					  // Interpolate stagnation values from file
					  xValue = fraction_through_pulse*p0FromFile[p0Datapoints-1][0];
					  i = 1;
					  while(p0FromFile[i][0] < xValue && i < p0Datapoints - 1) ++i;	// xValue lies between [i] and [i-1], or lies outside range
  					
					  P0 = p0FromFile[i-1][1] + 
							  (((p0FromFile[i][1] - p0FromFile[i-1][1])/(p0FromFile[i][0] - p0FromFile[i-1][0])) // Gradient
							  *(xValue - p0FromFile[i-1][0]));

					  xValue = fraction_through_pulse*T0FromFile[T0Datapoints-1][0];
					  i = 1;
					  while(T0FromFile[i][0] < xValue && i < T0Datapoints - 1) ++i;		// xValue lies between [i] and [i-1], or lies outside range
					  T0 = T0FromFile[i-1][1] + 
							  (((T0FromFile[i][1] - T0FromFile[i-1][1])/(T0FromFile[i][0] - T0FromFile[i-1][0])) // Gradient
							  *(xValue - T0FromFile[i-1][0]));
        }
		else
        {
  /*
					  // Test for start of a new pulse (does not include the start of the new pulse at the start of the simulation)
					  if(whole_number_of_pulses != int(integer_part))
					  {
						  START_OF_NEW_PULSE = true;
					  }
  */
					  double exact_no_of_pulses =(time_in_period/(1/FREQ)) + (fmod(phaseAngle, 360)/360); // Adjust by phase angle first
					  double integer_part_pulse;
					  modf(exact_no_of_pulses, &integer_part_pulse);
					  double fraction_through_pulse = exact_no_of_pulses - integer_part_pulse; // Decimal part, i.e., fraction through pulse following phasing
					  time_in_period = fraction_through_pulse*(1/FREQ); // Recalculate time_in_period

					  P0 = P0_LOW + ( ((-1*cos(2*PI*time_in_period/(PHI*(1/FREQ))) + 1)/2) * (P0_HIGH - P0_LOW));	// Cosine waveform
					  T0 = fileT0;
  					
					  //if(exact_no_of_pulses > 1)
					  //{
					  //pPpt->Out("\n");
					  //pPpt->Out("\n");
					  //pPpt->Out("FREQ = "); pPpt->Out(FREQ); pPpt->Out("\n");
					  //pPpt->Out("period = "); pPpt->Out(1/FREQ); pPpt->Out("\n");
					  //pPpt->Out("time_in_period = "); pPpt->Out(time_in_period); pPpt->Out("\n");
					  //pPpt->Out("phaseAngle = "); pPpt->Out(phaseAngle); pPpt->Out("\n");
					  //pPpt->Out("exact_no_of_pulses = "); pPpt->Out(exact_no_of_pulses); pPpt->Out("\n");
					  //pPpt->Out("integer_part_pulse = "); pPpt->Out(integer_part_pulse); pPpt->Out("\n");
					  //pPpt->Out("fraction_through_pulse = "); pPpt->Out(fraction_through_pulse); pPpt->Out("\n");
					  //pPpt->Out("P0 = "); pPpt->Out(P0); pPpt->Out("\n");
					  //pPpt->Out("\n");
					  //pPpt->Out("\n");
					  //}
          }
      }
      else // Step or linear P0 variation
      {	
        if(!SQUARE) // Linear variation
        {
          if(time < START_TIME)
          {
						  P0 = P0_START;
						  T0 = T0_START;
				  }
          else
          {
            if(time < START_TIME + VAR_TIME)
            {
							  P0 = P0_START + (P0_END - P0_START)*((time - START_TIME)/VAR_TIME);
							  T0 = T0_START + (T0_END - T0_START)*((time - START_TIME)/VAR_TIME);
            }
					  else
            {
							  P0 = P0_END;
							  T0 = T0_END;
					  }
				  }
			  }
        else // Step variation
        {
          if(time < START_TIME)
          {
						  P0 = P0_START;
						  T0 = T0_START;
				  }
				  else
          {
            if(time < START_TIME + VAR_TIME)
            {
							  P0 = P0_START;
							  T0 = T0_START;
					  }
					  else
            {
                P0 = P0_END;
							  T0 = T0_END;
					  }
				  }
        }
      }
    }
    else
    {
		  P0 = fileP0;
		  T0 = fileT0;
	  }

    // Calculate variable phi according to input parameters
    if(VAR_PHI)
    {    
	    if(!SQUARE) // Linear variation
      {
		    if(time < START_TIME) phi = PHI_START;
			  else
        {
					  if(time < START_TIME + VAR_TIME) phi = PHI_START + (PHI_END - PHI_START)*((time - START_TIME)/VAR_TIME);
					  else phi = PHI_END;
			  }
		  }
		  else // Step variation
      {
		    if(time < START_TIME) phi = PHI_START;
        else
        {
			    if(time < START_TIME + VAR_TIME) phi = PHI_END;
				  else phi = PHI_START;
			  }
      }
    }
    
	  // If provided desired waveforms are in static form rather than stagnation, convert to stagnation
	  if(STATIC)
	  {
  //cout << "Des p = " << P0 << endl;
  //cout << "Des T = " << T0 << endl;
		  P0 = TotalPressureBar(pPpt, P0, pBN[ONE_SIDE]->T, pBN[ONE_SIDE]->U*pPipe[ONE_SIDE]->AREF); 
		  T0 = TotalTemperature(pPpt, pBN[ONE_SIDE]->T, pBN[ONE_SIDE]->U*pPipe[ONE_SIDE]->AREF);
    }
  }
  //cout << "P0 final = " << P0 << endl << endl;

	// Call appropriate boundary method
	// ----------------------------------------------------------------------------------------------------
	pPpt->HOMENTROPIC ? HE(pPpt, timestep, time) : NHE(pPpt, timestep, time);
}

void CEndEnvironment::HE(CProperties* pPpt, int timestep, double time)
// ====================================================================================================
// End environment for flow to/from stagnation conditions; homentropic flow
// ====================================================================================================
{
	double lambda_in_n, lambda_out_n;
	double *lambda_in_c, *lambda_out_c;
	int *pipe_flow;
	lambda_in_c = new double [NPIPES]; lambda_out_c = new double [NPIPES]; pipe_flow = new int [NPIPES];
	bool* CHOKED; CHOKED = new bool [NPIPES]; 
	CHOKED[ONE_SIDE] = false; // Reset here for all cases

	lambda_in_n = (*(pCLIN[ONE_SIDE]))[R+1];
	lambda_out_n = (*(pCLOUT[ONE_SIDE]))[R+1];
	lambda_in_c[ONE_SIDE] = lambda_in_n;
	lambda_out_c[ONE_SIDE] = lambda_out_n;
	
	double A0 = pow(P0/pPpt->PREF, pPpt->Q); // Generalised flow direction test - identical to PIp in n-h method (Q = (gammaAir()-1.0)/(2.0*gammaAir()))
//if(!from_cylinder) rpipe_flow_old = rpipe_flow; // Store old value

	if(fabs(lambda_in_n - A0) < pPpt->ZERO_TOL || phi==0) // NOFLOW
	{
		pipe_flow[ONE_SIDE] = NOFLOW;
		lambda_in_c[ONE_SIDE] = lambda_in_n;
		lambda_out_c[ONE_SIDE] = lambda_in_n;
	}
	else
	{
		if(lambda_in_n < A0) // INFLOW
		{
			pipe_flow[ONE_SIDE] = INFLOW;
			common_HI_code(pPpt, lambda_in_c[ONE_SIDE], lambda_out_c[ONE_SIDE], CHOKED[ONE_SIDE], A0, T0);
		}
		else // OUTFLOW - need new values for lambda_out_c only
		{
			pipe_flow[ONE_SIDE] = OUTFLOW;
			lambda_in_c[ONE_SIDE] = lambda_in_n;
			if(phi >= 1) lambda_out_c[ONE_SIDE] =  2*A0 - lambda_in_n;	// Run open procedure
			else lambda_out_c[ONE_SIDE] = common_HN_code(pPpt, timestep, time, lambda_in_n, phi, P0, CHOKED[ONE_SIDE], pBN[ONE_SIDE]->T); // Run nozzle procedure
		}
	}
	common_UPDATE_H(pPpt, lambda_in_c, lambda_out_c, pipe_flow, CHOKED);
	delete [] lambda_in_c;
	delete [] lambda_out_c;
	delete [] pipe_flow;
	delete [] CHOKED;
}

void CEndEnvironment::NHE(CProperties* pPpt, int timestep, double time)
// ====================================================================================================
// End environment for flow to/from stagnation conditions; non-homentropic flow
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".NHE\n");}

	double lambda_in_n, lambda_out_n, AA_n, PIp, T, k;
	int pipe_flow_n;
	double *lambda_in_c, *lambda_out_c, *AA_c;
	int *pipe_flow;
	lambda_in_c = new double [NPIPES]; lambda_out_c = new double [NPIPES]; AA_c = new double [NPIPES];	pipe_flow = new int [NPIPES];
	bool* CHOKED; CHOKED = new bool [NPIPES]; 
	CHOKED[ONE_SIDE] = false; // Reset here for all cases

	lambda_in_n = (*(pCLIN[ONE_SIDE]))[R+1];
	lambda_out_n = (*(pCLOUT[ONE_SIDE]))[R+1];
	AA_n = pBN[ONE_SIDE]->AA[R+1];
	pipe_flow_n = *(pend_flow[ONE_SIDE]);

	// Iterative flow direction test
	// ----------------------------------------------------------------------------------------------------
	T = T0; // Assume inflow, use T0
	bool FINISHED = false;
	do{
		k = pPpt->gammaAir(T); 
		PIp = pow(P0/pPpt->PREF,(k-1)/(2*k));
        if(fabs(lambda_in_n/AA_n - PIp) < pPpt->ZERO_TOL || phi==0){ // ZERO FLOW
			T = T0;
			FINISHED = true;
		}
		else{
            if(lambda_in_n/AA_n < PIp){ // INFLOW - T should be T0	
                if(fabs(T - T0) < pPpt->ZERO_TOL) FINISHED = true;
                else{
					T = T0;
					//FINISHED = true; // Else can cause endless loop
				}
			}
            else{ // OUTFLOW - T should be pBN[ONE_SIDE]->T
				if(fabs(T - pBN[ONE_SIDE]->T) < pPpt->ZERO_TOL) FINISHED = true;
                else{
					T = pBN[ONE_SIDE]->T;
					FINISHED = true; // Else can cause endless loop
				}
			}
		}
	}while(!FINISHED);

	// Identify flow direction and call appropriate boundary method
	// ----------------------------------------------------------------------------------------------------
	if(fabs(lambda_in_n/AA_n - PIp) < pPpt->ZERO_TOL || phi==0){ // NOFLOW
		pipe_flow[ONE_SIDE] = NOFLOW;
		lambda_in_c[ONE_SIDE] = lambda_in_n;
		lambda_out_c[ONE_SIDE] = lambda_in_n;
		AA_c[ONE_SIDE] = AA_n;
	}
	else{
		if(lambda_in_n/AA_n > PIp){ // OUTFLOW
			pipe_flow[ONE_SIDE] = OUTFLOW;
			lambda_in_c[ONE_SIDE] = lambda_in_n; 
			AA_c[ONE_SIDE] = AA_n;
			lambda_out_c[ONE_SIDE] = common_NHN_code(pPpt, lambda_in_n, AA_n, phi, P0, CHOKED[ONE_SIDE], T, time);
		}
        else{ // INFLOW
            pipe_flow[ONE_SIDE] = INFLOW;
			common_NHI_code(pPpt, lambda_in_n, lambda_out_n, AA_n, lambda_in_c[ONE_SIDE], lambda_out_c[ONE_SIDE], AA_c[ONE_SIDE], (pPpt->USE_PHI ? phi : 1.0), P0, T, CHOKED[ONE_SIDE], SONIC, T, timestep, time, true/*true=constant pressure (valve) model or false=pressure loss (port) model*/, pPpt->NHI_TOL);
		}
	}

	// Update pipe conditions
	// ----------------------------------------------------------------------------------------------------
	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out_c, AA_c, pipe_flow, CHOKED);
    delete [] lambda_in_c;
    delete [] lambda_out_c;
    delete [] AA_c;
    delete [] pipe_flow;
    delete [] CHOKED;
}


void CEndEnvironment::Anechoic(CProperties* pPpt, int timestep, double time)
//--------------------------------------------------//
// Non-homentropic anechoic termination				//
// ------------------------------------				//
//													//
//--------------------------------------------------//
{
	double lambda_in_n, lambda_out_n, AA_n;
	double *lambda_in_c, *lambda_out_c, *AA_c;
	int *pipe_flow;
	lambda_in_c = new double [NPIPES]; lambda_out_c = new double [NPIPES]; AA_c = new double [NPIPES];	pipe_flow = new int [NPIPES];
	bool* CHOKED; CHOKED = new bool [NPIPES]; 
	CHOKED[ONE_SIDE] = false; // Reset here for all cases

	lambda_in_n = (*(pCLIN[ONE_SIDE]))[R+1];
	lambda_out_n = (*(pCLOUT[ONE_SIDE]))[R+1];
	AA_n = pBN[ONE_SIDE]->AA[R+1];

/*
//	lambda_in_c[ONE_SIDE] = lambda_in_n;
//	lambda_out_c[ONE_SIDE] = lambda_out_n;
//	AA_c[ONE_SIDE] = AA_n;


	// lambda_out is equivalent to the lambda_in arriving from the imaginary domain

	// OUTFLOW
	pipe_flow[ONE_SIDE] = OUTFLOW;
			
	// Need new values for lambda_out_c only
	lambda_in_c[ONE_SIDE] = lambda_in_n;
	AA_c[ONE_SIDE] = AA_n;

	// Calculate stagnation that would be achieved given the velocity at the boundary
	double U_temp = (lambda_in_n - lambda_out_n)/(pPpt->gammaAir() - 1);
	double A_temp = ((lambda_in_n + lambda_out_n)/2);
	double a_temp = A_temp*pPipe[ONE_SIDE]->AREF;
	double T_temp = pow(a_temp,2)/(pPpt->gammaAir()*pPpt->R_air);

	double p_temp = pow(A_temp/AA_n, (2*pPpt->gammaAir())/(pPpt->gammaAir()-1))*pPpt->PREF;	
	double rho_temp = p_temp/(pPpt->R_air*T_temp);


	//double p0_temp = p_anechoic - 0.5*rho_temp*(U_temp*pPipe[ONE_SIDE]->AREF,2);
	double p0_temp = p_temp;// + 0.5*rho_temp*(U_temp*pPipe[ONE_SIDE]->AREF,2);

	//cout << "anechoic p0_temp = " << p0_temp << endl;

// Run nozzle procedure; accounts for fully open ends since then phi = 1
//lambda_out_c[ONE_SIDE] = common_NHN_code(pPpt, lambda_in_n, AA_n, phi, P0, CHOKED[ONE_SIDE]);
//lambda_out_c[ONE_SIDE] = NHNozzle(pPpt, lambda_in_n, lambda_out_n, AA_n, 1, p0_temp, CHOKED[ONE_SIDE]);
//if(time>0.036) if(ID==1) cout << "OUTFLOW\n";

lambda_out_c[ONE_SIDE] = this->lambda_in_an;
this->lambda_out_an = lambda_in_n;

lambda_out_n = (*(pCLOUT[ONE_SIDE]))[R+1];
*/

	lambda_in_c[ONE_SIDE] = lambda_in_n;
	AA_c[ONE_SIDE] = AA_n;

	// Extrapolate a value for lambda_out from interior domain
	if(pBN[ONE_SIDE]->side == ODD)
	{
		lambda_out_c[ONE_SIDE] = 
			pPipe[ONE_SIDE]->Node[0+1].CL1[R+1];
/*
		// At odd (left-hand) end, CL1 is right-running, i.e., lambda_out
		lambda_out_c[ONE_SIDE] = 
			pPipe[ONE_SIDE]->Node[0+1].CL1[R+1] - 
				((pPipe[ONE_SIDE]->Node[0+2].CL1[R+1] - pPipe[ONE_SIDE]->Node[0+1].CL1[R+1])
				/
				(pPipe[ONE_SIDE]->Node[0+2].X - pPipe[ONE_SIDE]->Node[0+1].X))
				*(pPipe[ONE_SIDE]->Node[0+1].X - pPipe[ONE_SIDE]->Node[0+0].X);
//*/
	}
	else // EVEN
	{
		lambda_in_c[ONE_SIDE] = 
			pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1-1].CL1[R+1];

		AA_c[ONE_SIDE] = 
			pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1-1].AA[R+1];

		lambda_out_c[ONE_SIDE] = 
			pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1-1].CL2[R+1]; 
/*
		// At even (right-hand) end, CL2 is left-running, i.e., lambda_out
		lambda_out_c[ONE_SIDE] = 
			pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1-1].CL2[R+1] + 
				((pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1-1].CL2[R+1] - pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1-2].CL2[R+1])
				/
				(pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1-1].X - pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1-2].X))
				*(pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1].X - pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1-1].X);
//*/
	}

	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out_c, AA_c, pipe_flow, CHOKED);




/*
	// Generalised flow direction test
	double PIp = pow(P0/pPpt->PREF, (pPpt->gammaAir()-1)/(2*pPpt->gammaAir()));

	if(fabs(lambda_in_n/AA_n - PIp) < 1e-6 || phi==0) 
	{
		// NOFLOW
		pipe_flow[ONE_SIDE] = NOFLOW;

		lambda_in_c[ONE_SIDE] = lambda_in_n;
		lambda_out_c[ONE_SIDE] = lambda_in_n;
		AA_c[ONE_SIDE] = AA_n;
if(time>0.036) if(ID==1) cout << "NOFLOW\n";
	}
	else
	{
		if(lambda_in_n/AA_n < PIp)	
		{
			// INFLOW
			pipe_flow[ONE_SIDE] = INFLOW;

			common_NHI_code(pPpt, lambda_in_n, lambda_out_n, AA_n, 
						lambda_in_c[ONE_SIDE], lambda_out_c[ONE_SIDE], AA_c[ONE_SIDE], 
						(pPpt->USE_PHI ? phi : 1.0), P0, T0, 
						CHOKED[ONE_SIDE], SONIC);
if(time>0.036) if(ID==1) cout << "INFLOW\n";
		}
		else 
		{
			// OUTFLOW
			pipe_flow[ONE_SIDE] = OUTFLOW;
			
			// Need new values for lambda_out_c only
			lambda_in_c[ONE_SIDE] = lambda_in_n;
			AA_c[ONE_SIDE] = AA_n;

// Run nozzle procedure; accounts for fully open ends since then phi = 1
//lambda_out_c[ONE_SIDE] = common_NHN_code(pPpt, lambda_in_n, AA_n, phi, P0, CHOKED[ONE_SIDE]);
lambda_out_c[ONE_SIDE] = NHNozzle(pPpt, lambda_in_n, lambda_out_n, AA_n, phi, P0, CHOKED[ONE_SIDE]);
if(time>0.036) if(ID==1) cout << "OUTFLOW\n";
		}
	}

	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out_c, AA_c, pipe_flow, CHOKED);
*/
}


void CEndEnvironment::PrintToScreen(CProperties* pPpt)
{
	//if(phi!=0) // No print out if closed end
	{
		pPpt->Out(Underline(Identify(), "=", "", strDesc)); pPpt->Out("\n");
//			if(ANECHOIC)
//			{
//				pPpt->Out("Anechoic lambda_in_an, lambda_in_an\t\t\t=\t"); pPpt->Out(lambda_in_an); pPpt->Out("\n");
//				pPpt->Out("Anechoic lambda_out_an, lambda_out_an\t\t\t=\t"); pPpt->Out(lambda_out_an); pPpt->Out("\n");
//			}
//			else
/*
		{
			if(!pPpt->HOMENTROPIC && SONIC) pPpt->Out("\t\tCHOKED flow in end environment\n");
			if(VAR_P0)
			{
				pPpt->Out("Varying stagnation pressure, P0\t\t\t=\t"); pPpt->Out(P0); pPpt->Out(" bar\n");
				pPpt->Out("Varying stagnation temperature, T0\t\t\t=\t"); pPpt->Out(T0); pPpt->Out(" K\n");
				pPpt->Out("\n");
			}
			if(VAR_PHI)
			{
				pPpt->Out("Varying nozzle area ratio, phi\t\t\t\t=\t"); pPpt->Out(phi); pPpt->Out("\n");
				pPpt->Out("\n");
			}
		}
*/		
		if(VAR_P0) pPpt->Out("Variable stagnation conditions:\n");
		else pPpt->Out("Constant stagnation conditions:\n");
		pPpt->Out("Stagnation pressure, P0\t\t\t\t\t=\t"); pPpt->Out(P0); pPpt->Out(" bar\n");
		pPpt->Out("Stagnation temperature, T0\t\t\t\t=\t"); pPpt->Out(T0); pPpt->Out(" K\n");
		pPpt->Out("Nozzle area ratio, phi\t\t\t\t\t=\t"); pPpt->Out(phi); pPpt->Out("\n");
		pPpt->Out("\n");
		pPpt->Out("\n");
	}
}

void CEndEnvironment::PrintToFile(CProperties* pPpt, int timestep, double time, double ca_elapsed, double ca, CPipe* Pipe)
// ============================================================ //
// Prints data to file.											//
//																//	
// ============================================================ //
{	
	if(timestep == 0) {

		std::string temp_str = "Object results file for ";
		temp_str += Identify();
		fprintf(OUTPUT_FILE,"%s\n", Underline(StringToChar(temp_str), "-"));	

		fprintf(OUTPUT_FILE,"%s", "Time(s)");
		if(!pPpt->CONTINUOUS) fprintf(OUTPUT_FILE,"\t%s%c%s\t%s%c%s", "Elapsed (", Deg(), "CA)", "Crank angle (", Deg(), "CA)");
    	
		// Print environment stagnation conditions to file
		fprintf(OUTPUT_FILE,"\t%s\t%s", "Environment stagnation pressure (bar)", "Environment stagnation temperature (K)");

		// Print internal boundary conditions to file
		fprintf(OUTPUT_FILE, "\t%s\t%s\t%s\t%s", "Static pressure (bar)", "Static temperature (K)", "Density (kg/m^3)", "Velocity (m/s)");
		
		// Print flow direction and whether choked to file
		fprintf(OUTPUT_FILE, "\t%s\t%s\t%s\t%s", "Direction", "Direction No.", "Choked", "Choked No.");

		// Print lambda_in and lambda_out to file
		fprintf(OUTPUT_FILE, "\t%s\t%s", "lambda_in", "lambda_out");

		// Print lambda_in_star and lambda_out_star to file
		fprintf(OUTPUT_FILE, "\t%s\t%s", "lambda_in_star", "lambda_out_star");
		
		if(ROTOR_EXIT_TAPPING) fprintf(OUTPUT_FILE,"\t%s\t%s", "PR", "MFP (kg.K^0.5/bar)");
		fprintf(OUTPUT_FILE,"\n");
	}

	if(timestep % freq == 0 && time >= print_from_time) { // Print data at the specified sampling frequency
		
		fprintf(OUTPUT_FILE,"%f", time);
		if(!pPpt->CONTINUOUS) fprintf(OUTPUT_FILE,"\t%f\t%f", ca_elapsed, ca);

		// Print environment stagnation conditions to file
		fprintf(OUTPUT_FILE, "\t%f\t%f", P0, T0);

		// Print internal boundary conditions to file
		fprintf(OUTPUT_FILE, "\t%f\t%f\t%f\t%f", pBN[ONE_SIDE]->p_dash * pPpt->PREF, pBN[ONE_SIDE]->T, pBN[ONE_SIDE]->rho, pBN[ONE_SIDE]->U*pPipe[ONE_SIDE]->AREF);

		// Print flow direction and whether choked to file
		fprintf(OUTPUT_FILE, "\t%s\t%d\t%s\t%d", InflowOrOutflow(*(pend_flow[ONE_SIDE])), InflowOrOutflowNumber(*(pend_flow[ONE_SIDE])), Choked(pBN[ONE_SIDE]->CHOKED), ChokedNo(pBN[ONE_SIDE]->CHOKED));

		// Print lambda_in and lambda_out to file
		fprintf(OUTPUT_FILE, "\t%f\t%f", (*(pCLIN[0]))[R + 1], (*(pCLOUT[0]))[R + 1]);

		// Print lambda_in_star and lambda_out_star to file
		fprintf(OUTPUT_FILE, "\t%f\t%f", CLIN_STAR[0][R + 1], CLOUT_STAR[0][R + 1]);
		
		if(ROTOR_EXIT_TAPPING) {
			// Interpolate p2_pipe for tapping placed at p2_loc		
			// Locate the two nodes either side of the location of this tapping
			int S=0;
			CNode p2_loc_node;
			if(Pipe[p2_pipe].N > 1) {
				while(Pipe[p2_pipe].Node[S].x < p2_loc*Pipe[p2_pipe].length && S < Pipe[p2_pipe].N - 1) ++S; // p2_loc is between S and S-1
				if(S==0) S=1;
			
				p2_loc_node = Pipe[p2_pipe].Node[S-1] +
					(Pipe[p2_pipe].Node[S] - Pipe[p2_pipe].Node[S-1])*
					((p2_loc*Pipe[p2_pipe].length - Pipe[p2_pipe].Node[S-1].x)
					/(Pipe[p2_pipe].Node[S].x - Pipe[p2_pipe].Node[S-1].x));
			}
			else p2_loc_node = Pipe[p2_pipe].Node[S]; // For single node pipes
		
			double mdot_temp;
			if(end[ONE_SIDE]==ODD) {
				mdot_temp = pPipe[ONE_SIDE]->Node[0].mdot;
			}
			else {
				mdot_temp = pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1].mdot;
			}

			fprintf(OUTPUT_FILE,"\t%f\t%f", P0/(p2_loc_node.p_dash*pPpt->PREF), mdot_temp*sqrt(T0)/P0);
		}
		fprintf(OUTPUT_FILE,"\n");
	}

	if(PRINT_DEBUG_FILE){
		if(timestep==0){
			std::string temp_str = "Debug file for ";
			temp_str += Identify();
			fprintf(OUTPUT_FILE_DEBUG,"%s\n", Underline(StringToChar(temp_str), "-"));	
			fprintf(OUTPUT_FILE_DEBUG,"%s", "Time(s)");
			if(!pPpt->CONTINUOUS) fprintf(OUTPUT_FILE_DEBUG,"\t%s%c%s\t%s%c%s", "Elapsed (", Deg(), "CA)", "Crank angle (", Deg(), "CA)");
    		fprintf(OUTPUT_FILE_DEBUG,"\t%s\t%s\t%s\t%s\t%s\t%s", "pCLIN[ONE_SIDE]", "pCLOUT[ONE_SIDE]", "Direction", "Direction No.", "Choked", "Choked No.");
    		fprintf(OUTPUT_FILE_DEBUG,"\n");
		}

		if(timestep%freq==0){ // Print data at the specified sampling frequency
			fprintf(OUTPUT_FILE_DEBUG,"%f", time);
			if(!pPpt->CONTINUOUS) fprintf(OUTPUT_FILE_DEBUG,"\t%f\t%f", ca_elapsed, ca);// Periodic
			fprintf(OUTPUT_FILE_DEBUG,"\t%f\t%f\t%s\t%d\t%s\t%d", (*(pCLIN[ONE_SIDE]))[R+1], (*(pCLOUT[ONE_SIDE]))[R+1], InflowOrOutflow(*(pend_flow[ONE_SIDE])), InflowOrOutflowNumber(*(pend_flow[ONE_SIDE])), Choked(pBN[ONE_SIDE]->CHOKED), ChokedNo(pBN[ONE_SIDE]->CHOKED));
			fprintf(OUTPUT_FILE_DEBUG,"\n");
		}
	}
}
void CEndEnvironment::Set_P0_gradual(double P0_temp)
{
	P0_desired = P0_temp; // Record requested value
}