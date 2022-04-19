// Transmissive.cpp: implementation of the CTransmissive class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Tools.h"
#include "Transmissive.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTransmissive::CTransmissive()
{

}

CTransmissive::~CTransmissive()
{
	delete [] Results;
}

// Copy constructor
CTransmissive::CTransmissive(const CTransmissive& inTransmissive)
{
	ID = inTransmissive.ID;
	EX = inTransmissive.EX;
	A0test = inTransmissive.A0test;	
	Ts = inTransmissive.Ts;			
//	TREF = inTransmissive.TREF;

	end = inTransmissive.end;
	
	// Yes, copy pointers, this is what is desired
	pPipe = inTransmissive.pPipe;
	pBN = inTransmissive.pBN;
	pCLIN = inTransmissive.pCLIN;
	pCLOUT = inTransmissive.pCLOUT;
	pend_flow = inTransmissive.pend_flow;
	pPathLine = inTransmissive.pPathLine;
}

CTransmissive& CTransmissive::operator=(const CTransmissive& inTransmissive)
{
	if(this != &inTransmissive)
	{
		ID = inTransmissive.ID;
		EX = inTransmissive.EX;
		A0test = inTransmissive.A0test;	
		Ts = inTransmissive.Ts;			
//		TREF = inTransmissive.TREF;

		end = inTransmissive.end;
	
		// Yes, copy pointers, this is what is desired
		pPipe = inTransmissive.pPipe;
		pBN = inTransmissive.pBN;
		pCLIN = inTransmissive.pCLIN;
		pCLOUT = inTransmissive.pCLOUT;
		pend_flow = inTransmissive.pend_flow;
		pPathLine = inTransmissive.pPathLine;
	}
	return *this;
}

void CTransmissive::Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rTRANSPIPES, int** &rTRANSPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, std::string param_dir, int assyid, string parent_assy_res_dir, int ntransm, int npipesinassy, string calling_object_str)
{
	int m;
	InitialiseGen(pPpt, pPipes, rPipe, rTRANSPIPES, rTRANSPIPES_ENDS, rENDCORR, id, ex, npipes, assyid, "CTransmissive", parent_assy_res_dir);
	std::string bcname_str = "TRANSM";
	ReadInput(pPpt, ConstructString(pPpt, param_dir, bcname_str, EX, ID));

	// Boundary name
	NAME = new int [NPIPES];
	NAME[ONE_SIDE] = TRANSMISSIVE;	// Always required and used

	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Initialise\n");}

//	LEFT_TO_RIGHT = true;
	WAVEFORM_TO_INTERNAL = true;

//	pipe_flow[ONE_SIDE] = INFLOW;


	// Apply properties
	// ================

	if(pPpt->HOMENTROPIC) A0test = pow(this->press_low/pPpt->PREF, pPpt->Q);
	else A0test = sqrt(this->Ts/(EX ? pPpt->TREFe : pPpt->TREFi));

	// Initialise tappings
	nPipesInAssy = npipesinassy;
	if (p2_pipe >= nPipesInAssy) {	// Check is chosen p2_pipe is actually available
		// If not, exit safely and print reasonto screen
		pPpt->Out("In the input file for ");
		pPpt->Out(Identify());
		pPpt->Out(" p2_pipe is set to ");
		pPpt->Out(p2_pipe);
		pPpt->Out(" but this pipe does not exist in the assembly. There are only ");
		pPpt->Out(nPipesInAssy);
		pPpt->Out(" pipe(s) in the assembly. Exiting.\n");
		exit(1);					
	}

	Measurements = new double [nmeasurements];
	MeasureNode = new CNode [nmeasurements];
	for(m=0; m<nmeasurements; ++m) Measurements[m] = loc_measure[m];

	Results = new double ** [nmeasurements];
	for(m=0; m<nmeasurements; ++m)
	{ 
		Results[m] = new double* [max_pts];
		for(int d=0; d<max_pts; ++d)
		{
			Results[m][d] = new double [num_props_measured+1+1];
			for(int p=0; p<(num_props_measured+1+1); ++p) Results[m][d][p] = 0;
		}
	}

	rowcounter=0;
	sample_factor=1;

	// Setup and open results file
	SetupFiles(pPpt);

	// Pulse calculations
	time_in_period = 0;
	time_in_period_prev = 0;
//	p_wave = 1;
//	p_wave_prev = 1; // Set in loadwaveform
//	MFR_wave = 1;
//	MFR_wave_prev = 1; // Set in loadwaveform
	i_prev = -1;

	u_WAVEFORM = 0;
	T_WAVEFORM = 300;
	
	whole_number_of_pulses = 0;
	START_OF_NEW_PULSE = true;
	START_OF_SIMULATION = true;
//	PRINT_PULSE_INFO = false;
//	Pressure_Desired_SET = false;
	iterations_during_pulse = 0;
//	iterations_during_last_pulse = 0;
//	iterations_during_last_last_pulse = 0;
//	iterations_during_first_recorded_pulse = 0;
//	iterations_during_first_recorded_pulse_guess = 0;
//	Interpolated_Pressure_Desired = 1;
//	dim_Pressure_To_Apply_This_Pulse = 0;
	MATCH = false;
	MATCHED = false;
	nTransm = ntransm;
/*
	if(nTransm>0)
	{
		if(ID==0) ATTEMPT_MATCH = true; else ATTEMPT_MATCH = false; // Start by allowing the first transmissive to match if there are multiple ones
	}
	else 
//*/
	ATTEMPT_MATCH = true;

	p_SEEN_LESS_THAN_DESIRED = true;
	v_DESIRED_GREATER_THAN_SEEN = true;
	del_p = 0.1;
	del_v = 0.1;

	if(!CONSTANT) period = (1/f); // Otherwise already read in directly from file
	
	// Labels
	TIME = 0;
	TIME_IN_PERIOD = 1;

	DESIRED_P = 2;
	APPLIED_P = 3;
	SEEN_P = 4;
	NEXT_P = 5;
	THIS_P = 6;

	DESIRED_T = 10;
	DESIRED_T0 = 9;
	APPLIED_T = 26;
	SEEN_T = 23;
	NEXT_T = 25;
	THIS_T = 24;
	
	DESIRED_V = 11;
	APPLIED_V = 12;
	SEEN_V = 13;
	NEXT_V = 14;
	THIS_V = 15;
	
	DESIRED_MFR = 7;
	APPLIED_MFR = 27;
	SEEN_MFR = 8;
	NEXT_MFR = 28;
	THIS_MFR = 29;
		
	VEL_GRAD = 16;
	P_VEL = 17;
	P_VEL_GRAD = 18;
	UNDER_EST = 19;
	UNDER_EST_PREV = 20;
	UNDER_EST_FACT = 21;
	K_LOSS = 22;

	WAVE_TIME = 0;
	WAVE_PARAMETER = 1;
	WAVE_MFR = 1;

	p_des_app_seen = new double* [nsamples];
	for(int i=0; i<nsamples; ++i)
	{
		p_des_app_seen[i] = new double [30]; 
		// time, time_in_period, pressure desired, applied, seen, next

		// Set appropriate times and unit pressures in these arrays
		p_des_app_seen[i][TIME] = 0 + i*(period/nsamples);
		p_des_app_seen[i][TIME_IN_PERIOD] = 0 + i*(period/nsamples);

		p_des_app_seen[i][DESIRED_P] = 1;
		p_des_app_seen[i][APPLIED_P] = 1;
		p_des_app_seen[i][SEEN_P] = 1;
		p_des_app_seen[i][NEXT_P] = 1;
		p_des_app_seen[i][THIS_P] = 1;
		
		p_des_app_seen[i][DESIRED_T] = 0;
		p_des_app_seen[i][DESIRED_T0] = 0;
		p_des_app_seen[i][APPLIED_T] = 0;
		p_des_app_seen[i][SEEN_T] = 0;
		p_des_app_seen[i][NEXT_T] = 0;
		p_des_app_seen[i][THIS_T] = 0;

		p_des_app_seen[i][DESIRED_V] = 0;
		p_des_app_seen[i][APPLIED_V] = 0;
		p_des_app_seen[i][SEEN_V] = 0;
		p_des_app_seen[i][NEXT_V] = 0;
		p_des_app_seen[i][THIS_V] = 0;

		p_des_app_seen[i][DESIRED_MFR] = 0;
		p_des_app_seen[i][APPLIED_MFR] = 0;
		p_des_app_seen[i][SEEN_MFR] = 0;
		p_des_app_seen[i][NEXT_MFR] = 0;
		p_des_app_seen[i][THIS_MFR] = 0;

		p_des_app_seen[i][VEL_GRAD] = 0;
		p_des_app_seen[i][P_VEL] = 0;
		p_des_app_seen[i][P_VEL_GRAD] = 0;
		p_des_app_seen[i][UNDER_EST] = 0;
		p_des_app_seen[i][UNDER_EST_PREV] = 0;
		p_des_app_seen[i][UNDER_EST_FACT] = 1;
		p_des_app_seen[i][K_LOSS] = 1;
	}

////////////////////////////////////////////////////////////////

	if(CONSTANT) // Constant conditions
	{
		for(int i_fill=0; i_fill<nsamples; ++i_fill)
		{
			p_des_app_seen[i_fill][DESIRED_P] = constantp;
			p_des_app_seen[i_fill][DESIRED_T0] = constantT;
			p_des_app_seen[i_fill][DESIRED_T] = Ts;
			p_des_app_seen[i_fill][DESIRED_V] = v;
			//p_des_app_seen[i_fill][DESIRED_MFR] = ?;
		}

		// Initialise p_wave_prev to the first value in the array, else get funny values on first pulse
		p_wave = p_des_app_seen[0][DESIRED_P];
		p_wave_prev = p_des_app_seen[0][DESIRED_P];
		//cout << "p_wave_prev = " << p_wave_prev << endl;

		T_wave = p_des_app_seen[0][DESIRED_T];
		T_wave_prev = p_des_app_seen[0][DESIRED_T];

		v_wave = p_des_app_seen[0][DESIRED_V];
		v_wave_prev = p_des_app_seen[0][DESIRED_V];

		//mfr_wave = p_des_app_seen[0][DESIRED_MFR];
		//mfr_wave_prev = p_des_app_seen[0][DESIRED_MFR];

		p_wave_desired = p_wave;
		T_wave_desired = T_wave;
		v_wave_desired = v_wave;
		mfr_wave_desired = mfr_wave;
	}
	else
	{
		// Calculate or read the desired waveform
		if(shape==3)
		{
			// Load pressure waveform
			LoadWaveform(/*FILE_ps*/ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->TRANSM_DIR), FILE_ps),
							pPpt, TOTAL_PRESS, ADJUST_TIME, ADJUST_PRESSURE, 
							press_high, press_low, press_min, DESIRED_P, p_wave, p_wave_prev);

			if(USE_FILE_Ts)
			{
				LoadWaveform(/*FILE_Ts*/ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->TRANSM_DIR), FILE_Ts), 
							pPpt, false, ADJUST_TIME, ADJUST_TSM, 
								Tsm_high, Tsm_low, Tsm_min, DESIRED_T, T_wave, T_wave_prev);

				LoadWaveform(/*WAVE_FILE_T0m*/ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->TRANSM_DIR), WAVE_FILE_T0m), 
							pPpt, false, ADJUST_TIME, ADJUST_T0M, 
							T0m_high, T0m_low, T0m_min, DESIRED_T0, T0m_wave, T0m_wave_prev);

				LoadWaveform(/*WAVE_FILE_v*/ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->TRANSM_DIR), WAVE_FILE_v), 
							pPpt, false, ADJUST_TIME, false, 
							300, 0, 0, DESIRED_V, v_wave, v_wave_prev);

				LoadWaveform(/*WAVE_FILE_MFR*/ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->TRANSM_DIR), WAVE_FILE_MFR), 
							pPpt, false, ADJUST_TIME, ADJUST_MFR, 
							MFR_high, MFR_low, MFR_min, DESIRED_MFR, mfr_wave, mfr_wave_prev);
			}
			else // Use constant conditions
			{
				for(int i_fill=0; i_fill<nsamples; ++i_fill)
				{
					p_des_app_seen[i_fill][DESIRED_T0] = Ts;
					p_des_app_seen[i_fill][DESIRED_T] = Ts;
					p_des_app_seen[i_fill][DESIRED_V] = v;
				}
			}
			
			//LoadWaveform(/*WAVE_FILE_MFR*/ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->TRANSM_DIR), WAVE_FILE_MFR), 
			//				pPpt, false, ADJUST_TIME, ADJUST_MFR, 
			//				MFR_high, MFR_low, MFR_min, DESIRED_MFR, MFR_wave, MFR_wave_prev);

	//		LoadWaveformMFR(/*WAVE_FILE_MFR*/ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->TRANSM_DIR), WAVE_FILE_MFR), pPpt);

			// Initialise p_wave_prev to the first value in the array, else get funny values on first pulse
			p_wave = p_des_app_seen[0][DESIRED_P];
			p_wave_prev = p_des_app_seen[0][DESIRED_P];
			
			T_wave = p_des_app_seen[0][DESIRED_T];
			T_wave_prev = p_des_app_seen[0][DESIRED_T];

			v_wave = p_des_app_seen[0][DESIRED_V];
			v_wave_prev = p_des_app_seen[0][DESIRED_V];

			mfr_wave = p_des_app_seen[0][DESIRED_MFR];
			mfr_wave_prev = p_des_app_seen[0][DESIRED_MFR];

			p_wave_desired = p_wave;
			T_wave_desired = T_wave;
			v_wave_desired = v_wave;
			mfr_wave_desired = mfr_wave;

	/*
			for(int i=0; i<this->nsamples; ++i)
			{
				char* tab = "  ";
				cout << tab << p_des_app_seen[i][TIME]
					 << tab << p_des_app_seen[i][TIME_IN_PERIOD]

					 << tab << p_des_app_seen[i][DESIRED_P]
					 << tab << p_des_app_seen[i][APPLIED_P]
					 << tab << p_des_app_seen[i][SEEN_P]
					 << tab << p_des_app_seen[i][NEXT_P]
					 << tab << p_des_app_seen[i][THIS_P]

					 << tab << p_des_app_seen[i][DESIRED_T]
					 << tab << p_des_app_seen[i][DESIRED_T0]
					 << tab << p_des_app_seen[i][SEEN_T]
					 << tab << p_des_app_seen[i][NEXT_T]
					 << tab << p_des_app_seen[i][THIS_T]

					 << tab << p_des_app_seen[i][DESIRED_V]
					 << tab << p_des_app_seen[i][APPLIED_V]
					 << tab << p_des_app_seen[i][SEEN_V]
					 << tab << p_des_app_seen[i][NEXT_V]
					 << tab << p_des_app_seen[i][THIS_V]

	//				 << tab << p_des_app_seen[i][DESIRED_MFR]
					 << tab << p_des_app_seen[i][APPLIED_MFR]
	//				 << tab << p_des_app_seen[i][SEEN_MFR]
					 << tab << p_des_app_seen[i][NEXT_MFR]
					 << tab << p_des_app_seen[i][THIS_MFR]

	//				 << tab << p_des_app_seen[i][VEL_GRAD]
	//				 << tab << p_des_app_seen[i][P_VEL]
	//				 << tab << p_des_app_seen[i][P_VEL_GRAD]
	//				 << tab << p_des_app_seen[i][UNDER_EST]
	//				 << tab << p_des_app_seen[i][UNDER_EST_PREV]
	//				 << tab << p_des_app_seen[i][UNDER_EST_FACT]
	//				 << tab << p_des_app_seen[i][K_LOSS]
	//				 << tab << p_des_app_seen[i][SEEN_T]

					 << "\n";
			}
	//*/
		}
		else
		{
			if(shape==4) ReadFourierFile(pPpt, ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->TRANSM_DIR), FOURIER_FILE));
							
			for(int i_fill=0; i_fill<nsamples; ++i_fill)
			{
				if(p_des_app_seen[i_fill][TIME_IN_PERIOD] < phi*period) // We are within the pulse section
				{
					if(shape==0) // Triangular
					{
						double gradient;
						gradient = (press_high - press_low)/(phi*period*0.5);
						if(p_des_app_seen[i_fill][TIME_IN_PERIOD] < phi*period*0.5) p_des_app_seen[i_fill][DESIRED_P] = press_low + gradient*p_des_app_seen[i_fill][TIME_IN_PERIOD]; // Positive gradient
						else p_des_app_seen[i_fill][DESIRED_P] = press_high - gradient*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - (phi*period*0.5)); // Negative gradient
					}
					else
					{
						if(shape==1) p_des_app_seen[i_fill][DESIRED_P] = press_high;	// Square waveform 
						else
						{
							if(shape==2)
							{
								p_des_app_seen[i_fill][DESIRED_P] = press_low + ( ((-1*cos(2*PI*p_des_app_seen[i_fill][TIME_IN_PERIOD]/(phi*period)) + 1)/2) * (press_high - press_low));// Cosine waveform
								//cout << "p_des_app_seen[" << i_fill << "][" << DESIRED_P << "] = " << p_des_app_seen[i_fill][DESIRED_P] << endl;
							}
							else
							{
								if(shape==4)
								{
									//fr_An_an_bn[col=fr(0),amp(1),an(2),bn(3)][row=0->ncomponents]
									double a0, f, an, bn, ti, sigsum; 
									ti = p_des_app_seen[i_fill][TIME_IN_PERIOD];
									sigsum = 0;
									a0 = fr_An_an_bn[1][0];			// Steady (dc), or mean level a0
							
	//pPpt->Out("a0 = "); pPpt->Out(a0); pPpt->Out("\n");
									
									for(int comp = 1; comp<ncomponents; ++comp) // Start at comp==1 since the first (a0) component is added separately
									{
										f = fr_An_an_bn[0][comp];
										an = fr_An_an_bn[2][comp];
										bn = fr_An_an_bn[3][comp];
										sigsum = sigsum + (an*cos(2*PI*f*ti) + bn*sin(2*PI*f*ti));
	/*
										pPpt->Out("comp = "); pPpt->Out(comp); pPpt->Out("\n");
										pPpt->Out("f = "); pPpt->Out(f); pPpt->Out("\n");
										pPpt->Out("an = "); pPpt->Out(an); pPpt->Out("\n");
										pPpt->Out("bn = "); pPpt->Out(bn); pPpt->Out("\n");
										pPpt->Out("ti = "); pPpt->Out(ti); pPpt->Out("\n");
										pPpt->Out("new sigsum = "); pPpt->Out(sigsum); pPpt->Out("\n");
	*/
									}
								
									p_des_app_seen[i_fill][DESIRED_P] = a0 + 1*sigsum;
	//pPpt->Out("p_des_app_seen["); pPpt->Out(i_fill); pPpt->Out("]["); pPpt->Out(DESIRED_P); pPpt->Out("] = "); pPpt->Out(p_des_app_seen[i_fill][DESIRED_P]); pPpt->Out("\n");
	//exit(1);							
								}
								else
								{
									pPpt->Out("Unknown input waveform for transmissive boundary. Exiting...\n");
									exit(1);
								}
							}
						}
					}
				}
				else p_des_app_seen[i_fill][DESIRED_P] = press_low;

				p_des_app_seen[i_fill][DESIRED_V] = v;
			}

			// Initialise p_wave_prev to the first value in the array, else get funny values on first pulse
			p_wave = p_des_app_seen[0][DESIRED_P];
			p_wave_prev = p_des_app_seen[0][DESIRED_P];
			//cout << "p_wave_prev = " << p_wave_prev << endl; exit(1);

			T_wave = p_des_app_seen[0][DESIRED_T];
			T_wave_prev = p_des_app_seen[0][DESIRED_T];

			v_wave = p_des_app_seen[0][DESIRED_V];
			v_wave_prev = p_des_app_seen[0][DESIRED_V];

			mfr_wave = p_des_app_seen[0][DESIRED_MFR];
			mfr_wave_prev = p_des_app_seen[0][DESIRED_MFR];

			p_wave_desired = p_wave;
			T_wave_desired = T_wave;
			v_wave_desired = v_wave;
			mfr_wave_desired = mfr_wave;
		}
	}

////////////////////////////////////////////////////////////////
	
	sum_flow_velocity_over_pulse = new double [nmeasurements];
	average_flow_velocity_last_pulse = new double [nmeasurements];

	sum_press_velocity_over_pulse_uplusa = new double [nmeasurements];
	average_press_velocity_last_pulse_uplusa = new double [nmeasurements];

	sum_press_velocity_over_pulse_uminusa = new double [nmeasurements];
	average_press_velocity_last_pulse_uminusa = new double [nmeasurements];

	st = new double [nmeasurements];
	beta = new double [nmeasurements];
	mst = new double [nmeasurements];
	pmst_uplusa = new double [nmeasurements];
	pmst_uminusa = new double [nmeasurements];

	ps_loc = new double [nmeasurements];
	Ts_loc = new double [nmeasurements];
	Velocity = new double [nmeasurements];
	p0_loc = new double [nmeasurements];
	T0_loc = new double [nmeasurements];
	MassFlowRate = new double [nmeasurements];
	Re = new double [nmeasurements];
	PR = new double [nmeasurements];
	MFP = new double [nmeasurements];
	

	for(m=0; m<nmeasurements; ++m)
	{
		sum_flow_velocity_over_pulse[m] = 0;
		average_flow_velocity_last_pulse[m] = 0;

		sum_press_velocity_over_pulse_uplusa[m] = 0;
		average_press_velocity_last_pulse_uplusa[m] = 0;

		sum_press_velocity_over_pulse_uminusa[m] = 0;
		average_press_velocity_last_pulse_uminusa[m] = 0;

		st[m] = 0;
		beta[m] = 0;
		mst[m] = 0;
		pmst_uplusa[m] = 0;
		pmst_uminusa[m] = 0;

		ps_loc[m] = 1;
		Ts_loc[m] = 300;
		Velocity[m] = 0;
		p0_loc[m] = 1;
		T0_loc[m] = 300;
		MassFlowRate[m] = 0;
		Re[m] = 0;
		PR[m] = 1;
		MFP[m] = 0;
	}
	Pressure_0_prev = 0;
	Pressure_0 = 0;
	MFR_0_prev = 0;
	MFR_0 = 0;
	VEL_0_prev = 0;
	VEL_0 = 0;
	vel_grad = 0;
	mfr_under_est = 0;
	P_VEL_0_prev = 0;
	P_VEL_0 = 0;
	TEMP_0_prev = 300;
	TEMP_0 = 300;
	UPSTREAM = 0;
	DOWNSTREAM = 1;
	M = new double [2];
	M[DOWNSTREAM] = 0.0;
	M[UPSTREAM] = 0.0;

	PRINT_HEADERS = false; // Column headers have not been printed yet
}

void CTransmissive::SetupFiles(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".SetupFiles\n");}
	std::string res_str;
	if(this->EX) res_str = "res_ex_transm"; else res_str = "res_in_transm";
	res_str += ID + 48;
	res_str += pPpt->strFileExt;		// Add the file extension
	FILE_LOC = fopen(ConstructString(pPpt, RES_DIR, res_str), "w");

	if(this->EX) res_str = "res_ex_transm"; else res_str = "res_in_transm";
	res_str += ID + 48;
	res_str += "_p_des";
	res_str += pPpt->strFileExt;		// Add the file extension
	FILE_p_des_app_seen = fopen(ConstructString(pPpt, RES_DIR, res_str), "w");
}

void CTransmissive::RunBoundary(CProperties* pPpt, CTime &rMyTime, double DELZe, double DELZi, double ca, double del_ca, double TIMEe, double TIMEi, int timestep, bool STEADY)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".RunBoundary\n");}
	if(!ANECHOIC_ONLY) // Do nothing if all we want is an anechoic termination 
	{
//	if(pPpt->HOMENTROPIC)
//		HCyl(pPpt, DELZe, DELZi);
//	else
	{
//		cout << "Running Transmissive!" << endl;
//		if(!(CALIBRATE))// && TRANSM_PTR_EXISTS))
		{
//if(timestep>=3369) cout << "top Timestep = " << timestep << endl;	
			PreProcess(rMyTime, DELZe, DELZi, pPpt, ca, del_ca, TIMEe, TIMEi, timestep, STEADY);
//if(timestep>=3369) cout << "middle Timestep = " << timestep << endl;	
			Transmissive(rMyTime, DELZe, DELZi, pPpt, ca, del_ca, TIMEe, TIMEi, timestep, STEADY);
//if(timestep>=3369) cout << "bottom Timestep = " << timestep << endl;	
		}
//		TransmissiveOld(rMyTime, DELZe, DELZi, pPpt, ca, del_ca, TIMEe, TIMEi, timestep, STEADY);
	}
	}
}


void CTransmissive::PreProcess(CTime &rMyTime, double DELZe, double DELZi, CProperties* pPpt, double ca, double del_ca, double TIMEe, double TIMEi, int timestep, bool STEADY)
{	
//--------------------------------------------------//
//--------------------------------------------------//
	int i_fill, m;
	//bool converged;
	//int counter, TEMP;
	int i=0;
/*	
	double lambda_in_c2, AA_c2, lambda_out1, lambda_out2;
	double CLIN_WAVEFORM, AA_WAVEFORM, A_WAVEFORM, U_WAVEFORM;

	// Calculate AA and lambda_in_star_n for the waveform being applied
	double* lambda_in_star_n;
	lambda_in_star_n = new double [2];
	double* lambda_in_star_test;
	lambda_in_star_test = new double [2];
*/
	double time;
	if(EX) time = TIMEe;
	else time = TIMEi;

	time_in_period_prev = time_in_period;
	time_in_period = fmod(TIMEe, period);
//cout << "time_in_period = " << time_in_period << endl;
	if(time_in_period_prev>time_in_period) // e.g. at the start of a new pulse
		time_in_period_prev = time_in_period_prev - period; // Will give a -ve but that's ok

	double integer_part;
	modf(TIMEe/period, &integer_part);

  if(START_OF_SIMULATION)
  {
    if(nMatchCycles>0)
    {
      double switch_matching_integer;
      modf(whole_number_of_pulses/nMatchCycles, &switch_matching_integer);
      if(int(fmod(double(ID), 2))==0) // If transmissive has even ID
      {
        if(int(fmod(switch_matching_integer, 2))==0 && INTELLIGENT) // If pulse is also even
          ATTEMPT_MATCH = true; 
        else ATTEMPT_MATCH = false;
      }
      else // Transmissive has odd ID
      {
        if(int(fmod(switch_matching_integer, 2))!=0 && INTELLIGENT) // If pulse is also odd
          ATTEMPT_MATCH = true; 
        else ATTEMPT_MATCH = false;
      }
    }
    else 
    {
      if(INTELLIGENT) ATTEMPT_MATCH = true; else ATTEMPT_MATCH = false;
    }
    cout << "ID = " << ID << ": ATTEMPT_MATCH = " << TrueOrFalse(ATTEMPT_MATCH) << endl;
  }

	// Test for start of a new pulse (does not include the start of the new pulse at the start of the simulation)
	if(whole_number_of_pulses != int(integer_part))
	{
		START_OF_NEW_PULSE = true;

    if(nMatchCycles>0)
    {
      double switch_matching_integer;
      modf(whole_number_of_pulses/nMatchCycles, &switch_matching_integer);
      if(int(fmod(double(ID), 2))==0) // If transmissive has even ID
      {
        if(int(fmod(switch_matching_integer, 2))==0 && INTELLIGENT) // If pulse is also even
          ATTEMPT_MATCH = true; 
        else ATTEMPT_MATCH = false;
      }
      else // Transmissive has odd ID
      {
        if(int(fmod(switch_matching_integer, 2))!=0 && INTELLIGENT) // If pulse is also odd
          ATTEMPT_MATCH = true; 
        else ATTEMPT_MATCH = false;
      }
    }
    else 
    {
      if(INTELLIGENT) ATTEMPT_MATCH = true; else ATTEMPT_MATCH = false;
    }
    cout << "ID = " << ID << ": ATTEMPT_MATCH = " << TrueOrFalse(ATTEMPT_MATCH) << endl;

		if(CONSTANT && whole_number_of_pulses > 0)
		{
			double average_p_seen = 0;		// Can be static or total
			double average_p_applied = 0;	// Static always
			double p_next;					// Static always

			for(i_fill=0; i_fill<nsamples; ++i_fill)
			{
				if(TOTAL_PRESS) // Convert to total values
					average_p_seen/*bar*/ += TotalPressureBar(pPpt, p_des_app_seen[i_fill][SEEN_P]/*bar*/, p_des_app_seen[i_fill][SEEN_T], p_des_app_seen[i_fill][SEEN_V]);
				else average_p_seen += p_des_app_seen[i_fill][SEEN_P];

				average_p_applied += p_des_app_seen[i_fill][APPLIED_P];
			}
			average_p_seen /= nsamples;
			average_p_applied /= nsamples;

			if(average_p_seen < constantp) // Static or total comparison
			{
				p_next = average_p_applied + del_p;
				if(!p_SEEN_LESS_THAN_DESIRED) del_p/=2;	// Change of sign so halve difference
				p_SEEN_LESS_THAN_DESIRED = true;
			}
			else
			{
				p_next = average_p_applied - del_p;
				if(p_SEEN_LESS_THAN_DESIRED) del_p/=2;	// Change of sign so halve difference
				p_SEEN_LESS_THAN_DESIRED = false;
			}

			for(i_fill=0; i_fill<nsamples; ++i_fill) p_des_app_seen[i_fill][NEXT_P] = p_next;

			pPpt->Out("\n");
			pPpt->Out("constantp = ");  pPpt->Out(constantp); if(TOTAL_PRESS) pPpt->Out(" (total)"); else pPpt->Out(" (static)"); pPpt->Out("\n");
			pPpt->Out("average_p_seen = "); pPpt->Out(average_p_seen); if(TOTAL_PRESS) pPpt->Out(" (total)"); else pPpt->Out(" (static)"); pPpt->Out("\n");
			pPpt->Out("average_p_applied = "); pPpt->Out(average_p_applied); pPpt->Out(" (static)");pPpt->Out("\n");
			pPpt->Out("del_p = "); pPpt->Out(del_p); pPpt->Out(" (static)"); pPpt->Out("\n");
			pPpt->Out("p_next = "); pPpt->Out(p_next); pPpt->Out(" (static)"); pPpt->Out("\n"); pPpt->Out("\n");
		}


		// Calculate averages then reset summing variables
		for(m=0; m<nmeasurements; ++m)
		{
			// Gas velocities
			if(iterations_during_pulse!=0)
				average_flow_velocity_last_pulse[m] = sum_flow_velocity_over_pulse[m]/iterations_during_pulse;
			else average_flow_velocity_last_pulse[m] = 0; // For the very first pulse
				sum_flow_velocity_over_pulse[m] = 0; // Reset for new pulse

			// Pressure wave velocities
			// uplusa
			if(iterations_during_pulse!=0)
				average_press_velocity_last_pulse_uplusa[m] = sum_press_velocity_over_pulse_uplusa[m]/iterations_during_pulse;
			else average_press_velocity_last_pulse_uplusa[m] = 0; // For the very first pulse
	
			sum_press_velocity_over_pulse_uplusa[m] = 0; // Reset for new pulse

			//uminusa
			if(iterations_during_pulse!=0)
				average_press_velocity_last_pulse_uminusa[m] = sum_press_velocity_over_pulse_uminusa[m]/iterations_during_pulse;
			else average_press_velocity_last_pulse_uminusa[m] = 0; // For the very first pulse

			sum_press_velocity_over_pulse_uminusa[m] = 0; // Reset for new pulse

			// Now calculate Strouhal numbers MSt and PMst
			st[m] = (f*pPipe[ONE_SIDE]->eff_length)/average_flow_velocity_last_pulse[m];
			beta[m] = (2*PI*f*pPipe[ONE_SIDE]->eff_length)/average_flow_velocity_last_pulse[m];
			//mst[m] = ((f)*pPipe[ONE_SIDE]->eff_length/average_flow_velocity_last_pulse[m])*(1/(2*phi));
			mst[m] = ((f)*pPipe[ONE_SIDE]->eff_length/average_flow_velocity_last_pulse[m])*(1/(2*phi_orig));
			//pmst_uplusa[m] = ((f)*pPipe[ONE_SIDE]->eff_length/average_press_velocity_last_pulse_uplusa[m])*(1/(2*phi));
			pmst_uplusa[m] = ((f)*pPipe[ONE_SIDE]->eff_length/average_press_velocity_last_pulse_uplusa[m])*(1/(2*phi_orig));
			//pmst_uminusa[m] = ((f)*pPipe[ONE_SIDE]->eff_length/average_press_velocity_last_pulse_uminusa[m])*(1/(2*phi));
			pmst_uminusa[m] = ((f)*pPipe[ONE_SIDE]->eff_length/average_press_velocity_last_pulse_uminusa[m])*(1/(2*phi_orig));
		}

		// Compare i_prev and the number of datapoints available - i_fill in any empty space
		if(i_prev < nsamples-1) // Fill between (nsamples-1) and i_prev
		{
			for(i_fill=i_prev+1; i_fill<nsamples; ++i_fill)
			{
				p_des_app_seen[i_fill][APPLIED_P] = p_des_app_seen[i_prev][APPLIED_P] + 
							((p_des_app_seen[i_prev][APPLIED_P] - p_des_app_seen[i_prev-1][APPLIED_P])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][SEEN_P] = p_des_app_seen[i_prev][SEEN_P] + 
							((p_des_app_seen[i_prev][SEEN_P] - p_des_app_seen[i_prev-1][SEEN_P])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][NEXT_P] = p_des_app_seen[i_prev][NEXT_P] + 
							((p_des_app_seen[i_prev][NEXT_P] - p_des_app_seen[i_prev-1][NEXT_P])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][THIS_P] = p_des_app_seen[i_prev][THIS_P] + 
							((p_des_app_seen[i_prev][THIS_P] - p_des_app_seen[i_prev-1][THIS_P])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference
				
				p_des_app_seen[i_fill][APPLIED_T] = p_des_app_seen[i_prev][APPLIED_T] + 
							((p_des_app_seen[i_prev][APPLIED_T] - p_des_app_seen[i_prev-1][APPLIED_T])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][SEEN_T] = p_des_app_seen[i_prev][SEEN_T] + 
							((p_des_app_seen[i_prev][SEEN_T] - p_des_app_seen[i_prev-1][SEEN_T])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][NEXT_T] = p_des_app_seen[i_prev][NEXT_T] + 
							((p_des_app_seen[i_prev][NEXT_T] - p_des_app_seen[i_prev-1][NEXT_T])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][THIS_T] = p_des_app_seen[i_prev][THIS_T] + 
							((p_des_app_seen[i_prev][THIS_T] - p_des_app_seen[i_prev-1][THIS_T])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][APPLIED_V] = p_des_app_seen[i_prev][APPLIED_V] + 
							((p_des_app_seen[i_prev][APPLIED_V] - p_des_app_seen[i_prev-1][APPLIED_V])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][SEEN_V] = p_des_app_seen[i_prev][SEEN_V] + 
							((p_des_app_seen[i_prev][SEEN_V] - p_des_app_seen[i_prev-1][SEEN_V])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][NEXT_V] = p_des_app_seen[i_prev][NEXT_V] + 
							((p_des_app_seen[i_prev][NEXT_V] - p_des_app_seen[i_prev-1][NEXT_V])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][THIS_V] = p_des_app_seen[i_prev][THIS_V] + 
							((p_des_app_seen[i_prev][THIS_V] - p_des_app_seen[i_prev-1][THIS_V])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][APPLIED_MFR] = p_des_app_seen[i_prev][APPLIED_MFR] + 
							((p_des_app_seen[i_prev][APPLIED_MFR] - p_des_app_seen[i_prev-1][APPLIED_MFR])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][SEEN_MFR] = p_des_app_seen[i_prev][SEEN_MFR] + 
							((p_des_app_seen[i_prev][SEEN_MFR] - p_des_app_seen[i_prev-1][SEEN_MFR])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][NEXT_MFR] = p_des_app_seen[i_prev][NEXT_MFR] + 
							((p_des_app_seen[i_prev][NEXT_MFR] - p_des_app_seen[i_prev-1][NEXT_MFR])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][THIS_MFR] = p_des_app_seen[i_prev][THIS_MFR] + 
							((p_des_app_seen[i_prev][THIS_MFR] - p_des_app_seen[i_prev-1][THIS_MFR])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference

				p_des_app_seen[i_fill][P_VEL] = p_des_app_seen[i_prev][P_VEL] + 
							((p_des_app_seen[i_prev][P_VEL] - p_des_app_seen[i_prev-1][P_VEL])
							/(p_des_app_seen[i_prev][TIME_IN_PERIOD] - p_des_app_seen[i_prev-1][TIME_IN_PERIOD]))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - p_des_app_seen[i_prev][TIME_IN_PERIOD]);		// Time difference			
			}
			i_prev = i_fill;  // Record up to which data point has been filled during this iteration
		}

		double currentError = 0;
		double maxError = 0;
		if(whole_number_of_pulses > 0) 
		// Only try to compare after the first pulse (the desired pulse is recorded on the first pulse)
		{
			MATCH = true;

			// Compare seen and desired values at each data point in p_des_app_seen
			// to see if they match within the specified tolerance
			i=0;
//			cout << "Testing pulse sample ";
			do{
				currentError = fabs(((p_des_app_seen[i][SEEN_P] - p_des_app_seen[i][DESIRED_P])/p_des_app_seen[i][DESIRED_P])*100);
				//if(fabs(((p_des_app_seen[i][SEEN_P] - p_des_app_seen[i][DESIRED_P])
				//			/p_des_app_seen[i][DESIRED_P])*100)>tol) MATCH = false;
				if(currentError>maxError)	maxError = currentError;
				if(currentError>tol)		MATCH = false;
				// We can compare at [i] since time_in_period is identical for each
//cout << "p_des_app_seen[" << i << "][" << SEEN_P << "] = " << p_des_app_seen[i][SEEN_P] << endl;
				++i;
			//}while(MATCH && i<nsamples);
			}while(i<nsamples);

			if(!MATCHED)
			{
				pPpt->Out(Identify());
				pPpt->Out(": PreProcess: end of pulse no. "); pPpt->Out(whole_number_of_pulses + 1);
				if(MATCH)// && !MATCHED)
				{
					//MATCHED = true; // Turn this off when using twin entry model
					pPpt->Out(" - MATCH acheived; ");
					if(SUSPEND) pPpt->STOP = true; // If desired, suspend the simulation to allow user to end the simulation every time a match is achieved
				}
				else pPpt->Out(" - no match; ");
				
				pPpt->Out("maxError = "); pPpt->Out(maxError); pPpt->Out("%\n");
			}
		}

		

		for(i=0; i<nsamples; ++i)
		{
			// Calculate actual time and write to file
			p_des_app_seen[i][TIME] = p_des_app_seen[i][TIME_IN_PERIOD] + whole_number_of_pulses*period;

			// Calculate the velocity gradient
			if(i!=0)
				p_des_app_seen[i][VEL_GRAD] = (p_des_app_seen[i][SEEN_V] - p_des_app_seen[i-1][SEEN_V])
											/(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]);
			else p_des_app_seen[i][VEL_GRAD] = (p_des_app_seen[i+1][SEEN_V] - p_des_app_seen[i][SEEN_V])
											/(p_des_app_seen[i+1][TIME_IN_PERIOD] - p_des_app_seen[i][TIME_IN_PERIOD]); 

// Or:
//			if(i!=nsamples-1) p_des_app_seen[i][VEL_GRAD] = (p_des_app_seen[i+1][VEL] - p_des_app_seen[i][VEL])
//											/(p_des_app_seen[i+1][TIME_IN_PERIOD] - p_des_app_seen[i][TIME_IN_PERIOD]);
//			else p_des_app_seen[i][VEL_GRAD] = (p_des_app_seen[i][VEL] - p_des_app_seen[i-1][VEL])
//											/(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]);
				

			// Calculate the pressure velocity gradient
			if(i!=0)
				p_des_app_seen[i][P_VEL_GRAD] = (p_des_app_seen[i][P_VEL] - p_des_app_seen[i-1][P_VEL])
											/(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]);
			else p_des_app_seen[i][P_VEL_GRAD] = (p_des_app_seen[i+1][P_VEL] - p_des_app_seen[i][P_VEL])
											/(p_des_app_seen[i+1][TIME_IN_PERIOD] - p_des_app_seen[i][TIME_IN_PERIOD]); 

// Or:
//			if(i!=nsamples-1) p_des_app_seen[i][P_VEL_GRAD] = (p_des_app_seen[i+1][P_VEL] - p_des_app_seen[i][P_VEL])
//											/(p_des_app_seen[i+1][TIME_IN_PERIOD] - p_des_app_seen[i][TIME_IN_PERIOD]);
//			else p_des_app_seen[i][P_VEL_GRAD] = (p_des_app_seen[i][P_VEL] - p_des_app_seen[i-1][P_VEL])
//											/(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]);
			

//			cout << "p_des_app_seen[i=" << i << "][SEEN_P] = " << p_des_app_seen[i][SEEN_P] << endl;
//			cout << "p_des_app_seen[i=" << i << "][APPLIED_P] = " << p_des_app_seen[i][APPLIED_P] << endl;
			if(whole_number_of_pulses >= print_from_pulse-1 && p_des_app_seen[i][TIME] >= print_from_time && i%freq==0) 
			// Only print after specified pulse and/or time and at the specified sampling frequency
			{	
				// Print column headers to file once only
				if (!PRINT_HEADERS) {
					PRINT_HEADERS = true;
					fprintf(FILE_p_des_app_seen, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
						"TIME",
						"TIME_IN_PERIOD",

						"DESIRED_P",
						"APPLIED_P",
						"SEEN_P",
						"NEXT_P",
						"THIS_P",

						"DESIRED_T",
						"DESIRED_T0",
						"APPLIED_T",
						"NEXT_T",
						"SEEN_T",
						"THIS_T",

						"DESIRED_V",
						"APPLIED_V",
						"SEEN_V",
						"NEXT_V",
						"THIS_V",

						"DESIRED_MFR",
						"APPLIED_MFR",
						"SEEN_MFR",
						"NEXT_MFR",
						"THIS_MFR",

						"VEL_GRAD",
						"P_VEL",
						"P_VEL_GRAD",
						"UNDER_EST",
						"UNDER_EST_PREV",
						"UNDER_EST_FACT",
						"K_LOSS");
					fprintf(FILE_p_des_app_seen, "\n");
				}

				// Print pulse data to file
				fprintf(FILE_p_des_app_seen,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", 
				p_des_app_seen[i][TIME], 
				p_des_app_seen[i][TIME_IN_PERIOD],

				p_des_app_seen[i][DESIRED_P], 
				p_des_app_seen[i][APPLIED_P], 
				p_des_app_seen[i][SEEN_P], 
				p_des_app_seen[i][NEXT_P], 
				p_des_app_seen[i][THIS_P],
				
				p_des_app_seen[i][DESIRED_T], 
				p_des_app_seen[i][DESIRED_T0], 
				p_des_app_seen[i][APPLIED_T], 
				p_des_app_seen[i][NEXT_T], 
				p_des_app_seen[i][SEEN_T], 
				p_des_app_seen[i][THIS_T],
				
				p_des_app_seen[i][DESIRED_V], 
				p_des_app_seen[i][APPLIED_V], 
				p_des_app_seen[i][SEEN_V], 
				p_des_app_seen[i][NEXT_V], 
				p_des_app_seen[i][THIS_V],
				
				p_des_app_seen[i][DESIRED_MFR], 
				p_des_app_seen[i][APPLIED_MFR], 
				p_des_app_seen[i][SEEN_MFR], 
				p_des_app_seen[i][NEXT_MFR], 
				p_des_app_seen[i][THIS_MFR],
				
				p_des_app_seen[i][VEL_GRAD],
				p_des_app_seen[i][P_VEL], 
				p_des_app_seen[i][P_VEL_GRAD],
				p_des_app_seen[i][UNDER_EST], 
				p_des_app_seen[i][UNDER_EST_PREV], 
				p_des_app_seen[i][UNDER_EST_FACT], 	
				p_des_app_seen[i][K_LOSS]);
				fprintf(FILE_p_des_app_seen,"\n");
			}
			
			// Copy NEXT_P data to THIS_P (must be done after printing to file) if no match, and only if we want to attempt matching this pulse
			if(!MATCHED && ATTEMPT_MATCH)
			{
				p_des_app_seen[i][THIS_P] = p_des_app_seen[i][NEXT_P];
				p_des_app_seen[i][THIS_T] = p_des_app_seen[i][NEXT_T];
				p_des_app_seen[i][THIS_V] = p_des_app_seen[i][NEXT_V];
				p_des_app_seen[i][THIS_MFR] = p_des_app_seen[i][NEXT_MFR];
			}
			
/*
			// Copy UNDER_EST data to UNDER_EST_PREV
			p_des_app_seen[i][UNDER_EST_PREV] = p_des_app_seen[i][UNDER_EST];
*/
		}
		iterations_during_pulse = 0; // Reset for new pulse
//		time_in_period_prev = 0; // Must be set to 0 else would be > time_in_period
		i_prev = -1; // No data points have been filled at the start of the pulse
	}
	else START_OF_NEW_PULSE = false;

	// Pulse calculations
	whole_number_of_pulses = int(integer_part);
	++iterations_during_pulse;
	
	// Calculate pressure to apply based on waveform parameters
	p_wave_prev = p_wave; // Save old value
	T_wave_prev = T_wave; // Save old value
	v_wave_prev = v_wave; // Save old value
	mfr_wave_prev = mfr_wave; // Save old value

	// Set the current velocity gradient
//	if() 
	{
		// Find the relevant data point in p_des_app_seen for the current time_in_period
		i=0;
		do{++i;}
		while(p_des_app_seen[i][TIME_IN_PERIOD]<time_in_period && i<nsamples-1);
	
		// Interpolate for velocity gradient at time_in_period between [i-1] and [i]
		vel_grad = p_des_app_seen[i-1][VEL_GRAD]
						+	((p_des_app_seen[i][VEL_GRAD] - p_des_app_seen[i-1][VEL_GRAD])/
							(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
							(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

		//cout << "vel_grad = " << vel_grad << endl;
	}

	if(INTELLIGENT && whole_number_of_pulses > 0) 
	// Only start converging onto desired pulse after the desired pulse has been recorded!
	{
		// Find the relevant data point in p_des_app_seen for the current time_in_period
		i=0;
		do{++i;}
		while(p_des_app_seen[i][TIME_IN_PERIOD]<time_in_period && i<nsamples-1);

		// Interpolate for pressure at time_in_period between [i-1] and [i]
		p_wave = p_des_app_seen[i-1][THIS_P]
						+	((p_des_app_seen[i][THIS_P] - p_des_app_seen[i-1][THIS_P])/
							(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
							(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

		T_wave = p_des_app_seen[i-1][THIS_T]
						+	((p_des_app_seen[i][THIS_T] - p_des_app_seen[i-1][THIS_T])/
							(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
							(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

		v_wave = p_des_app_seen[i-1][THIS_V]
						+	((p_des_app_seen[i][THIS_V] - p_des_app_seen[i-1][THIS_V])/
							(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
							(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

		mfr_wave = p_des_app_seen[i-1][THIS_MFR]
						+	((p_des_app_seen[i][THIS_MFR] - p_des_app_seen[i-1][THIS_MFR])/
							(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
							(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

		p_wave_desired = p_des_app_seen[i-1][DESIRED_P]
						+	((p_des_app_seen[i][DESIRED_P] - p_des_app_seen[i-1][DESIRED_P])/
							(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
							(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

		T_wave_desired = p_des_app_seen[i-1][DESIRED_T]
						+	((p_des_app_seen[i][DESIRED_T] - p_des_app_seen[i-1][DESIRED_T])/
							(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
							(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);
		
		v_wave_desired = p_des_app_seen[i-1][DESIRED_V]
						+	((p_des_app_seen[i][DESIRED_V] - p_des_app_seen[i-1][DESIRED_V])/
							(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
							(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

		mfr_wave_desired = p_des_app_seen[i-1][DESIRED_MFR]
						+	((p_des_app_seen[i][DESIRED_MFR] - p_des_app_seen[i-1][DESIRED_MFR])/
							(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
							(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);
	}
	else // Apply a set waveform
	{
		// Interpolate from desired array
		i=0;
		do ++i;
		while(p_des_app_seen[i][TIME_IN_PERIOD]<time_in_period && i<nsamples-1);

		p_wave = p_des_app_seen[i-1][DESIRED_P] +
			 ((p_des_app_seen[i][DESIRED_P] - p_des_app_seen[i-1][DESIRED_P])/
					(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
					(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

		p_wave_desired = p_wave; // They are the same in this case

		T_wave = p_des_app_seen[i-1][DESIRED_T] +
			 ((p_des_app_seen[i][DESIRED_T] - p_des_app_seen[i-1][DESIRED_T])/
					(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
					(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

		T_wave_desired = T_wave; // They are the same in this case

		v_wave = p_des_app_seen[i-1][DESIRED_V] +
			 ((p_des_app_seen[i][DESIRED_V] - p_des_app_seen[i-1][DESIRED_V])/
					(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
					(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

		v_wave_desired = v_wave; // They are the same in this case

		mfr_wave = p_des_app_seen[i-1][DESIRED_MFR] +
			 ((p_des_app_seen[i][DESIRED_MFR] - p_des_app_seen[i-1][DESIRED_MFR])/
					(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
					(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

		mfr_wave_desired = mfr_wave; // They are the same in this case
	}
/*
cout << "p_wave = " << p_wave << endl;
cout << "T_wave = " << T_wave << endl;
cout << "v_wave = " << v_wave << endl;
cout << "mfr_wave = " << mfr_wave << endl;
//*/

/*
	// Have now calculated a static pressure to apply
	// Use this with the constant T value to get rho, then derive a velocity to apply
	// based on the mass flow profile.
	// p/rho = RT
	double rho = (p_wave*1e5)/(pPpt->R_air*T_WAVEFORM);
	double mdot=0;
	// mdot = rho*area*vel

	// Interpolate mdot from desired array
	i=0;
	do ++i;
	while(p_des_app_seen[i][TIME_IN_PERIOD]<time_in_period && i<nsamples-1);

	mdot = p_des_app_seen[i-1][DESIRED_MFR] +
		 ((p_des_app_seen[i][DESIRED_MFR] - p_des_app_seen[i-1][DESIRED_MFR])/
				(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
				(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

	u_WAVEFORM = mdot/(rho*this->pBN->f);

	if(u_WAVEFORM>150) u_WAVEFORM = 150;
*/

///*

	if(CONSTANT)
	{
		u_WAVEFORM = pBN[ONE_SIDE]->U*pPipe[ONE_SIDE]->AREF; // Lets u float as a result of the simulation. Otherwise problem is over-specified.
		
		if(TOTAL_PRESS) T_WAVEFORM = StaticTemperature(pPpt, constantT, u_WAVEFORM);		
		else T_WAVEFORM = constantT;
	}
	else
	{
		if(shape==2) // Cosine
		{
			//u_WAVEFORM = v;// m/s, read from file
			u_WAVEFORM = pBN[ONE_SIDE]->U*pPipe[ONE_SIDE]->AREF; // Lets u float as a result of the simulation. Otherwise problem is over-specified.
			T_WAVEFORM = Ts;
		}
		else
		{
			if(USE_FILE_Ts)
			{
				// Find position in pulse				
				i=0;
				do ++i;
				while(p_des_app_seen[i][TIME_IN_PERIOD]<time_in_period && i<nsamples-1);
			
				// Interpolate T0m and Tsm from their desired arrays
				double T0m_temp = p_des_app_seen[i-1][DESIRED_T0] +
					 ((p_des_app_seen[i][DESIRED_T0] - p_des_app_seen[i-1][DESIRED_T0])/
							(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
							(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

				double Tsm_temp = p_des_app_seen[i-1][DESIRED_T] +
					 ((p_des_app_seen[i][DESIRED_T] - p_des_app_seen[i-1][DESIRED_T])/
							(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
							(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

				T_wave = Tsm_temp;
				

				double C1_temp = sqrt( ((2*pPpt->gammaAir(Tsm_temp)*pPpt->R_air)/(pPpt->gammaAir(Tsm_temp)-1))*(T0m_temp - Tsm_temp) );
				
				//double v_des_temp = p_des_app_seen[i-1][DESIRED_V] +
				//	 ((p_des_app_seen[i][DESIRED_V] - p_des_app_seen[i-1][DESIRED_V])/
				//			(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
				//			(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

				//v_wave = C1_temp;
				//v_wave = v_des_temp;

				// Derive the velocity to use form the desired mass flow rate
				double tempRho, tempF;
				//tempRho = (p_wave*1e5)/(pPpt->R_air*T_wave);
				tempRho = (pBN[ONE_SIDE]->p_dash*pPpt->PREF*1e5)/(pPpt->R_air*pBN[ONE_SIDE]->T);
				tempF = pBN[ONE_SIDE]->f;
				v_wave = mfr_wave/(tempRho*tempF);
				//v_wave = pBN[ONE_SIDE]->U*pPipe[ONE_SIDE]->AREF; // Lets u float as a result of the simulation. Otherwise problem is over-specified.


				//mfr_wave = ((p_wave*1e5)/(pPpt->R_air*T_wave))*v_wave*pBN[ONE_SIDE]->f;

				p_wave = ((mfr_wave/(v_wave*pBN[ONE_SIDE]->f))*(pPpt->R_air*T_wave))/1e5;
			}
			else // Use constant values from file
			{
				v_wave = v;
				T_wave = Ts;
				mfr_wave = ((p_wave*1e5)/(pPpt->R_air*T_wave))*v_wave*pBN[ONE_SIDE]->f;
			}
			
			p_wave = p_wave;
			T_WAVEFORM = T_wave;
			u_WAVEFORM = v_wave;
			//mfr_wave = ((p_wave*1e5)/(pPpt->R_air*T_wave))*v_wave*pBN[ONE_SIDE]->f;
		}
	}

/*
cout << "p_wave = " << p_wave << endl;
cout << "T_WAVEFORM = " << T_WAVEFORM << endl;
cout << "u_WAVEFORM = " << u_WAVEFORM << endl;
cout << "mfr_wave = " << mfr_wave << endl;
//*/

//cout << "PreProcess: T_WAVEFORM = " << T_WAVEFORM << endl;

//	if(fabs(v_wave)<100) u_WAVEFORM = v_wave;
//	else u_WAVEFORM = 100*(fabs(v_wave)/v_wave);

//*/

//cout << "Calculated rho = " << rho << endl;
//cout << "Interpolated desired mdot = " << mdot << endl;
//cout << "Corresponding u_WAVEFORM = " << u_WAVEFORM << endl;
	
//	u_WAVEFORM = v;// m/s, read from file
//	T_WAVEFORM = Ts;


/*
	if(CALIBRATE)
	{
		if(STEADY) // Record operating point and increment p_wave
		{

			PR[0]
			MFP[0]


		}
	}
*/


/*	
	if(CALIBRATE)
	{
		char pause;
		if(
			//pPipe[ONE_SIDE]->Steady(pPpt) // If the main pipe has steady flow
			true
			//STEADY // If all the pipes have steady flow
			)
		{
			double step = double(match_counter)/double(cal_points);
			step += (cal_PR_high - cal_PR_low)/100;
			//			if(match_counter==1) cout << "step = " << step << endl;
			
			double p_cal_des = cal_PR_low + step*(cal_PR_high - cal_PR_low); // Assumes p2 is 1 bar
			//cout << "double(match_counter/cal_points) = " << double(match_counter/cal_points) << endl;
			//cin >> pause;
			double p_cal_seen = pBN[ONE_SIDE]->p_dash*pPpt->PREF;
			//if(!cal_MATCH) p_wave = p_wave_prev + (p_cal_des - p_cal_seen); // Original

			if(!cal_MATCH)
			{
				if(CHANGE_PR)
				{
					if(INCREASE_PR)
					{
						if(SWITCH)
						{
							//cout << "SWITCH" << endl;
							p_wave = p_wave_prev_prev + del_change_PR;
							cout << "Transmissive:" << endl;
							cout << "Increasing p_wave from " << p_wave_prev_prev << " to " << p_wave << endl;
							cout << endl;
						}
						else
						{
							p_wave = p_wave_prev + del_change_PR;
							cout << "Transmissive:" << endl;
							cout << "Increasing p_wave from " << p_wave_prev << " to " << p_wave << endl;
							cout << endl;
						}
					}
					else
					{
						if(SWITCH)
						{
							//cout << "SWITCH" << endl;
							p_wave = p_wave_prev_prev - del_change_PR;
							if(p_wave < 1) p_wave = 1; // Forward flow only
							cout << "Transmissive:" << endl;
							cout << "Decreasing p_wave from " << p_wave_prev_prev << " to " << p_wave << endl;
							cout << endl;
						}
						else
						{
							p_wave = p_wave_prev - del_change_PR;
							if(p_wave < 1) p_wave = 1; // Forward flow only
							cout << "Transmissive:" << endl;
							cout << "Decreasing p_wave from " << p_wave_prev << " to " << p_wave << endl;
							cout << endl;
						}
					}
					CHANGE_PR = false; // Reset needed here
				}
				else
				{
					p_wave = p_wave_prev;
					//p_wave = p_wave_prev + 0.1*(p_cal_des - p_cal_seen); // Relaxed
					//cout << "Attempting to attain required PR" << endl;	
				}
			}
			else
			{
				p_wave = p_wave_prev;
				cout << "Keeping p_wave the same" << endl;	
			}
			
//			cout << "Transm: STEADY!" << endl;
//			cout << "p_cal_des = " << p_cal_des << endl;
//			cout << "p_cal_seen = " << p_cal_seen << endl;
//			cout << "p_wave_prev = " << p_wave_prev << endl;
//			cout << "p_wave = " << p_wave << endl;
//			cin >> pause;

			// Original 1e-4
			if(fabs(p_cal_des - p_cal_seen) < 1e-2 && !cal_MATCH)
			{
				if(WAIT_COUNTER >= 100) cal_MATCH = true;  // WAIT_COUNTER to allow flow to settle
				else ++WAIT_COUNTER;
				
				cout << endl;
				cout << "Transm: MATCH!" << endl;
				cout << "p_cal_des = " << p_cal_des << endl;
				cout << "p_cal_seen = " << p_cal_seen << endl;
				cout << "p_wave = " << p_wave << endl;
				cin >> pause;		
			}
		}
		else p_wave = p_wave_prev;
		
		u_WAVEFORM = v;// m/s, read from file
		T_WAVEFORM = Ts;
	}
//*/

//	cout << "Transm: p_wave = " << p_wave << endl;
/*
	if(CALIBRATE)
	{
		char pause;
		if(pPipe[ONE_SIDE]->Steady()) // If the pipe has steady flow
		{
			cal_MATCH = true;
		}
	}
*/


	U_WAVEFORM = u_WAVEFORM/pPipe[ONE_SIDE]->AREF;


	// Write the applied p_wave to the p_des_app_seen[][APPLIED_P] array

	// Find the relevant data point in p_des_app_seen for the current time_in_period
	i=0;
	do{++i;}
	while(p_des_app_seen[i][TIME_IN_PERIOD]<time_in_period && i<nsamples-1);
	
  if(START_OF_SIMULATION) i_prev = -1;
//	cout << "time_in_period = " << time_in_period << endl;
//	cout << "i = " << i << endl;
//	cout << "i_prev = " << i_prev << endl;

	// Extrapolate for pressure at time[i] and fill in for any gaps in p_des_app_seen
	for(i_fill=i; i_fill>i_prev; --i_fill)
	{
//cout << "time_in_period = " << time_in_period << endl;
//cout << "p_des_app_seen[i_fill=" << i_fill << "][APPLIED_P] = " << p_des_app_seen[i_fill][APPLIED_P] << endl;
//cout << "HERE" << endl;
		
		p_des_app_seen[i_fill][APPLIED_P] = p_wave +
							((p_wave - p_wave_prev)/(time_in_period - time_in_period_prev))		// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - time_in_period);			// Time difference

		p_des_app_seen[i_fill][APPLIED_T] = T_wave +
							((T_wave - T_wave_prev)/(time_in_period - time_in_period_prev))		// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - time_in_period);			// Time difference

		p_des_app_seen[i_fill][APPLIED_V] = v_wave +
							((v_wave - v_wave_prev)/(time_in_period - time_in_period_prev))		// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - time_in_period);			// Time difference

		p_des_app_seen[i_fill][APPLIED_MFR] = mfr_wave +
							((mfr_wave - mfr_wave_prev)/(time_in_period - time_in_period_prev))	// Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - time_in_period);			// Time difference
	}
  /*
  for(i_fill=i; i_fill>-1; --i_fill)
	{
cout << "p_des_app_seen[i_fill=" << i_fill << "][APPLIED_P] = " << p_des_app_seen[i_fill][APPLIED_P] << endl;
  }
//*/

  //i_prev = i;  // Record up to which data point has been filled during this iteration

  if(START_OF_SIMULATION) START_OF_SIMULATION = false;
}

void CTransmissive::Transmissive(CTime &rMyTime, double DELZe, double DELZi, CProperties* pPpt, double ca, double del_ca, double TIMEe, double TIMEi, int timestep, bool STEADY)
//--------------------------------------------------//
// Non-homentropic Transmissive Boundary			//
// -------------------------------------			//
// Applies a desired pressure waveform				//
// at a non-reflecting transmissive boundary.		//
//--------------------------------------------------//
{
	// Start of normal boundary condition calculations
	// ===============================================
	//if(CALIBRATE) p_wave = p_wave_calibrate;

	bool converged;
	int counter, TEMP;
	double lambda_in_c2, AA_c2, lambda_out1, lambda_out2;
	double CLIN_WAVEFORM, AA_WAVEFORM, A_WAVEFORM;//, U_WAVEFORM;

	// Calculate AA and lambda_in_star_n for the waveform being applied
	double* lambda_in_star_n;
	lambda_in_star_n = new double [2];
	double* lambda_in_star_test;
	lambda_in_star_test = new double [2];

	int *pipe_flow; pipe_flow = new int [NPIPES];
	A_WAVEFORM = sqrt(pPpt->gammaAir(T_WAVEFORM)*pPpt->R_air*T_WAVEFORM)/pPipe[ONE_SIDE]->AREF;
	AA_WAVEFORM = A_WAVEFORM/pow(p_wave/pPpt->PREF, (pPpt->gammaAir(T_WAVEFORM)-1)/(2*pPpt->gammaAir(T_WAVEFORM)));
//cout << "p_wave = " << p_wave << endl;
//cout << "T_WAVEFORM = " << T_WAVEFORM << endl;
//cout << "A_WAVEFORM = " << A_WAVEFORM << endl;
//cout << "AA_WAVEFORM = " << AA_WAVEFORM << endl;
	if(end[ONE_SIDE] == ODD)
	{
		INTERNAL = RIGHT_SIDE;
		WAVEFORM = LEFT_SIDE;

		CLIN_WAVEFORM = A_WAVEFORM + ((pPpt->gammaAir(T_WAVEFORM)-1)/2)*U_WAVEFORM;
	}
	else
	{
		INTERNAL = LEFT_SIDE;
		WAVEFORM = RIGHT_SIDE;

		CLIN_WAVEFORM = A_WAVEFORM - ((pPpt->gammaAir(T_WAVEFORM)-1)/2)*U_WAVEFORM;
	}

	lambda_in_star_n[INTERNAL] = (*(pCLIN[ONE_SIDE]))[R+1] / pBN[ONE_SIDE]->AA[1];
//	lambda_in_star_n[INTERNAL] = (*(pCLIN[ONE_SIDE]))[R+1] / 1;

//	AA_WAVEFORM = 1;
	lambda_in_star_n[WAVEFORM] = CLIN_WAVEFORM / AA_WAVEFORM;

//	cout << "time\t=\t" << time << "\t(*(pCLIN[ONE_SIDE]))[R+1]\t=\t" << (*(pCLIN[ONE_SIDE]))[R+1] << endl;
//	cout << "time\t=\t" << time << "\tpBN[ONE_SIDE]->AA[1]\t=\t" << pBN[ONE_SIDE]->AA[1] << endl;

//	cout << "time\t=\t" << time << "\tlambda_in_star_n[WAVEFORM]\t=\t" << lambda_in_star_n[WAVEFORM] << endl;
//	cout << "time\t=\t" << time << "\tlambda_in_star_n[INTERNAL]\t=\t" << lambda_in_star_n[INTERNAL] << endl;
//	cout << endl;

	// Direction test
	// ==============
	if(lambda_in_star_n[WAVEFORM] - lambda_in_star_n[INTERNAL] > 1e-12)
	{
		// Flow from waveform to internal; waveform is upstream, internal is downstream
//		cout << "Flow is from waveform to internal\n";
		UPSTREAM = WAVEFORM;
		DOWNSTREAM = INTERNAL;
		if(!WAVEFORM_TO_INTERNAL)
		{
//			cout << "Flow through transmissive boundary is changing direction, to transmissive-to-internal\n";
//			cin >> pause;
		}
		WAVEFORM_TO_INTERNAL = true;
	}
	else
	{
		if(lambda_in_star_n[INTERNAL] - lambda_in_star_n[WAVEFORM] > 1e-12)
		{
			// Flow from internal to waveform; internal is upstream, waveform is downstream
//			cout << "Flow is from internal to waveform\n";
			UPSTREAM = INTERNAL;
			DOWNSTREAM = WAVEFORM;
			if(WAVEFORM_TO_INTERNAL)
			{
//				cout << "Flow through transmissive boundary is changing direction, to internal-to-transmissive\n";
//				cin >> pause;
			}
			WAVEFORM_TO_INTERNAL = false;
		}
		else
		{
///			cout << "No flow through transmissive boundary\n";
//			cin >> pause;

			// No flow - treat as closed end
			pipe_flow[ONE_SIDE] = NOFLOW;
		
			(*(pCLOUT))[R+1] = (*(pCLIN))[R+1];
		
			// No need to update lambda_in or AA
			return;
		}
	}
	
	// Main algorithm
	// ==============
	converged = false;
	counter = 0;
	do
	{
//if(timestep>=3369) cout << "top Timestep = " << timestep << endl;	
		Loss(pPpt, lambda_in_c2, AA_c2, lambda_out1, lambda_out2, CLIN_WAVEFORM, AA_WAVEFORM, pipe_flow[ONE_SIDE], timestep);
//if(timestep>=3369) cout << "bottom Timestep = " << timestep << endl;	
		// Test to see if new values cause change in flow direction
		lambda_in_star_test[UPSTREAM] = lambda_in_star_n[UPSTREAM];
		lambda_in_star_test[DOWNSTREAM] = lambda_in_c2/AA_c2;

///*
		if(lambda_in_star_test[UPSTREAM]>lambda_in_star_test[DOWNSTREAM])
		{
			// Then we still maintain flow from UPSTREAM to DOWNSTREAM
			converged = true;
		}
		else // Flow has switched direction with the new values calculated
		{
			cout << "Transmissive: flow has switched direction within algorithm\n";
			//cout << "I haven't tested this code yet\n";
			//cin >> pause;
			TEMP = UPSTREAM;
			UPSTREAM = DOWNSTREAM;
			DOWNSTREAM = TEMP;
			// Then loop back and re-run main algorithm
		}
		++ counter;
		if(counter>10)
		{
			cout << "counter = " << counter << endl;
			converged = true;
		}
//*/
//		converged = true;
	}while(!converged);

//////////////////////////////////////////
	double *lambda_in_c, *lambda_out_c, *AA_c;
	lambda_in_c = new double [NPIPES]; lambda_out_c = new double [NPIPES]; AA_c = new double [NPIPES];	
		
	lambda_in_c[ONE_SIDE] = lambda_in_c2;
	AA_c[ONE_SIDE] = AA_c2;
	
	if(INTERNAL==UPSTREAM) // OUTFLOW from pipe
	{
		lambda_out_c[ONE_SIDE] = lambda_out1;
	}
	else // INFLOW into pipe
	{
		lambda_out_c[ONE_SIDE] = lambda_out2;
	}

	bool* CHOKED; CHOKED = new bool [NPIPES]; CHOKED[ONE_SIDE] = false;

	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out_c, AA_c, pipe_flow, CHOKED);
//////////////////////////////////////////
/*
	// Update for inflow and outflow
	// =============================
//	(*(pCLOUT[UPSTREAM]))[1] = lambda_out1;	// Always update lambda_out
//	(*(pCLOUT[DOWNSTREAM]))[1] = lambda_out2;	// Always update lambda_out

//	cout << "lambda_out1 = " << lambda_out1 << endl;
//	cout << "lambda_out2 = " << lambda_out2 << endl;

	if(INTERNAL==UPSTREAM) (*(pCLOUT[ONE_SIDE]))[1] = lambda_out1;	// Always update lambda_out
	else (*(pCLOUT[ONE_SIDE]))[1] = lambda_out2;	// Always update lambda_out

	for(int p=0; p<2; ++p)
	{
//		pipe_flow_old[p] = *pend_flow[p];	// Record previous flow directions
//		*pend_flow[p] = pipe_flow[p];

		pipe_flow_old = *(pend_flow[ONE_SIDE]);	// Record previous flow directions
		*(pend_flow[ONE_SIDE]) = pipe_flow;

//		if(*pend_flow[p]==INFLOW) // Only INFLOW creates pathlines
		if(*(pend_flow[ONE_SIDE])==INFLOW) // Only INFLOW creates pathlines
		{
//			(*(pCLIN[p]))[1] = lambda_in_c2; // Only need to update on INFLOW, because OUTFLOWs will maintain the uncorrected value
//			(*(pPathLine[p])).AA = AA_c2; // Only INFLOW will have a new pathline at the end of the pipe ready to go

			(*(pCLIN[ONE_SIDE]))[1] = lambda_in_c2; // Only need to update on INFLOW, because OUTFLOWs will maintain the uncorrected value
			(*(pPathLine[ONE_SIDE])).AA = AA_c2; // Only INFLOW will have a new pathline at the end of the pipe ready to go

//			pBN[p]->AA[1] = AA_c2; // This seems to be necessary to get the correct pipe pressure
			pBN[ONE_SIDE]->AA[1] = AA_c2; // This seems to be necessary to get the correct pipe pressure

			if(pipe_flow_old != *(pend_flow[ONE_SIDE]))
			// If old flow direction was not INFLOW then there was no new path line created
			// so adjust existing one by setting its XK value to the appropriate end
			{
				//cout << "(*(pPathLine[p])).XK = " << (*(pPathLine[p])).XK << endl;
				//cout << "this->end[p] = " << this->end[p] << endl; 
//				if(this->end[p]==ODD){(*(pPathLine[p])).XK = 0; (*(pPathLine[p])).XK_old = (*(pPathLine[p])).XK;}
//				else{(*(pPathLine[p])).XK = pPipe[p]->XPIPE; (*(pPathLine[p])).XK_old = (*(pPathLine[p])).XK;}

				if(end[ONE_SIDE]==ODD){(*(pPathLine[ONE_SIDE])).XK = 0; (*(pPathLine[ONE_SIDE])).XK_old = (*(pPathLine[ONE_SIDE])).XK;}
				else{(*(pPathLine[ONE_SIDE])).XK = pPipe[ONE_SIDE]->XPIPE; (*(pPathLine[ONE_SIDE])).XK_old = (*(pPathLine[ONE_SIDE])).XK;}
			}
			else 
			// rpipe_flow_old == *pend_flow
			// Check that a new path line HAS been created correctly if INFLOW
			{
//				if(*(pend_flow[p])==INFLOW && ((end[p]==ODD && (*(pPathLine[p])).XK!=0) || (end[p]==EVEN && (*(pPathLine[p])).XK!=pPipe[p]->XPIPE)))
				if(*(pend_flow[ONE_SIDE])==INFLOW && ((end[ONE_SIDE]==ODD && (*(pPathLine[ONE_SIDE])).XK!=0) || (end[ONE_SIDE]==EVEN && (*(pPathLine[ONE_SIDE])).XK!=pPipe[ONE_SIDE]->XPIPE)))
				{
					cout << "Transmissive: should've created new path line but its XK is not right!!" << endl;
//					if(this->end[p]==ODD)
					if(this->end[ONE_SIDE]==ODD)
//					cout << "(*(pPathLine[p])).XK = " << (*(pPathLine[p])).XK << ", but it should be 0" << endl;
					cout << "(*(pPathLine[ONE_SIDE])).XK = " << (*(pPathLine[ONE_SIDE])).XK << ", but it should be 0" << endl;
//					else cout << "(*(pPathLine[p])).XK = " << (*(pPathLine[p])).XK << ", but it should be " << pPipe[p]->XPIPE << endl;
					else cout << "(*(pPathLine[ONE_SIDE])).XK = " << (*(pPathLine[ONE_SIDE])).XK << ", but it should be " << pPipe[ONE_SIDE]->XPIPE << endl;
				}
			}
		}
	}
*/
	return;
}

void CTransmissive::Loss(CProperties* pPpt, double &rlambda_in_c2, double &rAA_c2, double &rlambda_out1, double &rlambda_out2, double &rCLIN_WAVEFORM, double &rAA_WAVEFORM, int &rpipe_flow, int timestep)
//--------------------------------------------------//
// Non-isentropic 				//
// ---------------------------------				//
// Developed from Benson's  model.//
//--------------------------------------------------//
{
	double lambda_in_n[2];
	double AA_n[2];
	double lambda_in[2];
	double lambda_out[2];
	double A_star[2];
	double U_star[2];
	//double M[2]; // Set as member variable instead

	double AA_1, AA2_over_AA1, AA_c2, lambda_in_star_c2, lambda_in_c2;
	double lambda_in_d1;
	double limit;
	double del_M2, temp_M2;
	double area_ratio, phi, a, b, alpha, K;

	bool M2_converged, M2_positive, GREATER;
	int counter;
	
	M2_converged = false;
	limit = 1e-6;

	// Set flow directions
	// ===================
	if(INTERNAL==UPSTREAM) rpipe_flow = OUTFLOW;
	else rpipe_flow = INFLOW;

	// Calculate area ratio parameter for this flow direction (<1)
	// ===========================================================
//	area_ratio = (pBN[UPSTREAM]->f_dash*pPpt->fref)/(pBN[DOWNSTREAM]->f_dash*pPpt->fref);
	area_ratio = 1;

	// Enter other system parameters
	// ================================
	phi = 1;

	// Enter inital characteristic values
	// ==================================
	if(INTERNAL==UPSTREAM)
	{
		lambda_in_n[UPSTREAM] = (*(pCLIN[ONE_SIDE]))[1];
		lambda_in_n[DOWNSTREAM] = rCLIN_WAVEFORM;
		AA_n[UPSTREAM] = pBN[ONE_SIDE]->AA[1];
//////
//		AA_n[UPSTREAM] = 1;

		AA_n[DOWNSTREAM] = rAA_WAVEFORM;
	}
	else
	{
		lambda_in_n[DOWNSTREAM] = (*(pCLIN[ONE_SIDE]))[1];
		lambda_in_n[UPSTREAM] = rCLIN_WAVEFORM;
		AA_n[DOWNSTREAM] = pBN[ONE_SIDE]->AA[1];
//////
//		AA_n[DOWNSTREAM] = 1;
		AA_n[UPSTREAM] = rAA_WAVEFORM;
	}

	lambda_in[UPSTREAM] = lambda_in_n[UPSTREAM];
	AA_1 = AA_n[UPSTREAM]; // Since no correction needed on this side

	// Second part - loop to converge on T2
	// ====================================
	GREATER = false;
	//M2_max = 1;
	//M[DOWNSTREAM] = M2_max/2;
	//del_M2 = M2_max/4.5;
	//M[DOWNSTREAM] = 0.0001;
	//M[DOWNSTREAM] = 0.00001;
	M[DOWNSTREAM] = 1.0;
	//del_M2 = 1/4.5;
	del_M2 = 0.1;

	counter = 0;
	double error, error_prev;
	error_prev = 1e6;
	error = 1e6;
	lambda_in_c2 = lambda_in_n[DOWNSTREAM];	// An initial estimate
	double minError = 1e6; // A large value
	double sign = 1.0;

//if(timestep>=3369) cout << "top Timestep = " << timestep << endl;
	do
	{
//if(timestep>=3369) cout << "counter = " << counter << endl;	
		++counter;
/*
		if(DEVICE==GAUZE)
		{
			PARAM_R = pBN[UPSTREAM]->Re; // Upstream Re no. upon entering this method, so not as accurate
//			cout << "Re = " << PARAM_R << endl;
		}
		else
		{
			if(DEVICE==THROTTLE || DEVICE==EGR) PARAM_R = M[DOWNSTREAM]; // Uses the inside loop value!
			else PARAM_R = 1; // Some arbitrary value
		}
*/
//		K = Interpolate(pPpt, DEVICE, PARAM_S, PARAM_R);
//		if(DEVICE==THROTTLE) K = X*K;	// Dynamic factor allows for unsteady flow
		K = 0; // Transmissive boundary, no loss.

/*
		// Find the relevant data point in p_des_app_seen for the current time_in_period
		int i=0;
		do{++i;}
		while(p_des_app_seen[i][TIME_IN_PERIOD]<time_in_period && i<nsamples-1);
	
		// Interpolate for K_LOSS at time_in_period between [i-1] and [i]
		this->mfr_under_est // Really K_LOSS
		= p_des_app_seen[i-1][K_LOSS]
						+	((p_des_app_seen[i][K_LOSS] - p_des_app_seen[i-1][K_LOSS])/
							(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
							(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

//		K = this->mfr_under_est;
//cout << "K = " << K << endl;
*/

		a = 2/(pPpt->gammaAir(pBN[ONE_SIDE]->T) - 1);
		b = a*pow(M[DOWNSTREAM],2) + pow(M[DOWNSTREAM],4);
		alpha = pow(area_ratio,2);

		M[UPSTREAM] = sqrt( ( (2*b*K + alpha*a) - sqrt( pow((2*b*K + alpha*a),2) - 4*b*(b*pow(K,2) - alpha) ) )
								/(2*(b*pow(K,2) - alpha)) ); // Eq. 8.168

		AA2_over_AA1 = pow( (1/(1 - K*pow(M[UPSTREAM],2))), (pPpt->gammaAir(pBN[ONE_SIDE]->T)-1)/(2*pPpt->gammaAir(pBN[ONE_SIDE]->T)) )
					   * sqrt( ( (2/(pPpt->gammaAir(pBN[ONE_SIDE]->T)-1)) + pow(M[UPSTREAM],2) )/( (2/(pPpt->gammaAir(pBN[ONE_SIDE]->T)-1)) + pow(M[DOWNSTREAM],2) ) ); // Eq. 8.172

		AA_c2 = (AA2_over_AA1)*AA_1;	// Eq. 8.173
		lambda_in_star_c2 = lambda_in_c2/AA_c2;	// Eq. 8.174

		A_star[UPSTREAM] = pow( (1/(1 - K*pow(M[UPSTREAM],2))), (pPpt->gammaAir(pBN[ONE_SIDE]->T)-1)/(2*pPpt->gammaAir(pBN[ONE_SIDE]->T)) )
						   * (lambda_in_star_c2/(1 - ((pPpt->gammaAir(pBN[ONE_SIDE]->T)-1)/2)*M[DOWNSTREAM]));	// Eq. 8.156	
		U_star[UPSTREAM] = M[UPSTREAM]*A_star[UPSTREAM];	// Eq. 8.175

		A_star[DOWNSTREAM] = lambda_in_star_c2/(1 - ((pPpt->gammaAir(pBN[ONE_SIDE]->T)-1)/2)*M[DOWNSTREAM]);	// Eq. 8.176
		U_star[DOWNSTREAM] = M[DOWNSTREAM]*A_star[DOWNSTREAM];	// Eq. 8.177
		
		lambda_in_d1 = AA_1*(A_star[UPSTREAM] + ((pPpt->gammaAir(pBN[ONE_SIDE]->T)-1)/2)*U_star[UPSTREAM]);	// Eq. 8.178

		M2_positive = false;

		//error_prev_prev = error_prev;
		error_prev = error;
		error = lambda_in[UPSTREAM] - lambda_in_d1;
		if(fabs(error)<minError) minError = error;	// Record lowest error so far
///*
		if(lambda_in_d1>lambda_in[UPSTREAM])
		{
			if(!GREATER) del_M2 = del_M2/2;	// Halve step size only if error has changed sign
	
			temp_M2 = M[DOWNSTREAM] - del_M2;
///*
			do
			{
				temp_M2 = M[DOWNSTREAM] - del_M2;
				if(temp_M2 <= 0) del_M2 = del_M2/2;
				else M2_positive = true;
			}while(!M2_positive);
//*/

			M[DOWNSTREAM] = temp_M2;

			//M[DOWNSTREAM] -= del_M2;
			GREATER = true;
		}
		else
		{
			if(GREATER) del_M2 = del_M2/2; // Halve step size only if error has changed sign

			temp_M2 = M[DOWNSTREAM] + del_M2;
///*
			do
			{
				temp_M2 = M[DOWNSTREAM] + del_M2;
				if(temp_M2 <= 0) del_M2 = del_M2/2;
				else M2_positive = true;
			}while(!M2_positive);
//*/

			M[DOWNSTREAM] = temp_M2;

			//M[DOWNSTREAM] += del_M2;
			GREATER = false;
		}
		lambda_in_c2 = lambda_in_n[DOWNSTREAM] + A_star[DOWNSTREAM]*(AA_c2 - AA_n[DOWNSTREAM]);	// Eq. 8.179


		if(counter>1000) // Loop not ending
		{
			pPpt->Out("\n");
			pPpt->Out(Identify()); pPpt->Out(":Loss: counter = "); pPpt->Out(counter); pPpt->Out("\n");
			pPpt->Out("timestep = "); pPpt->Out(timestep); pPpt->Out("\n");
			pPpt->Out("error = "); pPpt->Out(error); pPpt->Out("\n");
			pPpt->Out("error_prev = "); pPpt->Out(error_prev); pPpt->Out("\n");
			pPpt->Out("\n");
			exit(1);
		}


		if(
			//fabs(lambda_in[UPSTREAM] - lambda_in_d1)<limit
			fabs(error)<limit
			//|| fabs(del_M2)<1e-20
			) 
			M2_converged = true;
//*/	

/*
//if(timestep>=3369)
//if(timestep>=3362)
//if(timestep>=1264)
//if(timestep>=3515)
{
	cout << "counter = " << counter << endl;
	cout << "lambda_in[UPSTREAM] = " << lambda_in[UPSTREAM] << endl;
	cout << "lambda_in_d1 = " << lambda_in_d1 << endl;
	//cout << "(lambda_in[UPSTREAM] - lambda_in_d1) = " << (lambda_in[UPSTREAM] - lambda_in_d1) << endl;
	//cout << "fabs(lambda_in[UPSTREAM] - lambda_in_d1) = " << fabs(lambda_in[UPSTREAM] - lambda_in_d1) << endl;
	cout << "error_prev = " << error_prev << endl;
	cout << "error = " << error << endl;
	//cout << "limit = " << limit << endl;
	cout << "M[DOWNSTREAM] = " << M[DOWNSTREAM] << endl;
	cout << "del_M2 = " << del_M2 << endl;
	cout << endl;
}
//*/

	}while(!M2_converged);

/*
{
	cout << endl;
	cout << "***CONVERGED***" << endl;
	cout << "counter = " << counter << endl;
	cout << "lambda_in[UPSTREAM] = " << lambda_in[UPSTREAM] << endl;
	cout << "lambda_in_d1 = " << lambda_in_d1 << endl;
	//cout << "(lambda_in[UPSTREAM] - lambda_in_d1) = " << (lambda_in[UPSTREAM] - lambda_in_d1) << endl;
	//cout << "fabs(lambda_in[UPSTREAM] - lambda_in_d1) = " << fabs(lambda_in[UPSTREAM] - lambda_in_d1) << endl;
	cout << "error_prev = " << error_prev << endl;
	cout << "error = " << error << endl;
	//cout << "limit = " << limit << endl;
	cout << "M[DOWNSTREAM] = " << M[DOWNSTREAM] << endl;
	cout << "del_M2 = " << del_M2 << endl;
	cout << "***CONVERGED***" << endl;
	cout << endl;
}
//*/


	lambda_out[UPSTREAM] = AA_1*(A_star[UPSTREAM] - ((pPpt->gammaAir(pBN[ONE_SIDE]->T)-1)/2)*U_star[UPSTREAM]);	// Eq. 8.180
	lambda_out[DOWNSTREAM] = AA_c2*(A_star[DOWNSTREAM] + ((pPpt->gammaAir(pBN[ONE_SIDE]->T)-1)/2)*U_star[DOWNSTREAM]);	// Eq. 8.181
/*
	cout << "CONVERGED:\n";
	cout << "lambda_out[UPSTREAM] = " << lambda_out[UPSTREAM] << endl;
	cout << "lambda_out[DOWNSTREAM] = " << lambda_out[DOWNSTREAM] << endl;
*/
	// Update lambda_in2 = lambda_in_c2, AA_2 = AA_c2, lambda_out2, lambda_out1, by reference
	// ======================================================================================
	rlambda_in_c2 = lambda_in_c2;
	rAA_c2 = AA_c2;
	rlambda_out1 = lambda_out[UPSTREAM];
	rlambda_out2 = lambda_out[DOWNSTREAM];

	// Update validation variables
	// ===========================
//	V1 = U_star[UPSTREAM]/lambda_in_star2;
//	del1 = A_star[UPSTREAM]/lambda_in_star2;
//	del2 = A_star[DOWNSTREAM]/lambda_in_star2;

	return;
}

void CTransmissive::PrintToFile(CProperties* pPpt, int timestep, double time, double ca)
// ============================================================ //
// Prints instantaneous transmissive boundary data to file.		//
// This function is called from the main function.				//
// ============================================================ //
{
	int f, m;
	f = 0;
	
	if(timestep==0)
	{
		fprintf(FILE_LOC,"%s", "RESULTS FILE FOR");
		if(EX) fprintf(FILE_LOC,"%s", " EXHAUST");
		else fprintf(FILE_LOC,"%s", " INTAKE");
		fprintf(FILE_LOC,"%s%i%s\n", " TRANSMISSIVE BOUNDARY [", ID, "]");
		//fprintf(FILE_LOC,"%s", "--------------------------------------------------");

		if(PRINT_PULSE_INFO)
		{
			// Print new pulse headings
			fprintf(FILE_LOC,"\n\t\t%s\n\n", "***START OF NEW PULSE***");		
			fprintf(FILE_LOC,"%s%i\n", "Pulse No. ", whole_number_of_pulses + 1);
			//fprintf(FILE_LOC,"%s\n", "-------------");
		}

		fprintf(FILE_LOC,"\n%s\t\t\t\t%s\t\t\t\t%s\t\t\t\t%s\n", "Timing", "Desired values", "Transmissive b.c. values applied", "Resultant values recorded");
		fprintf(FILE_LOC,"\t\t\t\t\t\t\t\t\t\t\t\t");
		for(m=0; m<nmeasurements; ++m) fprintf(FILE_LOC,"%s%i%s%.1f%s\t\t\t\t\t\t\t\t\t\t", "Loc. [", m, "] at ", loc_measure[m]*100, "% length of attached pipe");
		fprintf(FILE_LOC,"\n");
		
		fprintf(FILE_LOC,"%s\t%s\t%s\t\t%s\t%s\t%s\t\t%s\t%s\t%s", "Simulation time(s)", "Time in pulse (s)", "Rotation in pulse (deg.)", "Desired P (bar)", "Desired T (K)", "Desired u (m/s)", "App. P (bar)", "App. T (K)", "App. u (m/s)");
		for(m=0; m<nmeasurements; ++m)
			fprintf(FILE_LOC,"\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
			"ps (bar)", "Ts (K)", "u (m/s)", "p0 (bar)", "T0 (K)", "MFR (kg/s)", "PR", "MFP", "Re");
		fprintf(FILE_LOC,"\n");
	}
	
	if(time >= print_from_time && whole_number_of_pulses >= print_from_pulse-1)// Only start recording data as specified
	{
		if(START_OF_NEW_PULSE && PRINT_PULSE_INFO && whole_number_of_pulses>0)
		{
			fprintf(FILE_LOC,"\n");
			fprintf(FILE_LOC,"\t\t%s\n", "***END OF PULSE***");
			fprintf(FILE_LOC,"\n");
			if(MATCH) fprintf(FILE_LOC,"\t\t%s\n", "***MATCHED PULSE***"); else fprintf(FILE_LOC,"\t\t%s\n", "***UNMATCHED PULSE***");
			fprintf(FILE_LOC,"\n");
			//fprintf(FILE_LOC,"%s%i%s\n", "Completed pulse no. " , whole_number_of_pulses, " parameters");
			fprintf(FILE_LOC,"%s%i%s\n", "Completed pulse no. " , whole_number_of_pulses, " parameters");
			fprintf(FILE_LOC,"%s\n", "------------------------------------------------------------------------------------------------");

			// Print headings
			fprintf(FILE_LOC,"%s\t\t\t", "Parameter");
			for(m=0; m<nmeasurements; ++m)
			{
				fprintf(FILE_LOC,"%s%i%s\t", "Loc. [", m, "]");
			}
			fprintf(FILE_LOC,"\n");
			
			fprintf(FILE_LOC,"%s\t\t\t", "Loc. along pipe (%)");
			for(m=0; m<nmeasurements; ++m) fprintf(FILE_LOC,"%f\t", loc_measure[m]*100);
			fprintf(FILE_LOC,"\n");
	
			// Print data - average_flow_velocity_last_pulse
			fprintf(FILE_LOC,"%s\t\t\t", "Mean flow vel., u (m/s)");
			for(m=0; m<nmeasurements; ++m) fprintf(FILE_LOC,"%f\t", average_flow_velocity_last_pulse[m]);
			fprintf(FILE_LOC,"\n");

			// Print data - average_press_velocity_last_pulse_uplusa
			fprintf(FILE_LOC,"%s\t\t\t", "Mean pres. vel., u+a (m/s)");
			for(m=0; m<nmeasurements; ++m) fprintf(FILE_LOC,"%f\t", average_press_velocity_last_pulse_uplusa[m]);
			fprintf(FILE_LOC,"\n");

			// Print data - average_press_velocity_last_pulse_uminusa
			fprintf(FILE_LOC,"%s\t\t\t", "Mean pres. vel., u-a (m/s)");
			for(m=0; m<nmeasurements; ++m) fprintf(FILE_LOC,"%f\t", average_press_velocity_last_pulse_uminusa[m]);
			fprintf(FILE_LOC,"\n");

			// Print data - St
			fprintf(FILE_LOC,"%s\t\t\t", "St");
			for(m=0; m<nmeasurements; ++m) fprintf(FILE_LOC,"%f\t", st[m]);
			fprintf(FILE_LOC,"\n");

			// Print data - beta
			fprintf(FILE_LOC,"%s\t\t\t", "beta");
			for(m=0; m<nmeasurements; ++m) fprintf(FILE_LOC,"%f\t", beta[m]);
			fprintf(FILE_LOC,"\n");

			// Print data - MSt
			fprintf(FILE_LOC,"%s\t\t\t", "MSt");
			for(m=0; m<nmeasurements; ++m) fprintf(FILE_LOC,"%f\t", mst[m]);
			fprintf(FILE_LOC,"\n");
	
			// Print data - PMSt_uplusa
			fprintf(FILE_LOC,"%s\t\t\t", "PMSt u+a");
			for(m=0; m<nmeasurements; ++m) fprintf(FILE_LOC,"%f\t", pmst_uplusa[m]);
			fprintf(FILE_LOC,"\n");
		
			// Print data - PMSt_uminusa
			fprintf(FILE_LOC,"%s\t\t\t", "PMSt u-a");
			for(m=0; m<nmeasurements; ++m) fprintf(FILE_LOC,"%f\t", pmst_uminusa[m]);
			fprintf(FILE_LOC,"\n");
			fprintf(FILE_LOC,"%s\n\n", "------------------------------------------------------------------------------------------------");
	
			// Print new pulse headings
			fprintf(FILE_LOC,"\n\t\t%s\n\n", "***START OF NEW PULSE***");		
			fprintf(FILE_LOC,"%s%i\n", "Pulse No. ", whole_number_of_pulses + 1);
			//fprintf(FILE_LOC,"%s\n", "-------------");
	
			fprintf(FILE_LOC,"\n%s\t\t\t\t%s\t\t\t\t%s\t\t\t\t%s\n", "Timing", "Desired values", "Transmissive b.c. values applied", "Resultant values recorded");
			fprintf(FILE_LOC,"\t\t\t\t\t\t\t\t\t\t\t\t");
			for(m=0; m<nmeasurements; ++m) fprintf(FILE_LOC,"%s%i%s%.1f%s\t\t\t\t\t\t\t\t\t\t", "Loc. [", m, "] at ", loc_measure[m]*100, "% length of attached pipe");
			fprintf(FILE_LOC,"\n");									// "%s%i%s%.1#f%s\t\t\t\t\t\t\t\t"
					
			fprintf(FILE_LOC,"%s\t%s\t%s\t\t%s\t%s\t%s\t\t%s\t%s\t%s", "Simulation time(s)", "Time in pulse (s)", "Rotation in pulse (deg.)", "Desired P (bar)", "Desired T (K)", "Desired u (m/s)", "App. P (bar)", "App. T (K)", "App. u (m/s)");

			for(m=0; m<nmeasurements; ++m)
				fprintf(FILE_LOC,"\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
				"ps (bar)", "Ts (K)", "u (m/s)", "p0 (bar)", "T0 (K)", "MFR (kg/s)", "PR", "MFP", "Re");
			fprintf(FILE_LOC,"\n");
		}

		if(timestep%freq==0) // Print data at the specified sampling frequency
		{
			fprintf(FILE_LOC,"%f\t%f\t%f\t\t%f\t%f\t%f\t\t%f\t%f\t%f", time, time_in_period, (time_in_period/period)*360, 
				p_wave_desired, T_wave_desired, v_wave_desired, p_wave, T_WAVEFORM, u_WAVEFORM);
			for(m=0; m<nmeasurements; ++m)
				fprintf(FILE_LOC,"\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", 
				ps_loc[m], Ts_loc[m], Velocity[m], p0_loc[m], T0_loc[m], MassFlowRate[m], PR[m], MFP[m], Re[m]);
			fprintf(FILE_LOC,"\n");
		}
	}
}

void CTransmissive::PrintToScreen(CProperties* pPpt)
{
	int m;
	pPpt->Out(Underline(Identify(), "=", "", strDesc)); pPpt->Out("\n");

	pPpt->Out("\tp_wave = "); pPpt->Out(p_wave); pPpt->Out("\n\n");

	if(START_OF_NEW_PULSE)
	{
		// Print out details for last pulse completed
		pPpt->Out("\t***START OF NEW PULSE***\n");
		pPpt->Out("\n");
		pPpt->Out("\tCompleted pulse parameters\n");
		pPpt->Out("\t--------------------------\n");
		for(m=0; m<nmeasurements; ++m)
		{
			pPpt->Out("\tMean gas flow velocity, u, location["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(average_flow_velocity_last_pulse[m]); pPpt->Out(" m/s\n");
			pPpt->Out("\tMean pressure wave velocity, u + a, location["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(average_press_velocity_last_pulse_uplusa[m]); pPpt->Out(" m/s\n");
			pPpt->Out("\tMean pressure wave velocity, u - a, location["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(average_press_velocity_last_pulse_uminusa[m]); pPpt->Out(" m/s\n");
			pPpt->Out("\n");
			pPpt->Out("\tSt, location["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(st[m]); pPpt->Out(" \n");
			pPpt->Out("\tbeta, location["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(beta[m]); pPpt->Out(" \n");
			pPpt->Out("\tMSt, location["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(mst[m]); pPpt->Out(" \n");
			pPpt->Out("\tPMSt_uplusa, location["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(pmst_uplusa[m]); pPpt->Out(" \n");
			pPpt->Out("\tPMSt_uminusa, location["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(pmst_uminusa[m]); pPpt->Out(" \n");
		}
		pPpt->Out("\n");
	}
	else
	{
		pPpt->Out("\t***CONTINUING PULSE***\n\n");
	}

	pPpt->Out("\tCurrent pulse parameters\n");
	pPpt->Out("\t------------------------\n");
//		pPpt->Out("\tperiod = "); pPpt->Out(period); pPpt->Out("s \n");
	pPpt->Out("\ttime_in_period = "); pPpt->Out(time_in_period); pPpt->Out(" s\n");
	pPpt->Out("\twhole_number_of_pulses = "); pPpt->Out(whole_number_of_pulses); pPpt->Out("\n");
	pPpt->Out("\titerations_during_pulse = "); pPpt->Out(iterations_during_pulse); pPpt->Out("\n");
	for(m=0; m<nmeasurements; ++m)
	{
		pPpt->Out("\tsum_flow_velocity_over_pulse["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(sum_flow_velocity_over_pulse[m]); pPpt->Out(" m/s\n");
		pPpt->Out("\tsum_press_velocity_over_pulse_uplusa["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(sum_press_velocity_over_pulse_uplusa[m]); pPpt->Out(" m/s\n");
		pPpt->Out("\tsum_press_velocity_over_pulse_uminusa["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(sum_press_velocity_over_pulse_uminusa[m]); pPpt->Out(" m/s\n");
		pPpt->Out("\n");
		pPpt->Out("\n");
		pPpt->Out("\tPR[");  pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(PR[m]); pPpt->Out(" \n");
		pPpt->Out("\tMFP["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(MFP[m]); pPpt->Out(" \n");
		pPpt->Out("\tRe["); pPpt->Out(m); pPpt->Out("] = "); pPpt->Out(Re[m]); pPpt->Out(" \n");
	}
	pPpt->Out("\n");
	//if(CALIBRATE)
	//{
	//	pPpt->Out(Underline("Calibration", "-", "\t"));
	//	pPpt->Out("\tp_wave = "); pPpt->Out(p_wave); pPpt->Out(" bar\n");
	//	pPpt->Out("\tActual PR = "); pPpt->Out(PR[0]); pPpt->Out("\n");
	//	pPpt->Out("\tActual MFP = "); pPpt->Out(MFP[0]); pPpt->Out("\n");
	//	pPpt->Out("\tMFP_des = "); pPpt->Out(MFP_des); pPpt->Out("\n");
	//}
	pPpt->Out("\n");
	//if(ATTEMPT_MATCH)
	{
		pPpt->Out("\tATTEMPT_MATCH\t=\t"); pPpt->Out(TrueOrFalse(ATTEMPT_MATCH)); pPpt->Out("\n");
	}
	pPpt->Out("\n");
	pPpt->Out("\n");
}

void CTransmissive::WriteFiles()
{
/*
	int n, columncounter, r, d;
	if(num_props_measured>0)
	{
		for(n=0; n<nmeasurements; ++n)
		{
			columncounter=0;

			fprintf(this->FILE_M[n], "%s\t", "Location (m)");
			++columncounter;

			fprintf(this->FILE_M[n], "%s\t", "Time (s)");
			++columncounter;

			if(FLOW_VELOCITY && columncounter<=num_props_measured+1)
			{
				fprintf(this->FILE_M[n], "%s\t", "Flow Velocity (m/s)");
				++columncounter;
			}
			if(PRESSURE_VELOCITY && columncounter<=num_props_measured+1)
			{
				fprintf(this->FILE_M[n], "%s\t", "Pressure Wave Velocity (m/s)");
				++columncounter;
			}

			for(r=0; r<rowcounter; ++r)
			{
				for(d=0; d<num_props_measured+1; ++d)
				{
					fprintf(this->FILE_M[n], "%f\t", Results[n][r][d]);
				}
				fprintf(this->FILE_M[n],"%f\n", Results[n][r][num_props_measured+1]);
			}
			fclose(this->FILE_M[n]);
		}
	}
*/
}

void CTransmissive::CloseFiles()
{
	fclose(FILE_LOC);
	fclose(FILE_p_des_app_seen);
}

void CTransmissive::RecordTappings(CProperties* pPpt, int timestep, double time, CPipe* Pipe)
{
		int n, r, d, S, columncounter;
		bool located;

//	if((timestep % sample_factor) == 0)
//	{
		for(n=0; n<this->nmeasurements; ++n)
		{
			// Interpolate at each tapping location
			if(Measurements[n]<0 || Measurements[n]>1) {cout << "Measurement location outside of pipe\n"; break;}

			// Locate the two nodes either side of the location of this tapping
			S=0;
			if(this->pPipe[ONE_SIDE]->N > 1)
			{
				located = false;
				while(!located)
				{
					if(this->pPipe[ONE_SIDE]->Node[S+1].x > (Measurements[n]*this->pPipe[ONE_SIDE]->length) ||
						fabs(this->pPipe[ONE_SIDE]->Node[S+1].x - (Measurements[n]*this->pPipe[ONE_SIDE]->length)) < 1e-6)
						located = true;
					else ++S;
				};
				
				// Then the location is between nodes S and S+1
				MeasureNode[n] =  this->pPipe[ONE_SIDE]->Node[S] +
					(this->pPipe[ONE_SIDE]->Node[S+1] - this->pPipe[ONE_SIDE]->Node[S])*
					(((Measurements[n]*this->pPipe[ONE_SIDE]->length) - this->pPipe[ONE_SIDE]->Node[S].x)
					/(this->pPipe[ONE_SIDE]->Node[S+1].x - this->pPipe[ONE_SIDE]->Node[S].x));
			}
			else MeasureNode[n] =  this->pPipe[ONE_SIDE]->Node[S]; // For single node pipes

			// Now update derived properties for this node
			// Reynold's number
			MeasureNode[n].Re = ( ( (MeasureNode[n].p_dash*pPpt->PREF*1e5)/(pPpt->R_air*MeasureNode[n].T) ) * ( fabs(MeasureNode[n].U)*pPipe[ONE_SIDE]->AREF ) * MeasureNode[n].d )
							//	/ pPpt->mu_air;
								/ pPpt->ViscosityAir(MeasureNode[n].T);
		}
//	}

	if((timestep % sample_factor) == 0)
	{
//		int n, r, d, S, columncounter;
//		bool located;

		// Check whether we are above the limit of number of points to record
		if(rowcounter>=this->max_pts)
		{
			// Go through results array and remove every second row, shift up
			for(n=0; n<this->nmeasurements; ++n)
				for(r=1; r <= int((this->max_pts-1)/2); ++r)
					for(d=0; d<= this->num_props_measured+1; ++d) Results[n][r][d] = Results[n][2*r][d];

			rowcounter = int((this->max_pts-1)/2) + 1;
			sample_factor *= 2; // Reduce sampling rate by half
//			cout << "Reducing sample rate to every " << sample_factor << " timesteps\n";
		}

		for(n=0; n<this->nmeasurements; ++n)
		{
			// Then extract the desired properties from the interpolated node
			columncounter=0;
			if(this->num_props_measured>0)
			{
//				Results[n][rowcounter][columncounter] = MeasureNode[n].X_REAL*pPpt->xref;
				Results[n][rowcounter][columncounter] = MeasureNode[n].x;
//cout << "x = " << Results[n][rowcounter][columncounter] << endl;
				++columncounter;
			}
			if(this->num_props_measured>0)
			{
				Results[n][rowcounter][columncounter] = time;
				++columncounter;
			}

//			if(STATIC_PRESSURE && columncounter<=this->num_props_measured+1)
			if(FLOW_VELOCITY && columncounter<=this->num_props_measured+1)
			{
				Results[n][rowcounter][columncounter] = MeasureNode[n].U*pPipe[ONE_SIDE]->AREF;
//cout << "u = " << Results[n][rowcounter][columncounter] << endl;
				++columncounter;
			}
//			if(TEMPERATURE && columncounter<=this->num_props_measured+1)
			if(PRESSURE_VELOCITY && columncounter<=this->num_props_measured+1)
			{
				Results[n][rowcounter][columncounter] = MeasureNode[n].U*pPipe[ONE_SIDE]->AREF;
				++columncounter;
			}
//			if(REYNOLDS_NO && columncounter<=this->num_props_measured+1)
//			{
//				Results[n][rowcounter][columncounter] = MeasureNode[n].Re;
//				++columncounter;
//			}
		}
		++ rowcounter;
	}

	// Pulse calculations - sums
	Pressure_0_prev = Pressure_0; // Saved
//	Pressure_0 = pBN->p_dash*pPpt->PREF;
	
	MFR_0_prev = MFR_0;
	MFR_0 = pBN[ONE_SIDE]->mdot;
//	MFR_0 = this->pPipe->Node[1].mdot;
	
	VEL_0_prev = VEL_0;
//	VEL_0 = pBN->U*pPipe->AREF;
	VEL_0 = pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1].U*pPipe[ONE_SIDE]->AREF; // Other end

	TEMP_0_prev = TEMP_0;
	TEMP_0 = pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1].T; // Other end

	P_VEL_0_prev = P_VEL_0;
//	P_VEL_0 = pBN->U*pPipe->AREF + pBN->A*pPipe->AREF; // U +PLUS A
	P_VEL_0 = pBN[ONE_SIDE]->U*pPipe[ONE_SIDE]->AREF - pBN[ONE_SIDE]->A*pPipe[ONE_SIDE]->AREF; // U -MINUS A
//	P_VEL_0 = fabs(pBN->U*pPipe->AREF) + pBN->A*pPipe->AREF; // fabs(U) +PLUS A
//	P_VEL_0 = fabs(pBN->U*pPipe->AREF) - pBN->A*pPipe->AREF; // fabs(U) -MINUS A


	for(int m=0; m<this->nmeasurements; ++m)
	{
//		sum_flow_velocity_over_pulse[m] += MeasureNode[m].U*pPipe->AREF;
		sum_flow_velocity_over_pulse[m] += fabs(MeasureNode[m].U)*pPipe[ONE_SIDE]->AREF;
//		sum_press_velocity_over_pulse_uplusa[m] += MeasureNode[m].U*pPipe->AREF + MeasureNode[m].A*pPipe->AREF;
		sum_press_velocity_over_pulse_uplusa[m] += fabs(MeasureNode[m].U)*pPipe[ONE_SIDE]->AREF + MeasureNode[m].A*pPipe[ONE_SIDE]->AREF;
//		sum_press_velocity_over_pulse_uminusa[m] += MeasureNode[m].U*pPipe->AREF - MeasureNode[m].A*pPipe->AREF;
		sum_press_velocity_over_pulse_uminusa[m] += fabs(MeasureNode[m].U)*pPipe[ONE_SIDE]->AREF - MeasureNode[m].A*pPipe[ONE_SIDE]->AREF;

		ps_loc[m]/*bar*/ = MeasureNode[m].p_dash*pPpt->PREF;
		Ts_loc[m] = MeasureNode[m].T;
		Velocity[m] = MeasureNode[m].U*pPipe[ONE_SIDE]->AREF;
		p0_loc[m]/*bar*/ = TotalPressureBar(pPpt, ps_loc[m]/*bar*/, Ts_loc[m], Velocity[m]);
		T0_loc[m] = TotalTemperature(pPpt, Ts_loc[m], Velocity[m]);
		MassFlowRate[m] = MeasureNode[m].mdot;

		double k, mdot_temp, T1_temp, T01_temp, C1_temp, p1_temp, p01_temp, p2_temp; 
		if(end[ONE_SIDE]==ODD){
			mdot_temp = pPipe[ONE_SIDE]->Node[0].mdot;
			T1_temp = pPipe[ONE_SIDE]->Node[0].T;
			C1_temp = pPipe[ONE_SIDE]->Node[0].U*pPipe[ONE_SIDE]->AREF;
			p1_temp/*bar*/ = pPipe[ONE_SIDE]->Node[0].p_dash*pPpt->PREF;
		}
		else{
			mdot_temp = pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1].mdot;
			T1_temp = pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1].T;
			C1_temp = pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1].U*pPipe[ONE_SIDE]->AREF;
			p1_temp/*bar*/ = pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1].p_dash*pPpt->PREF;
		}

		// Calculate PR, MFP etc
		// Interpolate p2_pipe for tapping placed at p2_loc; locate the two nodes either side of the location of this tapping
		S=0;
		CNode p2_loc_node;
		if(Pipe[p2_pipe].N > 1)
		{
			while(Pipe[p2_pipe].Node[S].x < p2_loc*Pipe[p2_pipe].length && S < Pipe[p2_pipe].N - 1) ++S; // p2_loc is between S and S-1
			if(S==0) S=1;
			
			p2_loc_node = Pipe[p2_pipe].Node[S-1] +
				(Pipe[p2_pipe].Node[S] - Pipe[p2_pipe].Node[S-1])*
				((p2_loc*Pipe[p2_pipe].length - Pipe[p2_pipe].Node[S-1].x)
				/(Pipe[p2_pipe].Node[S].x - Pipe[p2_pipe].Node[S-1].x));
		}
		else p2_loc_node = Pipe[p2_pipe].Node[S]; // For single node pipes
		p2_temp/*bar*/ = p2_loc_node.p_dash*pPpt->PREF;

		k = pPpt->gammaAir(T1_temp);
		T01_temp = T1_temp*(1 + ((k-1)/2)*(pow(C1_temp,2)/(k*pPpt->R_air*T1_temp)));
		p01_temp/*bar*/ = p1_temp*pow(T01_temp/T1_temp, k/(k-1));
	
/*
		// If filling and emptying recalculate
		if(pPpt->METHOD != 3) 
		{

			if(end[ONE_SIDE]==ODD)
			{
				mdot_temp = pPipe[ONE_SIDE]->Node[0].mdot;
				T1_temp = pPipe[ONE_SIDE]->Node[0].T;
				C1_temp = pPipe[ONE_SIDE]->Node[0].U*pPipe[ONE_SIDE]->AREF;
				p1_temp = pPipe[ONE_SIDE]->Node[0].p_dash*pPpt->PREF;
			}
			else
			{
				mdot_temp = pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1].mdot;
				T1_temp = pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1].T;
				C1_temp = pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1].U*pPipe[ONE_SIDE]->AREF;
				p1_temp = pPipe[ONE_SIDE]->Node[pPipe[ONE_SIDE]->N-1].p_dash*pPpt->PREF;
			}

			mdot_temp = pPipe[ONE_SIDE]->Node[0].mdot;
			T1_temp = pPipe[ONE_SIDE]->FV[0].T;
			C1_temp = pPipe[ONE_SIDE]->Node[0].U*pPipe[ONE_SIDE]->AREF;
			p1_temp = pPipe[ONE_SIDE]->FV[0].p;


			T01_temp = T1_temp*(1 + ((pPpt->gammaAir()-1)/2)*(pow(C1_temp,2)/(pPpt->gammaAir()*pPpt->R_air*T1_temp)));
			p01_temp = p1_temp*( pow(T01_temp/T1_temp, pPpt->gammaAir()/(pPpt->gammaAir()-1)) );

		}
*/

		// Form of results:
		// ================
		// PR
		// ==
//		PR[m] = p1_temp/p2_temp;
		PR[m] = p01_temp/p2_temp;					// // My standard MFP, also MHI the same
//		PR[m] = p01_temp/1;
//		PR[m] = p1_temp/p2_temp;
//		PR[m] = p1_temp/1;
		

		// MFP - for calibration the units of MFP must the same as on the steady curve provided
		// ====================================================================================
//		MFP[m] = mdot_temp*sqrt(T01_temp)/p1_temp;
		MFP[m] = mdot_temp*sqrt(T01_temp)/p01_temp;  // (kg/s)K^0.5/bar
//		MFP[m] = mdot_temp*sqrt(T1_temp)/p01_temp;   // MHI GRTP with P in kPa
//		MFP[m] = mdot_temp*sqrt(T01_temp)/p1_temp;
//		MFP[m] = mdot_temp*sqrt(T1_temp)/p1_temp;
//		MFP[m] = mdot_temp*sqrt(T1_temp)/p01_temp;

		Re[m] = MeasureNode[m].Re;

		//// Select whether desired pressure is a total or static pressure
		//if(TOTAL_PRESS) Pressure_0 = p01_temp;
		//else 
			Pressure_0 = p1_temp;

//cout << "Pressure_0 = " << Pressure_0 << endl;

//		MFR_0 = ((p1_temp*1e5)/(pPpt->R_air*T1_temp))*this->pBN->f*C1_temp;
//		MFR_0 = ((p01_temp*1e5)/(pPpt->R_air*T1_temp))*this->pBN->f*C1_temp;
//		MFR_0 = ((p01_temp*1e5)/(pPpt->R_air*T01_temp))*this->pBN->f*C1_temp;
	}

	// Write the seen pressure to the p_des_app_seen array
	// Find the relevant data point in p_des_app_seen for the current time_in_period
	int i=0;
	do{++i;}
	while(p_des_app_seen[i][TIME_IN_PERIOD]<time_in_period && i<nsamples-1);
	// time_in_period lies between p_des_app_seen[i-1][TIME_IN_PERIOD] and p_des_app_seen[i][TIME_IN_PERIOD]

	// Extrapolate for pressure at time [i] and fill in for any gaps in 
	// p_des_app_seen that were jumped between iterations
	
	if(whole_number_of_pulses == 0 && iterations_during_pulse ==1) // On the very first time through here
	{
		Pressure_0_prev = Pressure_0;
		MFR_0_prev = MFR_0;
		VEL_0_prev = VEL_0;
		P_VEL_0_prev = P_VEL_0;
		TEMP_0_prev = TEMP_0;
	}

/*
	for(int filly=0; filly<nsamples; ++filly)
	{
		cout << "p_des_app_seen[filly=" << filly << "][SEEN_P=" << SEEN_P << "] = " << p_des_app_seen[filly][SEEN_P] << endl;
	}
*/
	for(int i_fill=i; i_fill>i_prev; --i_fill)
	{
		if(fabs(p_des_app_seen[i_fill][TIME_IN_PERIOD] - time_in_period) < 1e-9
			|| fabs(time_in_period - time_in_period_prev) < 1e-9)
		{
			p_des_app_seen[i_fill][SEEN_P] = Pressure_0;
		}
		else
			p_des_app_seen[i_fill][SEEN_P] = Pressure_0 +
							((Pressure_0 - Pressure_0_prev)/(time_in_period - time_in_period_prev)) // Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - time_in_period); // Time difference
/*
		if(i_fill==0 || i_fill==1)
		{
			cout << "i = " << i << endl;
			cout << "i_prev = " << i_prev << endl;
			cout << "Pressure_0 = " << Pressure_0 << endl;
			cout << "Pressure_0_prev = " << Pressure_0_prev << endl;
			cout << "time_in_period = " << time_in_period << endl;
			cout << "time_in_period_prev = " << time_in_period_prev << endl;
			cout << "p_des_app_seen[i_fill=" << i_fill << "][TIME_IN_PERIOD=" << TIME_IN_PERIOD << "] = " << p_des_app_seen[i_fill][TIME_IN_PERIOD] << endl;
			cout << "p_des_app_seen[i_fill=" << i_fill << "][SEEN_P=" << SEEN_P << "] = " << p_des_app_seen[i_fill][SEEN_P] << endl;
			cout << endl;
		}
//*/
		p_des_app_seen[i_fill][SEEN_MFR] = MFR_0 +
							((MFR_0 - MFR_0_prev)/(time_in_period - time_in_period_prev)) // Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - time_in_period); // Time difference

		p_des_app_seen[i_fill][SEEN_V] = VEL_0 +
							((VEL_0 - VEL_0_prev)/(time_in_period - time_in_period_prev)) // Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - time_in_period); // Time difference

		p_des_app_seen[i_fill][SEEN_T] = TEMP_0 +
							((TEMP_0 - TEMP_0_prev)/(time_in_period - time_in_period_prev)) // Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - time_in_period); // Time difference

		p_des_app_seen[i_fill][P_VEL] = P_VEL_0 +
							((P_VEL_0 - P_VEL_0_prev)/(time_in_period - time_in_period_prev)) // Gradient
							*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - time_in_period); // Time difference


		if(!MATCHED && ATTEMPT_MATCH) // Leave the NEXT_P column with previous pulse values when a match has been found, or when not attempting to match
		{
			double relax_factor = 1.0;
			p_des_app_seen[i_fill][NEXT_P] = p_des_app_seen[i_fill][APPLIED_P] +
									relax_factor*
									(p_des_app_seen[i_fill][DESIRED_P] - p_des_app_seen[i_fill][SEEN_P]);

			p_des_app_seen[i_fill][NEXT_T] = p_des_app_seen[i_fill][APPLIED_T] +
									relax_factor*
									(p_des_app_seen[i_fill][DESIRED_T] - p_des_app_seen[i_fill][SEEN_T]);

			p_des_app_seen[i_fill][NEXT_V] = p_des_app_seen[i_fill][APPLIED_V] +
									relax_factor*
									(p_des_app_seen[i_fill][DESIRED_V] - p_des_app_seen[i_fill][SEEN_V]);

			p_des_app_seen[i_fill][NEXT_MFR] = p_des_app_seen[i_fill][APPLIED_MFR] +
									relax_factor*
									(p_des_app_seen[i_fill][DESIRED_MFR] - p_des_app_seen[i_fill][SEEN_MFR]);
		}
/*
		p_des_app_seen[i_fill][UNDER_EST] = (p_des_app_seen[i_fill][DESIRED_MFR] - p_des_app_seen[i_fill][MFR]);

		if( (p_des_app_seen[i_fill][UNDER_EST] < 0 && p_des_app_seen[i_fill][UNDER_EST_PREV] > 0)
			|| (p_des_app_seen[i_fill][UNDER_EST] > 0 && p_des_app_seen[i_fill][UNDER_EST_PREV] < 0) )
		// If change of sign between cycles reduce the factor
		{
			p_des_app_seen[i_fill][UNDER_EST_FACT] *= 0.5;
		}

		if( fabs((p_des_app_seen[i_fill][UNDER_EST] - p_des_app_seen[i_fill][UNDER_EST_PREV])
					/p_des_app_seen[i_fill][UNDER_EST_PREV])*100 < 1 
			&& fabs(p_des_app_seen[i_fill][UNDER_EST]/p_des_app_seen[i_fill][DESIRED_MFR])*100 > 10 )
		// If not much change in the under estimate between cycles and still not close to the desired mfr
		{
			// Reset loss value
			p_des_app_seen[i_fill][K_LOSS] = 1;
			// Reset loss factor, maintain sign
			p_des_app_seen[i_fill][UNDER_EST_FACT] = p_des_app_seen[i_fill][UNDER_EST_FACT]/fabs(p_des_app_seen[i_fill][UNDER_EST_FACT]);
			// Change direction of convergence
			p_des_app_seen[i_fill][UNDER_EST_FACT] *= -1;
		}
		
		p_des_app_seen[i_fill][K_LOSS] -= p_des_app_seen[i_fill][UNDER_EST_FACT]*p_des_app_seen[i_fill][UNDER_EST];

		if(p_des_app_seen[i_fill][K_LOSS]<0) p_des_app_seen[i_fill][K_LOSS]=0;
*/
/*
		double desired_mfr = p_des_app_seen[i-1][DESIRED_MFR] +
				((p_des_app_seen[i][DESIRED_MFR] - p_des_app_seen[i-1][DESIRED_MFR])/
				(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
				(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);

		double actual_mfr = p_des_app_seen[i-1][MFR] +
				((p_des_app_seen[i][MFR] - p_des_app_seen[i-1][MFR])/
				(p_des_app_seen[i][TIME_IN_PERIOD] - p_des_app_seen[i-1][TIME_IN_PERIOD]))*
				(time_in_period - p_des_app_seen[i-1][TIME_IN_PERIOD]);
		
		mfr_under_est = desired_mfr - actual_mfr;
		
//		K = 1 - 1e-5*mfr_under_est; // Increase K when actual_mfr > desired_mfr etc.
//		if(K<0) K=0;
*/

	}
	i_prev = i;  // Record up to which data point has been filled during this iteration




//	for(int q=0; q<nsamples; ++q) cout << "p_des_app_seen[q=" << q << "][SEEN_P] = " << p_des_app_seen[q][SEEN_P] << endl;

}

void CTransmissive::ReadInput(CProperties* pPpt, char *InputFile)
{
	int last_entry;
	CommonReadInput(InputFile, pPpt->NUM_PARAMS, this->labels, this->values, this->strings, last_entry);	
	int r;
	for(r=0; r<last_entry+1; ++r)
	{
		// ====================================================================================================
		// Parameter file for Exhaust Transmissive [0]
		// ====================================================================================================

		// Optional object description (max. 500 characters - use underscores for spaces)
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "strDesc") == 0) strDesc = strings[r];

		if(strcmp(labels[r], "ANECHOIC_ONLY") == 0) ANECHOIC_ONLY = DoubleToBool(values[r]);
		if(strcmp(labels[r], "CONSTANT") == 0) CONSTANT = DoubleToBool(values[r]);
		if(strcmp(labels[r], "TOTAL_PRESS") == 0) TOTAL_PRESS = DoubleToBool(values[r]);

			// If constant conditions, i.e., CONSTANT == 1 == true
			// ----------------------------------------------------------------------------------------------------
			if(strcmp(labels[r], "constantp") == 0) constantp = values[r];
			if(strcmp(labels[r], "constantT") == 0) constantT = values[r];
			if(strcmp(labels[r], "period") == 0) period = values[r];

			// Else varying conditions, i.e., CONSTANT == 0 == false
			// ----------------------------------------------------------------------------------------------------

			// Waveform setup
			// ----------------------------------------------------------------------------------------------------
			if(strcmp(labels[r], "shape") == 0) shape = int(values[r]);
			if(strcmp(labels[r], "f") == 0) f = values[r];
			if(strcmp(labels[r], "phi") == 0) phi = values[r];
			if(strcmp(labels[r], "phi_orig") == 0) phi_orig = values[r];
			if(strcmp(labels[r], "v") == 0) v = values[r];
			if(strcmp(labels[r], "Ts") == 0) Ts = values[r];
			
			// Amplitudes
			// ----------------------------------------------------------------------------------------------------
			if(strcmp(labels[r], "press_low") == 0) press_low = values[r];
			if(strcmp(labels[r], "press_high") == 0) press_high = values[r];
			if(strcmp(labels[r], "press_min") == 0) press_min = values[r];
			if(strcmp(labels[r], "Tsm_low") == 0) Tsm_low = values[r];
			if(strcmp(labels[r], "Tsm_high") == 0) Tsm_high = values[r];
			if(strcmp(labels[r], "Tsm_min") == 0) Tsm_min = values[r];
			if(strcmp(labels[r], "T0m_low") == 0) T0m_low = values[r];
			if(strcmp(labels[r], "T0m_high") == 0) T0m_high = values[r];
			if(strcmp(labels[r], "T0m_min") == 0) T0m_min = values[r];
			if(strcmp(labels[r], "MFR_low") == 0) MFR_low = values[r];
			if(strcmp(labels[r], "MFR_high") == 0) MFR_high = values[r];
			if(strcmp(labels[r], "MFR_min") == 0) MFR_min = values[r];
			
				// Waveform from file (shape==3)
				// ----------------------------------------------------------------------------------------------------
				if(strcmp(labels[r], "FILE_ps") == 0) FILE_ps = strings[r];
				if(strcmp(labels[r], "USE_FILE_Ts") == 0) USE_FILE_Ts = DoubleToBool(values[r]);
						
					// Temperature data from file (USE_FILE_Ts==1==true)
					// ----------------------------------------------------------------------------------------------------
					if(strcmp(labels[r], "FILE_Ts") == 0) FILE_Ts = strings[r];
					if(strcmp(labels[r], "WAVE_FILE_T0m") == 0) WAVE_FILE_T0m = strings[r];
					if(strcmp(labels[r], "WAVE_FILE_v") == 0) WAVE_FILE_v = strings[r];
					if(strcmp(labels[r], "WAVE_FILE_MFR") == 0) WAVE_FILE_MFR = strings[r];
			
				if(strcmp(labels[r], "ADJUST_TIME") == 0) ADJUST_TIME = DoubleToBool(values[r]);
				if(strcmp(labels[r], "ADJUST_PRESSURE") == 0) ADJUST_PRESSURE = DoubleToBool(values[r]);
				if(strcmp(labels[r], "ADJUST_TSM") == 0) ADJUST_TSM = DoubleToBool(values[r]);
				if(strcmp(labels[r], "ADJUST_T0M") == 0) ADJUST_T0M = DoubleToBool(values[r]);			
				if(strcmp(labels[r], "ADJUST_MFR") == 0) ADJUST_MFR = DoubleToBool(values[r]);

				// Fourier series waveform (shape==4). Set USE_FILE_Ts = 0.
				// ----------------------------------------------------------------------------------------------------
				if(strcmp(labels[r], "FOURIER_FILE") == 0) FOURIER_FILE = strings[r];
				
		// Pipe configuration
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "p2_pipe") == 0) p2_pipe = int(values[r]);
		if(strcmp(labels[r], "p2_loc") == 0) p2_loc = values[r];

		// Matching
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "INTELLIGENT") == 0) INTELLIGENT = DoubleToBool(values[r]);
		if(strcmp(labels[r], "nsamples") == 0) nsamples = int(values[r]);
		if(strcmp(labels[r], "tol") == 0) tol = double(values[r]);
		if(strcmp(labels[r], "SUSPEND") == 0) SUSPEND = DoubleToBool(values[r]);
    if(strcmp(labels[r], "nMatchCycles") == 0) nMatchCycles = int(values[r]);
//		if(strcmp(labels[r], "nmatches") == 0) nmatches = int(values[r]);

		// Measurements
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "USE_DEF_FREQ") == 0) USE_DEF_FREQ = DoubleToBool(values[r]);
		if(strcmp(labels[r], "freq") == 0) freq = int(values[r]);
		if(strcmp(labels[r], "print_from_pulse") == 0) print_from_pulse = int(values[r]);
		if(strcmp(labels[r], "print_from_time") == 0) print_from_time = values[r];
		if(strcmp(labels[r], "PRINT_PULSE_INFO") == 0) PRINT_PULSE_INFO = DoubleToBool(values[r]);
		if(strcmp(labels[r], "PRINT_MOVIE_FILE") == 0) PRINT_MOVIE_FILE = DoubleToBool(values[r]);
		if(strcmp(labels[r], "PRINT_TO_SCREEN") == 0) PRINT_TO_SCREEN = DoubleToBool(values[r]);

		// Locations (as a fraction of the attached pipe's physical length, measured from the odd end)
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "nmeasurements") == 0)
		{
			nmeasurements = int(values[r]);
			loc_measure = new double [nmeasurements];
		}
		if(strcmp(labels[r], "loc_measure0") == 0 && nmeasurements>0) loc_measure[0] = values[r];
		if(strcmp(labels[r], "loc_measure1") == 0 && nmeasurements>1) loc_measure[1] = values[r];
		if(strcmp(labels[r], "loc_measure2") == 0 && nmeasurements>2) loc_measure[2] = values[r];
		if(strcmp(labels[r], "loc_measure3") == 0 && nmeasurements>3) loc_measure[3] = values[r];
		if(strcmp(labels[r], "loc_measure4") == 0 && nmeasurements>4) loc_measure[4] = values[r];
		if(strcmp(labels[r], "loc_measure5") == 0 && nmeasurements>5) loc_measure[5] = values[r];
		if(strcmp(labels[r], "loc_measure6") == 0 && nmeasurements>6) loc_measure[6] = values[r];
		if(strcmp(labels[r], "loc_measure7") == 0 && nmeasurements>7) loc_measure[7] = values[r];
		if(strcmp(labels[r], "loc_measure8") == 0 && nmeasurements>8) loc_measure[8] = values[r];
		if(strcmp(labels[r], "loc_measure9") == 0 && nmeasurements>9) loc_measure[9] = values[r];
		if(strcmp(labels[r], "loc_measure10") == 0 && nmeasurements>10) loc_measure[10] = values[r];

		if(strcmp(labels[r], "max_pts") == 0) max_pts = int(values[r]);

		// Measurements to record
//		if(strcmp(labels[r], "num_props_measured") == 0) num_props_measured = int(values[r]);
		if(strcmp(labels[r], "FLOW_VELOCITY") == 0) FLOW_VELOCITY = DoubleToBool(values[r]);
		if(strcmp(labels[r], "PRESSURE_VELOCITY") == 0) PRESSURE_VELOCITY = DoubleToBool(values[r]);
//		if(strcmp(labels[r], "STATIC_PRESSURE") == 0) STATIC_PRESSURE = DoubleToBool(values[r]);
//		if(strcmp(labels[r], "TEMPERATURE") == 0) TEMPERATURE = DoubleToBool(values[r]);
		if(strcmp(labels[r], "MASS_FLOW_RATE") == 0) MASS_FLOW_RATE = DoubleToBool(values[r]);
		if(strcmp(labels[r], "MASS_FLOW_PARAMETER") == 0) MASS_FLOW_PARAMETER = DoubleToBool(values[r]);
		if(strcmp(labels[r], "PRESSURE_RATIO") == 0) PRESSURE_RATIO = DoubleToBool(values[r]);
//		if(strcmp(labels[r], "REYNOLDS_NO") == 0) REYNOLDS_NO = DoubleToBool(values[r]);

		//// Calibration of APLDev loss curve; uses v and T0 above
		//// ----------------------------------------------------------------------------------------------------
		//if(strcmp(labels[r], "CALIBRATE") == 0) CALIBRATE = DoubleToBool(values[r]);
		//if(strcmp(labels[r], "STEADY_FILE") == 0) STEADY_FILE = strings[r];
		//if(strcmp(labels[r], "cal_PR_low") == 0) cal_PR_low = values[r];
		//if(strcmp(labels[r], "cal_PR_high") == 0) cal_PR_high = values[r];
		//if(strcmp(labels[r], "cal_points") == 0) cal_points = int(values[r]);
		//if(strcmp(labels[r], "tol_MFP") == 0) tol_MFP = values[r];

		// ====================================================================================================
		// End of file
		// ====================================================================================================
	}
	// Set some derived parameters
	if(USE_DEF_FREQ) freq = pPpt->freq;		// Use the default sampling rate
	num_props_measured = int(FLOW_VELOCITY) + int(PRESSURE_VELOCITY) + int(MASS_FLOW_RATE) + int(MASS_FLOW_PARAMETER) + int(PRESSURE_RATIO);
}

void CTransmissive::LoadWaveform(char* InputFile, CProperties* pPpt, 
	bool TOTAL, bool ADJUST_TIME, bool ADJUST_PARAMETER, 
	double parameter_high, double parameter_low, double parameter_min, 
	int PARAMETER_DESIRED, 
	double &rparam_wave, double &rparam_wave_prev)
//--------------------------------------------------//
// Transmissive boundary user waveform				//
// -----------------------------------				//
// Loads parameter vs. time during pulse 			//
//													//
//--------------------------------------------------//
{
	int c = 0;
	//float temp;
	double temp;
	int num_columns = 2;
	int col, row;
	datapoints_parameter = 0;
	
	FILE* stream;
	stream = fopen(InputFile, "r");

	if(stream == NULL)
	{
		//printf("Error opening waveform parameter file\n");
		pPpt->Out("Error opening waveform parameter file: "); pPpt->Out(InputFile); pPpt->Out("\n");
		exit(1);
	}
	else
	{
		fseek(stream, 0L, SEEK_SET);	// Set pointer to beginning of file
		do
		{
			fscanf(stream, "%d", &c);	// Scans for the line number (ints), puts it in c
			//cout << "c = " << c << endl;
			if(c>datapoints_parameter) datapoints_parameter = c;
			fscanf(stream, "%lf", &temp);		// Runs over the time value
			//cout << temp << "\t";
			fscanf(stream, "%lf", &temp);		// Runs over the parameter value
			//cout << temp << "\n";
			//cin >> pause;
		}while(fscanf(stream, "%l")!=EOF);
		//cout << "Parameter raw file datapoints = " << datapoints_parameter << endl;
		//cin >> pause;
		// Can now dimension timing arrays
		parameter_raw = new double* [num_columns]; 
		for(col=0; col<num_columns; ++col) parameter_raw[col] = new double [datapoints_parameter];

		fseek(stream, 0L, SEEK_SET);	// Reset pointer to beginning of file

		for (row=0; row<datapoints_parameter; row++)
		{
			fscanf(stream, "%d", &c);	// Scans past the line number
			for (col=0; col<num_columns; col++)
			{
				// Finds the time and parameter values
				fscanf(stream, "%lf", &parameter_raw[col][row]);	
				//cout << "parameter_raw[col][row] = " << parameter_raw[col][row] << endl;
			}
//cout << "parameter_raw[time][" << row << "] = " << parameter_raw[0][row] << "\t";
//cout << "parameter_raw[parameter][" << row << "] = " << parameter_raw[1][row] << endl;

			if(TOTAL)
			{

//cout << "parameter_raw[time][" << row << "] = " << parameter_raw[0][row] << "\t";
//cout << "parameter_raw[pres][" << row << "] = " << parameter_raw[1][row] << "\t";

				// Now convert the parameter_raw to static parameter, if necessary
				// e.g. ps = p0 * ( 1 + ((gamma-1)/2) * (c^2/(gamma*R*T)) )^(-gamma/(gamma-1))
				parameter_raw[1][row] = parameter_raw[1][row]*
											pow(( 1 + ((pPpt->gammaAir(this->Ts)-1)/2)*(pow(this->v,2)/(pPpt->gammaAir(this->Ts)*pPpt->R_air*this->Ts)) ),
											-pPpt->gammaAir(this->Ts)/(pPpt->gammaAir(this->Ts)-1));

//cout << "parameter_raw[pres][" << row << "] = " << parameter_raw[1][row] << endl;
			}
		}
	}
	fclose(stream);

	// Now adjust raw values of time and parameter based on the required waveform parameters
	parameter_wave = new double* [num_columns];
	for(col=0; col<num_columns; ++col) parameter_wave[col] = new double [datapoints_parameter];

	// Find the read-in waveform high and low pressure
	double raw_low = parameter_raw[WAVE_PARAMETER][0];
	double raw_high = parameter_raw[WAVE_PARAMETER][0];
	if(ADJUST_PARAMETER) 
	{
		for(row=0; row<datapoints_parameter; ++row)
		{
			if(parameter_raw[WAVE_PARAMETER][row] < raw_low) raw_low = parameter_raw[WAVE_PARAMETER][row];
			if(parameter_raw[WAVE_PARAMETER][row] > raw_high) raw_high = parameter_raw[WAVE_PARAMETER][row];
		}
	}

	//cout << "TIME\tPARAMETER\n";
	for(row=0; row<datapoints_parameter; ++row)
	{
		if(ADJUST_TIME)
		{
			// Make the read-in wave fit inside a time of phi*period
			parameter_wave[WAVE_TIME][row] = (parameter_raw[WAVE_TIME][row]
				/parameter_raw[WAVE_TIME][datapoints_parameter-1]) // The time at the last datapoint
				*(phi*period);
		}
		else parameter_wave[WAVE_TIME][row] = parameter_raw[WAVE_TIME][row];


		if(ADJUST_PARAMETER)
		{
			// Make the read-in wave conform to parameter_low and parameter_high
			parameter_wave[WAVE_PARAMETER][row] = parameter_low + 
								((parameter_raw[WAVE_PARAMETER][row] - raw_low)/
								(raw_high - raw_low))
								*(parameter_high - parameter_low);

			if(parameter_wave[WAVE_PARAMETER][row] < parameter_min) parameter_wave[WAVE_PARAMETER][row] = parameter_min;
		}
		else
		{
			//if(parameter_raw[WAVE_PARAMETER][row] < parameter_min) parameter_wave[WAVE_PARAMETER][row] = parameter_min;
			//else parameter_wave[WAVE_PARAMETER][row] = parameter_raw[WAVE_PARAMETER][row];
			parameter_wave[WAVE_PARAMETER][row] = parameter_raw[WAVE_PARAMETER][row];
		}

//if(ID==1)
//cout << parameter_wave[WAVE_TIME][row] << "\t" << parameter_wave[WAVE_PARAMETER][row] << endl;
		//if(row == 200) cin >> pause;
	}

	// Now populate p_des_app_seen[][PARAMETER_DESIRED]
	for(int i_fill=0; i_fill<nsamples; ++i_fill)
	{
		row=0;
		do ++row;
		while(parameter_wave[WAVE_TIME][row] < p_des_app_seen[i_fill][TIME_IN_PERIOD] && row<datapoints_parameter-1);

		if(p_des_app_seen[i_fill][TIME_IN_PERIOD] <= phi*period)
		{
			// Interpolate between [row-1] and [row]
			p_des_app_seen[i_fill][PARAMETER_DESIRED] = parameter_wave[WAVE_PARAMETER][row-1] + 
											((parameter_wave[WAVE_PARAMETER][row] - parameter_wave[WAVE_PARAMETER][row-1])/
											(parameter_wave[WAVE_TIME][row] - parameter_wave[WAVE_TIME][row-1]))
											*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - parameter_wave[WAVE_TIME][row-1]);
		}
		else
		{
//cout << "Outside of phi*period" << endl;
			p_des_app_seen[i_fill][PARAMETER_DESIRED] = raw_low;
		}

// USEFUL for debugging:
//if(ID==0)
//cout << p_des_app_seen[i_fill][TIME_IN_PERIOD] << "\t" << p_des_app_seen[i_fill][PARAMETER_DESIRED] << endl;

/*
		if(TOTAL_PRESS)
		{
			// Now convert the same p_des_app_seen[i_fill][DESIRED_P] to static pressure, if necessary
			// ps = p0 * ( 1 + ((gamma-1)/2) * (c^2/(gamma*R*T)) )^(-gamma/(gamma-1))

			p_des_app_seen[i_fill][DESIRED_P] = p_des_app_seen[i_fill][DESIRED_P]*
											pow(( 1 + ((pPpt->gammaAir()-1)/2)*(pow(this->v,2)/(pPpt->gammaAir()*pPpt->R_air*this->Ts)) ),
												-pPpt->gammaAir()/(pPpt->gammaAir()-1));
		}

//cout << p_des_app_seen[i_fill][TIME_IN_PERIOD] << "\t" << p_des_app_seen[i_fill][DESIRED_P] << endl;
*/
	}

	// Initialise parameter_wave_prev to the first value in the array, else get funny values on first pulse
	rparam_wave = p_des_app_seen[0][PARAMETER_DESIRED];
	rparam_wave_prev = p_des_app_seen[0][PARAMETER_DESIRED];
	//cout << "rparam_wave_prev = " << rparam_wave_prev << endl;

//cout << endl;
//cout << InputFile << endl;
//if(ID==1)
//exit(1);
}

/*
void CTransmissive::LoadWaveformMFR(char* InputFile, CProperties* pPpt)
//--------------------------------------------------//
// Transmissive boundary user waveform				//
// -----------------------------------				//
// Loads MFR vs. time during pulse 					//
//													//
//--------------------------------------------------//
{
	int c = 0;
	float temp;
	int num_columns = 2;
	int col, row;
	datapoints_MFR = 0;

	char pause;
	
	FILE* stream;
	stream = fopen(InputFile, "r");

	if(stream == NULL)
	{
		printf("Error opening waveform MFR file\n");
		exit(1);
	}
	else
	{
		fseek(stream, 0L, SEEK_SET);	// Set pointer to beginning of file
		do
		{
			fscanf(stream, "%d", &c);	// Scans for the line number (ints), puts it in c
			//cout << "c = " << c << endl;
			if(c>datapoints_MFR) datapoints_MFR = c;
			fscanf(stream, "%lf", &temp);		// Runs over the time value
			fscanf(stream, "%lf", &temp);		// Runs over the MFR value
			//cin >> pause;
		}while(fscanf(stream, "%l")!=EOF);
		//cout << "MFR raw file datapoints = " << datapoints_MFR << endl;
		//cin >> pause;
		// Can now dimension timing arrays
		mfr_raw = new double* [num_columns]; 
		for(col=0; col<num_columns; ++col) mfr_raw[col] = new double [datapoints_MFR];

		fseek(stream, 0L, SEEK_SET);	// Reset pointer to beginning of file

		for (row=0; row<datapoints_MFR; row++)
		{
			fscanf(stream, "%d", &c);	// Scans past the line number
			for (col=0; col<num_columns; col++)
			{
				// Finds the time and MFR values
				fscanf(stream, "%lf", &mfr_raw[col][row]);	
				//cout << "mfr_raw[col][row] = " << mfr_raw[col][row] << endl;
			}
			//cout << "mfr_raw[time][" << row << "] = " << mfr_raw[0][row] << "\t";
			//cout << "mfr_raw[pres][" << row << "] = " << mfr_raw[1][row] << endl;

//			if(TOTAL_MFR)
//			{
//
////cout << "mfr_raw[time][" << row << "] = " << mfr_raw[0][row] << "\t";
////cout << "mfr_raw[pres][" << row << "] = " << mfr_raw[1][row] << "\t";
//
//				// Now convert the mfr_raw to static MFR, if necessary
//				// ps = p0 * ( 1 + ((gamma-1)/2) * (c^2/(gamma*R*T)) )^(-gamma/(gamma-1))
//				mfr_raw[1][row] = mfr_raw[1][row]*
//											pow(( 1 + ((pPpt->gammaAir()-1)/2)*(pow(this->v,2)/(pPpt->gammaAir()*pPpt->R_air*this->Ts)) ),
//											-pPpt->gammaAir()/(pPpt->gammaAir()-1));
//
////cout << "mfr_raw[pres][" << row << "] = " << mfr_raw[1][row] << endl;
//			}

		}
	}
	fclose(stream);

	// Now adjust raw values of time and pressure based on the required waveform parameters
	mfr_wave = new double* [num_columns];
	for(col=0; col<num_columns; ++col) mfr_wave[col] = new double [datapoints_MFR];

	// Find the read-in waveform high and low MFR
	double raw_low = mfr_raw[WAVE_MFR][0];
	double raw_high = mfr_raw[WAVE_MFR][0];
	if(ADJUST_MFR) 
	{
		for(row=0; row<datapoints_MFR; ++row)
		{
			if(mfr_raw[WAVE_MFR][row] < raw_low) raw_low = mfr_raw[WAVE_MFR][row];
			if(mfr_raw[WAVE_MFR][row] > raw_high) raw_high = mfr_raw[WAVE_MFR][row];
		}
	}

	//cout << "TIME\tMFR\n";
	for(row=0; row<datapoints_MFR; ++row)
	{
		if(ADJUST_TIME)
		{
			// Make the read-in wave fit inside a time of phi*period
			mfr_wave[WAVE_TIME][row] = (mfr_raw[WAVE_TIME][row]
				/mfr_raw[WAVE_TIME][datapoints_MFR-1]) // The time at the last datapoint
				*(phi*period);
		}
		else mfr_wave[WAVE_TIME][row] = mfr_raw[WAVE_TIME][row];


		if(ADJUST_MFR)
		{
			// Make the read-in wave conform to press_low and press_high
			mfr_wave[WAVE_MFR][row] = MFR_low + 
								((mfr_raw[WAVE_MFR][row] - raw_low)/
								(raw_high - raw_low))
								*(MFR_high - MFR_low);
		}
		else
		{
			if(mfr_raw[WAVE_MFR][row] < 0) mfr_wave[WAVE_MFR][row] = 0;
			else mfr_wave[WAVE_MFR][row] = mfr_raw[WAVE_MFR][row];
		}
//cout << mfr_wave[WAVE_TIME][row] << "\t" << mfr_wave[WAVE_MFR][row] << endl;
		//if(row == 200) cin >> pause;
	}

//cout << this->v << endl;
//cout << this->Ts << endl;

	// Now populate p_des_app_seen[][DESIRED_MFR]
	for(int i_fill=0; i_fill<nsamples; ++i_fill)
	{
		row=0;
		do ++row;
		while(mfr_wave[WAVE_TIME][row] < p_des_app_seen[i_fill][TIME_IN_PERIOD] && row<datapoints_MFR-1);

		// Interpolate between [row-1] and [row]
		p_des_app_seen[i_fill][DESIRED_MFR] = mfr_wave[WAVE_MFR][row-1] + 
										((mfr_wave[WAVE_MFR][row] - mfr_wave[WAVE_MFR][row-1])/
										(mfr_wave[WAVE_TIME][row] - mfr_wave[WAVE_TIME][row-1]))
										*(p_des_app_seen[i_fill][TIME_IN_PERIOD] - mfr_wave[WAVE_TIME][row-1]);

//cout << p_des_app_seen[i_fill][TIME_IN_PERIOD] << "\t" << p_des_app_seen[i_fill][DESIRED_MFR] << endl;


//		if(TOTAL_MFR)
//		{
//			// Now convert the same p_des_app_seen[i_fill][DESIRED_P] to static pressure, if necessary
//			// ps = p0 * ( 1 + ((gamma-1)/2) * (c^2/(gamma*R*T)) )^(-gamma/(gamma-1))
//
//			p_des_app_seen[i_fill][DESIRED_P] = p_des_app_seen[i_fill][DESIRED_P]*
//											pow(( 1 + ((pPpt->gammaAir()-1)/2)*(pow(this->v,2)/(pPpt->gammaAir()*pPpt->R_air*this->Ts)) ),
//												-pPpt->gammaAir()/(pPpt->gammaAir()-1));
//		}
//
////cout << p_des_app_seen[i_fill][TIME_IN_PERIOD] << "\t" << p_des_app_seen[i_fill][DESIRED_P] << endl;

	}

	// Initialise p_wave_prev to the first value in the array, else get funny values on first pulse
	MFR_wave = p_des_app_seen[0][DESIRED_MFR];
	MFR_wave_prev = p_des_app_seen[0][DESIRED_MFR];
	//cout << "MFR_wave_prev = " << MFR_wave_prev << endl;
}
*/
void CTransmissive::LoadSteady(char* InputFile)
//--------------------------------------------------//
// Transmissive boundary							//
// -----------------------------------				//
// Loads steady PR-MFP characteristic	 			//
//													//
//--------------------------------------------------//
{
	int c = 0;
	double temp;
	int num_columns = 2;
	int col, row;
	int datapoints = 0;
	
	FILE* stream;
	stream = fopen(InputFile, "r");

	if(stream == NULL)
	{
		printf("Error opening steady characteristic file\n");
		exit(1);
	}
	else
	{
		fseek(stream, 0L, SEEK_SET);	// Set pointer to beginning of file
		do
		{
			fscanf(stream, "%d", &c);	// Scans for the line number (ints), puts it in c
			//cout << "c = " << c << endl;
			if(c>datapoints) datapoints = c;
			fscanf(stream, "%lf", &temp);		// Runs over the time value
			fscanf(stream, "%lf", &temp);		// Runs over the parameter value
			//cin >> pause;
		}while(fscanf(stream, "%l")!=EOF);
		numsteadypts = datapoints;
//cout << "Raw file datapoints = " << datapoints << endl;
//cin >> pause;
		// Can now dimension timing arrays
		steadyPRMFP = new double* [num_columns]; 
		for(col=0; col<num_columns; ++col) steadyPRMFP[col] = new double [datapoints];

		fseek(stream, 0L, SEEK_SET);	// Reset pointer to beginning of file

		for (row=0; row<datapoints; row++)
		{
			fscanf(stream, "%d", &c);	// Scans past the line number
			for (col=0; col<num_columns; col++)
			{
				// Finds the time and parameter values
				fscanf(stream, "%lf", &steadyPRMFP[col][row]);	
				//cout << "steadyPRMFP[col][row] = " << steadyPRMFP[col][row] << endl;
			}
//cout << "steadyPRMFP[PR][" << row << "] = " << steadyPRMFP[0][row] << "\t";
//cout << "steadyPRMFP[MFP][" << row << "] = " << steadyPRMFP[1][row] << endl;
		}
	}
	fclose(stream);
}

double CTransmissive::InterpolateSteady(double pr)
{
	int point = 0;
	int PRlabel = 0; int MFPlabel = 1;

	// Find the two points on the steady curve that pr lies between and interpolate
	while(steadyPRMFP[PRlabel][point] < pr && point<numsteadypts-1) ++point;
	// pr lies between point and point-1
/*
	// --------------------
	// Linear interpolation
	if(point==0) point=1; // If pr out of range, interpolate using first two points
	
	// Interpolate these points
	return
		steadyPRMFP[MFPlabel][point-1]
		+ 
		((pr - steadyPRMFP[PRlabel][point-1])/(steadyPRMFP[PRlabel][point] - steadyPRMFP[PRlabel][point-1])
		*(steadyPRMFP[MFPlabel][point] - steadyPRMFP[MFPlabel][point-1]));
	// --------------------
*/
	// --------------------
	// Quadratic interpolation
	// Select three most appropriate points and interpolate quadratically; y = Ax^2 + Bx + C
	int point1, point2, point3;
	if(point==0)
	{
		if(numsteadypts>=3){point1 = point; point2 = point + 1; point3 = point + 2;}
		else 
		{
			if(numsteadypts>=2){point1 = point; point2 = point + 1; point3 = point2;}
			else {point1 = point; point2 = point1; point3 = point1;}
		}
	}
	else
	{
		if(point==numsteadypts-1)
		{
			if(numsteadypts>=3){point1 = point; point2 = point - 1; point3 = point - 2;}
			else 
			{
				if(numsteadypts>=2){point1 = point; point2 = point - 1; point3 = point2;}
				else {point1 = point; point2 = point1; point3 = point1;}
			}
		}
		else {point1 = point; point2 = point - 1; point3 = point + 1;}
	}

	double A, B, C, x1, x2, x3, y1, y2, y3;
	x1 = steadyPRMFP[PRlabel][point1]; y1 = steadyPRMFP[MFPlabel][point1];
	x2 = steadyPRMFP[PRlabel][point2]; y2 = steadyPRMFP[MFPlabel][point2];
	x3 = steadyPRMFP[PRlabel][point3]; y3 = steadyPRMFP[MFPlabel][point3];
	A = ((y3-y2)/((x3-x2)*(x3-x1))) - ((y1-y2)/((x1-x2)*(x3-x1)));
	B = ((y1-y2) + A*(pow(x2,2)-pow(x1,2)))/(x1-x2);//	B = ((y3-y2) + A*(pow(x2,2)-pow(x3,2)))/(x2-x3);
	C = y2 - A*pow(x2,2) - B*x2;//	C = y1 - A*pow(x1,2) - B*x1;//	C = y3 - A*pow(x3,2) - B*x3;

	return A*pow(pr,2) + B*pr + C;
	// --------------------
}

void CTransmissive::ReadFourierFile(CProperties* pPpt, char* InputFile)
//--------------------------------------------------//
// Transmissive boundary							//
// -----------------------------------				//
// Loads steady PR-MFP characteristic	 			//
//													//
//--------------------------------------------------//
{
	int c = 0;
	//float temp;
	double temp;
	int num_columns = 4; // Frequency, Amplitude, an, bn
	int col, row;
	int datapoints = 0;
	
	FILE* stream;
	stream = fopen(InputFile, "r");

	if(stream == NULL)
	{
		printf("Error opening Fourier file\n");
		exit(1);
	}
	else
	{
		fseek(stream, 0L, SEEK_SET);	// Set pointer to beginning of file
		do
		{
			fscanf(stream, "%d", &c);	// Scans for the line number (ints), puts it in c
//cout << "c = " << c << endl;
			if(c>datapoints) datapoints = c;
			fscanf(stream, "%lf", &temp);		// Runs over the frequency value
			fscanf(stream, "%lf", &temp);		// Runs over the amplitude value
			fscanf(stream, "%lf", &temp);		// Runs over the an value
			fscanf(stream, "%lf", &temp);		// Runs over the bn value
//cin >> pause;
		}while(fscanf(stream, "%l")!=EOF);
		ncomponents = datapoints;
//cout << "Raw file datapoints = " << datapoints << endl;
//cin >> pause;

		// Can now dimension timing arrays
		fr_An_an_bn = new double* [num_columns]; 
		for(col=0; col<num_columns; ++col) fr_An_an_bn[col] = new double [datapoints];

		fseek(stream, 0L, SEEK_SET);	// Reset pointer to beginning of file

		for (row=0; row<datapoints; row++)
		{
			fscanf(stream, "%d", &c);	// Scans past the line number
			for (col=0; col<num_columns; col++)
			{
				// Finds the time and parameter values
				fscanf(stream, "%lf", &fr_An_an_bn[col][row]);	
				//cout << "fr_An_an_bn[col][row] = " << fr_An_an_bn[col][row] << endl;
			}
///*
			pPpt->Out("fr_An_an_bn[fr]["); pPpt->Out(row); pPpt->Out("] = "); pPpt->Out(fr_An_an_bn[0][row]); pPpt->Out("\t");
			pPpt->Out("fr_An_an_bn[An]["); pPpt->Out(row); pPpt->Out("] = "); pPpt->Out(fr_An_an_bn[1][row]); pPpt->Out("\t");
			pPpt->Out("fr_An_an_bn[an]["); pPpt->Out(row); pPpt->Out("] = "); pPpt->Out(fr_An_an_bn[2][row]); pPpt->Out("\t");
			pPpt->Out("fr_An_an_bn[bn]["); pPpt->Out(row); pPpt->Out("] = "); pPpt->Out(fr_An_an_bn[3][row]); pPpt->Out("\n");
//*/
		}
		pPpt->Out("\n");
	}
	fclose(stream);
//cin >> pause;
}

void CTransmissive::ListProperties(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ListProperties\n");}

	// ====================================================================================================
	// Parameter file for Exhaust Transmissive [0]
	// ====================================================================================================

	pPpt->Out(Underline(Identify(), "=", "\t", strDesc)); pPpt->Out("\n");

	if(ANECHOIC_ONLY)
	{
		pPpt->Out("\tFunction as anechoic termination, ANECHOIC_ONLY\t=\t");
		pPpt->Out(TrueOrFalse(ANECHOIC_ONLY));
		pPpt->Out("\n");
	}
	else
	{
		if(TOTAL_PRESS)
		{
			pPpt->Out("\tDesired pressure is total, TOTAL_PRESS\t\t=\t");
			pPpt->Out(TrueOrFalse(TOTAL_PRESS));
			pPpt->Out("\n");
		}
		else 
		{
			pPpt->Out("\tDesired pressure is static, TOTAL_PRESS\t\t=\t");
			pPpt->Out(TrueOrFalse(TOTAL_PRESS));
			pPpt->Out("\n");
		}
		pPpt->Out("\n");
		if(CONSTANT)
		{
			// If constant conditions, i.e., CONSTANT == 1 == true
			// ----------------------------------------------------------------------------------------------------
			pPpt->Out(Underline("Constant setup", "-", "\t"));
			pPpt->Out("\tDesired pressure, constantp\t\t\t=\t"); pPpt->Out(constantp); pPpt->Out(" bar\n");
			pPpt->Out("\tDesired temperature, constantT\t\t\t=\t"); pPpt->Out(constantT); pPpt->Out(" K\n");
			pPpt->Out("\tPeriod across which to evaluate a match, period\t=\t"); pPpt->Out(period); pPpt->Out(" s\n");
		}
		else
		{
			// Else varying conditions, i.e., CONSTANT == 0 == false
			// ----------------------------------------------------------------------------------------------------

			// Waveform setup
			// ----------------------------------------------------------------------------------------------------
			pPpt->Out(Underline("Waveform setup", "-", "\t"));
			pPpt->Out("\t");
			if(shape==0) pPpt->Out("Triangular waveform, shape\t\t\t=\t");
			else
			{
				if(shape==1) pPpt->Out("Square waveform, shape\t\t\t=\t");
				else
				{
					if(shape==2) pPpt->Out("Cosine waveform, shape\t\t\t\t=\t");
					else
					{
						if(shape==3)
						{ 
							pPpt->Out("Waveform from file, shape\t\t\t=\t");
						}
						else
						{
							if(shape==4)
							{
								//pPpt->Out("Fourier series; frequency components read from "); pPpt->Out(FOURIER_FILE); pPpt->Out("\n");
								pPpt->Out("Fourier series waveform, shape\t\t\t=\t");
							}
							else pPpt->Out("Unknown waveform, shape\t\t\t=\t");
						}
					}	
				}
			}
			pPpt->Out(shape); pPpt->Out("\n");
			pPpt->Out("\tOscillation frequency, f\t\t\t=\t"); pPpt->Out(f); pPpt->Out(" Hz\n");
			pPpt->Out("\tPulse length as fraction of wavelength, phi\t=\t"); pPpt->Out(phi); pPpt->Out("\n");
			pPpt->Out("\tOrig. pulse length fraction, phi_orig\t\t=\t"); pPpt->Out(phi_orig); pPpt->Out("\n");
			if(shape!=3)
			{
				pPpt->Out("\tStatic temperature at boundary, Ts\t\t=\t"); pPpt->Out(Ts); pPpt->Out(" K\n");
				pPpt->Out("\tVelocity at boundary, v\t\t\t\t=\t"); pPpt->Out(v); pPpt->Out(" m/s\n");
			}
			pPpt->Out("\tCharacteristic length, L\t\t\t=\t"); pPpt->Out(pPipe[ONE_SIDE]->eff_length); pPpt->Out(" m\n");
			pPpt->Out("\n");

			// Amplitudes
			// ----------------------------------------------------------------------------------------------------
			if(shape==1 || shape==2)
			{
				pPpt->Out(Underline("Amplitudes", "-", "\t"));
				
				pPpt->Out("\tpress_low\t=\t"); pPpt->Out(press_low); pPpt->Out(" bar\n");
				pPpt->Out("\tpress_high\t=\t"); pPpt->Out(press_high); pPpt->Out(" bar\n");
				pPpt->Out("\tpress_min\t=\t"); pPpt->Out(press_min); pPpt->Out(" bar\n");
				pPpt->Out("\n");
			
				pPpt->Out("\tTsm_low\t=\t\t"); pPpt->Out(Tsm_low); pPpt->Out(" K\n");
				pPpt->Out("\tTsm_high\t=\t"); pPpt->Out(Tsm_high); pPpt->Out(" K\n");
				pPpt->Out("\tTsm_min\t=\t\t"); pPpt->Out(Tsm_min); pPpt->Out(" K\n");
				pPpt->Out("\n");
			
				pPpt->Out("\tT0m_low\t=\t\t"); pPpt->Out(T0m_low); pPpt->Out(" K\n");
				pPpt->Out("\tT0m_high\t=\t"); pPpt->Out(T0m_high); pPpt->Out(" K\n");
				pPpt->Out("\tT0m_min\t=\t\t"); pPpt->Out(T0m_min); pPpt->Out(" K\n");
				pPpt->Out("\n");
		
				pPpt->Out("\tMFR_low\t=\t\t"); pPpt->Out(MFR_low); pPpt->Out(" Kg/s\n");
				pPpt->Out("\tMFR_high\t=\t"); pPpt->Out(MFR_high); pPpt->Out(" Kg/s\n");
				pPpt->Out("\tMFR_min\t=\t\t"); pPpt->Out(MFR_min); pPpt->Out(" Kg/s\n");
				pPpt->Out("\n");
			}

			if(shape==3)
			{
				// Waveform from file (shape==3)
				// ----------------------------------------------------------------------------------------------------
				pPpt->Out(Underline("Waveform from file (shape==3)", "-", "\t"));
				pPpt->Out("\tStatic pressure file, FILE_ps\t\t\t=\t"); pPpt->Out(FILE_ps); pPpt->Out("\n");
				if(USE_FILE_Ts)
				{
					//pPpt->Out("\twaveform total temperature read from "); pPpt->Out(WAVE_FILE_T0m); pPpt->Out("\n");
					pPpt->Out("\tStatic temperature file, FILE_Ts\t\t=\t"); pPpt->Out(FILE_Ts); pPpt->Out("\n");
				}
				else
				{
					pPpt->Out("\tUse constant static temperature, Ts\t\t=\t"); pPpt->Out(Ts); pPpt->Out(" K\n");
				}
				pPpt->Out("\n");

				if(ADJUST_TIME || ADJUST_PRESSURE || ADJUST_TSM || ADJUST_T0M || ADJUST_MFR) pPpt->Out(Underline("Adjustments", "-", "\t"));
				if(ADJUST_TIME)
				{
					pPpt->Out("\tAdjust timing to fit pulse period, ADJUST_TIME\t=\t"); pPpt->Out(TrueOrFalse(ADJUST_TIME)); pPpt->Out("\n");
				}
				if(ADJUST_PRESSURE)
				{
					pPpt->Out("\tAdjusting file pressures to:\n");
					pPpt->Out("\tpress_low\t=\t"); pPpt->Out(press_low); pPpt->Out(" bar\n");
					pPpt->Out("\tpress_high\t=\t"); pPpt->Out(press_high); pPpt->Out(" bar\n");
					pPpt->Out("\tpress_min\t=\t"); pPpt->Out(press_min); pPpt->Out(" bar\n");
					pPpt->Out("\n");
				}
				if(ADJUST_TSM) 
				{
					pPpt->Out("\tAdjusting file Ts values to:\n");
					pPpt->Out("\tTsm_low\t=\t"); pPpt->Out(Tsm_low); pPpt->Out(" K\n");
					pPpt->Out("\tTsm_high\t=\t"); pPpt->Out(Tsm_high); pPpt->Out(" K\n");
					pPpt->Out("\tTsm_min\t=\t"); pPpt->Out(Tsm_min); pPpt->Out(" K\n");
					pPpt->Out("\n");
				}
				if(ADJUST_T0M)
				{
					pPpt->Out("\tAdjusting file T0 values to:\n");
					pPpt->Out("\tT0m_low\t=\t"); pPpt->Out(T0m_low); pPpt->Out(" K\n");
					pPpt->Out("\tT0m_high\t=\t"); pPpt->Out(T0m_high); pPpt->Out(" K\n");
					pPpt->Out("\tT0m_min\t=\t"); pPpt->Out(T0m_min); pPpt->Out(" K\n");
					pPpt->Out("\n");
				}
				if(ADJUST_MFR)
				{
					pPpt->Out("\tAdjusting MFR to:\n");
					pPpt->Out("\tMFR_low\t=\t"); pPpt->Out(MFR_low); pPpt->Out(" Kg/s\n");
					pPpt->Out("\tMFR_high\t=\t"); pPpt->Out(MFR_high); pPpt->Out(" Kg/s\n");
					pPpt->Out("\tMFR_min\t=\t"); pPpt->Out(MFR_min); pPpt->Out(" Kg/s\n");
					pPpt->Out("\n");
				}
			}
/*
			pPpt->Out("\tshape\t\t=\twaveform MFR read from "); pPpt->Out(WAVE_FILE_MFR);
			if(ADJUST_TIME) pPpt->Out(", adjusted timings");
			else pPpt->Out(", unadjusted timings");
			if(ADJUST_MFR) pPpt->Out(", adjusted MFR\n");
			else pPpt->Out(", unadjusted MFR\n");
			
			if(TEST)
			{
				pPpt->Out("\tLet u float, TEST\t\t\t\t=\t"); pPpt->Out(TrueOrFalse(TEST)); pPpt->Out("\n");
			}
			else
			{
				pPpt->Out("\tu derived from Tsm, T0m, TEST\t\t\t=\t"); pPpt->Out(TrueOrFalse(TEST)); pPpt->Out("\n");
			}
			pPpt->Out("\twaveform velocity read from "); pPpt->Out(WAVE_FILE_v); pPpt->Out("\n");
*/			
		}
		//pPpt->Out("\n");

		// Pipe configuration
		// ----------------------------------------------------------------------------------------------------
		pPpt->Out(Underline("Pipe configuration", "-", "\t"));
		pPpt->Out("\tRotor exit tapping located in pipe, p2_pipe\t=\t");
		if(this->EX)
		{
			pPpt->Out("Exhaust["); pPpt->Out(p2_pipe); pPpt->Out("]\n");
		}
		else
		{
			pPpt->Out("Intake["); pPpt->Out(p2_pipe); pPpt->Out("]\n");
		}
		pPpt->Out("\t- at a length fraction of, p2_loc\t\t=\t"); pPpt->Out(p2_loc); pPpt->Out("\n");
		pPpt->Out("\n");
		
		// Matching
		// ----------------------------------------------------------------------------------------------------
		pPpt->Out(Underline("Matching", "-", "\t"));
		if(INTELLIGENT) pPpt->Out("\tConverging onto desired pulse, INTELLIGENT\t=\t");
		else pPpt->Out("\tNot converging onto desired pulse, INTELLIGENT\t=\t");
		pPpt->Out(TrueOrFalse(INTELLIGENT)); pPpt->Out("\n");
		pPpt->Out("\tNumber of points sampled per pulse, nsamples\t=\t"); pPpt->Out(nsamples); pPpt->Out("\n");
		pPpt->Out("\tPulse match tolerance, tol\t\t\t=\t"); pPpt->Out(tol); pPpt->Out("%\n");
		if(SUSPEND) pPpt->Out("\tSuspend simulation once matched, SUSPEND\t=\t");
		else		pPpt->Out("\tDo not suspend simulation once matched, SUSPEND\t=\t");
		pPpt->Out(TrueOrFalse(SUSPEND)); pPpt->Out("\n");
    pPpt->Out("\tCycles between turning on/off, nMatchCycles\t=\t"); pPpt->Out(nMatchCycles); pPpt->Out("\n");
		pPpt->Out("\n");

//	// Now calculipe->eff_length/this->v)*(1/(2*this->phi));
//	this->pmst = ((2*PI*this->f)*this->pPipe->eff_length/this->v)*(1/(2*this->phi));
//	cout << "\tmst\t\t=\t" << mst << "\n";
//	cout << "\tpmst\t\t=\t" << pmst << "\n";

		// Measurements
		// ----------------------------------------------------------------------------------------------------
		pPpt->Out(Underline("Measurements", "-", "\t"));
		if(USE_DEF_FREQ)
		{
			if(freq==1) pPpt->Out("\tUsing default sampling rate (pPpt->freq)\t=\tonce per timestep\n");
			else
			{
				pPpt->Out("\tUsing default sampling rate (pPpt->freq)\t=\tonce per "); pPpt->Out(freq); pPpt->Out(" timesteps\n");
			}
		}
		else
		{
			if(freq==1) pPpt->Out("\tUsing local sampling rate (freq)\t\t=\tonce per timestep\n");
			pPpt->Out("\tUsing local sampling rate (freq)\t\t=\tonce per "); pPpt->Out(freq); pPpt->Out(" timesteps\n");
		}
		pPpt->Out("\tRecord onwards from pulse no., print_from_pulse\t=\t"); pPpt->Out(print_from_pulse); pPpt->Out("\n");
		pPpt->Out("\tRecord onwards from, print_from_time\t\t=\t"); pPpt->Out(print_from_time); pPpt->Out(" s\n");
		
		if(PRINT_PULSE_INFO) pPpt->Out("\tIncluding pulse data in file, PRINT_PULSE_INFO\t=\t"); else pPpt->Out("\tNot printing pulse data to file, PRINT_PULSE_INFO\t=\t");
		pPpt->Out(TrueOrFalse(PRINT_PULSE_INFO)); pPpt->Out("\n");

		if(PRINT_MOVIE_FILE) pPpt->Out("\tPrinting movie file, PRINT_MOVIE_FILE\t\t=\t"); else pPpt->Out("\tNot printing movie file, PRINT_MOVIE_FILE\t=\t"); 
		pPpt->Out(TrueOrFalse(PRINT_MOVIE_FILE)); pPpt->Out("\n");

		if(PRINT_TO_SCREEN) pPpt->Out("\tPrinting data to screen, PRINT_TO_SCREEN\t=\t"); else pPpt->Out("\tNot printing data to screen, PRINT_TO_SCREEN\t=\t");
		pPpt->Out(TrueOrFalse(PRINT_TO_SCREEN)); pPpt->Out("\n");
		pPpt->Out("\n");
		
		// Locations (as a fraction of the attached pipe's physical length, measured from the odd end)
		// ----------------------------------------------------------------------------------------------------
		pPpt->Out(Underline("Locations", "-", "\t"));
		pPpt->Out("\tNo. measuring locations (nmeasurements)\t\t=\t"); pPpt->Out(nmeasurements); pPpt->Out("\n");
		if(nmeasurements>0)
		{
			pPpt->Out("\n\tAttached pipe measuring locations\n");
			pPpt->Out("\t---------------------------------\n");
			pPpt->Out("\tAs a fraction of attached pipe physical length, measured from the odd end:\n");
		}
		for(int m=0; m<nmeasurements; ++m)
		{
			pPpt->Out("\t Location "); pPpt->Out(m); pPpt->Out(" = "); pPpt->Out(loc_measure[m]); pPpt->Out("\n");
		}
		if(num_props_measured>0)
		{
			pPpt->Out("\n\tMax. data points per location, max_pts\t=\t"); pPpt->Out(max_pts); pPpt->Out("\n");
			pPpt->Out("\tRecording:\n");
			if(FLOW_VELOCITY) pPpt->Out("\t Flow velocity\n");
			if(PRESSURE_VELOCITY) pPpt->Out("\t Pressure wave velocity\n");
//			if(STATIC_PRESSURE) pPpt->Out("\t Static pressure\n");
//			if(TEMPERATURE) pPpt->Out("\t Temperature\n");
			if(MASS_FLOW_RATE) pPpt->Out("\t Mass flow rate\n");
			if(MASS_FLOW_PARAMETER) pPpt->Out("\t MASS_FLOW_PARAMETER\n");
			if(PRESSURE_RATIO) pPpt->Out("\t PRESSURE_RATIO\n");		
//			if(REYNOLDS_NO) pPpt->Out("\t Reynold's No.\n");
		}
		pPpt->Out("\n");

		//if(CALIBRATE)
		//{
		//	pPpt->Out(Underline("Calibration of associated APLDEv objects", "-", "\t"));
		//	pPpt->Out("\tSteady PR-MFP characteristic read from\t\t=\t"); pPpt->Out(STEADY_FILE); pPpt->Out("\n");
		//	pPpt->Out("\tLower end of calibration PR range (cal_PR_low)\t=\t"); pPpt->Out(cal_PR_low); pPpt->Out("\n");
		//	pPpt->Out("\tUpper end of calibration PR range (cal_PR_high)\t=\t"); pPpt->Out(cal_PR_high); pPpt->Out("\n");
		//	pPpt->Out("\tNo. of points to calibrate against (cal_points)\t=\t"); pPpt->Out(cal_points); pPpt->Out("\n");
		//	pPpt->Out("\tTolerance when matching desired MFP (tol_MFP)\t=\t"); pPpt->Out(tol_MFP); pPpt->Out(" (fraction)\n");
		//}
	}

	pPpt->Out("\n");

	// ====================================================================================================
	// End of file
	// ====================================================================================================
}