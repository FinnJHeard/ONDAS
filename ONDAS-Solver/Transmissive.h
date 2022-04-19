// Transmissive.h: interface for the CTransmissive class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TRANSMISSIVE_H__61F3FD35_3221_4F8D_9B30_CBF782333F2F__INCLUDED_)
#define AFX_TRANSMISSIVE_H__61F3FD35_3221_4F8D_9B30_CBF782333F2F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Boundary.h"
#include "PathLine.h"
#include "Pipe.h"
#include "Properties.h"

class CTransmissive : public CBoundary  
{
public:
	CTransmissive();
	virtual ~CTransmissive();
	CTransmissive(const CTransmissive& inTransmissive);	// Copy constructor
	CTransmissive& operator=(const CTransmissive& inTransmissive); // Overload operator=

	void Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rTRANSPIPES, int** &rTRANSPIPES_ENDS, double* &rENDCORR, 
					int id, bool ex, int npipes, std::string param_dir, int assyid, string parent_assy_res_dir, int ntransm, int nPipesInAssy, string calling_object_str);
	void ReadInput(CProperties* pPpt, char *InputFile);
	void ListProperties(CProperties* pPpt);
	
	void RunBoundary(CProperties* pPpt, CTime &rMyTime, double DELZe, double DELZi, double ca, double del_ca, double TIMEe, double TIMEi, int timestep, bool STEADY);	
	void PreProcess(CTime &rMyTime, double DELZe, double DELZi, CProperties* pPpt, double ca, double del_ca, double TIMEe, double TIMEi, int timestep, bool STEADY);
	void Transmissive(CTime &rMyTime, double DELZe, double DELZi, CProperties* pPpt, double ca, double del_ca, double TIMEe, double TIMEi, int timestep, bool STEADY);
	void Loss(CProperties* pPpt, double &rlambda_in_c2, double &rAA_c2, double &rlambda_out1, double &rlambda_out2, double &rCLIN_WAVEFORM, double &rAA_WAVEFORM, int &rpipe_flow, int timestep);

	void SetupFiles(CProperties* pPpt);
	void WriteFiles();
	void CloseFiles();
	void RecordTappings(CProperties* pPpt, int timestep, double time, CPipe* Pipe);
	void PrintToFile(CProperties* pPpt, int timestep, double time, double ca);
	void PrintToScreen(CProperties* pPpt);

	//void LoadWaveform(char* InputFile, CProperties* pPpt);
	void LoadWaveform(char* InputFile, CProperties* pPpt, 
						bool TOTAL, bool ADJUST_TIME, bool ADJUST_PARAMETER, 
						double parameter_high, double parameter_low, double parameter_min, 
						int PARAMETER_DESIRED, 
						double&rparam_wave, double &rparam_wave_prev);
	//void LoadWaveformMFR(char* InputFile, CProperties* pPpt);
	void LoadSteady(char* InputFile);
	double InterpolateSteady(double pr);

	void ReadFourierFile(CProperties* pPpt, char* InputFile);

public:
	// ====================================================================================================
	// Parameter file for Exhaust Transmissive [0]
	// ====================================================================================================
			
	bool ANECHOIC_ONLY;				// Use transmissive boundary as anechoic termination only (1=true) (0=false)
	bool CONSTANT;					// Desired boundary conditions are constant (1==true) or varying (0==false)

		// If constant conditions, i.e., CONSTANT == 1 == true
		// ----------------------------------------------------------------------------------------------------
		double constantp;				// Desired total pressure (bar)
		double constantT;				// Desired total temperature (K)
		double period;					// Time period across which to evaluate a match (s)

		// Else varying conditions, i.e., CONSTANT == 0 == false
		// ----------------------------------------------------------------------------------------------------

		// Waveform setup
		// ----------------------------------------------------------------------------------------------------
		int shape;						// Waveform: triangular(0), square(1), cosine(2), file(3), Fourier series(4)
		double f;						// Oscillation frequency (Hz)
		double phi;						// Pulse length as fraction of wavelength (0-1)
		double phi_orig;				// phi(0-1) based on original excitation
		double v;						// Velocity at boundary (m/s)
		double Ts;						// Static temperature at boundary (K)

		// Amplitudes
		// ----------------------------------------------------------------------------------------------------
		double press_low;				// Waveform low pressure (static, bar)
		double press_high;				// Waveform high pressure (static, bar)
		double press_min;				// Minimum waveform pressure allowed (static, bar)

		double Tsm_low;					// Waveform low Tsm (K)
		double Tsm_high;				// Waveform high Tsm (K)
		double Tsm_min;					// Minimum waveform Tsm allowed (K)

		double T0m_low;					// Waveform low T0m (K)
		double T0m_high;				// Waveform high T0m (K)
		double T0m_min;					// Minimum waveform T0m allowed (K)

		double MFR_low;					// Waveform low MFR (Kg/s)
		double MFR_high;				// Waveform high MFR (Kg/s)
		double MFR_min;					// Minimum waveform MFR allowed (Kg/s)

			// Pulse data from files (shape==3)
			// ----------------------------------------------------------------------------------------------------
			char* FILE_ps;					// File containing static pressure data (line number, angle, value)
			bool TOTAL_PRESS;				// Is file data total (1==true) or static (0==false) pressure? HAS NO EFFECT?!?!?!?
			bool USE_FILE_Ts;				// Use temperature data from file (1==true) or use constant Ts value above (0==false)

				// Temperature data from file (USE_FILE_Ts==1==true)
				// ----------------------------------------------------------------------------------------------------
				char* FILE_Ts;					// File containing static temperature data (line number, angle, value)
				char* WAVE_FILE_T0m;			// Total temperature
				char* WAVE_FILE_v;				// Velocity
				char* WAVE_FILE_MFR;			// Mass flow rate

			bool ADJUST_TIME;				// Adjust file timing values to conform to desired pulse period (1==true) or use file as is (0==false)
			bool ADJUST_PRESSURE;			// Adjust file pressure values to conform to the above (1==true) or use file as is (0==false)
			bool ADJUST_TSM;				// Adjust file Tsm values to conform to the above (1==true) or use file as is (0==false)
			bool ADJUST_T0M;				// Adjust file T0m values to conform to the above (1==true) or use file as is (0==false)
			bool ADJUST_MFR;				// Adjust file MFR values to conform to the above (1==true) or use file as is (0==false)

			// Fourier series (shape==4)
			// ----------------------------------------------------------------------------------------------------
			char* FOURIER_FILE;				// File containing list of frequency components, amplitudes and coefficients

	// Pipe configuration
	// ----------------------------------------------------------------------------------------------------
	int	nPipesInAssy;				// Number of pipes in corresponding assembly
	int p2_pipe;					// Exhaust or intake pipe number to place p2 when calculating PR
	double p2_loc;					// Location of p2 tapping along p2_pipe as a fraction of pipe length

	// Matching
	// ----------------------------------------------------------------------------------------------------
	bool INTELLIGENT;				// Converge (1==true) onto desired pulse or not (0==false) in general
	bool ATTEMPT_MATCH;				// Converge (1==true) onto desired pulse or not (0==false) this cycle
	int nsamples;					// The number of data points to record to represent a pulse
	double tol;						// The tolerance (%) within which a pulse match will occur
	bool SUSPEND;					// Suspend simulation when match is achieved (1==true) or not (0==false)
	int nMatchCycles;				// Number of cycles to wait between turning matching on/off
	int nTransm;					// Total number of transmissive boundaries
	
	// Measurements
	// ----------------------------------------------------------------------------------------------------
	bool USE_DEF_FREQ;				// Use default sampling rate (1==true) or 'freq' below (0==false)
	int freq;						// Print data to file every freq timesteps
	int print_from_pulse;			// Only start printing data from the start of this pulse onwards
	//NOW in boundary.h: double print_from_time;			// Only start printing data after this simulation time onwards (s)
	bool PRINT_PULSE_INFO;			// Include (true==1) end of pulse info in file or not (false==0)

	// Locations (as a fraction of the attached pipe's physical length, measured from the odd end)
	// ----------------------------------------------------------------------------------------------------
	int nmeasurements;				// Number of measuring locations per pipe; up to 11:
	int max_pts;					// Max. data points per location
	bool FLOW_VELOCITY;				// Record flow velocity of gas i.e., u (1=true) or not (0=false)
	bool PRESSURE_VELOCITY;			// Record pressure wave velocity i.e., u+/-a (1=true) or not (0=false)
	//bool STATIC_PRESSURE;			// Record static pressure (1=true) or not (0=false)
	//bool TEMPERATURE; 			// Record temperature (1=true) or not (0=false)
	bool MASS_FLOW_RATE;			// Record mass flow rate (1=true) or not (0=false)
	bool MASS_FLOW_PARAMETER;		// Record mass flow parameter, m_dot*sqrt(T01)/p01 (1=true) or not (0=false)
	bool PRESSURE_RATIO;			// Record pressure ratio, p01/p2 (1=true) or not (0=false)
	//bool REYNOLDS_NO;				// Record nodal Reynold's No. (1=true) or not (0=false)

	//// Calibration of APLDev loss curve; uses v and T0 above
	//// ----------------------------------------------------------------------------------------------------
	//bool CALIBRATE;					// Calibrate the loss curve during this simulaton (1=true) (0=false)
	//char* STEADY_FILE;				// File containing the steady PR-MFP characteristic
	//double cal_PR_low;				// Lower end of calibration pressure ratio range
	//double cal_PR_high;				// Upper end of calibration pressure ratio range
	//double cal_points;				// No. of even-spaced points in range at which to calibrate
	//double tol_MFP;					// Tolerance when attempting match the desired MFP (fraction)

// ====================================================================================================
// End of file
// ====================================================================================================

private:
	
	// To be sure of identifying the correct item in the two item configuration vectors
	int exhaust_side;
	int intake_side;

	double* M;

public:

	int INTERNAL;
	int WAVEFORM;
	int UPSTREAM;
	int DOWNSTREAM;
//	bool LEFT_TO_RIGHT;
	bool WAVEFORM_TO_INTERNAL;

	double A0test;
//	double TREF;		// Reference temperature (depends whether inflow is in exhaust or intake) (K)

	double u_WAVEFORM, T_WAVEFORM;	// Read-in velocity and temperature
	double U_WAVEFORM;

	//double period;					// The calculated time period of the pulse (s)
	double time_in_period;			// Time after start of current pulse (s)
	double time_in_period_prev;		// Previous iteration's time after start of current pulse (s)

	double p_wave;					// Calculated waveform pressure to be applied (bar)
	double p_wave_calibrate;		// The value of p_wave to be applied during calibration
	double p_wave_prev;				// Previous iteration's calculated waveform pressure (bar)
	double p_wave_desired;			// Instantaneous desired pressure (bar)
	
	double T0m_wave;				// Calculated waveform total temperature to be applied (K)
	double T0m_wave_prev;			// Previous iteration's calculated waveform total temperature (K)
	double T_wave;					// Calculated waveform static temperature to be applied (K)
	double T_wave_prev;				// Previous iteration's calculated waveform static temperature (K)
	double T_wave_desired;			// Instantaneous desired static temperature (K)
	
	double v_wave;					// Calculated waveform velocity to be applied (m/s)
	double v_wave_prev;				// Previous iteration's calculated waveform velocity (m/s)
	double v_wave_desired;			// Instantaneous desired velocity (m/s)

	double mfr_wave;				// Calculated waveform velocity to be applied (m/s)
	double mfr_wave_prev;			// Previous iteration's calculated waveform velocity (m/s)
	double mfr_wave_desired;		// Instantaneous desired velocity (m/s)
	//double MFR_wave;				// Calculated waveform pressure to be applied (bar)
	//double MFR_wave_prev;			// Previous iteration's calculated waveform pressure (bar)
	
	int i_prev;						// The index of the data point last filled in Pressure_Desired

	int whole_number_of_pulses;		// Counts completed pulses since start of simulation
	bool START_OF_NEW_PULSE;		// Set to true at start of new pulse
  bool START_OF_SIMULATION;   // Flag to indicate start of simulation
	int iterations_during_pulse;	// Counts number of iterations in order to calculate averages

	double* sum_flow_velocity_over_pulse; // Sum of velocities over pulse - one for each measuring location (m/s)
	double* average_flow_velocity_last_pulse; // Records aveage flow velocity over the last completed pulse (m/s)
	
	double* sum_press_velocity_over_pulse_uplusa; // Sum of velocities over pulse - one for each measuring location (m/s)
	double* average_press_velocity_last_pulse_uplusa; // Records aveage flow velocity over the last completed pulse (m/s)
	
	double* sum_press_velocity_over_pulse_uminusa; // Sum of velocities over pulse - one for each measuring location (m/s)
	double* average_press_velocity_last_pulse_uminusa; // Records aveage flow velocity over the last completed pulse (m/s)
	
	double* st;					// Strouhal number at each location
	double* beta;				// Reduced frequency at each location
	double* mst;					// Waveform modified Strouhal number at each location
	double* pmst_uplusa;			// Waveform pressure modified Strouhal number at each location
	double* pmst_uminusa;			// Waveform pressure modified Strouhal number at each location

	double Pressure_0;
	double Pressure_0_prev;
	double MFR_0;
	double MFR_0_prev;
	double VEL_0;
	double VEL_0_prev;
	double vel_grad; // Previously recorded velocity gradient for the current time in the pulse
	double mfr_under_est;
	double P_VEL_0;
	double P_VEL_0_prev;
	double TEMP_0;
	double TEMP_0_prev;

	// Applied and seen pressure - arrays contain data for one pulse period
	double** p_des_app_seen;				// The interpolated (pulse from file) or calculated (waveform pulse) 
									// pressure to be achieved, and the applied pressure
	// Labels
	int TIME;
	int TIME_IN_PERIOD;

	int DESIRED_P;
	int APPLIED_P;
	int SEEN_P;
	int NEXT_P;
	int THIS_P;
	
	int DESIRED_T;
	int DESIRED_T0;
	int APPLIED_T;
	int SEEN_T;
	int NEXT_T;
	int THIS_T;

	int DESIRED_V;
	int APPLIED_V;
	int SEEN_V;
	int NEXT_V;
	int THIS_V;

	int DESIRED_MFR;
	int APPLIED_MFR;
	int SEEN_MFR;
	int NEXT_MFR;
	int THIS_MFR;
	
	int VEL_GRAD;
	int P_VEL;
	int P_VEL_GRAD;
	int UNDER_EST;  // MFR underestimate
	int UNDER_EST_PREV;  // MFR underestimate from previous cycle
	int UNDER_EST_FACT;  // factor to multiply the underestimate by to get loss value
	int K_LOSS;  // loss value
	int WAVE_TIME;
	int WAVE_PARAMETER;
	int WAVE_MFR;

	

		// Fourier series (shape==4)
		// ----------------------------------------------------------------------------------------------------
		double** fr_An_an_bn;			// Array containing frequency component data
		int ncomponents;				// Number of frequency components (including dc) conatined in FOURIER_FILE and fr_An_an_bn
	
	// Matching
	// ----------------------------------------------------------------------------------------------------
	bool MATCH;						// Does the current cycle achieve a pulse or not
	bool MATCHED;					// Set to true once a full pulse match is achieved
	int datapoints_parameter;		// Number of rows in waveform pressure array
	int datapoints_MFR;				// Number of rows in waveform MFR array
	double** parameter_raw;			// Array containing pressure wave raw data; prior to adjustement calculations
	//double** mfr_raw;				// Array containing MFR wave raw data; prior to adjustement calculations
	double** parameter_wave;			// Array containing pressure timing; 2 columns (time then pressure)
	//double** mfr_wave;				// Array containing MFR timing; 2 columns (time then MFR) --- MOVED INTO GENERAL LOAD WAVEFORM FUNCTION?
	int rowcounter;
	int sample_factor;

	bool p_SEEN_LESS_THAN_DESIRED;
	bool v_DESIRED_GREATER_THAN_SEEN;
	double del_p;
	double del_v;

	// Measurements
	// ----------------------------------------------------------------------------------------------------
	double* Measurements;
	CNode* MeasureNode;
	FILE *FILE_LOC;					// File containing data at specified measuring locations in attached pipe
	FILE *FILE_p_des_app_seen;
	bool PRINT_HEADERS;				// Have column headers been printed

	double*** Results;
	double* loc_measure;			// List of normalized positions along attached pipe of measuring location
	int num_props_measured;			// Number of properties to measure, not counting time
	
	// Locations
	// --------------------------------------------------------------------------------		
	double* ps_loc;					// Static pressure at each measuring location (bar)
	double* Ts_loc;					// Static temperature at each measuring location (K)
	double* Velocity;				// Velocity at each measuring location (m/s)
	double* p0_loc;					// Total pressure at each measuring location (bar)
	double* T0_loc;					// Total temperature at each measuring location (K)
	double* MassFlowRate;			// Mass flow rate at each location (kg/s)
	double* PR;						// PR at each location
	double* MFP;					// MFP at each location
	double* Re;						// Re at each location

	// Calibration of APLDev loss curve
	// --------------------------------------------------------------------------------
	double** steadyPRMFP;			// Array containing read-in steady PR-MFP data
	int numsteadypts;				// Number of data points in steady curve data
	bool cal_MATCH;					// Does the seen pressure match the desired
	int match_counter;				// Counts how many of cal_points have been matched
	int WAIT_COUNTER;				// Enforces a gap between steady tests
	double MFP_des;					// Desired MFP
	bool CHANGE_PR;					// Flag to instruct transmissive boundary to adjust PR
	bool INCREASE_PR;				// Flag to instruct tranmsissive boundary to increase (1=true) or decrease (0=false) PR
	double del_change_PR;			// Pressure increment by which to adjust p_wave (bar)
	bool SWITCH;
	double p_wave_prev_prev;
};

#endif // !defined(AFX_TRANSMISSIVE_H__61F3FD35_3221_4F8D_9B30_CBF782333F2F__INCLUDED_)
