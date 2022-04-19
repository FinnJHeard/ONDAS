// APLDev.cpp: implementation of the CAPLDev class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "APLDev.h"
#include "Tools.h"
#include "time.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CAPLDev::CAPLDev()
{

}

CAPLDev::~CAPLDev()
{

}

void CAPLDev::Initialise(CTime* pMyTime, CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rAPLDevPIPES, int** &rAPLDevPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, CTransmissive* TransmPtr, std::string param_dir, int assyid, bool transm_ptr_exists, string parent_assy_res_dir, CEndEnvironment* &rEndEnv, CTransmissive* &rTransm, CJunction* &rJunc, string calling_object_str)
{
	if (pPpt->SHOW_calls) { pPpt->Out("CAPLDev.Initialise ("); pPpt->Out(calling_object_str); pPpt->Out(")\n"); }

	InitialiseGen(pPpt, pPipes, rPipe, rAPLDevPIPES, rAPLDevPIPES_ENDS, rENDCORR, id, ex, npipes, assyid, calling_object_str, parent_assy_res_dir);

	// Join transmissive pointer
	pTransm = TransmPtr;
	TRANSM_PTR_EXISTS = transm_ptr_exists;
	REFINING = false;
	LAST_RUN_LESS = false;
	delK = 0.1;//1;

	// Boundary name
	NAME = new int [NPIPES];
	NAME[LEFT_SIDE] = APLDEV_LEFT;	
	NAME[RIGHT_SIDE] = APLDEV_RIGHT;

	// Read input file
	std::string bcname_str = "APLD";
	
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Initialise\n");}

	ReadInput(pPpt, ConstructString(pPpt, param_dir, bcname_str, EX, ID));

	// Join pointers to entry boundary conditions (either transmissive or end environment) and rotor inlet nodes for each entry
	int e;
	pRotorInletNode = new CNode* [nEntries]; // Pointers to pipe nodes closest to rotor boundary in each limb
	if(transmEnt)
	{
		pEntryTransm = new CTransmissive* [nEntries];
		for(e=0; e<nEntries; ++e)
		{
			// Pointer to transmissive boundary condition at entry [e]
			pEntryTransm[e] = &(rTransm[entry[e]]);

			// Pointer to node closest to rotor boundary in entry limb [e] - below assumes even end of entry pipe
			pRotorInletNode[e] = &(pEntryTransm[e]->pPipe[ONE_SIDE]->Node[pEntryTransm[e]->pPipe[ONE_SIDE]->N-1]);
		} 
	}
	else
	{
		pEntryEndEnv = new CEndEnvironment* [nEntries];
		for(e=0; e<nEntries; ++e)
		{
			// Pointer to end environment boundary condition at entry [e]
			pEntryEndEnv[e] = &(rEndEnv[entry[e]]);

			// Pointer to node closest to rotor boundary in entry limb [e] - below assumes even end of entry pipe
			pRotorInletNode[e] = &(pEntryEndEnv[e]->pPipe[ONE_SIDE]->Node[pEntryEndEnv[e]->pPipe[ONE_SIDE]->N-1]);

			// Check that the entry end environment is configured for calibration if so required
			if(CALIBRATE)
			{
				if(pEntryEndEnv[e]->Get_VAR_P0())
				{
					pPpt->Out("\nWarning: End Environment ["); pPpt->Out(entry[e]); pPpt->Out("] is not set for APLDev calibration since VAR_P0 = "); pPpt->Out(TrueOrFalse(pEntryEndEnv[e]->Get_VAR_P0())); 
					pPpt->Out(". Change VAR_PO to false and continue? (y/n):");
					char cont;
					cin >> cont;
					//pPpt->Out("\n");
					if(cont=='n' || cont=='N')
					{
						pPpt->Out("User answers NO. Exiting...\n\n");
						exit(1);
					}
					else
					{
						pPpt->Out("User answers YES\n\n");
					}
					pEntryEndEnv[e]->Set_VAR_P0(false);
				}
			}
		}
	}
	pExitEndEnv = &(rEndEnv[exitEndEnv]);		// Join pointer to exit end environment boundary condition
	pRotorInletJunction = &(rJunc[juncInlet]);	// Join pointer to rotor inlet junction

	// Labels for pEntryEndEnv/pEntryTransm
	SINGLE_ENTRY = 0;
	FIRST_ENTRY = 0;
	SECOND_ENTRY = 1;

	// Set up result files
	std::string res_str;
	if(EX) res_str = "res_ex_apld"; else res_str = "res_in_apld";
	SetupFiles(pPpt, res_str); // Open results file
	
	pipe_flow = new int [2];
	pipe_flow_old = new int [2];

	LEFT_TO_RIGHT = true;
	M2_TOO_LOW = false;
	K_COUNTER = 0;

	// Set arbitrary values
	UPSTREAM = LEFT_SIDE;
	DOWNSTREAM = RIGHT_SIDE;

	M[UPSTREAM] = 0;
	M[DOWNSTREAM] = 0;

	if(INTERP_M2) PARAM_R = M[DOWNSTREAM];
	else PARAM_R = pBN[UPSTREAM]->Re;

	// Record operating point variables
	nLocs = 3;		// Number of locations (per limb)
	LIMB_INLET = 0;	// Limb location at the housing inlet
	LIMB_EXIT = 1;	// Limb location closest to the rotor
	ROTOR = 2;		// Location just upstream of the rotor

	p = new double* [nEntries];
	T = new double* [nEntries];
	m = new double* [nEntries];
	u = new double* [nEntries];
	p0 = new double* [nEntries];
	T0 = new double* [nEntries];

	int limbID, loc;
	for(limbID=0; limbID<nEntries; ++limbID)
	{
		p[limbID] = new double [nLocs];
		T[limbID] = new double [nLocs];
		m[limbID] = new double [nLocs];
		u[limbID] = new double [nLocs];
		p0[limbID] = new double [nLocs];
		T0[limbID] = new double [nLocs];

		for(loc=0; loc<nLocs; ++loc)
		{
			p[limbID][loc] = 1;
			T[limbID][loc] = 300;
			m[limbID][loc] = 0;
			u[limbID][loc] = 0;
			p0[limbID][loc] = 1;
			T0[limbID][loc] = 300;
		}
	}

	m_stage = new double [nLocs];	// Stage mass flow rates summed over entries at limb inlets or exits
	m_abs_sum = new double [nLocs];
	p_sum = new double [nLocs];		// Static pressures summed over entries at LIMB_INLET and LIMB_EXIT
	lambda = new double [nLocs];	// Admission ratio, lambda, at [LIMB_INLET] and [LIMB_EXIT] though only [LIMB_EXIT] is really used, [LIMB_INLET] is for information
	PR = new double [nEntries];
	MFP = new double [nEntries];
	Ws = new double [nEntries];
	tipVel = new double [nEntries];
	UCs = new double [nEntries];
	tempEta = new double [nEntries];

	mRotorInletJunction = 0;

	for(limbID=0; limbID<nEntries; ++limbID)
	{
		PR[limbID] = 1;
		MFP[limbID] = 0;
		Ws[limbID] = 0;
		tipVel[limbID] = 0;
		UCs[limbID] = 0;
		tempEta[limbID] = 1;
	}

	// Initialize starting values
	if(CONSTANT_K)
	{
		K_OVERALL = K_VALUE;
		eta_OVERALL = etaConst;
		D_OVERALL = D_VALUE;
	}
	else
	{
		K_OVERALL = K_baseline;
	}
	eta_OVERALL = 1.0;
	D_OVERALL = 0.0;
	OUTSIDE_RANGE_OVERALL = true;

	// Zero sums
	m_stage[LIMB_INLET] = 0; m_stage[LIMB_EXIT] = 0; m_stage[ROTOR] = 0;
	m_abs_sum[LIMB_INLET] = 0; m_abs_sum[LIMB_EXIT] = 0; m_abs_sum[ROTOR] = 0;
	p_sum[LIMB_INLET] = 0; p_sum[LIMB_EXIT] = 0; p_sum[ROTOR] = 0;
	lambda[LIMB_INLET] = 0.5; lambda[LIMB_EXIT] = 0.5;	// Both limbs zero flow, so equal admission
	PR_sum = 0;
	UCs_sum = 0;
	PR_stage_mean = 1;
	MFP_stage = 0;
	UCs_stage_mean = 0;
	Ws_stage_inlet = 0;
	UCs_stage_shaft = 0;
	Ws_stage_shaft = 0;
	Wact_stage_shaft = 0;

	// Labels
	// ======

	// Labels for vector MFP_vs_PR
	// ---------------------------
	lblMFPvsPR_p0_applied = 0;
	lblMFPvsPR_P0 = 1;
	lblMFPvsPR_PR = 2;
	lblMFPvsPR_MFP = 3;
	lblMFPvsPR_M2 = 4;
	lblMFPvsPR_UCs = 5;
	lblMFPvsPR_eta = 6;

	// Labels for steadyPRMFP
	// ----------------------
	lblSpeed = 0;	// Label for speed in steadyPRMFP and K_loss_file
	lblPR = 1;		// Label for PR in steadyPRMFP
	lblMFP = 2;		// Label for MFP in steadyPRMFP
	lblUCs = 3;	// Label for U/Cis in steadyPRMFP
	lblEta = 4;		// Label for eta in steadyPRMFP
	
	// Labels for K_loss_file (loss file or restart file when continuing calibration)
	// ------------------------------------------------------------------------------
	lblSpeedFile = 0;			// Label for speed in K_loss_file
	lblM2File = 1;				// Label for M2 in K_loss_file
	lblKFile = 2;				// Label for K in K_loss_file
	lblKPrevFile = 3;			// Label for KPrev in K_loss_file
	lblKPrevPrevFile = 4;		// Label for KPrevPrev in K_loss_file
	lblErrorFile = 5;			// Label for Error in K_loss_file
	lblErrorPrevFile = 6;		// Label for ErrorPrev in K_loss_file
	lblErrorPrevPrevFile = 7;	// Label for ErrorPrevPrev in K_loss_file
	lblFactorFile = 8;			// Label for Factor in K_loss_file
	lblKUCsFile = 9;			// Label for U/Cis in K_loss_file
	lblKEtaFile = 10;			// Label for Eta in K_loss_file
	lblKDFile = 11;				// Label for D in K_loss_file

	// Calibration data
	// ================
	if(CALIBRATE)
	{
		// Setup result files
		if(EX) res_str = "res_ex_apld_cal_prog_mfp_pr"; else res_str = "res_in_apld_cal_prog_mfp_pr";
		pCalibrationProgressMFPPR = fopen(ConstructString(pPpt, RES_DIR, res_str, ID, pPpt->strFileExt), "w");

		if(EX) res_str = "res_ex_apld_cal_final_mfp_pr"; else res_str = "res_in_apld_cal_final_mfp_pr";
		pCalibrationFinalMFPPR = fopen(ConstructString(pPpt, RES_DIR, res_str, ID, pPpt->strFileExt), "w");

		if(EX) res_str = "res_ex_apld_cal_prog_k"; else res_str = "res_in_apld_cal_prog_k";
		pCalibrationProgressK = fopen(ConstructString(pPpt, RES_DIR, res_str, ID, pPpt->strFileExt), "w");

		if(EX) res_str = "res_ex_apld_cal_final_k"; else res_str = "res_in_apld_cal_final_k";
		pCalibrationFinalK = fopen(ConstructString(pPpt, RES_DIR, res_str, ID, pPpt->strFileExt), "w");

    if(EX) res_str = "res_ex_apld_cal_final_k_dynasty"; else res_str = "res_in_apld_cal_final_k_dynasty";
		pCalibrationFinalKDynasty = fopen(ConstructString(pPpt, RES_DIR, res_str, ID, pPpt->strFileExt), "w");

		// Load steady mass flow characteristic (into steadyPRMFP[0,numSteadySpeeds-1][0,numSteadyPoints-1][0=speed,1=PR,2=MFP]) against which to calibrate
		LoadSteadyCharacteristics(pPpt, ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->APLD_DIR), STEADY_FILE));

		// Print processed/trimmed map to file
		if(EX) res_str = "res_ex_apld_cal_trimmed_map"; else res_str = "res_in_apld_cal_trimmed_map";
		pProcessedMap = fopen(ConstructString(pPpt, RES_DIR, res_str, ID, pPpt->strFileExt), "w");
		for(int tempSteadyCurve=0; tempSteadyCurve<numSteadySpeeds; ++tempSteadyCurve)
		{
			for(int tempRow=0; tempRow<numSteadyPoints[tempSteadyCurve]; ++tempRow)
			{
				fprintf(pProcessedMap,"%f\t", steadyPRMFP[tempSteadyCurve][tempRow][lblSpeed]);
				fprintf(pProcessedMap,"%f\t", steadyPRMFP[tempSteadyCurve][tempRow][lblPR]);
				fprintf(pProcessedMap,"%f\t", steadyPRMFP[tempSteadyCurve][tempRow][lblMFP]);
				fprintf(pProcessedMap,"%f\t", steadyPRMFP[tempSteadyCurve][tempRow][lblUCs]);
				if(tempSteadyCurve==numSteadySpeeds-1 && tempRow==numSteadyPoints[tempSteadyCurve]-1) fprintf(pProcessedMap,"%f", steadyPRMFP[tempSteadyCurve][tempRow][lblEta]);
				else fprintf(pProcessedMap,"%f\n", steadyPRMFP[tempSteadyCurve][tempRow][lblEta]);
			}
		}
		fclose(pProcessedMap);

		if(EX) res_str = "res_ex_apld_cal_final_k_restart"; else res_str = "res_in_apld_cal_final_k_restart";
		pFinalKRestart = fopen(ConstructString(pPpt, RES_DIR, res_str, ID, pPpt->strFileExt), "w");

		// Labels for K_loss_coeffs array
		int numKColumns = 11/*10*/;
		lblM2 = 0; lblK = 1; lblKPrev = 2; lblKPrevPrev = 3; lblError = 4; lblErrorPrev = 5; lblErrorPrevPrev = 6; lblFactor = 7; lblKUCs = 8; lblKEta = 9; lblKD = 10;

		// Create baseline K_loss_coeffs array using regularly-spaced M2 values and constant K values, one data point for each one in steadyPRMFP
		K_loss_coeffs = new double** [numSteadySpeeds];
		int k=0;
		for(steadyCurve=0; steadyCurve<numSteadySpeeds; ++steadyCurve)
		{
			K_loss_coeffs[steadyCurve] = new double* [numSteadyPoints[steadyCurve]];
			for(k=0; k<numSteadyPoints[steadyCurve]; ++k)
			{
				K_loss_coeffs[steadyCurve][k] = new double [numKColumns];
				K_loss_coeffs[steadyCurve][k][lblM2] = 0.0 + k*((1.0 - 0.0)/(numSteadyPoints[steadyCurve] - 1));
				K_loss_coeffs[steadyCurve][k][lblK] = K_baseline;													// Constant K as specified in input file
				K_loss_coeffs[steadyCurve][k][lblKPrev] = K_loss_coeffs[steadyCurve][k][lblK];		// Stores previous K value
				K_loss_coeffs[steadyCurve][k][lblKPrevPrev] = K_loss_coeffs[steadyCurve][k][lblK];	// Stores previous K value
				K_loss_coeffs[steadyCurve][k][lblError] = 0.0;														// Stores error in MFP (desired - measured) at this K value
				K_loss_coeffs[steadyCurve][k][lblErrorPrev] = 0.0;													// Stores previous error
				K_loss_coeffs[steadyCurve][k][lblErrorPrevPrev] = 0.0;												// Stores previous error
				K_loss_coeffs[steadyCurve][k][lblFactor] = 1.0;														// Initial multiplication factor
				K_loss_coeffs[steadyCurve][k][lblKUCs] = 0.7;														// Initial U/Cis value at this M2
				K_loss_coeffs[steadyCurve][k][lblKEta] = 0.7;														// Initial eta value to apply at this M2
				K_loss_coeffs[steadyCurve][k][lblKD] = 0.0;															// Initial D value to apply at this M2
			}
		}

		// Create baseline D_temp array
		D_temp = new double** [numSteadySpeeds];
		k=0;
		for(steadyCurve=0; steadyCurve<numSteadySpeeds; ++steadyCurve)
		{
			D_temp[steadyCurve] = new double* [numSteadyPoints[steadyCurve]];
			for(k=0; k<numSteadyPoints[steadyCurve]; ++k)
			{
				D_temp[steadyCurve][k] = new double [1];
				D_temp[steadyCurve][k][0] = 0.0;	// Initialize the temp D value
			}
		}
	
		if(!CONST_BASELINE) // Load an existing/restart K profile into K_loss_file, to be used as a starting point for calibration, but choose which one in LOSS_FILE[admissionMapNo]
		{
			// Load the selected maps
			K_loss_file = new double*** [1];	// Only ever need one map when using it as the restart file for calibration
			numLossSpeeds = new int [1];		// Since there will only be one map
			numLossPoints = new int* [1];
			int ONLY_ONE_MAP = 0;
			LoadLossCoefficients(pMyTime, pPpt, ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->APLD_DIR), RESTART_LOSS_FILE), ONLY_ONE_MAP);
			
			// Now check whether the chosen loss file be loaded into the runtime array by comparing array dimensions.
			bool tempMatch = false;
			if(numLossSpeeds[ONLY_ONE_MAP]==numSteadySpeeds)
			{
				pPpt->Out("\nNumber of speed curves in restart K file (numLossSpeeds[ONLY_ONE_MAP="); pPpt->Out(ONLY_ONE_MAP); pPpt->Out("]=");
				pPpt->Out(numLossSpeeds[ONLY_ONE_MAP]);
				pPpt->Out(") matches that in steady mass flow characteristic file (numSteadySpeeds=");
				pPpt->Out(numSteadySpeeds);
				pPpt->Out(")\n");
				tempMatch = true;
				int tempCurve=0;
				do
				{
					if(numLossPoints[ONLY_ONE_MAP][tempCurve]!=numSteadyPoints[tempCurve])
					{
						pPpt->Out("\nbut number of data points on loss curve (numLossPoints[ONLY_ONE_MAP="); 
						pPpt->Out(ONLY_ONE_MAP); 
						pPpt->Out("][tempCurve=");
						pPpt->Out(tempCurve);
						pPpt->Out("]=");
						pPpt->Out(numLossPoints[ONLY_ONE_MAP][tempCurve]);
						pPpt->Out(") does not match that on steady mass flow characteristic curve (numSteadyPoints[tempCurve=");
						pPpt->Out(tempCurve);
						pPpt->Out("]=");
						pPpt->Out(numSteadyPoints[tempCurve]);
						pPpt->Out(")\n");
						tempMatch = false;
					}
					++tempCurve;
				}while(tempMatch && tempCurve<numLossSpeeds[ONLY_ONE_MAP]/*or numSteadySpeeds*/);
			}
			else
			{
				pPpt->Out("\nNumber of speed curves in restart K file (numLossSpeeds[ONLY_ONE_MAP="); pPpt->Out(ONLY_ONE_MAP); pPpt->Out("]=");
				pPpt->Out(numLossSpeeds[ONLY_ONE_MAP]);
				pPpt->Out(") does not match that in steady mass flow characteristic file (numSteadySpeeds=");
				pPpt->Out(numSteadySpeeds);
				pPpt->Out(")\n");
			}

			if(tempMatch) // Yes we can load the file into the runtime array
			{
				pPpt->Out("Loading into K_loss_coeffs from file values\n\n");
				for(steadyCurve=0; steadyCurve<numSteadySpeeds; ++steadyCurve)
				{
					for(k=0; k<numSteadyPoints[steadyCurve]; ++k)
					{
						K_loss_coeffs[steadyCurve][k][lblM2] = K_loss_file[ONLY_ONE_MAP][steadyCurve][k][lblM2File];						// Restores M2 value
						K_loss_coeffs[steadyCurve][k][lblK] = K_loss_file[ONLY_ONE_MAP][steadyCurve][k][lblKFile];							// Restores K value
						K_loss_coeffs[steadyCurve][k][lblKPrev] = K_loss_file[ONLY_ONE_MAP][steadyCurve][k][lblKPrevFile];					// Restores previous K value
						K_loss_coeffs[steadyCurve][k][lblKPrevPrev] = K_loss_file[ONLY_ONE_MAP][steadyCurve][k][lblKPrevPrevFile];			// Restores previous previous K value
						K_loss_coeffs[steadyCurve][k][lblError] = K_loss_file[ONLY_ONE_MAP][steadyCurve][k][lblErrorFile];					// Restores error in MFP (desired - measured) at this K value
						K_loss_coeffs[steadyCurve][k][lblErrorPrev] = K_loss_file[ONLY_ONE_MAP][steadyCurve][k][lblErrorPrevFile];			// Restores previous error
						K_loss_coeffs[steadyCurve][k][lblErrorPrevPrev] = K_loss_file[ONLY_ONE_MAP][steadyCurve][k][lblErrorPrevPrevFile];	// Restores previous previous error
						K_loss_coeffs[steadyCurve][k][lblFactor] = K_loss_file[ONLY_ONE_MAP][steadyCurve][k][lblFactorFile];				// Restores multiplication factor
						K_loss_coeffs[steadyCurve][k][lblKUCs] = K_loss_file[ONLY_ONE_MAP][steadyCurve][k][lblKUCsFile];					// Restores U/Cis value
						K_loss_coeffs[steadyCurve][k][lblKEta] = K_loss_file[ONLY_ONE_MAP][steadyCurve][k][lblKEtaFile];					// Restores Eta value
						K_loss_coeffs[steadyCurve][k][lblKD] = K_loss_file[ONLY_ONE_MAP][steadyCurve][k][lblKDFile];						// Restores D value
					}
				}
			}
			else
			{
				pPpt->Out("\nCannot load into K_loss_coeffs from file values due to different dimensions. Intstead using a constant initial K = ");
				pPpt->Out(K_loss_coeffs[0/*steadyCurve*/][0/*k*/][lblK]);
				pPpt->Out("\n\n");
				CONST_BASELINE = true;
			}
		}
///*
		pPpt->Out("CAPLDev::Initialise\n"); 
		pPpt->Out("K_loss_coeffs[steadyCurve][k][j]:\n"); 
		pPpt->Out("-------------------------------------------------\n");
		pPpt->Out("M2\t\tK\tKPrev\tKPrevPrev\tError\t\tErrorPrev\tErrorPrevPrev\tFactor\tU/Cis\tEta\tD\n");
		for(steadyCurve=0; steadyCurve<numSteadySpeeds; ++steadyCurve)
		{
			for(k=0; k<numSteadyPoints[steadyCurve]; ++k)
			{
				pPpt->Out(K_loss_coeffs[steadyCurve][k][lblM2]); pPpt->Out("\t");
				pPpt->Out(K_loss_coeffs[steadyCurve][k][lblK]); pPpt->Out("\t");
				pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKPrev]); pPpt->Out("\t");
				pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKPrevPrev]); pPpt->Out("\t\t");
				pPpt->Out(K_loss_coeffs[steadyCurve][k][lblError]); pPpt->Out("\t");
				pPpt->Out(K_loss_coeffs[steadyCurve][k][lblErrorPrev]); pPpt->Out("\t\t");
				pPpt->Out(K_loss_coeffs[steadyCurve][k][lblErrorPrevPrev]); pPpt->Out("\t\t");
				pPpt->Out(K_loss_coeffs[steadyCurve][k][lblFactor]); pPpt->Out("\t");
				pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKUCs]); pPpt->Out("\t");
				pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKEta]); pPpt->Out("\t");
				pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKD]); pPpt->Out("\n");
			}
			pPpt->Out("\n");
		}
		//exit(1);
//*/
		// Set up calibration for first non-zero point of first speed curve
		steadyCurve = 0;
		row = 0; // steadyPRMFP conatins only non-zero MFP data, so row==0 is the first non-zero data point
		PR_point = steadyPRMFP[steadyCurve][row][lblPR];
		MFP_point = steadyPRMFP[steadyCurve][row][lblMFP];
		UCs_point = steadyPRMFP[steadyCurve][row][lblUCs];
		p0_applied = PR_point;//*(pPipe[1]->Node[pPipe[1]->N-1].p_dash*pPpt->PREF); // p01=PR*p2
		
		for(int e=0; e<nEntries; ++e)
		{
      pEntryEndEnv[e]->Set_CALIBRATE(true);  // Set CALIBRATE flag in each attached entry end environment
			pEntryEndEnv[e]->Set_P0_gradual(p0_applied); // Set first pressure level to be applied for calibration purposes
//cout << "Setting " << pEntryEndEnv[e]->Identify() << " to p0_applied = " << p0_applied << ", pEntryEndEnv[e]->Get_P0 = " << pEntryEndEnv[e]->Get_P0() << endl;
//exit(1);
		}

		pass_counter = 0;
		counter_between_steady = 0;
		PASS_INCREASING = true; // Start with an increasing calibration direction
		RESET_INITIAL_CONDITIONS = false;
		//FAKE_STEADY = false;

		tCurrent = 0;
		tRecord = 0;
		
		currentError = 0; 
		error_max = 0; 
		error_max_percent = 0; 
		max_error_element = 0;
	}
	else // Read resistance coefficients from file(s)
	{
		// If single entry need only one full admission map, else need the full admission map plus one partial admission map per entry
		K = new double [nAdmissionMaps];				
		eta = new double [nAdmissionMaps];
		D = new double [nAdmissionMaps];
		OUTSIDE_RANGE = new bool [nAdmissionMaps];
		K_loss_file = new double*** [nAdmissionMaps];
		numLossSpeeds = new int [nAdmissionMaps];		
		numLossPoints = new int* [nAdmissionMaps];
		
		int admissionMapNo;
		for(admissionMapNo=0; admissionMapNo<(nAdmissionMaps); ++admissionMapNo)
		{
			K[admissionMapNo] = K_baseline;
			eta[admissionMapNo] = 1.0;
			D[admissionMapNo] = 0.0;
			OUTSIDE_RANGE[admissionMapNo] = true;
			LoadLossCoefficients(pMyTime, pPpt, ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->APLD_DIR), LOSS_FILE[admissionMapNo]), admissionMapNo);
		}
		RESET_INITIAL_CONDITIONS = false; // Used only for calibration purposes

		error = 0; 
		minError = 0;
    counter = 0;
	}

	timeInCyclePrev = 0;
	timeInCycle = 0;

	SKIP = false;

  nMa = 10;//11;//3;//6;//11;
  lblMa = 0;
  lblMaError = 1;
  MaLevel;
  errorMap = new double* [nMa];
  for(MaLevel=0; MaLevel<nMa; ++MaLevel) errorMap[MaLevel] = new double [2];

  T0_us =300; //a_ds=347;	D_OVERALL = 0; 
  COMPUTETIME = false;
  elapsedtime = 0; melapsedtime = 0;
  
}

void CAPLDev::InitialisePost(CTime* pMyTime, CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rAPLDevPIPES, int** &rAPLDevPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, CTransmissive* TransmPtr, std::string param_dir, int assyid, bool transm_ptr_exists, string parent_assy_res_dir, CEndEnvironment* &rEndEnv, CTransmissive* &rTransm, CJunction* &rJunc, string calling_object_str)
{
	pInnerInlet = &(rPipe[innerInlet]);			// Join pointer to pipe immediately upstream of rotor inlet junction, inner limb
	pInnerFollowing = &(rPipe[innerFollowing]);	// Join pointer to pipe immediately downstream of rotor inlet junction, inner limb
	pOuterInlet = &(rPipe[outerInlet]);			// Join pointer to pipe immediately upstream of rotor inlet junction, outer limb
	pOuterFollowing = &(rPipe[outerFollowing]);	// Join pointer to pipe immediately downstream of rotor inlet junction, outer limb
/*
cout << pInnerInlet->ID << endl;
cout << pInnerFollowing->ID << endl;
cout << pOuterInlet->ID << endl;
cout << pOuterFollowing->ID << endl;
exit(1);
//*/
}

void CAPLDev::InitialisePostConfig(CProperties* pPpt)
{
  // General APLD variables
	// ====================================================================================================
  velHead = 0;
  delP = 0;
  pUs = pBN[UPSTREAM]->p_dash*pPpt->PREF;
}

void CAPLDev::RunBoundary(CProperties* pPpt, double time, double vel_grad, double mfr_under_est, CTime &rMyTime, double DELZe, double DELZi, double ca, double del_ca, double TIMEe, double TIMEi, int timestep, bool STEADY)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".RunBoundary\n");}
  //if(time>0.117){cout << "timestep = " << timestep << endl; cout << "HERE1\n";}
	RecordTurbineOperatingPoint(pPpt);
  //if(time>0.117){cout << "timestep = " << timestep << endl; cout << "HERE2\n";}
	General(pPpt, time, vel_grad, mfr_under_est, timestep, STEADY);
  //if(time>0.117){cout << "timestep = " << timestep << endl; cout << "HERE3\n";}
	if(time>CALIBRATION_DELAY && CALIBRATE) Calibrate(pPpt, rMyTime, DELZe, DELZi, ca, del_ca, TIMEe, TIMEi, timestep, STEADY);
  //if(time>0.117){cout << "timestep = " << timestep << endl; cout << "HERE4\n";}
}

void CAPLDev::General(CProperties* pPpt, double time, double vel_grad, double mfr_under_est, int timestep, bool STEADY)
//--------------------------------------------------//
// Non-homentropic flow direction test				//
// -----------------------------------				//
// Determines flow direction through				//
// adiabatic pressure loss device.					//
//--------------------------------------------------//
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".General\n");}

	//bool converged;
	//int counter;

	double lambda_in_c2, AA_c2, lambda_out1, lambda_out2;
	double* lambda_in_star_n;
	lambda_in_star_n = new double [2];
	
  double* lambda_in_star_test;
	lambda_in_star_test = new double [2];

	lambda_in_star_n[LEFT_SIDE] = (*(pCLIN[LEFT_SIDE]))[R+1] / pBN[LEFT_SIDE]->AA[1];
	lambda_in_star_n[RIGHT_SIDE] = (*(pCLIN[RIGHT_SIDE]))[R+1] / pBN[RIGHT_SIDE]->AA[1];
//	pPpt->Out("lambda_in_star_n[LEFT_SIDE] = "); pPpt->Out(lambda_in_star_n[LEFT_SIDE]); pPpt->Out("\n");
//	pPpt->Out("lambda_in_star_n[RIGHT_SIDE] = "); pPpt->Out(lambda_in_star_n[RIGHT_SIDE]); pPpt->Out("\n");

  // Direction test
  // ==============
  if(fabs(lambda_in_star_n[RIGHT_SIDE] - lambda_in_star_n[LEFT_SIDE]) <= pPpt->ZERO_TOL)
  {
//	  	pPpt->Out("No flow in adiabatic pressure loss device\n");
//	  	cin >> pause;

	  // Set arbitrary values
	  UPSTREAM = LEFT_SIDE;
	  DOWNSTREAM = RIGHT_SIDE;

	  // No flow - treat both as closed ends
	  //pipe_flow[LEFT_SIDE] = NOFLOW;
	  //pipe_flow[RIGHT_SIDE] = NOFLOW;
	  *(pend_flow[LEFT_SIDE]) = NOFLOW;
	  *(pend_flow[RIGHT_SIDE]) = NOFLOW;

	  (*(pCLOUT[LEFT_SIDE]))[R+1] = (*(pCLIN[LEFT_SIDE]))[R+1];
	  (*(pCLOUT[RIGHT_SIDE]))[R+1] = (*(pCLIN[RIGHT_SIDE]))[R+1];

	  if(CONSTANT_K)
	  {
		  K_OVERALL = K_VALUE;
		  eta_OVERALL = etaConst;					// Associated constant eta value
		  D_OVERALL = D_VALUE;						// Associated constant D value
	  }
	  else
	  {
		  K_OVERALL = 0;	// Actually undefined for no flow as it stands at the moment
		  etaConst = 0;
		  D_OVERALL = 0;
	  }
	  return; // No need to update lambda_in or AA
  }
  else
  {
	  if(lambda_in_star_n[LEFT_SIDE] - lambda_in_star_n[RIGHT_SIDE] > pPpt->ZERO_TOL)
	  {
		  // Flow from left to right; left is upstream, right is downstream
//			pPpt->Out("Flow is left to right\n");
		  UPSTREAM = LEFT_SIDE;
		  DOWNSTREAM = RIGHT_SIDE;
		  if(!LEFT_TO_RIGHT)
		  {
//				pPpt->Out("Flow in adiabatic pressure loss device is changing direction, to left-to-right\n");
//				cin >> pause;
		  }
		  LEFT_TO_RIGHT = true;
	  }
	  else
	  {
		  // Flow from right to left; left is downstream, right is upstream
//			pPpt->Out("Flow is right to left\n");
		  UPSTREAM = RIGHT_SIDE;
		  DOWNSTREAM = LEFT_SIDE;
		  if(LEFT_TO_RIGHT)
		  {
//				pPpt->Out("Flow in adiabatic pressure loss device is changing direction, to right-to-left\n");
//				cin >> pause;
		  }
		  LEFT_TO_RIGHT = false;
	  }
  }	

  // Main algorithm
	// ==============
  Loss(pPpt, time, lambda_in_c2, AA_c2, lambda_out1, lambda_out2, vel_grad, mfr_under_est, timestep, STEADY);

	//////////////////////////////////////////
	double *lambda_in_c, *lambda_out_c, *AA_c;
	lambda_in_c = new double [NPIPES]; lambda_out_c = new double [NPIPES]; AA_c = new double [NPIPES];
				
	lambda_out_c[UPSTREAM] = lambda_out1;
	lambda_out_c[DOWNSTREAM] = lambda_out2;
	lambda_in_c[DOWNSTREAM] = lambda_in_c2; // lambda_in correction needed for DOWNSTREAM only
	AA_c[DOWNSTREAM] = AA_c2;

	bool* CHOKED; CHOKED = new bool [NPIPES];
	CHOKED[UPSTREAM] = false;
	CHOKED[DOWNSTREAM] = false;
	
	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out_c, AA_c, pipe_flow, CHOKED);
	//////////////////////////////////////////

	// Delete local arrays
	delete [] lambda_in_star_n;
	delete [] lambda_in_star_test;
	delete [] lambda_in_c;
	delete [] lambda_out_c;
	delete [] AA_c;
	delete [] CHOKED;
}


void CAPLDev::Loss(CProperties* pPpt, double time, double &rlambda_in_c2, double &rAA_c2, double &rlambda_out1, double &rlambda_out2, double vel_grad, double mfr_under_est, int timestep, bool STEADY)
// ====================================================================================================
// Developed from Benson's model
// ====================================================================================================
{
  bool SWITCH = false;

  //do{

	// Set flow directions
	// ----------------------------------------------------------------------------------------------------
	pipe_flow[UPSTREAM] = OUTFLOW;
	pipe_flow[DOWNSTREAM] = INFLOW;

	// Calculate area ratio parameter for this flow direction (<1)
	// ----------------------------------------------------------------------------------------------------
	alpha = pow(pBN[UPSTREAM]->f_dash/pBN[DOWNSTREAM]->f_dash, 2);
	
	// Enter other system parameters
	// ----------------------------------------------------------------------------------------------------
	phi = 1;

	// Enter inital characteristic values
	// ----------------------------------------------------------------------------------------------------
	lambda_in_n[UPSTREAM] = (*(pCLIN[UPSTREAM]))[1];
	lambda_in_n[DOWNSTREAM] = (*(pCLIN[DOWNSTREAM]))[1];
	AA_n[UPSTREAM] = (pBN[UPSTREAM])->AA[1];
	AA_n[DOWNSTREAM] = (pBN[DOWNSTREAM])->AA[1];
	lambda_in[UPSTREAM] = lambda_in_n[UPSTREAM];
	AA_1 = AA_n[UPSTREAM];                      // Since no correction needed on this side
  lambda_in_c2 = lambda_in_n[DOWNSTREAM];     // Initial guess

	// Loop to converge on M[DOWNSTREAM]
	// ----------------------------------------------------------------------------------------------------
  //M2_max = 0.1;//0.5;//1.0;
	//if(timestep<=1)
  //  M[DOWNSTREAM] = M2_max/2;	// 1/2 of M2_max as in Benson p449
	//del_M2 = M2_max/4;			  // 1/4 of M2_max as in Benson p449

  lowMa = -1.0;
  highMa = 1.0;

  counter = 0;
  do
  {
    ++counter;
    for(MaLevel=0; MaLevel<nMa; ++MaLevel)
    {
      errorMap[MaLevel][lblMa] = lowMa + MaLevel*((highMa - lowMa)/(nMa - 1));
      errorMap[MaLevel][lblMaError] = 1e6;
    }
    for(MaLevel=0; MaLevel<nMa; ++MaLevel) errorMap[MaLevel][lblMaError] = LossEquations(pPpt, time, timestep, errorMap[MaLevel][lblMa]);
 
    //for(MaLevel=0; MaLevel<nMa; ++MaLevel){cout << "errorMap[MaLevel=" << MaLevel << "]:\t" << errorMap[MaLevel][lblMa] << "\t\t" << errorMap[MaLevel][lblMaError] << endl;}

    FINISHED = false;
    MaLevel = 0;
    do
    {
      if(errorMap[MaLevel][lblMaError] < 0 && errorMap[MaLevel+1][lblMaError] > 0)
      {
        FINISHED = true;
        lblNeg = MaLevel;
        lblPos = MaLevel + 1;
      }
      ++MaLevel;
    }while(MaLevel+1<nMa && !FINISHED);

    lowMa = errorMap[lblNeg][lblMa];
    highMa = errorMap[lblPos][lblMa];

    if(!FINISHED)
    {
      cout << Identify() << ":Loss: lblNeg & lblPos not found\n";
      cout << Identify() << ":Loss: probably due to supersonic flow\n";
      cout << "timestep = " << timestep << endl;
      cout << "counter = " << counter << endl;

      int temp_precision = cout.precision(); cout << setprecision(12);
    for(MaLevel=0; MaLevel<nMa; ++MaLevel)
		{
			cout << "errorMap[MaLevel=" << MaLevel << "]:\t" << errorMap[MaLevel][lblMa] << "\t\t" << errorMap[MaLevel][lblMaError] << endl;
		}
 
    cout << endl;
    cout << "errorMap[lblNeg=" << lblNeg << "][lblMa, lblMaError] = " << errorMap[lblNeg][lblMa] << ", " << errorMap[lblNeg][lblMaError] << endl;
    cout << "errorMap[lblPos=" << lblPos << "][lblMa, lblMaError] = " << errorMap[lblPos][lblMa] << ", " << errorMap[lblPos][lblMaError] << endl;
    cout << endl << endl;
    cout << setprecision(temp_precision);

	cout << "Warning, error here. Continue? (y/n) \n";
	char cont;
	cin >> cont;
	cout << "\n";
	if(cont=='n' || cont=='N')
		exit(1);
	}
/*
    if(timestep>=452)
    //if(counter>6)
    {
    cout << "timestep = " << timestep << endl;
    cout << "counter = " << counter << endl;
    int temp_precision = cout.precision(); cout << setprecision(12);
    for(MaLevel=0; MaLevel<nMa; ++MaLevel)
    {
      cout << "errorMap[MaLevel=" << MaLevel << "]:\t" << errorMap[MaLevel][lblMa] << "\t\t" << errorMap[MaLevel][lblMaError] << endl;
    }
 
    cout << endl;
    cout << "errorMap[lblNeg=" << lblNeg << "][lblMa, lblMaError] = " << errorMap[lblNeg][lblMa] << ", " << errorMap[lblNeg][lblMaError] << endl;
    cout << "errorMap[lblPos=" << lblPos << "][lblMa, lblMaError] = " << errorMap[lblPos][lblMa] << ", " << errorMap[lblPos][lblMaError] << endl;
    cout << endl << endl;
    cout << setprecision(temp_precision);
    }
//*/
  }while( (fabs(errorMap[lblNeg][lblMaError]) > pPpt->LOSS_TOL || fabs(errorMap[lblPos][lblMaError]) > pPpt->LOSS_TOL) );// && !SWITCH);

  // Run equations one final time using mean of two values to set variables using converged M2
  LossEquations(pPpt, time, timestep, 0.5*(errorMap[lblNeg][lblMa] + errorMap[lblPos][lblMa]));


/*
  counter = 0;
	error = 1e6;
  bool GREATER = true;    // Arbitrary setting; indicates sign of error
	do
	{
    ++counter;
    error = LossEquations(pPpt, time, timestep, M[DOWNSTREAM]);

    if(lambda_in_d1>lambda_in[UPSTREAM])
    {
	    if(!GREATER)
      {
        del_M2 /= 2;
      }
      GREATER = true;

	    do
	    {
	  	  temp_M2 = M[DOWNSTREAM] - del_M2;
        if(temp_M2 <= pPpt->ZERO_TOL) 
        {
          del_M2 /= 2;
        }
	    }while(temp_M2 <= pPpt->ZERO_TOL);
	    M[DOWNSTREAM] = temp_M2;
    }
    else
    {
	    if(GREATER)
      {
        del_M2 /= 2;
      }
      GREATER = false;

	    do
	    {
	  	  temp_M2 = M[DOWNSTREAM] + del_M2;
	  	  if(temp_M2 <= pPpt->ZERO_TOL) 
        {
          del_M2 /= 2;
        }
	    }while(temp_M2 <= pPpt->ZERO_TOL);
	    M[DOWNSTREAM] = temp_M2;
    }
    
    if(counter>1000)
	  {
		  pPpt->Out("\n");
		  pPpt->Out(Identify()); pPpt->Out(": Loss: counter = "); pPpt->Out(counter); pPpt->Out("\n");
      pPpt->Out("time = "); pPpt->Out(time); pPpt->Out("\n");
      pPpt->Out("timestep = "); pPpt->Out(timestep); pPpt->Out("\n");
		  pPpt->Out("M[DOWNSTREAM] = "); pPpt->Out(M[DOWNSTREAM]); pPpt->Out("\n");
      pPpt->Out("K_OVERALL = "); pPpt->Out(K_OVERALL); pPpt->Out("\n");
		  //pPpt->Out("FIXING_K = "); pPpt->Out(TrueOrFalse(FIXING_K)); pPpt->Out("\n");
      pPpt->Out("lambda_in_d1 = "); pPpt->Out(lambda_in_d1); pPpt->Out("\n");
      pPpt->Out("lambda_in[UPSTREAM] = "); pPpt->Out(lambda_in[UPSTREAM]); pPpt->Out("\n");
		  pPpt->Out("error = "); pPpt->Out(error); pPpt->Out("\n");
		  //pPpt->Out("error_prev = "); pPpt->Out(error_prev); pPpt->Out("\n");
      pPpt->Out("del_M2 = "); pPpt->Out(del_M2); pPpt->Out("\n");
		  pPpt->Out("\n");
      exit(1);
    }
  }while(fabs(error) > pPpt->LOSS_TOL);
*/

  if(counter>50)
  {
    pPpt->Out(Identify()); pPpt->Out(": LossEquations: final counter is > 50, = "); pPpt->Out(counter); pPpt->Out("\n");
  }

  // Loop output
  // ----------------------------------------------------------------------------------------------------
	lambda_out[UPSTREAM] = AA_1*(A_star[UPSTREAM] - ((k_us-1)/2)*U_star[UPSTREAM]);			  // Eq. 8.180
	lambda_out[DOWNSTREAM] = AA_c2*(A_star[DOWNSTREAM] + ((k_ds-1)/2)*U_star[DOWNSTREAM]);// Eq. 8.181

  double u_us = U_star[UPSTREAM]*AA_1*pPipe[UPSTREAM]->AREF;
  double A_us = A_star[UPSTREAM]*AA_1;
	double T_us = pow(A_us,2)*(EX ? pPpt->TREFe : pPpt->TREFi);
  double p_dash_us = pow(A_us/AA_1, (2*k_us)/(k_us-1));		
	double rho_us = ((p_dash_us*pPpt->PREF)*1e5)/(pPpt->R_air*T_us);
  velHead = (0.5*rho_us*pow(u_us,2))/1e5;
  delP = ((2*K_OVERALL/k_us)*(0.5*rho_us*pow(u_us,2)))/1e5;
  pUs = p_dash_us*pPpt->PREF;

  a_ds = A_star[DOWNSTREAM]*AA_c2*pPipe[DOWNSTREAM]->AREF;
  T0_us = TotalTemperature(pPpt, T_us, u_us);
  d_temp = 2*dh0_turb*pow(M[DOWNSTREAM],2)/(pow(a_ds,2)*T0_us);
  Wact_turb = mRotorInletJunction*D_OVERALL*pow(a_ds,2)/(2*pow(M[DOWNSTREAM],2));

  // Update lambda_in2 = lambda_in_c2, AA_2 = AA_c2, lambda_out2, lambda_out1, by reference
  // ----------------------------------------------------------------------------------------------------
  rlambda_in_c2 = lambda_in_c2;
  rAA_c2 = AA_c2;
  rlambda_out1 = lambda_out[UPSTREAM];
  rlambda_out2 = lambda_out[DOWNSTREAM];
	return;
}



double CAPLDev::LossEquations(CProperties* pPpt, double time, int timestep, double Ma)
{ 
  M[DOWNSTREAM] = Ma;
	double localError;
  localError = 1e6;

  if(!CALIBRATE)
  {
	  if(time<0.01)  {CONSTANT_K = true;}		// To deal with transition region before flow field established
	  else {CONSTANT_K = false;}
  }

  bool INTERPOLATE = true;
 	if(INTERPOLATE)
	{
		//if(CALIBRATE)
		if(!CONSTANT_K && CALIBRATE)
		{
			// Interpolate dynamic values contained in K_loss_coeffs vector
			int k=1;
			while(K_loss_coeffs[steadyCurve][k][0] <= M[DOWNSTREAM] && k < numSteadyPoints[steadyCurve] - 1) ++k;
			// Interpolate between [k] and [k-1]

			// Interpolate for new K
			double K_new, D_new;
			if(fabs(K_loss_coeffs[steadyCurve][k][0/*M2*/] - K_loss_coeffs[steadyCurve][k-1][0/*M2*/]) < pPpt->ZERO_TOL) // Avoid DIV0
			{
				//pPpt->Out("Avoiding division by zero!\n");
				// Select a value of K based on whether M[DOWNSTREAM] is closer to [k] or [k-1]
				if(fabs(K_loss_coeffs[steadyCurve][k][0/*M2*/] - M[DOWNSTREAM]) 
				 < fabs(K_loss_coeffs[steadyCurve][k-1][0/*M2*/] - M[DOWNSTREAM])) 
				{
					K_new = K_loss_coeffs[steadyCurve][k][1/*K*/];
					D_new = K_loss_coeffs[steadyCurve][k][10/*D*/];
				}
				else
				{
					K_new = K_loss_coeffs[steadyCurve][k-1][1/*K*/];
					D_new = K_loss_coeffs[steadyCurve][k-1][10/*D*/];
				}
			}
			else
			{
				K_new = K_loss_coeffs[steadyCurve][k-1][1/*K*/] + 
						((K_loss_coeffs[steadyCurve][k][1/*K*/] - K_loss_coeffs[steadyCurve][k-1][1/*K*/])
						/(K_loss_coeffs[steadyCurve][k][0/*M2*/] - K_loss_coeffs[steadyCurve][k-1][0/*M2*/]))
						*(M[DOWNSTREAM] - K_loss_coeffs[steadyCurve][k-1][0/*M2*/]);
				D_new = K_loss_coeffs[steadyCurve][k-1][10/*D*/] + 
						((K_loss_coeffs[steadyCurve][k][10/*D*/] - K_loss_coeffs[steadyCurve][k-1][10/*D*/])
						/(K_loss_coeffs[steadyCurve][k][0/*M2*/] - K_loss_coeffs[steadyCurve][k-1][0/*M2*/]))
						*(M[DOWNSTREAM] - K_loss_coeffs[steadyCurve][k-1][0/*M2*/]);
			}

			// Only interpolate else can get negative values for K
			FIXING_K = false;
			if(K_loss_coeffs[steadyCurve][k-1][0/*M2*/] <= K_loss_coeffs[steadyCurve][k][0/*M2*/]) // As you would expect - M2 increasing with each data point
			{
				if(M[DOWNSTREAM] < K_loss_coeffs[steadyCurve][k-1][0/*M2*/])	// M[DOWNSTREAM] smaller than smallest M2 in the map... taken as the smallest K or D.
				{
					FIXING_K = true;
					K_OVERALL = K_loss_coeffs[steadyCurve][k-1][1/*K*/];
					D_OVERALL = K_loss_coeffs[steadyCurve][k-1][10/*D*/];
				}
				else
				{
					if(M[DOWNSTREAM] > K_loss_coeffs[steadyCurve][k][0/*M2*/])	// M[DOWNSTREAM] larger than largest M2 in the map... taken as the largest K or D.
					{
						FIXING_K = true;
						K_OVERALL = K_loss_coeffs[steadyCurve][k][1/*K*/];
						D_OVERALL = K_loss_coeffs[steadyCurve][k][10/*D*/];
					}
					else	// The actual interpolated value
					{
						K_OVERALL = K_new;
						D_OVERALL = D_new;
					}
				}
			}
			else // Also allow for this situation - M2 decreasing with each data point
			{
				if(M[DOWNSTREAM] > K_loss_coeffs[steadyCurve][k-1][0/*M2*/])
				{
					FIXING_K = true;
					K_OVERALL = K_loss_coeffs[steadyCurve][k-1][1/*K*/];
					D_OVERALL = K_loss_coeffs[steadyCurve][k-1][10/*D*/];
				}
				else
				{
					if(M[DOWNSTREAM] < K_loss_coeffs[steadyCurve][k][0/*M2*/])
					{
						FIXING_K = true;
						K_OVERALL = K_loss_coeffs[steadyCurve][k][1/*K*/];
						D_OVERALL = K_loss_coeffs[steadyCurve][k][10/*D*/];
					}
					else	// The actual interpolated value
					{
						K_OVERALL = K_new;
						D_OVERALL = D_new;
					}
				}
			}

			// Calibrating, so just use the K and eta values for the current element
			K_OVERALL = K_loss_coeffs[steadyCurve][row][lblK/*1*/];
			eta_OVERALL = K_loss_coeffs[steadyCurve][row][lblKEta];
			//D_OVERALL = K_loss_coeffs[steadyCurve][row][lblKD];
		}
		else
		{
			if(CONSTANT_K)
			{
				K_OVERALL = K_VALUE;
				eta_OVERALL = etaConst;
				D_OVERALL = D_VALUE;
			}
			else
			{
				if(INTERP_M2)	
				{
					PARAM_R = M[DOWNSTREAM];		  // Interpolate based on downstream Mach number, M2
				}
				else
				{
					PARAM_R = pBN[UPSTREAM]->Re;	// Interpolate based on upstream Re (recorded upon entering this method)
				}
				/*if(time>0.5)
				{
					starttime = clock();	cout<<"Start = "<<starttime<<"\t";	// This is the interpolation loop start time.
					COMPUTETIME = true;
				}*/
				int admissionMapNo;
				OUTSIDE_RANGE_OVERALL = false;					// Reset
				for(admissionMapNo=0; admissionMapNo<(nAdmissionMaps); ++admissionMapNo)
				{
					InterpolateLossCoefficients(pPpt, time, PARAM_S, PARAM_R, K[admissionMapNo], eta[admissionMapNo], OUTSIDE_RANGE[admissionMapNo],
												lblSpeed, lblM2File, lblKFile, lblKEtaFile, false/*Show extrapolation warnings*/, admissionMapNo);// Interpolate value of K, eta
					InterpolateLossCoefficients(pPpt, time, PARAM_S, PARAM_R, K[admissionMapNo], D[admissionMapNo], OUTSIDE_RANGE[admissionMapNo],
												lblSpeed, lblM2File, lblKFile, lblKDFile, false/*Show extrapolation warnings*/, admissionMapNo);// Interpolate value of K, D
					if(K[admissionMapNo]<0) K[admissionMapNo] = 0;					// Must not have negative K else -1.#IND
					if(OUTSIDE_RANGE[admissionMapNo]) OUTSIDE_RANGE_OVERALL = true; // If operation of any map lies outside calibrated range, overall operation is also outside range
				}

				// Determine K_OVERALL
				if(nAdmissionMaps==1 || lambda[LIMB_EXIT]==0.5 || !USE_LAMBDA)
				{
					K_OVERALL = K[0/*FULL ADMISSION*/];		// If single entry or equal admission just use full admission map value
					eta_OVERALL = eta[0/*FULL ADMISSION*/];	// If single entry or equal admission just use full admission map value
					D_OVERALL = D[0/*FULL ADMISSION*/];	// If single entry or equal admission just use full admission map value
				}
				else										// Linearly interpolate between two K values (full admission and one partial admission value)
				{
					if(LINEAR_INTERP)	// linearly interpolate
					{
						if(lambda[LIMB_EXIT] > 0.5)						// Admission is weighted to flow through entry [0]
						{
							K_OVERALL = weightPartial*
										(((K[1/*PARTIAL ADMISSION [0]*/] - K[0/*FULL ADMISSION*/])/0.5)			// Gradient
										*(lambda[LIMB_EXIT] - 0.5))												// Multiplier 
										+ 
										(1 - weightPartial)*
										K[0/*FULL ADMISSION*/];												// Addition
							eta_OVERALL = weightPartial*
										(((eta[1/*PARTIAL ADMISSION [0]*/] - eta[0/*FULL ADMISSION*/])/0.5)	// Gradient
										*(lambda[LIMB_EXIT] - 0.5))												// Multiplier 
										+
										(1 - weightPartial)*
										eta[0/*FULL ADMISSION*/];												// Addition
							D_OVERALL = weightPartial*
										(((D[1/*PARTIAL ADMISSION [0]*/] - D[0/*FULL ADMISSION*/])/0.5)			// Gradient
										*(lambda[LIMB_EXIT] - 0.5))												// Multiplier 
										+
										(1 - weightPartial)*
										D[0/*FULL ADMISSION*/];													// Addition
						}
						else							// Admission is weighted to flow through entry [1]
						{
							K_OVERALL = weightPartial*
										(((K[2/*PARTIAL ADMISSION [1]*/] - K[0/*FULL ADMISSION*/])/0.5)			// Gradient
										*(0.5 - lambda[LIMB_EXIT]))												// Multiplier 
										+ 
										(1 - weightPartial)*
										K[0/*FULL ADMISSION*/];												// Addition
							eta_OVERALL = weightPartial*
										(((eta[2/*PARTIAL ADMISSION [1]*/] - eta[0/*FULL ADMISSION*/])/0.5)	// Gradient
										*(0.5 - lambda[LIMB_EXIT]))												// Multiplier 
										+ 
										(1 - weightPartial)*
										eta[0/*FULL ADMISSION*/];												// Addition
							D_OVERALL = weightPartial*
										(((D[2/*PARTIAL ADMISSION [1]*/] - D[0/*FULL ADMISSION*/])/0.5)			// Gradient
										*(0.5 - lambda[LIMB_EXIT]))												// Multiplier 
										+	
										(1 - weightPartial)*
										D[0/*FULL ADMISSION*/];												// Addition
						}
					}
					else	// polynomially interpolate
					{
						if(lambda[LIMB_EXIT] > 0.5)						// Admission is weighted to flow through entry [0]
						{
							K_OVERALL = ((K[1]-K[0]) / (pow(0.5,POLY_DEGREE))			// Gradient
										*pow((lambda[LIMB_EXIT]-0.5),POLY_DEGREE) )		// Multiplier
										+ 
										K[0];											// Addition

							eta_OVERALL = ((eta[1]-eta[0]) / (pow(0.5,POLY_DEGREE))		// Gradient
										*pow((lambda[LIMB_EXIT]-0.5),POLY_DEGREE) )		// Multiplier
										+ 
										eta[0];											// Addition

							D_OVERALL = ((D[1]-D[0]) / (pow(0.5,POLY_DEGREE))			// Gradient
										*pow((lambda[LIMB_EXIT]-0.5),POLY_DEGREE) )		// Multiplier
										+ 
										D[0];											// Addition
						}
						else
						{
							K_OVERALL = ((K[2]-K[0]) / (pow(0.5,POLY_DEGREE))			// Gradient
										*pow((0.5-lambda[LIMB_EXIT]),POLY_DEGREE) )		// Multiplier
										+ 
										K[0];											// Addition

							eta_OVERALL = ((eta[2]-eta[0]) / (pow(0.5,POLY_DEGREE))		// Gradient
										*pow((0.5-lambda[LIMB_EXIT]),POLY_DEGREE) )		// Multiplier
										+ 
										eta[0];											// Addition

							D_OVERALL = ((D[2]-D[0]) / (pow(0.5,POLY_DEGREE))			// Gradient
										*pow((0.5-lambda[LIMB_EXIT]),POLY_DEGREE) )		// Multiplier
										+ 
										D[0];											// Addition
						}
					}
				}
			}
		}
	}
	if(CALIBRATE && PR_stage_mean<PR_THRESHOLD) D_OVERALL=0;		// apply PR threshold
	D_OVERALL = D_OVERALL*T0_us/1000;
	
	
	
/*	if(COMPUTETIME)
	{
		endtime = clock();	cout<<"End = "<<endtime<<"\t";	// This is the interpolation loop end time.
		elapsedtime = elapsedtime+(endtime-starttime);	pPpt->Out("Elapsed time = "); pPpt->Out(elapsedtime); pPpt->Out("\n");	
		MeanlineCalc(pPpt, time, timestep);
	}*/

  // Boundary method equations using M[DOWNSTREAM] as the independent variable
  // ----------------------------------------------------------------------------------------------------
	k_ds = pPpt->gammaAir(pBN[DOWNSTREAM]->T);
	k_us = pPpt->gammaAir(pBN[UPSTREAM]->T);
 	a = 2/(k_ds - 1);		// Intermediate value as in Benson p448
	b = a*pow(M[DOWNSTREAM],2) + pow(M[DOWNSTREAM],4);	// Intermediate value as in Benson p448

/*	M[UPSTREAM] = sqrt( ( (2*b*K_OVERALL + alpha*a) - sqrt( pow((2*b*K_OVERALL + alpha*a),2) - 4*b*(b*pow(K_OVERALL,2) - alpha) ) ) 
							/(2*(b*pow(K_OVERALL,2) - alpha)) );*/ // Eq. 8.168

	M[UPSTREAM] = sqrt( ( (2*(b+D_OVERALL)*K_OVERALL + alpha*a) - sqrt( pow((2*(b+D_OVERALL)*K_OVERALL + alpha*a),2) - 4*(b+D_OVERALL)*((b+D_OVERALL)*pow(K_OVERALL,2) - alpha) ) ) 
							/(2*((b+D_OVERALL)*pow(K_OVERALL,2) - alpha)) ); // Revised Eq. 8.168

/*	AA2_over_AA1 = pow( (1/(1 - K_OVERALL*pow(M[UPSTREAM],2))), (k_us-1)/(2*k_us) )
				   * sqrt( ( (2/(k_us-1)) + pow(M[UPSTREAM],2) )/( (2/(k_ds-1)) + pow(M[DOWNSTREAM],2) ) );*/ // Eq. 8.172

	AA2_over_AA1 = pow( (1/(1 - K_OVERALL*pow(M[UPSTREAM],2))), (k_us-1)/(2*k_us) )
				   * sqrt( ( (2/(k_us-1)) + pow(M[UPSTREAM],2) )/
				   ( (2/(k_ds-1)) + pow(M[DOWNSTREAM],2) + (D_OVERALL/pow(M[DOWNSTREAM],2)) ) ); // Revised Eq. 8.172

	AA_c2 = (AA2_over_AA1)*AA_1;						// Eq. 8.173

/*
  lambda_in_star_c2 = lambda_in_c2/AA_c2;				// Eq. 8.174
	A_star[UPSTREAM] = pow( (1/(1 - K_OVERALL*pow(M[UPSTREAM],2))), (k_us-1)/(2*k_us) )
					   * (lambda_in_star_c2/(1 - ((k_ds-1)/2)*M[DOWNSTREAM]));	// Eq. 8.156	
	A_star[DOWNSTREAM] = lambda_in_star_c2/(1 - ((k_ds-1)/2)*M[DOWNSTREAM]);	// Eq. 8.176
//*/    
///*
  lambda_in_star_n2 = lambda_in_n[DOWNSTREAM]/AA_c2;				// Eq. 8.174

  A_star[UPSTREAM] = pow( (1/(1 - K_OVERALL*pow(M[UPSTREAM],2))), (k_us-1)/(2*k_us) )
					   * (lambda_in_star_n2/(1 - ((k_ds-1)/2)*M[DOWNSTREAM]));	// Eq. 8.156
  A_star[DOWNSTREAM] = lambda_in_star_n2/(1 - ((k_ds-1)/2)*M[DOWNSTREAM]);	// Eq. 8.176
//*/

  U_star[UPSTREAM] = M[UPSTREAM]*A_star[UPSTREAM];	// Eq. 8.175     
  U_star[DOWNSTREAM] = M[DOWNSTREAM]*A_star[DOWNSTREAM];	// Eq. 8.177
  
  lambda_in_d1 = AA_1*(A_star[UPSTREAM] + ((k_us - 1)/2)*U_star[UPSTREAM]);	// Eq. 8.178
  localError/*error*/ = lambda_in_d1 - lambda_in[UPSTREAM];
	
	lambda_in_c2 = lambda_in_n[DOWNSTREAM] + A_star[DOWNSTREAM]*(AA_c2 - AA_n[DOWNSTREAM]);	// Eq. 8.179
	

	return localError;
}



void CAPLDev::Calibrate(CProperties* pPpt, CTime &rMyTime, double DELZe, double DELZi, double ca, double del_ca, double TIMEe, double TIMEi, int timestep, bool STEADY)
//--------------------------------------------------//
//--------------------------------------------------//
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Calibrate\n");}
	++counter_between_steady;
	if(EX) tCurrent = TIMEe; else tCurrent = TIMEi; // Note the current simulation time

  bool ENDENV_DESIRED_P0_READY = true;
  for(int e=0; e<nEntries; ++e) if(!pEntryEndEnv[e]->Get_CALIBRATE_READY()) ENDENV_DESIRED_P0_READY = false;

  //if(ENDENV_DESIRED_P0_READY) cout << "ENDENV_DESIRED_P0_READY" << endl;

	if(ENDENV_DESIRED_P0_READY && (STEADY || counter_between_steady >= waitCount) && ((tCurrent - tRecord) >= tInterval))	// Record operating point
  //if(ENDENV_DESIRED_P0_READY && (STEADY) && ((tCurrent - tRecord) >= tInterval))	// Record operating point
	{
		tRecord = tCurrent; // Note the simulation time at which this point is being recorded

    /*
		if(counter_between_steady >= waitCount)
		{
			//pPpt->Out("counter_between_steady = "); pPpt->Out(counter_between_steady); pPpt->Out("\n");
			pPpt->Out("\n"); pPpt->Out("CAPLDev::Calibrate: assumed steady conditions ");
		}
		else
    */
      pPpt->Out("CAPLDev::Calibrate: steady conditions ");
		counter_between_steady = 0;
		pPpt->tol_steady_multiplier = 1;
		pPpt->Out("at operating point ["); pPpt->Out(row); pPpt->Out("] of [0-"); pPpt->Out(numSteadyPoints[steadyCurve]-1); 
		pPpt->Out("] on speed curve ["); pPpt->Out(steadyCurve); pPpt->Out("] of [0-"); pPpt->Out(numSteadySpeeds-1); pPpt->Out("]\n");
		
		PARAM_S = steadyPRMFP[steadyCurve][0/*tempRow*/][lblSpeed];// Change PARAM_S from the user-specified value in the input file, to that of the current speed curve being calibrated
		p_exit/*bar*/ = pExitEndEnv->pBN[ONE_SIDE]->p_dash*pPpt->PREF;														

		for(int entryID=0; entryID<nEntries; ++entryID) // For each entry
		{
			p[entryID][LIMB_INLET]/*bar*/ = pEntryEndEnv[entryID]->pBN[ONE_SIDE]->p_dash*pPpt->PREF;
			T[entryID][LIMB_INLET] = pEntryEndEnv[entryID]->pBN[ONE_SIDE]->T;		
			m[entryID][LIMB_INLET] = pEntryEndEnv[entryID]->pBN[ONE_SIDE]->mdot;											
			u[entryID][LIMB_INLET] = pEntryEndEnv[entryID]->pBN[ONE_SIDE]->U*pEntryEndEnv[entryID]->pPipe[ONE_SIDE]->AREF;	
			p0[entryID][LIMB_INLET]/*bar*/ = TotalPressureBar(pPpt, p[entryID][LIMB_INLET], T[entryID][LIMB_INLET], u[entryID][LIMB_INLET]);
			T0[entryID][LIMB_INLET] = TotalTemperature(pPpt, T[entryID][LIMB_INLET], u[entryID][LIMB_INLET]);
			tipVel[entryID] = PI*(PARAM_S*sqrt(T0[entryID][LIMB_INLET]))*(tipDiameter*1e-3);

			if(1-pow(p_exit/p0[entryID][LIMB_INLET],(pPpt->gammaAir(T[entryID][LIMB_INLET])-1)/pPpt->gammaAir(T[entryID][LIMB_INLET]))<=0) // Prevent division by zero
				UCs[entryID] = 0;
			else
				UCs[entryID] = tipVel[entryID]/sqrt(2*pPpt->cpAir(T[entryID][LIMB_INLET])*T0[entryID][LIMB_INLET]*(1-pow(p_exit/p0[entryID][LIMB_INLET],(pPpt->gammaAir(T[entryID][LIMB_INLET])-1)/pPpt->gammaAir(T[entryID][LIMB_INLET]))));

			// Find UCs element closest to tempUCs, and interpolate/extrapolate using it and the next closest
			int tempRow, rowNextClosest, rowClosest;
			double deltaUCs;

			rowNextClosest = 0;
			rowClosest = 0;
			deltaUCs = 1e6; // A large number
			for(tempRow = 0; tempRow < numSteadyPoints[steadyCurve]; ++tempRow)
			{
				if(fabs(steadyPRMFP[steadyCurve][tempRow][lblUCs] - UCs[entryID]) < deltaUCs)
				{
					rowClosest = tempRow;
					deltaUCs = fabs(steadyPRMFP[steadyCurve][tempRow][lblUCs] - UCs[entryID]);
				}
			}

			// Found row closest to tempUCs, now find appropriate other row for interpolation
			if(rowClosest==numSteadyPoints[steadyCurve]-1) rowNextClosest = rowClosest - 1; // Must do this
			else
			{
				if(rowClosest==0) rowNextClosest = rowClosest + 1; // Must do this
				else
				{
					if(UCs[entryID] > steadyPRMFP[steadyCurve][rowClosest][lblUCs]) // Look for row having a value greater than tempUCs
					{
						if(steadyPRMFP[steadyCurve][rowClosest+1][lblUCs] > UCs[entryID]) rowNextClosest = rowClosest + 1;
						else rowNextClosest = rowClosest - 1;
					}
					else // Look for row having a value less than tempUCs
					{
						if(steadyPRMFP[steadyCurve][rowClosest+1][lblUCs] < UCs[entryID]) rowNextClosest = rowClosest + 1;
						else rowNextClosest = rowClosest - 1;
					}
				}
			}

			if(rowClosest <= rowNextClosest)
			{
				tempEta[entryID] = steadyPRMFP[steadyCurve][rowClosest][lblEta] + 
							(((steadyPRMFP[steadyCurve][rowNextClosest][lblEta] - steadyPRMFP[steadyCurve][rowClosest][lblEta])
							/(steadyPRMFP[steadyCurve][rowNextClosest][lblUCs] - steadyPRMFP[steadyCurve][rowClosest][lblUCs]))	// Gradient
							*(UCs[entryID] - steadyPRMFP[steadyCurve][rowClosest][lblUCs]));
			}
			else
			{
				tempEta[entryID] = steadyPRMFP[steadyCurve][rowNextClosest][lblEta] + 
							(((steadyPRMFP[steadyCurve][rowClosest][lblEta] - steadyPRMFP[steadyCurve][rowNextClosest][lblEta])
							/(steadyPRMFP[steadyCurve][rowClosest][lblUCs] - steadyPRMFP[steadyCurve][rowNextClosest][lblUCs]))	// Gradient
							*(UCs[entryID] - steadyPRMFP[steadyCurve][rowNextClosest][lblUCs]));
			}
		}
		
		double* temp;
		temp = new double [7];
		temp[lblMFPvsPR_p0_applied]	= p0_applied;
		temp[lblMFPvsPR_P0]	= pEntryEndEnv[0]->Get_P0();
		temp[lblMFPvsPR_M2]	= M[DOWNSTREAM];

		K_loss_coeffs[steadyCurve][row][lblM2] = M[DOWNSTREAM]; // Record M2 for this row/datapoint in K_loss_coeffs

		if(nEntries==1) // Single-entry turbine
		{
			temp[lblMFPvsPR_MFP] = m[SINGLE_ENTRY][LIMB_INLET]*sqrt(T0[SINGLE_ENTRY][LIMB_INLET])/p0[SINGLE_ENTRY][LIMB_INLET]; // MFP = mdot_temp*sqrt(T01_temp)/p01_temp;
//cout << "temp[lblMFPvsPR_MFP] SINGLE_ENTRY = " << temp[lblMFPvsPR_MFP] << endl;
//cout << endl;

			temp[lblMFPvsPR_PR] = p0[SINGLE_ENTRY][LIMB_INLET]/p_exit; // PR = p01_temp/p_exit
//cout << "temp[lblMFPvsPR_PR] SINGLE_ENTRY = " << temp[lblMFPvsPR_PR] << endl;
//cout << endl;

			temp[lblMFPvsPR_UCs] = UCs[SINGLE_ENTRY];	
			temp[lblMFPvsPR_eta] = tempEta[SINGLE_ENTRY]; // Need to associate this interpolated value of eta with the K value at this operating point

			K_loss_coeffs[steadyCurve][row][lblKUCs] = UCs[SINGLE_ENTRY];	// Record U/Cis for this row/datapoint in K_loss_coeffs, but just for info
			K_loss_coeffs[steadyCurve][row][lblKEta] = tempEta[SINGLE_ENTRY];	// Record eta for this row/datapoint in K_loss_coeffs
		}
		else
		{
			if(nEntries==2) // Twin-entry turbine
			{
				// Must sum MFPs from both entries
				temp[lblMFPvsPR_MFP] = (m[FIRST_ENTRY][LIMB_INLET]*sqrt(T0[FIRST_ENTRY][LIMB_INLET])/p0[FIRST_ENTRY][LIMB_INLET]) 
										+ (m[SECOND_ENTRY][LIMB_INLET]*sqrt(T0[SECOND_ENTRY][LIMB_INLET])/p0[SECOND_ENTRY][LIMB_INLET]); // MFP = mdot_temp*sqrt(T01_temp)/p01_temp;
//pPpt->Out("temp[lblMFPvsPR_MFP] FIRST_ENTRY = "); pPpt->Out((mdot_temp[FIRST_ENTRY]*sqrt(T01_temp[FIRST_ENTRY])/p01_temp[FIRST_ENTRY])); pPpt->Out("\n");
//pPpt->Out("temp[lblMFPvsPR_MFP] SECOND_ENTRY = "); pPpt->Out((mdot_temp[SECOND_ENTRY]*sqrt(T01_temp[SECOND_ENTRY])/p01_temp[SECOND_ENTRY])); pPpt->Out("\n");
//pPpt->Out("temp[lblMFPvsPR_MFP] COMBINED = "); pPpt->Out(temp[lblMFPvsPR_MFP]); pPpt->Out("\n");
//pPpt->Out("\n");

				// Average PR across entries
				temp[lblMFPvsPR_PR] = (0.5*(p0[FIRST_ENTRY][LIMB_INLET] + p0[SECOND_ENTRY][LIMB_INLET]))/p_exit; // PR = p01_temp/p_exit
//pPpt->Out("temp[lblMFPvsPR_PR] FIRST_ENTRY = "); pPpt->Out(p01_temp[FIRST_ENTRY]/p_exit); pPpt->Out("\n");
//pPpt->Out("temp[lblMFPvsPR_PR] SECOND_ENTRY = "); pPpt->Out(p01_temp[SECOND_ENTRY]/p_exit); pPpt->Out("\n");
//pPpt->Out("temp[lblMFPvsPR_PR] AVERAGE = "); pPpt->Out(temp[lblMFPvsPR_PR]); pPpt->Out("\n");
//pPpt->Out("\n");

				// Average other parameters
				temp[lblMFPvsPR_UCs] = 0.5*(UCs[FIRST_ENTRY] + UCs[SECOND_ENTRY]);
				temp[lblMFPvsPR_eta] = 0.5*(tempEta[FIRST_ENTRY] + tempEta[SECOND_ENTRY]);
//pPpt->Out("temp[lblMFPvsPR_UCs] FIRST_ENTRY = "); pPpt->Out(tempUCs[FIRST_ENTRY]); pPpt->Out("\n");
//pPpt->Out("temp[lblMFPvsPR_UCs] SECOND_ENTRY = "); pPpt->Out(tempUCs[SECOND_ENTRY]); pPpt->Out("\n");
//pPpt->Out("temp[lblMFPvsPR_UCs] AVERAGE = "); pPpt->Out(temp[lblMFPvsPR_UCs]); pPpt->Out("\n");
//pPpt->Out("\n");
//pPpt->Out("temp[lblMFPvsPR_eta] FIRST_ENTRY = "); pPpt->Out(tempEta[FIRST_ENTRY]); pPpt->Out("\n");
//pPpt->Out("temp[lblMFPvsPR_eta] SECOND_ENTRY = "); pPpt->Out(tempEta[SECOND_ENTRY]); pPpt->Out("\n");
//pPpt->Out("temp[lblMFPvsPR_eta] AVERAGE = "); pPpt->Out(temp[lblMFPvsPR_eta]); pPpt->Out("\n");
//pPpt->Out("\n");

				K_loss_coeffs[steadyCurve][row][lblKUCs] = 0.5*(UCs[FIRST_ENTRY] + UCs[SECOND_ENTRY]);	// Record U/Cis for this row/datapoint in K_loss_coeffs, but just for info (eta is linked to K, not U/Cis)
				K_loss_coeffs[steadyCurve][row][lblKEta] = 0.5*(tempEta[FIRST_ENTRY] + tempEta[SECOND_ENTRY]);		// Record eta for this row/datapoint in K_loss_coeffs
			}
			else
			{
				pPpt->Out("nEntries = "); pPpt->Out(nEntries); 
				pPpt->Out(" - not currently capable of dealing with more than two turbine entries. Exiting...\n");
				exit(1);
			}
		}

///*
temp[lblMFPvsPR_PR]	= PR_point; // Use map data point value of PR rather than the measured value
//pPpt->Out("but set PR = PR_point instead, temp[lblMFPvsPR_PR] = "); pPpt->Out(temp[lblMFPvsPR_PR]); pPpt->Out("\n");
//pPpt->Out("\n");
//temp[lblMFPvsPR_UCs] = UCs_point; // Use data point value of UCs rather than the measured value
//*/

//cout << "CAPLDev::Calibrate: p01_temp recorded, p01_temp = " << p01_temp << endl;
//cout << "CAPLDev::Calibrate: p_exit recorded, p_exit = " << p_exit << endl;
//cout << "CAPLDev::Calibrate: PR recorded, temp[lblMFPvsPR_PR] = " << temp[lblMFPvsPR_PR] << endl;
//pPpt->Out("\n");

		// Test MFP error of current point
		
		D_temp[steadyCurve][row][0] = d_temp*1000;
		//K_loss_coeffs[steadyCurve][row][lblKD] = D_temp[steadyCurve][row][0];
		
		double currentPR, currentMFP, interpolatedMFP;
		double interpolatedUCs, interpolatedEta;
		currentPR = temp[lblMFPvsPR_PR];
		currentMFP = temp[lblMFPvsPR_MFP];
		LinearlyInterpolateSteadyCharacteristics(steadyCurve, currentPR, interpolatedMFP, interpolatedUCs, interpolatedEta);
		currentError = interpolatedMFP - currentMFP;	// Calculate error

		pPpt->Out("current PARAM_S = "); pPpt->Out(PARAM_S); pPpt->Out("\n"); 
		pPpt->Out("currentPR = "); pPpt->Out(currentPR); pPpt->Out("\n");
		pPpt->Out("currentMFP = "); pPpt->Out(currentMFP); pPpt->Out("\n");
		pPpt->Out("interpolatedMFP = "); pPpt->Out(interpolatedMFP); pPpt->Out("\n");
		pPpt->Out("currentError = "); pPpt->Out(currentError); pPpt->Out(" or "); pPpt->Out(fabs(currentError/interpolatedMFP)*100); pPpt->Out("%\n");
		pPpt->Out("interpolatedUCs = "); pPpt->Out(interpolatedUCs); pPpt->Out("\n");
		pPpt->Out("interpolatedEta = "); pPpt->Out(interpolatedEta); pPpt->Out("\n");
		pPpt->Out("currentD (d_temp) = "); pPpt->Out(d_temp); pPpt->Out("\n");
		pPpt->Out("currentD_OVERALL = "); pPpt->Out(D_OVERALL); pPpt->Out("\n");

		if(fabs(currentError/interpolatedMFP)*100 >= tol_MFP && !SKIP) 
		{
			K_loss_coeffs[steadyCurve][row][lblM2] = M[DOWNSTREAM]; // Set the M2 value in K_loss_coeffs to that currently being experienced

			// Check if error has changed sign; if new error is smaller then continue, else restore
			if(currentError*K_loss_coeffs[steadyCurve][row][lblError] < 0)// && pass_counter > 1) 
			// Divide new error by previous error and test sign, but only from second pass onwards
			{
				// Error has switched - restore previous values and halve the predicted adjustment
				pPpt->Out("Error has changed sign");
	bool RESTORE = false;
				// If new error is larger than previously restore old values
				if(fabs(currentError) > fabs(K_loss_coeffs[steadyCurve][row][lblError]) && RESTORE)
				{
					pPpt->Out(", but is larger than before. Restore previous values.\n");
					K_loss_coeffs[steadyCurve][row][lblK] = K_loss_coeffs[steadyCurve][row][lblKPrev];
					K_loss_coeffs[steadyCurve][row][lblKPrev] = K_loss_coeffs[steadyCurve][row][lblKPrevPrev];
					K_loss_coeffs[steadyCurve][row][lblError] = K_loss_coeffs[steadyCurve][row][lblErrorPrev];
					currentError = K_loss_coeffs[steadyCurve][row][lblError];
					K_loss_coeffs[steadyCurve][row][lblErrorPrev] = K_loss_coeffs[steadyCurve][row][lblErrorPrevPrev];
				}
				else	// Continue to update values
				{
					if(RESTORE) pPpt->Out(", but is smaller than before. Halve factor but continue to update values.\n");
					else pPpt->Out(", so halve factor.\n");
					K_loss_coeffs[steadyCurve][row][lblErrorPrevPrev] = K_loss_coeffs[steadyCurve][row][lblErrorPrev];
					K_loss_coeffs[steadyCurve][row][lblErrorPrev] = K_loss_coeffs[steadyCurve][row][lblError];
					K_loss_coeffs[steadyCurve][row][lblError] = currentError;
				}
				//if(K_loss_coeffs[steadyCurve][row][lblFactor] > 1) K_loss_coeffs[steadyCurve][row][lblFactor] = 1; // Reset to the standard value
				K_loss_coeffs[steadyCurve][row][lblFactor] *= 0.5; // Halve the movement in both cases since changed sign
			}
			else	// Continue to update values
			{
				pPpt->Out("Error has not changed sign. Continue to update values.\n");
				K_loss_coeffs[steadyCurve][row][lblErrorPrevPrev] = K_loss_coeffs[steadyCurve][row][lblErrorPrev];
				K_loss_coeffs[steadyCurve][row][lblErrorPrev] = K_loss_coeffs[steadyCurve][row][lblError];
				K_loss_coeffs[steadyCurve][row][lblError] = currentError;
			}

			// Adjust K by false position method
			K_loss_coeffs[steadyCurve][row][lblKPrevPrev] = K_loss_coeffs[steadyCurve][row][lblKPrev];	// Save previous K value
			K_loss_coeffs[steadyCurve][row][lblKPrev] = K_loss_coeffs[steadyCurve][row][lblK];			// Save previous K value
			double del_K = 1.0;//0.5;
			//if(fabs(currentError) > tol_MFP) // Only adjust if this element has not met the convergence criteria
			{
				if(currentError<0) // Increase K
				{
					K_loss_coeffs[steadyCurve][row][lblK] += K_loss_coeffs[steadyCurve][row][lblFactor]*del_K;
				}
				else // Decrease K
				{
					K_loss_coeffs[steadyCurve][row][lblK] -= K_loss_coeffs[steadyCurve][row][lblFactor]*del_K;
				}
				
				while(K_loss_coeffs[steadyCurve][row][lblK]<=0) // && K_loss_coeffs[steadyCurve][row][lblFactor]>0.01}
				{	
					K_loss_coeffs[steadyCurve][row][lblK] = K_loss_coeffs[steadyCurve][row][lblKPrev]; // Return to previous value
					K_loss_coeffs[steadyCurve][row][lblFactor] /= 2;

					if(currentError<0) // Increase K
					{
						K_loss_coeffs[steadyCurve][row][lblK] += K_loss_coeffs[steadyCurve][row][lblFactor]*del_K;
					}
					else // Decrease K
					{
						K_loss_coeffs[steadyCurve][row][lblK] -= K_loss_coeffs[steadyCurve][row][lblFactor]*del_K;
					}

					if(K_loss_coeffs[steadyCurve][row][lblFactor]<=0.01) SKIP = true;
					//if(K_loss_coeffs[steadyCurve][row][lblK] < 0) K_loss_coeffs[steadyCurve][row][lblK] = 0;		// Prohibit K from going negative
				}
			}

			

			//pPpt->Out("After adjustment:\n");
			//pPpt->Out("row = "); pPpt->Out(row); pPpt->Out("\n");
			pPpt->Out("M2 = "); pPpt->Out(K_loss_coeffs[steadyCurve][row][lblM2]); 
			pPpt->Out("\tK = "); pPpt->Out(K_loss_coeffs[steadyCurve][row][lblK]); 
			pPpt->Out("\tK (prev) = "); pPpt->Out(K_loss_coeffs[steadyCurve][row][lblKPrev]); 
			pPpt->Out("\tK (prev prev) = "); pPpt->Out(K_loss_coeffs[steadyCurve][row][lblKPrevPrev]); 
			pPpt->Out("\tError = "); pPpt->Out(K_loss_coeffs[steadyCurve][row][lblError]); 
			pPpt->Out("\tError (previous) = "); pPpt->Out(K_loss_coeffs[steadyCurve][row][lblErrorPrev]); 
			pPpt->Out("\tError (prev prev) = "); pPpt->Out(K_loss_coeffs[steadyCurve][row][lblErrorPrevPrev]); 
			pPpt->Out("\tFactor = "); pPpt->Out(K_loss_coeffs[steadyCurve][row][lblFactor]);
			pPpt->Out("\tU/Cis = "); pPpt->Out(K_loss_coeffs[steadyCurve][row][lblKUCs]);
			pPpt->Out("\tEta = "); pPpt->Out(K_loss_coeffs[steadyCurve][row][lblKEta]);
			pPpt->Out("\tD = "); pPpt->Out(D_temp[steadyCurve][row][0]);
			pPpt->Out("\n");
			pPpt->Out("\n");
		}
		else // This point has converged on this pass
		{
			if(SKIP)
			{	
				pPpt->Out("SKIPPING\n");
				SKIP = false; // Reset if necessary
			}

			pPpt->Out("Point converged to within the user-specified tolerance of "); pPpt->Out(tol_MFP); pPpt->Out("%\n");
			pPpt->Out("For which M2 = "); pPpt->Out(K_loss_coeffs[steadyCurve][row][lblM2]); pPpt->Out("\n");
			pPpt->Out("Corresponding K = "); pPpt->Out(K_loss_coeffs[steadyCurve][row][lblK]); pPpt->Out("\n");
			pPpt->Out("Corresponding D = "); pPpt->Out(D_temp[steadyCurve][row][0]); pPpt->Out("\n");
			pPpt->Out("\n");

			
			
			// Add this temp to the vector in the right place (depending on pass direction)
			if(PASS_INCREASING) MFP_vs_PR.push_back(temp);	// If increasing pass, add to end of vector 
			else MFP_vs_PR.insert(MFP_vs_PR.begin(), temp);	// If decreasing pass, add to start of vector	

			// Record which point has the largest final error
			if(fabs(currentError) > fabs(error_max))
			{
				error_max = currentError;
				error_max_percent = fabs(currentError/interpolatedMFP)*100;//fabs(error/interpMFP)*100;
				max_error_element = row;	// Track the largest error
			}

			// Adjust applied pressure according to direction of pass
			if((PASS_INCREASING && row < numSteadyPoints[steadyCurve] - 1) || (!PASS_INCREASING && row > 0))
			{
				if(PASS_INCREASING) ++row;
				else --row;

				PR_point = steadyPRMFP[steadyCurve][row][lblPR];
				MFP_point = steadyPRMFP[steadyCurve][row][lblMFP];
				UCs_point = steadyPRMFP[steadyCurve][row][lblUCs];
				p0_applied = PR_point;//*p_exit; // p01=PR*p2
	
				for(int e=0; e<nEntries; ++e)
				{
					pEntryEndEnv[e]->Set_P0_gradual(p0_applied);	// Request the EndEnv to gradually adjust itself to requested stagnation pressure
					//cout << "Point has converged on this pass. Setting " << pEntryEndEnv[e]->Identify() << " to p0_applied = " << p0_applied << ", pEntryEndEnv[e]->Get_P0 = " << pEntryEndEnv[e]->Get_P0() << endl;
					//exit(1);
				}
			}
			else // Current pass has finished
			{
				++ pass_counter;

				// Dump steady characteristic to file
				if(pass_counter==1) // If this is the first pass for this speed curve calibration
				{
					fprintf(pCalibrationProgressMFPPR,"Starting calibration of speed curve [");
					fprintf(pCalibrationProgressMFPPR,"%d", steadyCurve);
					fprintf(pCalibrationProgressMFPPR,"]\n\n");

					fprintf(pCalibrationProgressK,"Starting calibration of speed curve [");
					fprintf(pCalibrationProgressK,"%d", steadyCurve);
					fprintf(pCalibrationProgressK,"]\n\n");
				}

				int i;
				for(i=0; i<int(MFP_vs_PR.size()); ++i) 
					fprintf(pCalibrationProgressMFPPR,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
						MFP_vs_PR[i][lblMFPvsPR_p0_applied], 
						MFP_vs_PR[i][lblMFPvsPR_P0],
						MFP_vs_PR[i][lblMFPvsPR_PR], 
						MFP_vs_PR[i][lblMFPvsPR_MFP], 
						MFP_vs_PR[i][lblMFPvsPR_M2], 
						MFP_vs_PR[i][lblMFPvsPR_UCs], 
						MFP_vs_PR[i][lblMFPvsPR_eta]);
				fprintf(pCalibrationProgressMFPPR,"\n");

				// Dump K loss coefficients to file
				int k;
				for(k=0; k<numSteadyPoints[steadyCurve]; ++k) 
					fprintf(pCalibrationProgressK,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
						K_loss_coeffs[steadyCurve][k][lblM2], 
						K_loss_coeffs[steadyCurve][k][lblK], 
						K_loss_coeffs[steadyCurve][k][lblKPrev], 
						K_loss_coeffs[steadyCurve][k][lblKPrevPrev], 
						K_loss_coeffs[steadyCurve][k][lblError], 
						K_loss_coeffs[steadyCurve][k][lblErrorPrev], 
						K_loss_coeffs[steadyCurve][k][lblErrorPrevPrev], 
						K_loss_coeffs[steadyCurve][k][lblFactor],
						K_loss_coeffs[steadyCurve][k][lblKUCs],
						K_loss_coeffs[steadyCurve][k][lblKEta],
						K_loss_coeffs[steadyCurve][k][lblKD]);
				fprintf(pCalibrationProgressK,"\n");
/*
				// Now compare measured and desired steady characteristics then adjust loss values appropriately
				int max_error_element;
				double measuredPR, measuredMFP, interpMFP, grad_temp, error, error_max, error_max_percent;
				double interpUCs, interpEta;
				double measuredUCs, measuredEta;
				error = 0;
				error_max = 0;
				error_max_percent = 0;
				max_error_element = 0;
				bool SWITCHED = false;

				//for(k=0; k<K_loss_coeffs.size(); ++k)
				for(k=0; k<numSteadyPoints[steadyCurve]; ++k)
				{
					// Interpolate actual steady characteristic for this loss coefficient's M2 value
					i=0;
					while((K_loss_coeffs[steadyCurve][k][lblM2] >= MFP_vs_PR[i][lblMFPvsPR_M2] || i==0) && i < MFP_vs_PR.size() - 1) ++i;
					// M2 for this K_loss_coeffs lies between i and i-1 of MFP_vs_PR[i][lblMFPvsPR_M2]
				
					// Do not permit extraploation, interpolation only, otherwise blows up
					if(K_loss_coeffs[steadyCurve][k][lblM2] < MFP_vs_PR[i-1][lblMFPvsPR_M2]) // Outside of interval - so would otherwise extrapolate
					{
						cout << "Warning: preventing extrapolation!!\n\n";

						// Choose values closest
						if(fabs(K_loss_coeffs[steadyCurve][k][lblM2] - MFP_vs_PR[i-1][lblMFPvsPR_M2]) < fabs(K_loss_coeffs[steadyCurve][k][lblM2] - MFP_vs_PR[i][lblMFPvsPR_M2]))
						{
							measuredPR = MFP_vs_PR[i-1][lblMFPvsPR_PR];
							measuredMFP = MFP_vs_PR[i-1][lblMFPvsPR_MFP];
							measuredUCs = MFP_vs_PR[i-1][lblMFPvsPR_UCs];
							measuredEta = MFP_vs_PR[i-1][lblMFPvsPR_eta];
						}
						else
						{
							measuredPR = MFP_vs_PR[i][lblMFPvsPR_PR];
							measuredMFP = MFP_vs_PR[i][lblMFPvsPR_MFP];
							measuredUCs = MFP_vs_PR[i][lblMFPvsPR_UCs];
							measuredEta = MFP_vs_PR[i][lblMFPvsPR_eta];
						}
					}
					else // Normal interpolation expressions should be ok
					{
						measuredPR = MFP_vs_PR[i-1][lblMFPvsPR_PR] + 
											(((MFP_vs_PR[i][lblMFPvsPR_PR] - MFP_vs_PR[i-1][lblMFPvsPR_PR])/(MFP_vs_PR[i][lblMFPvsPR_M2] - MFP_vs_PR[i-1][lblMFPvsPR_M2]))		// Gradient
											*(K_loss_coeffs[steadyCurve][k][lblM2] - MFP_vs_PR[i-1][lblMFPvsPR_M2]));
						measuredMFP = MFP_vs_PR[i-1][lblMFPvsPR_MFP] + 
											(((MFP_vs_PR[i][lblMFPvsPR_MFP] - MFP_vs_PR[i-1][lblMFPvsPR_MFP])/(MFP_vs_PR[i][lblMFPvsPR_M2] - MFP_vs_PR[i-1][lblMFPvsPR_M2]))	// Gradient
											*(K_loss_coeffs[steadyCurve][k][lblM2] - MFP_vs_PR[i-1][lblMFPvsPR_M2]));
						measuredUCs = MFP_vs_PR[i-1][lblMFPvsPR_UCs] + 
											(((MFP_vs_PR[i][lblMFPvsPR_UCs] - MFP_vs_PR[i-1][lblMFPvsPR_UCs])/(MFP_vs_PR[i][lblMFPvsPR_M2] - MFP_vs_PR[i-1][lblMFPvsPR_M2]))	// Gradient
											*(K_loss_coeffs[steadyCurve][k][lblM2] - MFP_vs_PR[i-1][lblMFPvsPR_M2]));
						measuredEta = MFP_vs_PR[i-1][lblMFPvsPR_eta] + 
											(((MFP_vs_PR[i][lblMFPvsPR_eta] - MFP_vs_PR[i-1][lblMFPvsPR_eta])/(MFP_vs_PR[i][lblMFPvsPR_M2] - MFP_vs_PR[i-1][lblMFPvsPR_M2]))	// Gradient
											*(K_loss_coeffs[steadyCurve][k][lblM2] - MFP_vs_PR[i-1][lblMFPvsPR_M2]));
					}
					LinearlyInterpolateSteadyCharacteristics(steadyCurve, measuredPR, interpMFP, interpUCs, interpEta);
					error = interpMFP - measuredMFP;	// Calculate error

					//pPpt->Out("\n\n\n");
					//pPpt->Out("Point ["); pPpt->Out(k); pPpt->Out("] of [0-"); pPpt->Out(numSteadyPoints[steadyCurve]-1); pPpt->Out("]:\n");
					////pPpt->Out("k = "); pPpt->Out(k); pPpt->Out("\n");
					////pPpt->Out("Check:\n");
					//pPpt->Out("MFP_vs_PR[i-1="); pPpt->Out(i-1); pPpt->Out("][lblMFPvsPR_M2] = "); pPpt->Out(MFP_vs_PR[i-1][lblMFPvsPR_M2]); pPpt->Out("\n");
					//pPpt->Out("K_loss_coeffs[steadyCurve="); pPpt->Out(steadyCurve); pPpt->Out("][k="); pPpt->Out(k); pPpt->Out("][lblM2] = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblM2]); pPpt->Out("\n");
					//pPpt->Out("MFP_vs_PR[i="); pPpt->Out(i); pPpt->Out("][lblMFPvsPR_M2] = "); pPpt->Out(MFP_vs_PR[i][lblMFPvsPR_M2]); pPpt->Out("\n");
					//pPpt->Out("\n");
					//pPpt->Out("Measured PR, measuredPR = "); pPpt->Out(measuredPR); pPpt->Out("\n");
					//pPpt->Out("Corresponding MFP interpolated from map, interpMFP = "); pPpt->Out(interpMFP); pPpt->Out("\n");
					//pPpt->Out("Actual measured MFP, measuredMFP = "); pPpt->Out(measuredMFP); pPpt->Out("\n");
					//pPpt->Out("This is an error of "); pPpt->Out(error); pPpt->Out(" ("); pPpt->Out((fabs(error)/interpMFP)*100); pPpt->Out("%)\n");
					//pPpt->Out("\n");
					////pPpt->Out("measuredPR = "); pPpt->Out(measuredPR); pPpt->Out("\n");
					////pPpt->Out("Corresponding U/Cis interpolated from map, interpUCs = "); pPpt->Out(interpUCs); pPpt->Out("\n");
					////pPpt->Out("But actual measured U/Cis, measuredUCs = "); pPpt->Out(measuredUCs); pPpt->Out("\n");
					////pPpt->Out("\n");
					////pPpt->Out("measuredEta = "); pPpt->Out(measuredEta); pPpt->Out("\n");
					////pPpt->Out("interpEta = "); pPpt->Out(interpEta); pPpt->Out("\n");
					////pPpt->Out("\n");
	

					// Check if error has changed sign; if new error is smaller then continue, else restore
					if(error*K_loss_coeffs[steadyCurve][k][lblError] < 0)// && pass_counter > 1) 
					// Divide new error by previous error and test sign, but only from second pass onwards
					{
						// Error has switched - restore previous values and halve the predicted adjustment
						pPpt->Out("Error has changed sign");
						SWITCHED = true;

	bool RESTORE = true;//false;

						// If new error is larger than previously restore old values
						if(fabs(error) > fabs(K_loss_coeffs[steadyCurve][k][lblError]) && RESTORE)
						{
							pPpt->Out(", but is larger than before. Restore previous values.\n");
							K_loss_coeffs[steadyCurve][k][lblK] = K_loss_coeffs[steadyCurve][k][lblKPrev];
							K_loss_coeffs[steadyCurve][k][lblKPrev] = K_loss_coeffs[steadyCurve][k][lblKPrevPrev];
							K_loss_coeffs[steadyCurve][k][lblError] = K_loss_coeffs[steadyCurve][k][lblErrorPrev];
							error = K_loss_coeffs[steadyCurve][k][lblError];
							K_loss_coeffs[steadyCurve][k][lblErrorPrev] = K_loss_coeffs[steadyCurve][k][lblErrorPrevPrev];
						}
						else	// Continue to update values
						{
							if(RESTORE) pPpt->Out(", but is smaller than before. Halve factor but continue to update values.\n");
							else pPpt->Out(", so halve factor.\n");
							K_loss_coeffs[steadyCurve][k][lblErrorPrevPrev] = K_loss_coeffs[steadyCurve][k][lblErrorPrev];
							K_loss_coeffs[steadyCurve][k][lblErrorPrev] = K_loss_coeffs[steadyCurve][k][lblError];
							K_loss_coeffs[steadyCurve][k][lblError] = error; //interpMFP - measuredMFP;
						}

						//if(K_loss_coeffs[steadyCurve][k][lblFactor] > 1) K_loss_coeffs[steadyCurve][k][lblFactor] = 1; // Reset to the standard value
						
						K_loss_coeffs[steadyCurve][k][lblFactor] *= 0.5; // Halve the movement in both cases since changed sign
					}
					else	// Continue to update values
					{
						pPpt->Out("Error has not changed sign. Continue to update values.\n");
						K_loss_coeffs[steadyCurve][k][lblErrorPrevPrev] = K_loss_coeffs[steadyCurve][k][lblErrorPrev];
						K_loss_coeffs[steadyCurve][k][lblErrorPrev] = K_loss_coeffs[steadyCurve][k][lblError];
						K_loss_coeffs[steadyCurve][k][lblError] = error;
					}

					if(fabs(error) > fabs(error_max))
					{
						//pPpt->Out("k = "); pPpt->Out(k); pPpt->Out(": error("); pPpt->Out(error); pPpt->Out(") > error_max ("); pPpt->Out(error_max); pPpt->Out(")\n");
						error_max = error;
						error_max_percent = fabs(error/interpMFP)*100;
						max_error_element = k;	// Track the largest error
					}

					// Adjust K by false position method
					K_loss_coeffs[steadyCurve][k][lblKPrevPrev] = K_loss_coeffs[steadyCurve][k][lblKPrev];	// Save previous K value
					K_loss_coeffs[steadyCurve][k][lblKPrev] = K_loss_coeffs[steadyCurve][k][lblK];			// Save previous K value
					double del_K = 1.0;//0.5;
					//if(fabs(error) > tol_MFP) // Only adjust if this element has not met the convergence criteria
					{
						if(error<0) K_loss_coeffs[steadyCurve][k][lblK] += K_loss_coeffs[steadyCurve][k][lblFactor]*del_K;	// Increase K
						else K_loss_coeffs[steadyCurve][k][lblK] -= K_loss_coeffs[steadyCurve][k][lblFactor]*del_K;			// Decrease K
						if(K_loss_coeffs[steadyCurve][k][lblK] < 0) K_loss_coeffs[steadyCurve][k][lblK] = 0;				// Prohibit K from going negative
					}

					//pPpt->Out("After adjustment:\n");
					pPpt->Out("M2 = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblM2]); 
					pPpt->Out("\tK = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblK]); 
					pPpt->Out("\tK (prev) = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKPrev]); 
					pPpt->Out("\tK (prev prev) = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKPrevPrev]); 
					pPpt->Out("\tError = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblError]); 
					pPpt->Out("\tError (previous) = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblErrorPrev]); 
					pPpt->Out("\tError (prev prev) = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblErrorPrevPrev]); 
					pPpt->Out("\tFactor = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblFactor]);
					//pPpt->Out("\tU/Cis = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKUCs]);
					//pPpt->Out("\tEta = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKEta]);
					pPpt->Out("\n");
					pPpt->Out("\n");
				}

				pPpt->Out("\n");
				pPpt->Out("Finished pass "); pPpt->Out(pass_counter); pPpt->Out(" for speed curve ["); pPpt->Out(steadyCurve); pPpt->Out("] of [0-"); pPpt->Out(numSteadySpeeds - 1); pPpt->Out("]\n");
				pPpt->Out("Element ["); pPpt->Out(max_error_element); pPpt->Out("] has greatest error "); pPpt->Out(error_max); pPpt->Out(", or "); pPpt->Out(error_max_percent); pPpt->Out("% on pass "); pPpt->Out(pass_counter); pPpt->Out(":\n");
				pPpt->Out("M2 = "); pPpt->Out(K_loss_coeffs[steadyCurve][max_error_element][lblM2]); 
				pPpt->Out("\tK = "); pPpt->Out(K_loss_coeffs[steadyCurve][max_error_element][lblK]); 
				pPpt->Out("\tK (prev) = "); pPpt->Out(K_loss_coeffs[steadyCurve][max_error_element][lblKPrev]); 
				pPpt->Out("\tK (prev prev) = "); pPpt->Out(K_loss_coeffs[steadyCurve][max_error_element][lblKPrevPrev]); 
				pPpt->Out("\tError = "); pPpt->Out(K_loss_coeffs[steadyCurve][max_error_element][lblError]); 
				pPpt->Out("\tError (previous) = "); pPpt->Out(K_loss_coeffs[steadyCurve][max_error_element][lblErrorPrev]); 
				pPpt->Out("\tError (prev prev) = "); pPpt->Out(K_loss_coeffs[steadyCurve][max_error_element][lblErrorPrevPrev]); 
				pPpt->Out("\tFactor = "); pPpt->Out(K_loss_coeffs[steadyCurve][max_error_element][lblFactor]);
				pPpt->Out("\tUCs = "); pPpt->Out(K_loss_coeffs[steadyCurve][max_error_element][lblKUCs]);
				pPpt->Out("\tEta = "); pPpt->Out(K_loss_coeffs[steadyCurve][max_error_element][lblKEta]);
				pPpt->Out("\n");
				pPpt->Out("\n");
//*/

//				pPpt->Out("\n");
//				pPpt->Out("Finished pass of speed curve ["); pPpt->Out(steadyCurve); pPpt->Out("] of [0-"); pPpt->Out(numSteadySpeeds - 1); pPpt->Out("]\n");
//				pPpt->Out("Element ["); pPpt->Out(max_error_element); pPpt->Out("] has greatest error "); pPpt->Out(error_max); pPpt->Out(", or "); pPpt->Out(error_max_percent); pPpt->Out("% on pass "); pPpt->Out(pass_counter); pPpt->Out(":\n");
				
				if(error_max_percent < tol_MFP)
				{
					pPpt->Out("Finished matching steady characteristic to within specified tolerance, tol_MFP = "); pPpt->Out(tol_MFP); pPpt->Out("%\n");
					pPpt->Out("Element ["); pPpt->Out(max_error_element); pPpt->Out("] has the error_max_percent = "); 
					pPpt->Out(error_max_percent);
					pPpt->Out("%\n\n");

					//for(k=0; k<K_loss_coeffs.size(); ++k)
					for(k=0; k<numSteadyPoints[steadyCurve]; ++k)
					{
						pPpt->Out("M2 = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblM2]); 
						pPpt->Out("\tK = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblK]); 
						pPpt->Out("\tK (prev) = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKPrev]); 
						pPpt->Out("\tK (prev prev) = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKPrevPrev]); 
						pPpt->Out("\tError = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblError]); 
						pPpt->Out("\tError (previous) = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblErrorPrev]); 
						pPpt->Out("\tError (prev prev) = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblErrorPrevPrev]); 
						pPpt->Out("\tFactor = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblFactor]);
						pPpt->Out("\tU/Cis = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKUCs]);
						pPpt->Out("\teta = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKEta]);
						pPpt->Out("\tD = "); pPpt->Out(K_loss_coeffs[steadyCurve][k][lblKD]);
						pPpt->Out("\n");
					}

					// Put some space in output files to signify completion of calibration for this curve
					fprintf(pCalibrationProgressMFPPR,"Finished calibration of speed curve [");
					fprintf(pCalibrationProgressMFPPR,"%d", steadyCurve);
					fprintf(pCalibrationProgressMFPPR,"]\n\n");

					fprintf(pCalibrationProgressK,"Finished calibration of speed curve [");
					fprintf(pCalibrationProgressK,"%d", steadyCurve);
					fprintf(pCalibrationProgressK,"]\n\n");

					// Output final steady mass flow characteristic, calibrated K loss curve, and K loss restart data, for the current speed to separate files
					for(k=0; k<numSteadyPoints[steadyCurve]; ++k)
					{
						fprintf(pCalibrationFinalMFPPR,"%f\t%f\t%f\t%f\t%f", 
							steadyPRMFP[steadyCurve][k][lblSpeed], 
								MFP_vs_PR[k][lblMFPvsPR_PR], 
									MFP_vs_PR[k][lblMFPvsPR_MFP], 
										MFP_vs_PR[k][lblMFPvsPR_UCs],
											MFP_vs_PR[k][lblMFPvsPR_eta]);

						fprintf(pCalibrationFinalK,"%f\t%f\t%f\t%f\t%f\t%f", 
							steadyPRMFP[steadyCurve][k][lblSpeed], 
								K_loss_coeffs[steadyCurve][k][lblM2], 
									K_loss_coeffs[steadyCurve][k][lblK],
										K_loss_coeffs[steadyCurve][k][lblKUCs],
											K_loss_coeffs[steadyCurve][k][lblKEta],
												D_temp[steadyCurve][k][0]);

						fprintf(pCalibrationFinalKDynasty,"%f\t%f\t%f\t%f\t%f", 
							steadyPRMFP[steadyCurve][k][lblSpeed], 
								K_loss_coeffs[steadyCurve][k][lblM2], 
									K_loss_coeffs[steadyCurve][k][lblK]*(2/pPpt->gamma_air),
										K_loss_coeffs[steadyCurve][k][lblKUCs],
											K_loss_coeffs[steadyCurve][k][lblKEta]);
										
						fprintf(pFinalKRestart,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", 
							steadyPRMFP[steadyCurve][k][lblSpeed], 
								K_loss_coeffs[steadyCurve][k][lblM2], 
									K_loss_coeffs[steadyCurve][k][lblK],
										K_loss_coeffs[steadyCurve][k][lblKPrev],
											K_loss_coeffs[steadyCurve][k][lblKPrevPrev],
												K_loss_coeffs[steadyCurve][k][lblError],
													K_loss_coeffs[steadyCurve][k][lblErrorPrev],
														K_loss_coeffs[steadyCurve][k][lblErrorPrevPrev],
															K_loss_coeffs[steadyCurve][k][lblFactor],
																K_loss_coeffs[steadyCurve][k][lblKUCs],
																	K_loss_coeffs[steadyCurve][k][lblKEta],
																		D_temp[steadyCurve][k][0]);

						if(k < numSteadyPoints[steadyCurve] - 1)
						{
							fprintf(pCalibrationFinalMFPPR,"\n");
							fprintf(pCalibrationFinalK,"\n");
              fprintf(pCalibrationFinalKDynasty,"\n");
							fprintf(pFinalKRestart,"\n");
						}
					}
					if(steadyCurve < numSteadySpeeds - 1)
					{
						fprintf(pCalibrationFinalMFPPR,"\n");
						fprintf(pCalibrationFinalK,"\n");
            fprintf(pCalibrationFinalKDynasty,"\n");
						fprintf(pFinalKRestart,"\n");
					}

					if(steadyCurve < numSteadySpeeds - 1) // There are still speed curves left to calibrate
					{
						++steadyCurve;			// Move onto the next speed curve
						
						// Carry out resets
						pass_counter = 0;		
						counter_between_steady = 0;
						MFP_vs_PR.erase(MFP_vs_PR.begin(), MFP_vs_PR.end()); // Delete contents of MFP_vs_PR

						// Carry on from the position at the end of the previous curve calibration, but check:
						if(PASS_INCREASING) row = 0;
						else row = numSteadyPoints[steadyCurve] - 1;

						PR_point = steadyPRMFP[steadyCurve][row][lblPR];
						MFP_point = steadyPRMFP[steadyCurve][row][lblMFP];
						UCs_point = steadyPRMFP[steadyCurve][row][lblUCs];
						p0_applied = PR_point;//*(pPipe[1]->Node[pPipe[1]->N-1].p_dash*pPpt->PREF); // p01=PR*p2
						for(int e=0; e<nEntries; ++e)
            {
              pEntryEndEnv[e]->Set_P0_gradual(p0_applied);
//cout << "Point has converged on this pass, last point on line. Setting " << pEntryEndEnv[e]->Identify() << " to p0_applied = " << p0_applied << ", pEntryEndEnv[e]->Get_P0 = " << pEntryEndEnv[e]->Get_P0() << endl;
//exit(1);
            }
						
						pPpt->Out("\n");
						pPpt->Out("Moving on to calibrate speed curve ["); pPpt->Out(steadyCurve); pPpt->Out("] of [0-"); pPpt->Out(numSteadySpeeds-1); pPpt->Out("]\n");
						pPpt->Out("\n");
					}
					else
					{
						// Close output files
						fclose(pCalibrationProgressMFPPR);
						fclose(pCalibrationFinalMFPPR);
						fclose(pCalibrationProgressK);
						fclose(pCalibrationFinalK);
            fclose(pCalibrationFinalKDynasty);
						fclose(pFinalKRestart);
						pPpt->Out("\nFinished calibrating all speed curves; closing output files.\nExiting...\n");
						exit(1);
					}
				}
				else
				{					
					pPpt->Out("Continuing onto a new pass, since error_max_percent = "); 
					pPpt->Out(error_max_percent);
					pPpt->Out("%, but tol_MFP = ");
					pPpt->Out(tol_MFP);
					pPpt->Out("\n\n");

					// Continue with another pass - first delete contents of MFP_vs_PR
					MFP_vs_PR.erase(MFP_vs_PR.begin(), MFP_vs_PR.end());

					// Switch direction upon completion of each pass
					if(PASS_INCREASING) PASS_INCREASING = false; 
					else PASS_INCREASING = true;
				}
			}
		}
		//delete [] temp; // If you delete temp the final MFP PR file doesn't write properly
	}
	return;
}

void CAPLDev::PrintToFile(CProperties* pPpt, int timestep, double time, double ca)
// ============================================================ //
// Prints instantaneous boundary data to file.					//
// This function is called from the main function.				//
// ============================================================ //
{
	if(CALIBRATE)
	{
		if(timestep==0)
		{
			if(EX)
			{
				if(DEVICE==1)
				{
					//fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR EXHAUST GAUZE BOUNDARY", "-", ID));
					if(!NAPLDev)
						fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR EXHAUST APL DEVICE BOUNDARY", "-", ID));
					else
						fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR EXHAUST NAPL DEVICE BOUNDARY", "-", ID));
				}
				else
				{
					if(DEVICE==2) fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR EXHAUST THROTTLE BOUNDARY", "-", ID));
					else
					{
						if(DEVICE==3) if(DEVICE==2) fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR EXHAUST EGR VALVE BOUNDARY", "-", ID));
						else fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR EXHAUST UNKNOWN APL DEVICE BOUNDARY", "-", ID));
					}
				}
			}
			else
			{
				if(DEVICE==1) fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR INTAKE GAUZE BOUNDARY", "-", ID));
				else
				{
					if(DEVICE==2) fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR INTAKE THROTTLE BOUNDARY", "-", ID));
					else
					{
						if(DEVICE==3) if(DEVICE==2) fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR INTAKE EGR VALVE BOUNDARY", "-", ID));
						else fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR INTAKE UNKNOWN APL DEVICE BOUNDARY", "-", ID));
					}
				}
			}
			fprintf(OUTPUT_FILE,"\n%s\n", "CALIBRATION --- NO RESULTS");
		}
	}
	else
	{
		int entryID, admissionMapNo;

		if(timestep==0)
		{
			if(EX)
			{
				if(DEVICE==1)
				{
					//fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR EXHAUST GAUZE BOUNDARY", "-", ID));
					if(!NAPLDev)
						fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR EXHAUST APL DEVICE BOUNDARY", "-", ID));
					else
						fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR EXHAUST NAPL DEVICE BOUNDARY", "-", ID));
				}
				else
				{
					if(DEVICE==2) fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR EXHAUST THROTTLE BOUNDARY", "-", ID));
					else
					{
						if(DEVICE==3) if(DEVICE==2) fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR EXHAUST EGR VALVE BOUNDARY", "-", ID));
						else fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR EXHAUST UNKNOWN APL DEVICE BOUNDARY", "-", ID));
					}
				}
			}
			else
			{
				if(DEVICE==1) fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR INTAKE GAUZE BOUNDARY", "-", ID));
				else
				{
					if(DEVICE==2) fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR INTAKE THROTTLE BOUNDARY", "-", ID));
					else
					{
						if(DEVICE==3) if(DEVICE==2) fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR INTAKE EGR VALVE BOUNDARY", "-", ID));
						else fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR INTAKE UNKNOWN APL DEVICE BOUNDARY", "-", ID));
					}
				}
			}
			fprintf(OUTPUT_FILE,"\n");

			if(DEVICE==1)
			{
				if(!NAPLDev)
				{
				fprintf(OUTPUT_FILE,"%s\t\t", "       ");
			  //fprintf(OUTPUT_FILE,"%s\t\t", "       ");
				fprintf(OUTPUT_FILE,"%s\t\t\t", "                ");
			  //fprintf(OUTPUT_FILE,"%s\t\t\t", "                ");
				fprintf(OUTPUT_FILE,"\t\t%s\t\t\t\t\t\t\t\t\t", "                    ");
			  //fprintf(OUTPUT_FILE,"\t\t%s\t\t\t\t\t\t\t\t\t\t\t", "K LOSS MAP OPERATION");
				fprintf(OUTPUT_FILE,"%s\n", "PERFORMANCE");
			  //fprintf(OUTPUT_FILE,"%s\n", "PERFORMANCE");


				fprintf(OUTPUT_FILE,"%s\t\t", "       ");
			  //fprintf(OUTPUT_FILE,"%s\t\t", "TIMING ");
				fprintf(OUTPUT_FILE,"%s\t\t\t", "                ");
			  //fprintf(OUTPUT_FILE,"%s\t\t\t", "APLDev OPERATION");
				fprintf(OUTPUT_FILE,"\t\t%s\t\t\t\t\t\t\t\t\t", "K LOSS MAP OPERATION");
				for(entryID=0; entryID<nEntries; ++entryID)
				{
					fprintf(OUTPUT_FILE,"%s%d%s%d%s\t\t\t\t\t\t\t\t\t\t\t", "ENTRY[", entryID, "] OF [0-", nEntries-1, "]");
					fprintf(OUTPUT_FILE,"%s%s%s\t\t\t\t\t\t\t\t", "       ", "         ", " ");
				}
				fprintf(OUTPUT_FILE,"%s\t\t\t\t\t\t\t\t\t\t%s\n", "ENTRIES COMBINED", "POWER PARAMETERS");

				fprintf(OUTPUT_FILE,"%s\t\t\t\t", "TIMING ");
	  //fprintf(OUTPUT_FILE,"%s\t%s\t%s\t\t", "Time(s)", "Time in cycle (s)", "Position in cycle (degrees)");
				fprintf(OUTPUT_FILE,"%s\t\t\t", "APLDev OPERATION");
			  //fprintf(OUTPUT_FILE,"%s\t%s\t\t", "PARAM_S", "PARAM_R");
				for(admissionMapNo=0; admissionMapNo<(nAdmissionMaps); ++admissionMapNo)
				{
					if(admissionMapNo==0) fprintf(OUTPUT_FILE,"%s\t\t\t", "FULL ADMISSION");
					else
					{
						if(admissionMapNo==1) fprintf(OUTPUT_FILE,"%s", "1ST");
						else fprintf(OUTPUT_FILE,"%s", "2ND");
						fprintf(OUTPUT_FILE,"%s\t\t", " PARTIAL");
					}
				}
				fprintf(OUTPUT_FILE,"%s\t\t\t\t\t\t", "COMBINED");
				for(entryID=0; entryID<nEntries; ++entryID)
				{
					fprintf(OUTPUT_FILE,"%s\t\t\t\t\t\t\t\t\t\t\t", "LIMB INLET");
					fprintf(OUTPUT_FILE,"%s\t\t\t\t\t\t\t\t", "LIMB EXIT");
				}
				//fprintf(OUTPUT_FILE,"%s\t\t\t\t\t\t\t\t\t\t\t", "LIMB INLET");
				//fprintf(OUTPUT_FILE,"%s\n", "LIMB EXIT");
				fprintf(OUTPUT_FILE,"\n");
				
				fprintf(OUTPUT_FILE,"%s\t%s\t%s\t\t", "Time(s)", "Time in cycle (s)", "Position in cycle (degrees)");
				fprintf(OUTPUT_FILE,"%s\t%s\t\t", "PARAM_S", "PARAM_R");
				for(admissionMapNo=0; admissionMapNo<(nAdmissionMaps); ++admissionMapNo)
				{
					fprintf(OUTPUT_FILE,"%s%d%s\t%s%d%s\t%s%d%s\t", "K[", admissionMapNo, "]", "OUTSIDE_RANGE[", admissionMapNo, "]", "OUTSIDE_RANGE[", admissionMapNo, "]");
				}
				fprintf(OUTPUT_FILE,"%s\t%s\t%s\t%s\t\t", "lambda[EXIT] (-)", "K_OVERALL", "OUTSIDE_RANGE_OVERALL", "OUTSIDE_RANGE_OVERALL");
				fprintf(OUTPUT_FILE,"%s", "\t");

				for(entryID=0; entryID<nEntries; ++entryID)
				{
					// LIMB INLET (flow plus performance parameters PR, MFP, as these are measured from the LIMB INLET)
					fprintf(OUTPUT_FILE,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t", "p (bar)", "T (K)", "m (kg/s)", "u (m/s)", "p0 (bar)", "T0 (K)", "PR (-)", "MFP (kg/s.K^0.5/bar)", "U/Cs (-)", "Ws (kW)");
					// LIMB EXIT (flow parameters only)
					fprintf(OUTPUT_FILE,"%s\t%s\t%s\t%s\t%s\t%s\t\t\t", "p (bar)", "T (K)", "m (kg/s)", "u (m/s)", "p0 (bar)", "T0 (bar)");
				}
				fprintf(OUTPUT_FILE,"%s\t%s\t%s\t", "m_stage[INLETS] (kg/s)", "m_stage[EXITS] (kg/s)", "m_stage[ROTOR] (kg/s)");
				fprintf(OUTPUT_FILE,"%s\t%s\t", "lambda[INLET] (-)", "lambda[EXIT] (-)");
				if(averagePR) fprintf(OUTPUT_FILE,"%s\t", "PR_stage_mean (-)");
				else fprintf(OUTPUT_FILE,"%s%d%s\t", "PR, ENTRY[", entryPR, "] (-)");
				fprintf(OUTPUT_FILE,"%s\t%s\t%s\t\t", "MFP_stage (kg/s.K^0.5/bar)", "UCs_stage_mean (-)", "Ws_stage_inlet (kW)");
				fprintf(OUTPUT_FILE,"%s\t%s\t%s\t%s", "UCs_stage_shaft (-)", "Ws_stage_shaft (kW)", "eta_OVERALL (-)", "Wact_stage_shaft (kW)");
			
				}
				else	// NAPLDev formulation
				{
					fprintf(OUTPUT_FILE,"%s\t\t", "       ");
					fprintf(OUTPUT_FILE,"%s\t\t\t", "                ");
					fprintf(OUTPUT_FILE,"\t\t%s\t\t\t\t\t\t\t\t\t\t\t", "                    ");
					fprintf(OUTPUT_FILE,"%s\n", "PERFORMANCE");
					
					fprintf(OUTPUT_FILE,"%s\t\t", "       ");
					fprintf(OUTPUT_FILE,"%s\t\t\t", "                ");
					fprintf(OUTPUT_FILE,"\t\t%s\t\t\t\t\t\t\t\t\t\t\t", "K LOSS MAP OPERATION");
					for(entryID=0; entryID<nEntries; ++entryID)
					{
						fprintf(OUTPUT_FILE,"%s%d%s%d%s\t\t\t\t\t\t\t\t\t\t\t", "ENTRY[", entryID, "] OF [0-", nEntries-1, "]");
						fprintf(OUTPUT_FILE,"%s%s%s\t\t\t\t\t\t\t\t", "       ", "         ", " ");
					}
					fprintf(OUTPUT_FILE,"%s\t\t\t\t\t\t\t\t\t\t\t%s\n", "ENTRIES COMBINED", "POWER PARAMETERS");

					fprintf(OUTPUT_FILE,"%s\t\t\t\t", "TIMING ");
					fprintf(OUTPUT_FILE,"%s\t\t\t", "APLDev OPERATION");
					for(admissionMapNo=0; admissionMapNo<(nAdmissionMaps); ++admissionMapNo)
					{
						if(admissionMapNo==0) fprintf(OUTPUT_FILE,"%s\t\t\t\t", "FULL ADMISSION");
						else
						{
							if(admissionMapNo==1) fprintf(OUTPUT_FILE,"%s", "1ST");
							else fprintf(OUTPUT_FILE,"%s", "2ND");
							fprintf(OUTPUT_FILE,"%s\t\t", " PARTIAL");
						}
					}
					fprintf(OUTPUT_FILE,"%s\t\t\t\t\t\t\t", "COMBINED");
					for(entryID=0; entryID<nEntries; ++entryID)
					{
						fprintf(OUTPUT_FILE,"%s\t\t\t\t\t\t\t\t\t\t\t", "LIMB INLET");
						fprintf(OUTPUT_FILE,"%s\t\t\t\t\t\t\t\t", "LIMB EXIT");
					}
					fprintf(OUTPUT_FILE,"\n");
				
					fprintf(OUTPUT_FILE,"%s\t%s\t%s\t\t", "Time(s)", "Time in cycle (s)", "Position in cycle (degrees)");
					fprintf(OUTPUT_FILE,"%s\t%s\t\t", "PARAM_S", "PARAM_R");
					for(admissionMapNo=0; admissionMapNo<(nAdmissionMaps); ++admissionMapNo)
					{
						fprintf(OUTPUT_FILE,"%s%d%s\t%s%d%s\t%s%d%s\t%s%d%s\t", "K[", admissionMapNo, "]", "D[", admissionMapNo, "]", "OUTSIDE_RANGE[", admissionMapNo, "]", "OUTSIDE_RANGE[", admissionMapNo, "]");
					}
					fprintf(OUTPUT_FILE,"%s\t%s\t%s\t%s\t%s\t\t", "lambda[EXIT] (-)", "K_OVERALL", "D_OVERALL", "OUTSIDE_RANGE_OVERALL", "OUTSIDE_RANGE_OVERALL");
					fprintf(OUTPUT_FILE,"%s", "\t");

					for(entryID=0; entryID<nEntries; ++entryID)
					{
						// LIMB INLET (flow plus performance parameters PR, MFP, as these are measured from the LIMB INLET)
						fprintf(OUTPUT_FILE,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t", "p (bar)", "T (K)", "m (kg/s)", "u (m/s)", "p0 (bar)", "T0 (K)", "PR (-)", "MFP (kg/s.K^0.5/bar)", "U/Cs (-)", "Ws (kW)");
						// LIMB EXIT (flow parameters only)
						fprintf(OUTPUT_FILE,"%s\t%s\t%s\t%s\t%s\t%s\t\t\t", "p (bar)", "T (K)", "m (kg/s)", "u (m/s)", "p0 (bar)", "T0 (bar)");
					}
					fprintf(OUTPUT_FILE,"%s\t%s\t%s\t", "m_stage[INLETS] (kg/s)", "m_stage[EXITS] (kg/s)", "m_stage[ROTOR] (kg/s)");
					fprintf(OUTPUT_FILE,"%s\t%s\t", "lambda[INLET] (-)", "lambda[EXIT] (-)");
					if(averagePR) fprintf(OUTPUT_FILE,"%s\t", "PR_stage_mean (-)");
					else fprintf(OUTPUT_FILE,"%s%d%s\t", "PR, ENTRY[", entryPR, "] (-)");
					fprintf(OUTPUT_FILE,"%s\t%s\t%s\t%s\t\t", "MFP_stage (kg/s.K^0.5/bar)", "UCs_stage_mean (-)", "Ws_stage_inlet (kW)", "Wact_turb (kW)");
					fprintf(OUTPUT_FILE,"%s\t%s\t%s\t%s", "UCs_stage_shaft (-)", "Ws_stage_shaft (kW)", "eta_OVERALL (-)", "Wact_stage_shaft (kW)");
				}
			}
			else
			{
				fprintf(OUTPUT_FILE,"%s\t%s\t%s\t\t", "Time(s)", "Time in cycle (s)", "Position in cycle (degrees)");
			}
			fprintf(OUTPUT_FILE,"\n");
		}
		
		if(timestep%freq==0 && time>=print_from_time) // Only print at the specified sampling frequency
		{
			timeInCyclePrev = timeInCycle; // Save previous time in cycle for comparsion
			timeInCycle = fmod(time,1/freqPulse); // Current time in cycle
			if(timeInCycle < timeInCyclePrev) fprintf(OUTPUT_FILE,"\n\n"); // Separate cycles

			fprintf(OUTPUT_FILE,"%f\t", time); // Absolute simulation time
			fprintf(OUTPUT_FILE,"%f\t", fmod(time,1/freqPulse)); // Time in pulse of frequency freqPulse
			fprintf(OUTPUT_FILE,"%f\t\t", (fmod(time,1/freqPulse)/(1/freqPulse))*degCycle); // Position in pulse cycle (degrees)
			fprintf(OUTPUT_FILE,"%f\t%f\t\t", PARAM_S, PARAM_R);
			
			for(admissionMapNo=0; admissionMapNo<(nAdmissionMaps); ++admissionMapNo)
			{
				fprintf(OUTPUT_FILE,"%f\t%f\t%s\t%d\t", K[admissionMapNo], D[admissionMapNo], TrueOrFalse(OUTSIDE_RANGE[admissionMapNo]), int(OUTSIDE_RANGE[admissionMapNo]));
			}
			fprintf(OUTPUT_FILE,"%f\t%f\t%f\t%s\t%d\t\t\t", lambda[LIMB_EXIT], K_OVERALL, D_OVERALL, TrueOrFalse(OUTSIDE_RANGE_OVERALL), int(OUTSIDE_RANGE_OVERALL));
			
			if(DEVICE==1)
			{
				int entryID;
				for(entryID=0; entryID<nEntries; ++entryID)
				{
					// LIMB INLET (flow plus performance parameters PR, MFP, as these are measured from the LIMB INLET)
					fprintf(OUTPUT_FILE,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t\t", p[entryID][LIMB_INLET], T[entryID][LIMB_INLET], m[entryID][LIMB_INLET], u[entryID][LIMB_INLET], p0[entryID][LIMB_INLET], T0[entryID][LIMB_INLET], PR[entryID], MFP[entryID], UCs[entryID], Ws[entryID]/1000);
					// LIMB EXIT (flow parameters only)
					fprintf(OUTPUT_FILE,"%f\t%f\t%f\t%f\t%f\t%f\t\t\t", p[entryID][LIMB_INLET], T[entryID][LIMB_INLET], m[entryID][LIMB_INLET], u[entryID][LIMB_INLET], p0[entryID][LIMB_EXIT], T0[entryID][LIMB_EXIT]);
				}
				fprintf(OUTPUT_FILE,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t\t", m_stage[LIMB_INLET], m_stage[LIMB_EXIT], m_stage[ROTOR], lambda[LIMB_INLET], lambda[LIMB_EXIT], PR_stage_mean, MFP_stage, UCs_stage_mean, Ws_stage_inlet/1000, Wact_turb/1000);
				fprintf(OUTPUT_FILE,"%f\t%f\t%f\t%f\t", UCs_stage_shaft, Ws_stage_shaft/1000, eta_OVERALL, Wact_stage_shaft/1000);
				fprintf(OUTPUT_FILE,"%f\t%f\t%f\t", m_dot_inner, m_dot_outer, lambda[LIMB_EXIT]);
			}
			fprintf(OUTPUT_FILE,"\n");
		}
	}


  if(PRINT_DEBUG_FILE)
  {
    if(timestep==0)
    {
			std::string temp_str = "Debug file for ";
			temp_str += Identify();
			fprintf(OUTPUT_FILE_DEBUG,"%s\n", Underline(StringToChar(temp_str), "-"));	
			fprintf(OUTPUT_FILE_DEBUG,"%s", "Time(s)");
/*
			if(!pPpt->CONTINUOUS) 
      {
        fprintf(OUTPUT_FILE_DEBUG,"\t%s%c%s\t%s%c%s", "Elapsed (", Deg(), "CA)", "Crank angle (", Deg(), "CA)");
      }
//*/
      //fprintf(OUTPUT_FILE_DEBUG,"\t%s\t%s\t%s\t%s\t%s\t%s", "pCLIN[ONE_SIDE]", "pCLOUT[ONE_SIDE]", "Direction", "Direction No.", "Choked", "Choked No.");
      fprintf(OUTPUT_FILE_DEBUG,"\t%s\t%s\t%s", "counter", "error", "minError");
    	fprintf(OUTPUT_FILE_DEBUG,"\n");
		}

    if(timestep%freq==0)
    { 
      fprintf(OUTPUT_FILE_DEBUG,"%f", time);
/*      
      if(!pPpt->CONTINUOUS) 
      {
        fprintf(OUTPUT_FILE_DEBUG,"\t%f\t%f", ca_elapsed, ca);// Periodic
      }
//*/
			//fprintf(OUTPUT_FILE_DEBUG,"\t%f\t%f\t%s\t%d\t%s\t%d", (*(pCLIN[ONE_SIDE]))[R+1], (*(pCLOUT[ONE_SIDE]))[R+1], InflowOrOutflow(*(pend_flow[ONE_SIDE])), InflowOrOutflowNumber(*(pend_flow[ONE_SIDE])), Choked(pBN[ONE_SIDE]->CHOKED), ChokedNo(pBN[ONE_SIDE]->CHOKED));
      fprintf(OUTPUT_FILE_DEBUG,"\t%d\t%f\t%f", counter, error, minError);
			fprintf(OUTPUT_FILE_DEBUG,"\n");
    }
  }
}

void CAPLDev::PrintToScreen(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".PrintToScreen\n");}

	pPpt->Out(Underline(Identify(), "=", "", strDesc)); pPpt->Out("\n");
	//pPpt->Out("\n");
	//pPpt->Out("Device type = "); pPpt->Out(Get_DEVICE()); pPpt->Out("\n");
	//pPpt->Out("Curve parameter value (e.g., solidity), PARAM_S = "); pPpt->Out(PARAM_S); pPpt->Out("\n");

	//if(INTERP_M2) pPpt->Out("Based on a downstream Mach number M2, PARAM_R = "); 
	//else pPpt->Out("Based on an upstream Reynolds number, PARAM_R = ");
	//pPpt->Out(PARAM_R); pPpt->Out("\n");
/*	
	if(INTERP_M2)
	{
		pPpt->Out("Based on a downstream Mach number M2, M[DOWNSTREAM] = ");
		pPpt->Out(M[DOWNSTREAM]); 
	}
	else
	{
		pPpt->Out("Based on an upstream Reynolds number, M[UPSTREAM] = ");
		pPpt->Out(M[UPSTREAM]); 
	}
	pPpt->Out("\n");
*/
	//pPpt->Out("Dynamic throttle factor, X = "); pPpt->Out(X); pPpt->Out("\n");
	//pPpt->Out("Interpolated loss coefficient, K = "); pPpt->Out(K); pPpt->Out("\n");
	
	if(CALIBRATE) pPpt->Out("Undergoing CALIBRATION\n");
	else
	{
		if(INTERP_M2) pPpt->Out("Downstream Mach number, M2\t\t\t\t=\t");
		else pPpt->Out("Upstream Reynolds number, Re1\t\t\t=\t");
		pPpt->Out(PARAM_R); pPpt->Out("\n");
    pPpt->Out("Instantaneous loss coefficient, K_OVERALL\t\t=\t"); pPpt->Out(K_OVERALL); pPpt->Out("\n");
    pPpt->Out("Velocity head (upstream), velHead\t\t\t=\t"); pPpt->Out(velHead); pPpt->Out(" bar\n");
    pPpt->Out("Static pressure drop, delP\t\t\t\t=\t"); pPpt->Out(delP); pPpt->Out(" bar\n");
    pPpt->Out("Upstream static pressure, pUs\t\t\t\t=\t"); pPpt->Out(pUs); pPpt->Out(" bar\n");
    pPpt->Out("Final loop counter, counter\t\t\t\t=\t"); pPpt->Out(counter); pPpt->Out("\n");
    pPpt->Out("\n");
    pPpt->Out(Underline("Turbine specific parameters", "-", ""));
		pPpt->Out("Rotor admission ratio, lambda[LIMB_EXIT]\t\t=\t"); pPpt->Out(lambda[LIMB_EXIT]); pPpt->Out("\n");
		pPpt->Out("Instantaneous efficiency, eta_OVERALL\t\t\t=\t"); pPpt->Out(eta_OVERALL); pPpt->Out("\n");
		//pPpt->Out("Main loop iterations, counter = "); pPpt->Out(counter); pPpt->Out("\n");
	}
	//pPpt->Out("\n");
	pPpt->Out("Rotor inlet mass flow, mRotorInletJunction\t\t=\t"); pPpt->Out(mRotorInletJunction); pPpt->Out(" kg/s\n");
	pPpt->Out("\n");
	pPpt->Out("\n");
}
/*
void CAPLDev::InterpolateGraph(CProperties* pPpt, int device_graph, double S, double R)
{
	int graph, curve, curve_above, curve_below, point_curve_above, point_curve_below, 
		point_curve_above_above, point_curve_above_below, 
		point_curve_below_above, point_curve_below_below;
	
	double value_curve_above, value_curve_below;

	// Find the appropriate loss coefficient graph
	graph=0;
	while(Coeff[graph][0][0][0] != device_graph) ++graph;
//cout << "Using graph[" << graph << "]\n";

	// Find the right curve
	curve=0;
	while(Coeff[graph][curve][0][1] < S && curve < (num_curves[graph] - 1)) ++curve;
	// The S required is  between curve (curve_above) and curve-1 (curve_below)
	curve_above = curve;
	curve_below = curve-1;
	if(curve==0) curve_below = curve_above;
//cout << "Interploating between curves [" << curve_above << "] and [" << curve_below << "]\n";

	// Find the two points on each of the curves that Re lies between and interpolate each
	// 'curve above'
	point_curve_above=0;
	while(Coeff[graph][curve_above][point_curve_above][2] < R && point_curve_above < (num_points[graph][curve_above] - 1)) ++point_curve_above;
	// Re lies between point_curve_above and point_curve_above-1
	point_curve_above_above = point_curve_above;
	point_curve_above_below = point_curve_above-1;
	if(point_curve_above==0) point_curve_above_below = point_curve_above_above;
	// Interpolate these points
	if(fabs(Coeff[graph][curve_above][point_curve_above_above][2] - Coeff[graph][curve_above][point_curve_above_below][2]) < 1e-6)
	{
		value_curve_above 
		= Coeff[graph][curve_above][point_curve_above_below][3] + 0;
	}
	else
	{
		value_curve_above 
			= Coeff[graph][curve_above][point_curve_above_below][3]
			+ 
			( ((R - Coeff[graph][curve_above][point_curve_above_below][2])
			 /(Coeff[graph][curve_above][point_curve_above_above][2] -Coeff[graph][curve_above][point_curve_above_below][2]))
			*(Coeff[graph][curve_above][point_curve_above_above][3] - Coeff[graph][curve_above][point_curve_above_below][3]) );
	}

	// 'curve below'
	point_curve_below=0;
	while(Coeff[graph][curve_below][point_curve_below][2] < R && point_curve_below < (num_points[graph][curve_below] - 1)) ++point_curve_below;
	// Re lies between point_curve_below and point_curve_below-1
	point_curve_below_above = point_curve_below;
	point_curve_below_below = point_curve_below-1;
	if(point_curve_below==0) point_curve_below_below = point_curve_below_above;
	// Interpolate these points
	if(fabs(Coeff[graph][curve_below][point_curve_below_above][2] - Coeff[graph][curve_below][point_curve_below_below][2]) < 1e-6)
	{
		value_curve_below 
		= Coeff[graph][curve_below][point_curve_below_below][3] + 0;
	}
	else
	{
		value_curve_below 
			= Coeff[graph][curve_below][point_curve_below_below][3]
			+ 
			( ((R - Coeff[graph][curve_below][point_curve_below_below][2])
			/(Coeff[graph][curve_below][point_curve_below_above][2] - Coeff[graph][curve_below][point_curve_below_below][2]))
			*(Coeff[graph][curve_below][point_curve_below_above][3] -Coeff[graph][curve_below][point_curve_below_below][3]) );
	}
	
	// Now interpolate the values based on the curves
	if(fabs(Coeff[graph][curve_above][0][1] - Coeff[graph][curve_below][0][1]) < 1e-6)
	{
		K
		= value_curve_below + 0;
	}
	else
	{
		K
		= value_curve_below
		+
		( ((S - Coeff[graph][curve_below][0][1])
		  /(Coeff[graph][curve_above][0][1] -Coeff[graph][curve_below][0][1]))
		 *(value_curve_above - value_curve_below) );
	}	
	//cout << "CAPLDev::Interpolate:K = " << K << endl;
}
*/
/*
void CAPLDev::QuadInterpolate(CProperties* pPpt, int device_graph, double S, double R)
{
	double A, B, C; 
	double x1, y1, x2, y2, x3, y3;

	int graph, curve, curve_above, curve_below, point_curve_above, point_curve_below, 
		point_curve_above_above, point_curve_above_below, point_curve_above_third, 
		point_curve_below_above, point_curve_below_below, point_curve_below_third;
	
	double value_curve_above, value_curve_below;

	// Find the appropriate loss coefficient graph
	graph=0;
	while(Coeff[graph][0][0][0] != device_graph) ++graph;
//cout << "Using graph[" << graph << "]\n";

	// Find the right curve
	curve=0;
	while(Coeff[graph][curve][0][1] < S) ++curve;
	// The S required is  between curve (curve_above) and curve-1 (curve_below)
	curve_above = curve;
	curve_below = curve-1;
	if(curve==0) curve_below = curve_above;
//cout << "Interploating between curves [" << curve_above << "] and [" << curve_below << "]\n";

	// Find the three points on each of the curves that Re lies between and interpolate each
	// 'curve above'
	point_curve_above=0;
	while(Coeff[graph][curve_above][point_curve_above][2] < R) ++point_curve_above;
	// Re lies between point_curve_above and point_curve_above-1
	point_curve_above_above = point_curve_above;
	point_curve_above_below = point_curve_above-1;

	// Select third point based on position through array (i.e. near begining or near end)
	if(point_curve_above-2 >= 0) point_curve_above_third = point_curve_above-2;
	else point_curve_above_third = point_curve_above + 1;

	if(point_curve_above==0) point_curve_above_below = point_curve_above_above;
	// Interpolate these points quadratically
	// y = Ax^2 + Bx + C
	x1 = Coeff[graph][curve_above][point_curve_above_above][2]; y1 = Coeff[graph][curve_above][point_curve_above_above][3];
	x2 = Coeff[graph][curve_above][point_curve_above_below][2]; y2 = Coeff[graph][curve_above][point_curve_above_below][3];
	x3 = Coeff[graph][curve_above][point_curve_above_third][2]; y3 = Coeff[graph][curve_above][point_curve_above_third][3];
	A = ((y3-y2)/((x3-x2)*(x3-x1))) - ((y1-y2)/((x1-x2)*(x3-x1)));
//	B = ((y3-y2) + A*(pow(x2,2)-pow(x3,2)))/(x2-x3);
	B = ((y1-y2) + A*(pow(x2,2)-pow(x1,2)))/(x1-x2);
//	C = y1 - A*pow(x1,2) - B*x1;
	C = y2 - A*pow(x2,2) - B*x2;
//	C = y3 - A*pow(x3,2) - B*x3;
	value_curve_above = A*pow(R,2) + B*R + C;

//cout << "Re = " << R << endl;
//cout << "x1 = " << x1 << "; y1 = " << y1 << endl;
//cout << "x2 = " << x2 << "; y2 = " << y2 << endl;
//cout << "x3 = " << x3 << "; y3 = " << y3 << endl;	
//cout << "quad curve above: " << A << "x^2 + " << B << "x + " << C << endl;
//cout << "quad value_curve_above = " << value_curve_above << endl;


//	if(fabs(Coeff[graph][curve_above][point_curve_above_above][2] - Coeff[graph][curve_above][point_curve_above_below][2]) < 1e-6)
//	{
//		value_curve_above 
//		= Coeff[graph][curve_above][point_curve_above_below][3] + 0;
//	}
//	else
//	{
//		value_curve_above 
//			= Coeff[graph][curve_above][point_curve_above_below][3]
//			+ 
//			( ((R - Coeff[graph][curve_above][point_curve_above_below][2])
//			 /(Coeff[graph][curve_above][point_curve_above_above][2] -Coeff[graph][curve_above][point_curve_above_below][2]))
//			*(Coeff[graph][curve_above][point_curve_above_above][3] - Coeff[graph][curve_above][point_curve_above_below][3]) );
//	}
//
//cout << "lin value_curve_above = " << value_curve_above << endl;


	// 'curve below'
	point_curve_below=0;
	while(Coeff[graph][curve_below][point_curve_below][2] < R) ++point_curve_below;
	// Re lies between point_curve_below and point_curve_below-1
	point_curve_below_above = point_curve_below;
	point_curve_below_below = point_curve_below-1;

	// Select third point based on position through array (i.e. near begining or near end)
	if(point_curve_below-2 >= 0) point_curve_below_third = point_curve_below-2;
	else point_curve_below_third = point_curve_below + 1;

	if(point_curve_below==0) point_curve_below_below = point_curve_below_above;
	// Interpolate these points quadratically
	// y = Ax^2 + Bx + C
	x1 = Coeff[graph][curve_below][point_curve_below_above][2]; y1 = Coeff[graph][curve_below][point_curve_below_above][3];
	x2 = Coeff[graph][curve_below][point_curve_below_below][2]; y2 = Coeff[graph][curve_below][point_curve_below_below][3];
	x3 = Coeff[graph][curve_below][point_curve_below_third][2]; y3 = Coeff[graph][curve_below][point_curve_below_third][3];
	A = ((y3-y2)/((x3-x2)*(x3-x1))) - ((y1-y2)/((x1-x2)*(x3-x1)));
//	B = ((y3-y2) + A*(pow(x2,2)-pow(x3,2)))/(x2-x3);
	B = ((y1-y2) + A*(pow(x2,2)-pow(x1,2)))/(x1-x2);
//	C = y1 - A*pow(x1,2) - B*x1;
	C = y2 - A*pow(x2,2) - B*x2;
//	C = y3 - A*pow(x3,2) - B*x3;
	value_curve_below = A*pow(R,2) + B*R + C;

//cout << "Re = " << R << endl;
//cout << "x1 = " << x1 << "; y1 = " << y1 << endl;
//cout << "x2 = " << x2 << "; y2 = " << y2 << endl;
//cout << "x3 = " << x3 << "; y3 = " << y3 << endl;	
//cout << "quad curve below: " << A << "x^2 + " << B << "x + " << C << endl;
//cout << "quad value_curve_below = " << value_curve_below << endl;

//	if(fabs(Coeff[graph][curve_below][point_curve_below_above][2] - Coeff[graph][curve_below][point_curve_below_below][2]) < 1e-6)
//	{
//		value_curve_below 
//		= Coeff[graph][curve_below][point_curve_below_below][3] + 0;
//	}
//	else
//	{
//		value_curve_below 
//			= Coeff[graph][curve_below][point_curve_below_below][3]
//			+ 
//			( ((R - Coeff[graph][curve_below][point_curve_below_below][2])
//			/(Coeff[graph][curve_below][point_curve_below_above][2] - Coeff[graph][curve_below][point_curve_below_below][2]))
//			*(Coeff[graph][curve_below][point_curve_below_above][3] -Coeff[graph][curve_below][point_curve_below_below][3]) );
//	}
//
//cout << "lin value_curve_below = " << value_curve_below << endl;

	// Now interpolate the values based on the curves
	if(fabs(Coeff[graph][curve_above][0][1] - Coeff[graph][curve_below][0][1]) < 1e-6)
	{
		K
		= value_curve_below + 0;
	}
	else
	{
		K
		= value_curve_below
		+
		( ((S - Coeff[graph][curve_below][0][1])
		  /(Coeff[graph][curve_above][0][1] -Coeff[graph][curve_below][0][1]))
		 *(value_curve_above - value_curve_below) );
	}	

//cout << "K = " << K << endl;
}
*/
void CAPLDev::ReadResistanceCoefficients(char *InputFile)
{
	float C, C_old, S, S_old, Re, K;
	// C, the type of device?
	// S, a parameter, e.g. solidity of a gauze
	// Re, x-axis, e.g. Reynolds number
	// K, y-axis, e.g. resistance coefficient
	//int graph, num_graphs, curve, point;
	//int *num_curves;
	//int **num_points;

	int max_graphs = 3;
	int max_curves_per_graph = 5;
	num_curves = new int [max_graphs];
	num_points = new int* [max_graphs];

	for(graph=0; graph<max_graphs; ++graph)
		num_points[graph] = new int [max_curves_per_graph];

	FILE *stream;
	stream = fopen(InputFile, "r");
	//cout << InputFile << endl;
		
	if(stream == NULL)
	{
		if(EX) cout << "Exhaust APLDev [" << this->ID << "]: Error opening input file...\n";
		else cout << "Intake APLDev [" << this->ID << "]: Error opening input file...\n";
		exit(1);
	}
	else
	{
		// Set pointer to beginning of file
		fseek(stream, 0L, SEEK_SET);
		// Go through the file to count points per curve, etc.
		graph=0;
		curve=0;
		point=0;
		fscanf(stream, "%f", &C );// The first value of C
		fscanf(stream, "%f", &S );// The first value of S
//		cout << "First C = " << C << endl;
//		cout << "First S = " << S << endl;
		fseek(stream, 0L, SEEK_SET);// Reset pointer to beginning of file
		do
		{
			C_old = C;
			S_old = S;
			fscanf(stream, "%f", &C );
			fscanf(stream, "%f", &S );
			fscanf(stream, "%f", &Re );
			fscanf(stream, "%f", &K );
		
			if((fabs(C_old - C)>1e-4))
			{
				// Changing graphs
				num_points[graph][curve] = point;
//cout << num_points[graph][curve] << " points on curve " << curve << endl;
				++curve;
				num_curves[graph] = curve;
//cout << num_curves[graph] << " curves on graph " << graph << endl;
				point=0;
				curve=0;
				++graph;
			}
			else
			{
				// Changing curves
				if((fabs(S_old - S)>1e-4))
				{ 
					num_points[graph][curve] = point;
					//num_points = point;
					//cout << num_points[graph][curve] << " points on curve " << curve << endl;
					point=0;
					++curve;					
				}
			}
		
//			cout << "graph = " << graph << endl;
//			cout << "curve = " << curve << endl;
//			cout << "Type" << C << ": Curve S = " << S << ", point " << point << ": Re = " << Re << ", K = " << K << endl;
			++point;			
		}while(fscanf(stream, "%l")!=EOF);

		num_points[graph][curve] = point;
		//cout << num_points[graph][curve] << " points on curve " << curve << endl;
		++curve;
		num_curves[graph] = curve;
		//cout << num_curves[graph] << " curves on graph " << graph << endl;
		++graph;
		num_graphs = graph;
		//cout << "num_graphs = " << num_graphs << endl;	
/*
		for(graph=0; graph<num_graphs; ++graph)
		{
			cout << "num_curves[" << graph << "] = " << num_curves[graph] << endl;
			for(curve=0; curve<num_curves[graph]; ++curve)
				cout << "num_points[" << graph << "][" << curve << "] = " << num_points[graph][curve] << endl;
		}
*/
		// Can now dimension the matrix
		Coeff = new double*** [num_graphs];
		for(graph=0; graph<num_graphs; ++graph) // For each graph
		{
			Coeff[graph] = new double** [num_curves[graph]];
			for(curve=0; curve<num_curves[graph]; ++curve) // For each curve
			{
				Coeff[graph][curve] = new double* [num_points[graph][curve]];
				for(point=0; point<num_points[graph][curve]; ++point)
					Coeff[graph][curve][point] = new double [4]; // 4 pieces of data per point (Type, S, Re, K)
			}
		}

		// Now go through and copy the data to the relevant place
		// Set pointer to beginning of file
		fseek(stream, 0L, SEEK_SET);

		graph=0;
		curve=0;
		point=0;
		fscanf(stream, "%f", &C );// The first value of C
		fscanf(stream, "%f", &S );// The first value of S
//		cout << "First C = " << C << endl;
//		cout << "First S = " << S << endl;
		fseek(stream, 0L, SEEK_SET);// Reset pointer to beginning of file
//		C=1;
//		S=0.3267;
		do
		{
			C_old = C;
			S_old = S;
			fscanf(stream, "%f", &C );
			fscanf(stream, "%f", &S );
			fscanf(stream, "%f", &Re );
			fscanf(stream, "%f", &K );
		
			if((fabs(C_old - C)>1e-4))
			{ 
				point=0;
				curve=0;
				++graph;
			}
			else
			{
				if((fabs(S_old - S)>1e-4))
				{ 
					point=0;
					++curve;
				}
			}
		
			//cout << "curve = " << curve << endl;
			//cout << "Type" << C << ": Curve S = " << S << ", point " << point << ": Re = " << Re << ", K = " << K << endl;
			Coeff[graph][curve][point][0] = C;
			Coeff[graph][curve][point][1] = S;
			Coeff[graph][curve][point][2] = Re;
			Coeff[graph][curve][point][3] = K;
			++point;
		}while(fscanf(stream, "%l")!=EOF);
			
		fclose(stream);
	}
/*
	int temp_precision = cout.precision(); cout << setprecision(6);
	for(graph=0; graph<num_graphs; ++graph) // For each graph
	{
		for(curve=0; curve<num_curves[graph]; ++curve) // For each curve
		{
			for(point=0; point<num_points[graph][curve]; ++point) // For each point
			{
				//cout << "Graph " << graph << ", curve " << curve << ", point " << point << ":\n";
				cout << Coeff[graph][curve][point][0] << "\t" 
						<< Coeff[graph][curve][point][1] << "\t" 
							<< Coeff[graph][curve][point][2] << "\t" 
								<< Coeff[graph][curve][point][3] << "\t" << endl; 
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << setprecision(temp_precision);
//*/
}
void CAPLDev::ReadInput(CProperties* pPpt, char *InputFile)
{
	int last_entry;
	CommonReadInput(InputFile, pPpt->NUM_PARAMS, this->labels, this->values, this->strings, last_entry);
	for(int r=0; r<last_entry+1; ++r)
	{
		// ====================================================================================================
		// Parameter file for Exhaust APLD[0]
		// ====================================================================================================

		// Optional object description (max. 500 characters - use underscores for spaces)
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "strDesc") == 0) strDesc = strings[r];

		if(strcmp(labels[r], "DEVICE") == 0) DEVICE = int(values[r]);
		if(strcmp(labels[r], "PARAM_S") == 0) PARAM_S = values[r];
		if(strcmp(labels[r], "X") == 0) X = values[r];
		if(strcmp(labels[r], "NAPLDev") == 0) NAPLDev = DoubleToBool(values[r]);

		// Configuration
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "nEntries") == 0) 
		{
			nEntries = int(values[r]); if(nEntries<1) nEntries = 1; // Require at least one
			nAdmissionMaps = nEntries==1?nEntries:nEntries+1;
			entry = new int [nEntries];
			LOSS_FILE = new char* [nAdmissionMaps]; // If single entry need only one full admission map, else need the full admission map plus one partial admission map per entry
		}
		if(strcmp(labels[r], "transmEnt") == 0) transmEnt = DoubleToBool(values[r]);
		if(strcmp(labels[r], "entry0") == 0 && nEntries>0) entry[0] = int(values[r]);
		if(strcmp(labels[r], "entry1") == 0 && nEntries>1) entry[1] = int(values[r]);
		if(strcmp(labels[r], "entry2") == 0 && nEntries>2) entry[2] = int(values[r]);
		if(strcmp(labels[r], "entry3") == 0 && nEntries>3) entry[3] = int(values[r]);
		if(strcmp(labels[r], "exitEndEnv") == 0) exitEndEnv = int(values[r]);
		if(strcmp(labels[r], "juncInlet") == 0) juncInlet = int(values[r]);
		if(strcmp(labels[r], "innerInlet") == 0) innerInlet = int(values[r]);
		if(strcmp(labels[r], "innerFollowing") == 0) innerFollowing = int(values[r]);
		if(strcmp(labels[r], "outerInlet") == 0) outerInlet = int(values[r]);
		if(strcmp(labels[r], "outerFollowing") == 0) outerFollowing = int(values[r]);
		if(strcmp(labels[r], "averagePR") == 0) averagePR = DoubleToBool(values[r]);
		if(strcmp(labels[r], "massAverageMFP") == 0) massAverageMFP = DoubleToBool(values[r]);
		
			// If not averaging to calculate PR, i.e., averagePR == 0 == false
			// ----------------------------------------------------------------------------------------------------
			if(strcmp(labels[r], "entryPR") == 0) entryPR = int(values[r]);

		// Loss coefficient
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "CONSTANT_K") == 0) CONSTANT_K = DoubleToBool(values[r]);
	
			// If applying constant loss coefficient K, i.e., CONSTANT_K == 1 == true
			// ----------------------------------------------------------------------------------------------------
			if(strcmp(labels[r], "K_VALUE") == 0) K_VALUE = values[r];
			if(strcmp(labels[r], "etaConst") == 0) etaConst = values[r];
			if(strcmp(labels[r], "D_VALUE") == 0) D_VALUE = values[r];
		
			// Else variable loss coefficient K, i.e., CONSTANT_K == 0 == false
			// ----------------------------------------------------------------------------------------------------	
			if(strcmp(labels[r], "INTERP_M2") == 0) INTERP_M2 = DoubleToBool(values[r]);
			if(strcmp(labels[r], "FULL_LOSS_FILE") == 0) LOSS_FILE[0] = strings[r];
			if(strcmp(labels[r], "PARTIAL_LOSS_FILE0") == 0 && nEntries>0) LOSS_FILE[1] = strings[r];
			if(strcmp(labels[r], "PARTIAL_LOSS_FILE1") == 0 && nEntries>1) LOSS_FILE[2] = strings[r];
			if(strcmp(labels[r], "PARTIAL_LOSS_FILE2") == 0 && nEntries>2) LOSS_FILE[3] = strings[r];
			if(strcmp(labels[r], "PARTIAL_LOSS_FILE3") == 0 && nEntries>3) LOSS_FILE[4] = strings[r];
			if(strcmp(labels[r], "USE_LAMBDA") == 0) USE_LAMBDA = DoubleToBool(values[r]);
			
				// If interpolating K between full and partial admission values, i.e., USE_LAMBDA == 1 == true
				// ----------------------------------------------------------------------------------------------------
				if(strcmp(labels[r], "ABS_LAMBDA") == 0) ABS_LAMBDA = DoubleToBool(values[r]);
				if(strcmp(labels[r], "LAMBDA_THRESHOLD") == 0) LAMBDA_THRESHOLD = DoubleToBool(values[r]);
				if(strcmp(labels[r], "LINEAR_INTERP") == 0) LINEAR_INTERP = DoubleToBool(values[r]);
				if(strcmp(labels[r], "weightPartial") == 0) weightPartial = values[r];
				if(strcmp(labels[r], "POLY_DEGREE") == 0) POLY_DEGREE = values[r];

							
			if(strcmp(labels[r], "PRINT_LOAD_LOSS") == 0) PRINT_LOAD_LOSS = DoubleToBool(values[r]);
			if(strcmp(labels[r], "PRINT_LOSS_FILES") == 0) PRINT_LOSS_FILES = DoubleToBool(values[r]);
			if(strcmp(labels[r], "CALIBRATE") == 0) CALIBRATE = DoubleToBool(values[r]);
		
				// If calibrating, i.e., CALIBRATE == 1 == true
				// ----------------------------------------------------------------------------------------------------
				if(strcmp(labels[r], "STEADY_FILE") == 0) STEADY_FILE = strings[r];
				if(strcmp(labels[r], "tolPR") == 0) tolPR = values[r];
				if(strcmp(labels[r], "tol_MFP") == 0) tol_MFP = values[r];
				if(strcmp(labels[r], "waitCount") == 0) waitCount = int(values[r]);
				if(strcmp(labels[r], "tInterval") == 0) tInterval = values[r];
				if(strcmp(labels[r], "CONST_BASELINE") == 0) CONST_BASELINE = DoubleToBool(values[r]);
				if(strcmp(labels[r], "K_baseline") == 0) K_baseline = values[r];
				if(strcmp(labels[r], "RESTART_LOSS_FILE") == 0) RESTART_LOSS_FILE = strings[r];
				if(strcmp(labels[r], "PRINT_STEADY_SCREEN") == 0) PRINT_STEADY_SCREEN = DoubleToBool(values[r]);
				if(strcmp(labels[r], "PRINT_STEADY_FILE") == 0) PRINT_STEADY_FILE = DoubleToBool(values[r]);
				if(strcmp(labels[r], "CALIBRATION_DELAY") == 0) CALIBRATION_DELAY = values[r];
				if(strcmp(labels[r], "PR_THRESHOLD") == 0) PR_THRESHOLD = values[r];

		// Measurements
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "USE_DEF_FREQ") == 0) USE_DEF_FREQ = DoubleToBool(values[r]);
		if(strcmp(labels[r], "freq") == 0) freq = int(values[r]);
		if(strcmp(labels[r], "freqPulse") == 0) freqPulse = values[r];
		if(strcmp(labels[r], "degCycle") == 0) degCycle = int(values[r]);
		if(strcmp(labels[r], "print_from_time") == 0) print_from_time = values[r];
    if(strcmp(labels[r], "PRINT_DEBUG_FILE") == 0) PRINT_DEBUG_FILE = DoubleToBool(values[r]);
		if(strcmp(labels[r], "PRINT_MOVIE_FILE") == 0) PRINT_MOVIE_FILE = DoubleToBool(values[r]);

		// Turbine geometry
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "tipDiameter") == 0) tipDiameter = values[r];

		// ====================================================================================================
		// End of file
		// ====================================================================================================
	}
	// Set some derived parameters
	if(USE_DEF_FREQ) freq = pPpt->freq;		// Use the default sampling rate

	if(!averagePR)
	{
		bool ENTRY_AVAILABLE = false;
		for(int e=0; e<nEntries; ++e) if(entry[e]==entryPR) ENTRY_AVAILABLE = true;
		if(!ENTRY_AVAILABLE)
		{
			pPpt->Out(Identify()); pPpt->Out("::ReadInput: ");
			pPpt->Out("ID of entry used for PR (entryPR="); pPpt->Out(entryPR); pPpt->Out(") is not listed as a valid entry. Exiting.\n\n");
			exit(1);
		}
	}
}

void CAPLDev::ListProperties(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ListProperties\n");}

	// ====================================================================================================
	// Parameter file for Exhaust APLD[0]
	// ====================================================================================================

	pPpt->Out(Underline(Identify(), "=", "\t", strDesc));
	pPpt->Out("\n");

	pPpt->Out("\tBase device type, DEVICE\t\t\t=\t"); pPpt->Out(Get_DEVICE()); pPpt->Out("\n");
	X = 1;	
	if(DEVICE==GAUZE)
	{
		pPpt->Out("\tSolidity, PARAM_S\t\t\t\t=\t"); pPpt->Out(PARAM_S); pPpt->Out("\n");
		pPpt->Out("\tX\t\t\t\t\t\t=\t"); pPpt->Out(X); pPpt->Out("\n");
	}
	else
	{
		if(DEVICE==THROTTLE)
		{
			pPpt->Out("\tSolidity, PARAM_S\t\t\t\t=\t"); pPpt->Out(PARAM_S); pPpt->Out("\n");
			pPpt->Out("\tX\t\t\t\t\t\t=\t"); pPpt->Out(X); pPpt->Out("\n");
		}
		else
		{
			if(DEVICE==EGR)
			{
				pPpt->Out("\tSolidity, PARAM_S\t\t\t\t=\t"); pPpt->Out(PARAM_S); pPpt->Out("\n");
				pPpt->Out("\tX\t\t\t\t\t\t=\t"); pPpt->Out(X); pPpt->Out("\n");
			}
			else
			{
				pPpt->Out("\tSolidity, PARAM_S\t\t\t\t=\tUNKNOWN"); pPpt->Out("\n");
				pPpt->Out("\tX\t\t\t\t\t\t=\tUNKNOWN"); pPpt->Out("\n");
			}
		}
	}
	pPpt->Out("\n");

	pPpt->Out("\tNumber of entries, nEntries\t\t\t=\t"); pPpt->Out(nEntries); pPpt->Out("\n");
	if(transmEnt)
	{	
		pPpt->Out("\tControlled by transm. boundaries, transmEnt\t=\t"); pPpt->Out(TrueOrFalse(transmEnt)); pPpt->Out("\n");
	}
	else
	{
		pPpt->Out("\tControlled by end env. boundaries, transmEnt\t=\t"); pPpt->Out(TrueOrFalse(transmEnt)); pPpt->Out("\n");
	}

	if(nEntries>1)
	{
		pPpt->Out("\tID of 1st entry ");
		if(transmEnt) pPpt->Out("(transmissive), "); else pPpt->Out("(end environment), ");
		pPpt->Out("entry[0]\t=\t");
		pPpt->Out(entry[0]); pPpt->Out("\n");

		pPpt->Out("\tID of 2nd entry ");
		if(transmEnt) pPpt->Out("(transmissive), "); else pPpt->Out("(end environment), ");
		pPpt->Out("entry[1]\t=\t");
		pPpt->Out(entry[1]); pPpt->Out("\n");
	}
	else
	{
		pPpt->Out("\tID of entry ");
		if(transmEnt) pPpt->Out("(transmissive), "); else pPpt->Out("(end env.), ");
		pPpt->Out("entry[0]\t\t=\t");
		pPpt->Out(entry[0]); pPpt->Out("\n");
	}
	pPpt->Out("\tID of exit, exitEndEnv\t\t\t\t=\t"); pPpt->Out(exitEndEnv); pPpt->Out("\n");
	pPpt->Out("\tID of rotor inlet junction, juncInlet\t\t=\t"); pPpt->Out(juncInlet); pPpt->Out("\n");
	pPpt->Out("\tID of pipe upstream of rotor inlet, inner limb\t=\t"); pPpt->Out(innerInlet); pPpt->Out("\n");
	pPpt->Out("\tID of pipe downstream of rotor inlet, inner limb\t=\t"); pPpt->Out(innerFollowing); pPpt->Out("\n");
	pPpt->Out("\tID of pipe upstream of rotor inlet, outer limb\t=\t"); pPpt->Out(outerInlet); pPpt->Out("\n");
	pPpt->Out("\tID of pipe downstream of rotor inlet, outer limb\t=\t"); pPpt->Out(outerFollowing); pPpt->Out("\n");
	if(averagePR)
	{
		pPpt->Out("\tAveraging PR across entries, averagePR\t\t=\t"); pPpt->Out(TrueOrFalse(averagePR)); pPpt->Out("\n");
	}
	else
	{
		pPpt->Out("\tNot averaging PR across entries, averagePR\t=\t"); pPpt->Out(TrueOrFalse(averagePR)); pPpt->Out("\n");
		pPpt->Out("\tID of entry used for PR, entryPR\t\t=\t"); pPpt->Out(entryPR); pPpt->Out("\n");
	}
	pPpt->Out("\n");

	// Loss coefficient
	// ----------------------------------------------------------------------------------------------------
	pPpt->Out(Underline("Loss coefficient", "-", "\t"));
	if(CONSTANT_K)
	{
		pPpt->Out("\tApplying constant loss coefficient, CONSTANT_K\t=\t"); pPpt->Out(TrueOrFalse(CONSTANT_K)); pPpt->Out("\n");

		// If applying constant loss coefficient K, i.e., CONSTANT_K == 1 == true
		// ----------------------------------------------------------------------------------------------------
		pPpt->Out("\tConstant loss coefficient K value, K_VALUE\t=\t"); pPpt->Out(K_VALUE); pPpt->Out("\n");
		pPpt->Out("\tAssociated constant efficiency value, etaConst\t=\t"); pPpt->Out(etaConst); pPpt->Out("\n");
		pPpt->Out("\tAssociated constant non-adiabatic term value, D\t=\t"); pPpt->Out(D_VALUE); pPpt->Out("\n");
		pPpt->Out("\n");
	}
	else
	{
		pPpt->Out("\tApplying variable loss coefficient, CONSTANT_K\t=\t"); pPpt->Out(TrueOrFalse(CONSTANT_K)); pPpt->Out("\n");

		// Else variable loss coefficient K, i.e., CONSTANT_K == 0 == false
		// ----------------------------------------------------------------------------------------------------
		int tempAdmissionMap = 0;
		//pPpt->Out(Underline("Loss curve", "-", "\t"));
		if(CALIBRATE)
		{
			pPpt->Out("\tCalibrating this run, CALIBRATE\t\t\t=\t"); pPpt->Out(TrueOrFalse(CALIBRATE)); pPpt->Out("\n");
			pPpt->Out("\tSteady PR-MFP characteristic, STEADY_FILE\t=\t"); pPpt->Out(STEADY_FILE); pPpt->Out("\n");
			pPpt->Out("\tMax. diff. between consecutive PR data, tolPR\t=\t"); pPpt->Out(tolPR); pPpt->Out("%\n");
			pPpt->Out("\tTolerance when matching desired MFP, tol_MFP\t=\t"); pPpt->Out(tol_MFP); pPpt->Out("%\n");
			pPpt->Out("\tMax. time steps to wait for steady, waitCount\t=\t"); pPpt->Out(waitCount); pPpt->Out("\n");
			pPpt->Out("\tInterval between recording points, tInterval\t=\t"); pPpt->Out(tInterval); pPpt->Out(" s\n");
			if(CONST_BASELINE) 
			{
				pPpt->Out("\tUsing constant initial K value, CONST_BASELINE\t=\t"); pPpt->Out(TrueOrFalse(CONST_BASELINE)); pPpt->Out("\n");
				pPpt->Out("\tInitial (constant) loss coeff. K, K_baseline\t=\t"); pPpt->Out(K_baseline); pPpt->Out("\n");
			}
			else
			{
				pPpt->Out("\tNot using constant initial K, CONST_BASELINE\t=\t"); pPpt->Out(TrueOrFalse(CONST_BASELINE)); pPpt->Out("\n");
				pPpt->Out("\tLoad initial K from file, RESTART_LOSS_FILE\t=\t"); pPpt->Out(RESTART_LOSS_FILE); pPpt->Out("\n");
				pPpt->Out("\tPrint loss coeffs loading, PRINT_LOAD_LOSS\t=\t"); pPpt->Out(TrueOrFalse(PRINT_LOAD_LOSS)); pPpt->Out("\n");
				if(PRINT_LOSS_FILES)
				{
					pPpt->Out("\tPrinting loss file, PRINT_LOSS_FILES\t\t=\t"); pPpt->Out(TrueOrFalse(PRINT_LOSS_FILES)); pPpt->Out("\n");
					pPpt->Out("\n");
					int desired_precision = 3;
					int temp_precision = cout.precision(); cout << setprecision(desired_precision);
					pPpt->Out("\t=====================================\n");
					pPpt->Out("\t************ RESTART DATA ***********\n");
					pPpt->Out("\t=====================================\n");
					pPpt->Out("\t"); pPpt->Out("Speed"); pPpt->Out("\t"); pPpt->Out("M2"); pPpt->Out("\t"); pPpt->Out("K"); pPpt->Out("\t"); pPpt->Out("U/Cis"); pPpt->Out("\t"); pPpt->Out("eta"); pPpt->Out("\t"); pPpt->Out("D"); pPpt->Out("\n");
					pPpt->Out("\t=====================================\n");
					tempAdmissionMap = 0; // Only one map when using a restart loss file for calibration
					for(int tempLossCurve=0; tempLossCurve<numLossSpeeds[tempAdmissionMap]; ++tempLossCurve)
					{
						for(int tempRow=0; tempRow<numLossPoints[tempAdmissionMap][tempLossCurve]; ++tempRow)
						{
							pPpt->Out("\t");
							pPpt->Out(K_loss_file[tempAdmissionMap][tempLossCurve][tempRow][lblSpeed]); pPpt->Out("\t");
							pPpt->Out(K_loss_file[tempAdmissionMap][tempLossCurve][tempRow][lblM2File]); pPpt->Out("\t");
							pPpt->Out(K_loss_file[tempAdmissionMap][tempLossCurve][tempRow][lblKFile]); pPpt->Out("\t");
							pPpt->Out(K_loss_file[tempAdmissionMap][tempLossCurve][tempRow][lblKUCsFile]); pPpt->Out("\t");
							pPpt->Out(K_loss_file[tempAdmissionMap][tempLossCurve][tempRow][lblKEtaFile]); pPpt->Out("\t");
							pPpt->Out(K_loss_file[tempAdmissionMap][tempLossCurve][tempRow][lblKDFile]); pPpt->Out("\n");
						}
						pPpt->Out("\t-------------------------------------\n");
					}
					pPpt->Out("\n");
					cout << setprecision(temp_precision);
				}
				else
				{
					pPpt->Out("\tNot printing loss file, PRINT_LOSS_FILES\t\t=\t"); pPpt->Out(TrueOrFalse(PRINT_LOSS_FILES)); pPpt->Out("\n");
				}
			}
			if(PRINT_STEADY_SCREEN)
			{
				pPpt->Out("\tProcessed map to screen, PRINT_STEADY_SCREEN\t=\t"); pPpt->Out(TrueOrFalse(PRINT_STEADY_SCREEN)); pPpt->Out("\n");
				pPpt->Out("\n");
				int desired_precision = 3;
				int temp_precision = cout.precision(); cout << setprecision(desired_precision);
				pPpt->Out("\t=====================================\n");
				pPpt->Out("\t"); pPpt->Out("Speed"); pPpt->Out("\t"); pPpt->Out("PR"); pPpt->Out("\t"); pPpt->Out("MFP"); pPpt->Out("\t"); pPpt->Out("U/Cis"); pPpt->Out("\t"); pPpt->Out("eta"); pPpt->Out("\n");
				pPpt->Out("\t=====================================\n");
				for(int tempSteadyCurve=0; tempSteadyCurve<numSteadySpeeds; ++tempSteadyCurve)
				{
					for(int tempRow=0; tempRow<numSteadyPoints[tempSteadyCurve]; ++tempRow)
					{
						pPpt->Out("\t");
						pPpt->Out(steadyPRMFP[tempSteadyCurve][tempRow][lblSpeed]); pPpt->Out("\t");
						pPpt->Out(steadyPRMFP[tempSteadyCurve][tempRow][lblPR]); pPpt->Out("\t");
						pPpt->Out(steadyPRMFP[tempSteadyCurve][tempRow][lblMFP]); pPpt->Out("\t");
						pPpt->Out(steadyPRMFP[tempSteadyCurve][tempRow][lblUCs]); pPpt->Out("\t");
						pPpt->Out(steadyPRMFP[tempSteadyCurve][tempRow][lblEta]); pPpt->Out("\n");
					}
					pPpt->Out("\t-------------------------------------\n");
				}
				cout << setprecision(temp_precision);
				pPpt->Out("\n");
			}
			if(PRINT_STEADY_FILE) pPpt->Out("\tProcessed map to file, PRINT_STEADY_FILE\t=\t"); pPpt->Out(TrueOrFalse(PRINT_STEADY_FILE)); pPpt->Out("\n");
		}
		else
		{
			pPpt->Out("\tInterpolating K based on:\n");
			if(INTERP_M2) pPpt->Out("\t- downstream Mach number");
			else pPpt->Out("\t- upstream Reynolds number");
			pPpt->Out(", INTERP_M2\t\t=\t"); pPpt->Out(TrueOrFalse(INTERP_M2)); pPpt->Out("\n");
			for(tempAdmissionMap=0; tempAdmissionMap<(nAdmissionMaps); ++tempAdmissionMap)
			{
				if(tempAdmissionMap==0) pPpt->Out("\tFull admission loss file, FULL_LOSS_FILE\t=\t");
				else
				{
					pPpt->Out("\tPartial admission loss file, PARTIAL_LOSS_FILE");
					pPpt->Out(tempAdmissionMap-1);
					pPpt->Out("\t=\t");
				}
				pPpt->Out(LOSS_FILE[tempAdmissionMap]); pPpt->Out("\n");
			}

			if(nAdmissionMaps>1) // If partial data is available
			{
				if(USE_LAMBDA) pPpt->Out("\tInterpolating based on admission");
				else		   pPpt->Out("\tNot using admission to interpolate");
				pPpt->Out(", USE_LAMBDA\t=\t"); pPpt->Out(TrueOrFalse(USE_LAMBDA)); pPpt->Out("\n");
				if(USE_LAMBDA) pPpt->Out("\tWeighting factor on partial K, weightPartial\t=\t"); pPpt->Out(weightPartial); pPpt->Out("\n");
			}
			if(PRINT_LOSS_FILES)
			{
				pPpt->Out("\tPrinting loss file, PRINT_LOSS_FILES\t\t=\t"); pPpt->Out(TrueOrFalse(PRINT_LOSS_FILES)); pPpt->Out("\n");
				pPpt->Out("\n");
				int desired_precision = 3;
				int temp_precision = cout.precision(); cout << setprecision(desired_precision);
				for(tempAdmissionMap=0; tempAdmissionMap<(nAdmissionMaps); ++tempAdmissionMap)
				{
					pPpt->Out("\t=====================================\n");
					if(tempAdmissionMap==0) pPpt->Out("\t*********** FULL ADMISSION **********\n");
					else
					{
						pPpt->Out("\tPARTIAL ADMISSION: ENTRY["); pPpt->Out(tempAdmissionMap-1); pPpt->Out("] FULL OPEN\n");
					}
					pPpt->Out("\t=====================================\n");
					pPpt->Out("\t"); pPpt->Out("Speed"); pPpt->Out("\t"); pPpt->Out("M2"); pPpt->Out("\t"); pPpt->Out("K"); pPpt->Out("\t"); pPpt->Out("U/Cis"); pPpt->Out("\t"); pPpt->Out("eta"); pPpt->Out("\t"); pPpt->Out("D"); pPpt->Out("\n");
					pPpt->Out("\t=====================================\n");
					for(int tempLossCurve=0; tempLossCurve<numLossSpeeds[tempAdmissionMap]; ++tempLossCurve)
					{
						for(int tempRow=0; tempRow<numLossPoints[tempAdmissionMap][tempLossCurve]; ++tempRow)
						{
							pPpt->Out("\t");
							pPpt->Out(K_loss_file[tempAdmissionMap][tempLossCurve][tempRow][lblSpeed]); pPpt->Out("\t");
							pPpt->Out(K_loss_file[tempAdmissionMap][tempLossCurve][tempRow][lblM2File]); pPpt->Out("\t");
							pPpt->Out(K_loss_file[tempAdmissionMap][tempLossCurve][tempRow][lblKFile]); pPpt->Out("\t");
							pPpt->Out(K_loss_file[tempAdmissionMap][tempLossCurve][tempRow][lblKUCsFile]); pPpt->Out("\t");
							pPpt->Out(K_loss_file[tempAdmissionMap][tempLossCurve][tempRow][lblKEtaFile]); pPpt->Out("\t");
							pPpt->Out(K_loss_file[tempAdmissionMap][tempLossCurve][tempRow][lblKDFile]); pPpt->Out("\n");
						}
						pPpt->Out("\t-------------------------------------\n");
					}
					pPpt->Out("\n");
				}
				cout << setprecision(temp_precision);
			}
			else
			{
				pPpt->Out("\tNot printing loss file, PRINT_LOSS_FILES\t\t=\t"); pPpt->Out(TrueOrFalse(PRINT_LOSS_FILES)); pPpt->Out("\n");
			}
			pPpt->Out("\tNot calibrating this run, CALIBRATE\t\t=\t"); pPpt->Out(TrueOrFalse(CALIBRATE)); pPpt->Out("\n");	
		}
		pPpt->Out("\n");
	}

	// Measurements
	// ----------------------------------------------------------------------------------------------------
	pPpt->Out(Underline("Measurements", "-", "\t"));
	if(USE_DEF_FREQ)
	{
		if(freq==1){pPpt->Out("\tUsing default sampling rate, pPpt->freq\t\t=\tonce per timestep\n");}
		else{pPpt->Out("\tUsing default sampling rate, pPpt->freq\t\t=\tonce per "); pPpt->Out(freq); pPpt->Out(" timesteps\n");}
	}
	else
	{
		if(freq==1){pPpt->Out("\tUsing local sampling rate, freq\t\t=\tonce per timestep\n");}
		else {pPpt->Out("\tUsing local sampling rate, freq\t\t\t=\tonce per "); pPpt->Out(freq); pPpt->Out(" timesteps\n");}
	}
	pPpt->Out("\tFrequency to separate cycle data, freqPulse\t=\t"); pPpt->Out(freqPulse); pPpt->Out(" Hz\n");
	pPpt->Out("\tNumber of degrees in cycle, degCycle\t\t=\t"); pPpt->Out(degCycle); pPpt->Out("\n");
	pPpt->Out("\tRecord onwards from, print_from_time\t\t=\t"); pPpt->Out(print_from_time); pPpt->Out(" s\n");
  if(PRINT_DEBUG_FILE) pPpt->Out("\tPrinting debug file\n"); else pPpt->Out("\tNot printing debug file\n");
  if(PRINT_MOVIE_FILE) pPpt->Out("\tPrinting movie file\n"); else pPpt->Out("\tNot printing movie file\n");
	pPpt->Out("\n");

	// Turbine geometry
	// ----------------------------------------------------------------------------------------------------
	pPpt->Out(Underline("Turbine geometry", "-", "\t"));
	pPpt->Out("\tTurbine rotor tip diameter, tipDiameter\t\t=\t"); pPpt->Out(tipDiameter); pPpt->Out(" mm\n");
	pPpt->Out("\n");

	pPpt->Out("\n");

	// ====================================================================================================
	// End of file
	// ====================================================================================================
}
char* CAPLDev::Get_DEVICE()
{
	char* temp;
	switch(DEVICE)
	{
	case UNKNOWN:
		temp = "UNKNOWN";
		break;
	case GAUZE:
		temp =  "Gauze";
		break;
	case THROTTLE:
		temp =  "Throttle";
		break;
	case EGR:
		temp =  "EGR Valve";
		break;
	default:
		temp =  "UNDEFINED";
		break;
	}
	return temp;
}
/*
void CAPLDev::LoadSteadyOld(char* InputFile)
//--------------------------------------------------//
// Loads steady PR-MFP characteristic	 			//
//													//
//--------------------------------------------------//
{
	int c = 0;
	double temp;
	int num_columns = 2;
	int col, row;
	int datapoints = 0;

	char pause;
	
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
//exit(1);
	}
	fclose(stream);
}
*/
void CAPLDev::LoadSteadyCharacteristics(CProperties* pPpt, char* InputFile)
//--------------------------------------------------//
// Loads steady PR-MFP characteristic into			//
// steadyPRMFP array, ignoring data points with		//
// zero MFP and those which are too close together	//
// for the calibration routine to handle.			//
//--------------------------------------------------//
{
	double tempSpeed, tempPR, tempPRPrev, tempMFP, tempUCs, tempEta, speedParam;
	int row, dataPoints, totalDataPoints;
	int num_columns = 5; // speed, PR, MFP, U/Cis, eta values
	
	FILE* stream;
	stream = fopen(InputFile, "r");

	if(stream == NULL)
	{
		pPpt->Out(Identify());
		pPpt->Out(":CAPLDev::LoadSteadyCharacteristics: Error while trying to open the steady characteristics file "); pPpt->Out(InputFile); pPpt->Out("\n");
		cout<<"Press any key to exit.\n";
		char cont;
		cin >> cont;
		exit(1);
	}
	else
	{
		numSteadySpeeds = 0;
		speedParam = 0;
		fseek(stream, 0L, SEEK_SET);			// Set pointer to beginning of file
		do
		{	
			fscanf(stream, "%lf", &tempSpeed);	// Runs over speed parameter
			//cout << tempSpeed << "\t";
			fscanf(stream, "%lf", &tempPR);		// Runs over PR parameter
			//cout << tempPR << "\t";
			fscanf(stream, "%lf", &tempMFP);	// Runs over MFP parameter
			//cout << tempMFP << "\t";
			fscanf(stream, "%lf", &tempUCs);	// Runs over U/Cis parameter
			//cout << tempUCs << "\t";
			fscanf(stream, "%lf", &tempEta);	// Runs over eta parameter
			//cout << tempEta << "\n";
			
			if(tempSpeed!=speedParam)			// Found new speed curve
			{
				speedParam = tempSpeed;			// Note new speed parameter
				++numSteadySpeeds;					// Increment counter				
			}
		}while(fscanf(stream, "%l")!=EOF);

		if(numSteadySpeeds==0)
		{
			pPpt->Out(Identify());
			pPpt->Out(":CAPLDev::LoadSteadyCharacteristics:\n"); 
			pPpt->Out("Need at least one speed curve in "); pPpt->Out(InputFile); pPpt->Out(" to operate the APL rotor boundary condition.\n\n");
			exit(1);
		}
		else
		{
			numSteadyPoints = new int [numSteadySpeeds];	// Can now dimension numSteadyPoints
		
			steadyCurve = 0;
			dataPoints = 0;
			totalDataPoints = 0;
			fseek(stream, 0L, SEEK_SET);			// Set pointer to beginning of file
			fscanf(stream, "%lf", &speedParam);		// Record first speed parameter
			fseek(stream, 0L, SEEK_SET);			// Set pointer to beginning of file
			do
			{	
				fscanf(stream, "%lf", &tempSpeed);	// Runs over speed parameter
				//cout << tempSpeed << "\t";
				tempPRPrev = tempPR;				// Save previous PR parameter
				fscanf(stream, "%lf", &tempPR);		// Runs over PR parameter
				//cout << tempPR << "\t";
				fscanf(stream, "%lf", &tempMFP);	// Runs over MFP parameter
				//cout << tempMFP << "\t";
				fscanf(stream, "%lf", &tempUCs);	// Runs over U/Cis parameter
				//cout << tempUCs << "\t";
				fscanf(stream, "%lf", &tempEta);	// Runs over eta parameter
				//cout << tempEta << "\n";
			
				if(tempSpeed!=speedParam)			// Found new speed curve
				{
					speedParam = tempSpeed;			// Note new speed parameter
					numSteadyPoints[steadyCurve] = dataPoints;	// Record number of data points on speed curve
					dataPoints = 0;					// Reset for next speed curve
					++steadyCurve;					// Increment counter				
				}

				if(fabs(tempMFP) > pPpt->ZERO_TOL)	// Only count non-zero MFP data points
				{
					if((fabs(tempPR - tempPRPrev)/tempPR)*100 > tolPR || dataPoints==0)	
					{
						++dataPoints;
						++totalDataPoints;
					}
					else							// Discard any data point too close to the previous one based on PR
					{
						pPpt->Out(Identify()); 
						pPpt->Out(":CAPLDev::LoadSteadyCharacteristics:\n");
						pPpt->Out("Not counting this data point on speed curve ["); pPpt->Out(steadyCurve); pPpt->Out("] as tempPR = "); pPpt->Out(tempPR); 
						pPpt->Out(", tempPRPrev = "); pPpt->Out(tempPRPrev); pPpt->Out(", % diff. = "); pPpt->Out((fabs(tempPR - tempPRPrev)/tempPR)*100); 
						pPpt->Out("\n\n");
						tempPR = tempPRPrev;		// Reset parameter value
					}
				}
				else
				{
					pPpt->Out(Identify()); 
					pPpt->Out(":CAPLDev::LoadSteadyCharacteristics:\n");
					pPpt->Out("Not counting this data point on speed curve ["); pPpt->Out(steadyCurve); pPpt->Out("] as tempMFP = "); pPpt->Out(tempMFP); 
					pPpt->Out("\n\n");
				}
			}while(fscanf(stream, "%l")!=EOF);
			numSteadyPoints[steadyCurve] = dataPoints;	// Record number of data points on final speed curve

			//cout << "Number of speed curves, numSteadySpeeds = " << numSteadySpeeds << endl;
			//for(steadyCurve=0; steadyCurve<numSteadySpeeds; ++steadyCurve) cout << "No. of non-zero MFP data points, numSteadyPoints[steadyCurve=" << steadyCurve << "] = " << numSteadyPoints[steadyCurve] << endl;
	
			// Dimension arrays
			steadyPRMFP = new double** [numSteadySpeeds];
			for(steadyCurve=0; steadyCurve<numSteadySpeeds; ++steadyCurve)
			{
				steadyPRMFP[steadyCurve] = new double* [numSteadyPoints[steadyCurve]];
				for(row=0; row<numSteadyPoints[steadyCurve]; ++row) steadyPRMFP[steadyCurve][row] = new double [num_columns];
			}
			// => steadyPRMFP[speed curve][data point or row][speed/PR/MFP value]

			fseek(stream, 0L, SEEK_SET);		// Reset pointer to beginning of file
			
			bool OVERWRITE, DISCARD;
			for(steadyCurve=0; steadyCurve<numSteadySpeeds; ++steadyCurve)
			{
				OVERWRITE = false;
				for(row=0; row<numSteadyPoints[steadyCurve]; ++row)
				{
					if(OVERWRITE)
					{
						--row;
						OVERWRITE = false;
					}
					
					do
					{
						DISCARD = false;					// Assume data is good
						fscanf(stream, "%lf", &tempSpeed);	// Runs speed parameter
						fscanf(stream, "%lf", &tempPR);		// Runs over PR parameter
						fscanf(stream, "%lf", &tempMFP);	// Runs over MFP parameter
						fscanf(stream, "%lf", &tempUCs);	// Runs over U/Cis parameter
						fscanf(stream, "%lf", &tempEta);	// Runs over eta parameter

						if(steadyCurve > 0)
						{
							// If data is not on new speed curve
							if(tempSpeed == steadyPRMFP[steadyCurve-1][numSteadyPoints[steadyCurve-1]-1][lblSpeed])
							{
								DISCARD = true;
								pPpt->Out(Identify()); 
								pPpt->Out(":CAPLDev::LoadSteadyCharacteristics:\n");
								pPpt->Out("Ignoring this data point as its speed ("); pPpt->Out(tempSpeed); 
								pPpt->Out(") is identical to that of the previous speed curve [");
								pPpt->Out(steadyCurve-1); pPpt->Out("] = "); pPpt->Out(steadyPRMFP[steadyCurve-1][numSteadyPoints[steadyCurve-1]-1][lblSpeed]); 
								pPpt->Out("\n\n");
							}
						}
					}while(DISCARD);
					
					steadyPRMFP[steadyCurve][row][lblSpeed] = tempSpeed;// Move speed parameter into array
					steadyPRMFP[steadyCurve][row][lblPR] = tempPR;		// Move PR parameter into array
					steadyPRMFP[steadyCurve][row][lblMFP] = tempMFP;	// Move MFP parameter into array
					steadyPRMFP[steadyCurve][row][lblUCs] = tempUCs;	// Move U/Cis parameter into array
					steadyPRMFP[steadyCurve][row][lblEta] = tempEta;	// Move eta parameter into array

					if(fabs(steadyPRMFP[steadyCurve][row][lblMFP]) < pPpt->ZERO_TOL)
					{
						OVERWRITE = true;	// If data point has zero MFP, prepare to overwrite it
					}
					else
					{
						if(row>0)
						{
							if(!((fabs(steadyPRMFP[steadyCurve][row][lblPR] - steadyPRMFP[steadyCurve][row-1][lblPR])
									/steadyPRMFP[steadyCurve][row-1][lblPR])*100 > tolPR))	
								OVERWRITE = true;	// If data point is too close to previous, prepare to overwrite it
						}
					}
				}
			}
///*
			//pPpt->Out(Identify()); pPpt->Out(":CAPLDev::LoadSteadyCharacteristics"); pPpt->Out("\n");
			pPpt->Out("CAPLDev::LoadSteadyCharacteristics\n");
			pPpt->Out("steadyPRMFP[steadyCurve][row][lbl]:\n");
			pPpt->Out("-----------------------------------\n");
			pPpt->Out("Speed\tPR\tMFP\tU/Cis\t\tEta\n");
			for(steadyCurve=0; steadyCurve<numSteadySpeeds; ++steadyCurve)
			{
				for(row=0; row<numSteadyPoints[steadyCurve]; ++row)
				{
					pPpt->Out(steadyPRMFP[steadyCurve][row][lblSpeed]); pPpt->Out("\t");
					pPpt->Out(steadyPRMFP[steadyCurve][row][lblPR]); pPpt->Out("\t");
					pPpt->Out(steadyPRMFP[steadyCurve][row][lblMFP]); pPpt->Out("\t");
					pPpt->Out(steadyPRMFP[steadyCurve][row][lblUCs]); pPpt->Out("\t");
					pPpt->Out(steadyPRMFP[steadyCurve][row][lblEta]); pPpt->Out("\n");
				}
				pPpt->Out("\n");
			}
//exit(1);
//*/
		}
	}
	fclose(stream);
}
void CAPLDev::LoadLossCoefficients(CTime* pMyTime, CProperties* pPpt, char* InputFile, int admissionMapNo)
//--------------------------------------------------//
//--------------------------------------------------//
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".LoadLossCoefficients\n");}
	double speedParam;
	double tempSpeed, tempM2, tempK, tempKPrev, tempKPrevPrev, tempError, tempErrorPrev, tempErrorPrevPrev, tempFactor, tempUCs, tempEta, tempD;
	double tempM2Prev;

	int row, dataPoints, totalDataPoints;
	int num_columns = 11/*10*/;
	double tolK = tolPR;
	
	FILE* stream;
	stream = fopen(InputFile, "r");

	if(stream == NULL)
	{
		pPpt->Out(Identify());
		pPpt->Out(":CAPLDev::LoadLossCoefficients: Error while trying to open the loss coefficients file "); pPpt->Out(InputFile); pPpt->Out("\n");
		cout<<"Press any key to exit.\n";
		char cont;
		cin >> cont;
		exit(1);
	}
	else
	{
		numLossSpeeds[admissionMapNo] = 0;
		speedParam = 0;
		fseek(stream, 0L, SEEK_SET);			// Set pointer to beginning of file
		do
		{	
			fscanf(stream, "%lf", &tempSpeed);	// Runs over speed parameter
			//cout << tempSpeed << "\t";
			fscanf(stream, "%lf", &tempM2);		// Runs over M2 parameter
			//cout << tempM2 << "\t";
			fscanf(stream, "%lf", &tempK);		// Runs over K parameter
			//cout << tempK << "\t";
			fscanf(stream, "%lf", &tempKPrev);	// Runs over KPrev parameter
			//cout << tempKPrev << "\t";
			fscanf(stream, "%lf", &tempKPrevPrev);// Runs over KPrevPrev parameter
			//cout << tempKPrevPrev << "\t";
			fscanf(stream, "%lf", &tempError);	// Runs over Error parameter
			//cout << tempError << "\t";
			fscanf(stream, "%lf", &tempErrorPrev);// Runs over ErrorPrev parameter
			//cout << tempErrorPrev << "\t";
			fscanf(stream, "%lf", &tempErrorPrevPrev);// Runs over ErrorPrevPrev parameter
			//cout << tempErrorPrevPrev << "\t";
			fscanf(stream, "%lf", &tempFactor);	// Runs over Factor parameter
			//cout << tempFactor << "\t";
			fscanf(stream, "%lf", &tempUCs);	// Runs over U/Cis parameter
			//cout << tempUCs << "\t";
			fscanf(stream, "%lf", &tempEta);	// Runs over Eta parameter
			//cout << tempEta << "\n";
			fscanf(stream, "%lf", &tempD);		// Runs over D parameter
			
			if(tempSpeed!=speedParam)			// Found new speed curve
			{
				speedParam = tempSpeed;			// Note new speed parameter
				++numLossSpeeds[admissionMapNo];// Increment counter				
			}
		}while(fscanf(stream, "%l")!=EOF);

		if(numLossSpeeds[admissionMapNo]==0)
		{
			pPpt->Out(Identify());
			pPpt->Out(":CAPLDev::LoadLossCoefficients:\n"); 
			pPpt->Out("Need at least one loss speed curve in "); pPpt->Out(InputFile); pPpt->Out(" to operate the A-/NA-PL rotor boundary condition.\n\n");
			cout<<"Press any key to exit.\n";
			char cont;
			cin >> cont;
			exit(1);
		}
		else
		{
			// If not calibrating, we need at least two speed curves on each map for interpolation purposes
			bool COPY_SINGLE_CURVE = false;
			if(!CALIBRATE && numLossSpeeds[admissionMapNo]==1)
			{
				pPpt->Out(Identify()); pPpt->Out(":LoadLossCoefficients: numLossSpeeds[admissionMapNo="); 
				pPpt->Out(admissionMapNo); pPpt->Out("] = "); pPpt->Out(numLossSpeeds[admissionMapNo]); 
				pPpt->Out(" so duplicating (and adding 0.001% speed) this curve to allow interpolation\n");
				pMyTime->Pause(pPpt, 0.5);
				numLossSpeeds[admissionMapNo]=2;
				COPY_SINGLE_CURVE = true;
			}
			numLossPoints[admissionMapNo] = new int [numLossSpeeds[admissionMapNo]]; // Can now dimension numLossPoints[admissionMapNo]
			lossCurve = 0;
			dataPoints = 0;
			totalDataPoints = 0;
			fseek(stream, 0L, SEEK_SET);			// Set pointer to beginning of file
			fscanf(stream, "%lf", &speedParam);		// Record first speed parameter
			fseek(stream, 0L, SEEK_SET);			// Set pointer to beginning of file
			do
			{	
				fscanf(stream, "%lf", &tempSpeed);	// Runs over speed parameter
				//cout << tempSpeed << "\t";
				tempM2Prev = tempM2;				// Save previous M2 parameter
				fscanf(stream, "%lf", &tempM2);		// Runs over M2 parameter
				//cout << tempM2 << "\t";
				fscanf(stream, "%lf", &tempK);		// Runs over K parameter
				//cout << tempK << "\t";
				fscanf(stream, "%lf", &tempKPrev);	// Runs over KPrev parameter
				//cout << tempKPrev << "\t";
				fscanf(stream, "%lf", &tempKPrevPrev);// Runs over KPrevPrev parameter
				//cout << tempKPrevPrev << "\t";
				fscanf(stream, "%lf", &tempError);	// Runs over Error parameter
				//cout << tempError << "\t";
				fscanf(stream, "%lf", &tempErrorPrev);// Runs over ErrorPrev parameter
				//cout << tempErrorPrev << "\t";
				fscanf(stream, "%lf", &tempErrorPrevPrev);// Runs over ErrorPrevPrev parameter
				//cout << tempErrorPrevPrev << "\t";
				fscanf(stream, "%lf", &tempFactor);	// Runs over Factor parameter
				//cout << tempFactor << "\t";
				fscanf(stream, "%lf", &tempUCs);	// Runs over U/Cis parameter
				//cout << tempUCs << "\t";
				fscanf(stream, "%lf", &tempEta);	// Runs over Eta parameter
				//cout << tempEta << "\n";
				fscanf(stream, "%lf", &tempD);		// Runs over D parameter
							
				if(tempSpeed!=speedParam)			// Found new speed curve
				{
					speedParam = tempSpeed;			// Note new speed parameter
					numLossPoints[admissionMapNo][lossCurve] = dataPoints;	// Record number of data points on speed curve
					dataPoints = 0;					// Reset for next speed curve
					++lossCurve;					// Increment counter				
				}

				if(fabs(tempK) > pPpt->ZERO_TOL)	// Only count non-zero K data points
				{
					//if((fabs(tempM2 - tempM2Prev)/tempM2)*100 > tolK || dataPoints==0)
					if(true)
					{
						++dataPoints;
						++totalDataPoints;
					}
					else							// Discard any data point too close to the previous one based on M2
					{
						pPpt->Out(Identify()); 
						pPpt->Out(":CAPLDev::LoadLossCoefficients:\n");
						pPpt->Out("Not counting this loss value on speed curve ["); pPpt->Out(lossCurve); pPpt->Out("] as tempM2 = "); pPpt->Out(tempM2); 
						pPpt->Out(", tempM2Prev = "); pPpt->Out(tempM2Prev); pPpt->Out(", % diff. = "); pPpt->Out((fabs(tempM2 - tempM2Prev)/tempM2)*100); 
						pPpt->Out("\n\n");
						tempM2 = tempM2Prev;		// Reset parameter value
					}
				}
				else
				{
					pPpt->Out(Identify()); 
					pPpt->Out(":CAPLDev::LoadLossCoefficients:\n");
					pPpt->Out("Not counting this loss value on speed curve ["); pPpt->Out(lossCurve); pPpt->Out("] as tempK = "); pPpt->Out(tempK); 
					pPpt->Out("\n\n");
				}
			}while(fscanf(stream, "%l")!=EOF);
			numLossPoints[admissionMapNo][lossCurve] = dataPoints;	// Record number of loss values on final speed curve
			if(COPY_SINGLE_CURVE) numLossPoints[admissionMapNo][1] = dataPoints;

			//pPpt->Out("Number of loss curves, numLossSpeeds = "); pPpt->Out(numLossSpeeds[admissionMapNo]); pPpt->Out("\n");
			//for(lossCurve=0; lossCurve<numLossSpeeds[admissionMapNo]; ++lossCurve) 
			//{
			//	pPpt->Out("No. of non-zero K loss values, numLossPoints[lossCurve="); pPpt->Out(lossCurve); pPpt->Out("] = "); pPpt->Out(numLossPoints[admissionMapNo][lossCurve]); pPpt->Out("\n");
			//}
			//exit(1);

			// Dimension arrays
			K_loss_file[admissionMapNo] = new double** [numLossSpeeds[admissionMapNo]];
			for(lossCurve=0; lossCurve<numLossSpeeds[admissionMapNo]; ++lossCurve)
			{
				K_loss_file[admissionMapNo][lossCurve] = new double* [numLossPoints[admissionMapNo][lossCurve]];
				for(row=0; row<numLossPoints[admissionMapNo][lossCurve]; ++row) K_loss_file[admissionMapNo][lossCurve][row] = new double [num_columns];
			}
			// => K_loss_file[admissionMapNo][speed curve][data point or row][speed/M2/K/U/Cis value]

			bool OVERWRITE, DISCARD;
			for(lossCurve=0; lossCurve<numLossSpeeds[admissionMapNo]; ++lossCurve)
			{
				fseek(stream, 0L, SEEK_SET);		// Reset pointer to beginning of file
				OVERWRITE = false;
				for(row=0; row<numLossPoints[admissionMapNo][lossCurve]; ++row)
				{
					if(OVERWRITE)
					{
						--row;
						OVERWRITE = false;
					}
					
					do
					{
						DISCARD = false;					// Assume data is good
						fscanf(stream, "%lf", &tempSpeed);	// Runs speed parameter
						//cout << tempSpeed << "\t";
						fscanf(stream, "%lf", &tempM2);		// Runs over M2 parameter
						//cout << tempM2 << "\t";
						fscanf(stream, "%lf", &tempK);		// Runs over K parameter
						//cout << tempK << "\t";
						fscanf(stream, "%lf", &tempKPrev);	// Runs over KPrev parameter
						//cout << tempKPrev << "\t";
						fscanf(stream, "%lf", &tempKPrevPrev);// Runs over KPrevPrev parameter
						//cout << tempKPrevPrev << "\t";
						fscanf(stream, "%lf", &tempError);	// Runs over Error parameter
						//cout << tempError << "\t";
						fscanf(stream, "%lf", &tempErrorPrev);// Runs over ErrorPrev parameter
						//cout << tempErrorPrev << "\t";
						fscanf(stream, "%lf", &tempErrorPrevPrev);// Runs over ErrorPrevPrev parameter
						//cout << tempErrorPrevPrev << "\t";
						fscanf(stream, "%lf", &tempFactor);	// Runs over Factor parameter
						//cout << tempFactor << "\t";
						fscanf(stream, "%lf", &tempUCs);	// Runs over U/Cis parameter
						//cout << tempUCs << "\t";
						fscanf(stream, "%lf", &tempEta);	// Runs over Eta parameter
						//cout << tempEta << "\n";
						fscanf(stream, "%lf", &tempD);		// Runs over D parameter

						if(lossCurve > 0)
						{
							// If data is not on new speed curve
							if(tempSpeed == K_loss_file[admissionMapNo][lossCurve-1][numLossPoints[admissionMapNo][lossCurve-1]-1][lblSpeed]
							&& !COPY_SINGLE_CURVE)
							{
								DISCARD = true;
								pPpt->Out(Identify()); 
								pPpt->Out(":CAPLDev::LoadLossCoefficients:\n");
								pPpt->Out("Ignoring this data point as its speed ("); pPpt->Out(tempSpeed); 
								pPpt->Out(") is identical to that of the previous loss curve [");
								pPpt->Out(lossCurve-1); pPpt->Out("] = "); pPpt->Out(K_loss_file[admissionMapNo][lossCurve-1][numLossPoints[admissionMapNo][lossCurve-1]-1][lblSpeed]); 
								pPpt->Out("\n\n");
							}
						}
					}while(DISCARD);
					
					K_loss_file[admissionMapNo][lossCurve][row][lblSpeedFile] = tempSpeed;					// Move speed parameter into array
					if(COPY_SINGLE_CURVE && lossCurve==1) K_loss_file[admissionMapNo][lossCurve][row][lblSpeedFile] *= 1.00001; // Add 0.001% to the speed to differentiate between duplicate curves when interpolating
					K_loss_file[admissionMapNo][lossCurve][row][lblM2File] = tempM2;						// Move M2 parameter into array
					K_loss_file[admissionMapNo][lossCurve][row][lblKFile] = tempK;							// Move K parameter into array
					K_loss_file[admissionMapNo][lossCurve][row][lblKPrevFile] = tempKPrev;					// Move KPrev parameter into array
					K_loss_file[admissionMapNo][lossCurve][row][lblKPrevPrevFile] = tempKPrevPrev;			// Move KPrevPrev parameter into array
					K_loss_file[admissionMapNo][lossCurve][row][lblErrorFile] = tempError;					// Move Error parameter into array
					K_loss_file[admissionMapNo][lossCurve][row][lblErrorPrevFile] = tempErrorPrev;			// Move ErrorPrev parameter into array
					K_loss_file[admissionMapNo][lossCurve][row][lblErrorPrevPrevFile] = tempErrorPrevPrev;	// Move ErrorPrevPrev parameter into array
					//K_loss_file[admissionMapNo][lossCurve][row][lblFactorFile] = tempFactor;				// Move Factor parameter into array
					//K_loss_file[admissionMapNo][lossCurve][row][lblFactorFile] = 1;						// When reloading this, use a value of 1 in case the file contains 0
					if(tempFactor==0) K_loss_file[admissionMapNo][lossCurve][row][lblFactorFile] = 1;		// When reloading this, use a value of 1 in case the file contains 0
					else K_loss_file[admissionMapNo][lossCurve][row][lblFactorFile] = tempFactor;			// Move Factor parameter into array
					K_loss_file[admissionMapNo][lossCurve][row][lblKUCsFile] = tempUCs;					// Move U/Cis parameter into array
					K_loss_file[admissionMapNo][lossCurve][row][lblKEtaFile] = tempEta;						// Move Eta parameter into array
					K_loss_file[admissionMapNo][lossCurve][row][lblKDFile] = tempD;							// Move D parameter into array

					if(fabs(K_loss_file[admissionMapNo][lossCurve][row][lblKFile]) < pPpt->ZERO_TOL)
					{
						OVERWRITE = true;	// If data point has zero K, prepare to overwrite it
					}
					//else
					//{
					//	if(row>0)
					//	{
					//		if(!((fabs(K_loss_file[admissionMapNo][lossCurve][row][lblX] - K_loss_file[admissionMapNo][lossCurve][row-1][lblX])
					//				/K_loss_file[admissionMapNo][lossCurve][row-1][lblX])*100 > tolK))	
					//			OVERWRITE = true;	// If data point is too close to previous, prepare to overwrite it
					//	}
					//}
				}
			}
			if(PRINT_LOAD_LOSS)
			{
				pPpt->Out("K_loss_file[admissionMapNo][lossCurve][row][lbl]:"); pPpt->Out("\n");
				pPpt->Out("---------------------------------"); pPpt->Out("\n");
				pPpt->Out("SpeedFile\tM2File\t\tKFile\tKPrevFile\tKPrevPrevFile\tErrorFile\tErrorPrevFile\tErrorPrevPrevFile\tFactorFile\tUCsFile\tEtaFile\tDFile\n");
				for(lossCurve=0; lossCurve<numLossSpeeds[admissionMapNo]; ++lossCurve)
				{
					for(row=0; row<numLossPoints[admissionMapNo][lossCurve]; ++row)
					{
						pPpt->Out(K_loss_file[admissionMapNo][lossCurve][row][lblSpeedFile]); pPpt->Out("\t\t");
						pPpt->Out(K_loss_file[admissionMapNo][lossCurve][row][lblM2File]); pPpt->Out("\t");
						pPpt->Out(K_loss_file[admissionMapNo][lossCurve][row][lblKFile]); pPpt->Out("\t");
						pPpt->Out(K_loss_file[admissionMapNo][lossCurve][row][lblKPrevFile]); pPpt->Out("\t\t");
						pPpt->Out(K_loss_file[admissionMapNo][lossCurve][row][lblKPrevPrevFile]); pPpt->Out("\t\t");
						pPpt->Out(K_loss_file[admissionMapNo][lossCurve][row][lblErrorFile]); pPpt->Out("\t");
						pPpt->Out(K_loss_file[admissionMapNo][lossCurve][row][lblErrorPrevFile]); pPpt->Out("\t\t");
						pPpt->Out(K_loss_file[admissionMapNo][lossCurve][row][lblErrorPrevPrevFile]); pPpt->Out("\t\t\t");
						pPpt->Out(K_loss_file[admissionMapNo][lossCurve][row][lblFactorFile]); pPpt->Out("\t\t");
						pPpt->Out(K_loss_file[admissionMapNo][lossCurve][row][lblKUCsFile]); pPpt->Out("\t");
						pPpt->Out(K_loss_file[admissionMapNo][lossCurve][row][lblKEtaFile]); pPpt->Out("\t");
						pPpt->Out(K_loss_file[admissionMapNo][lossCurve][row][lblKDFile]); pPpt->Out("\n");
					}
					if(lossCurve!=numLossSpeeds[admissionMapNo]-1) pPpt->Out("\n");
				}
			}
		}
	}
	fclose(stream);
}
void CAPLDev::QuadraticallyInterpolateSteadyCharacteristics(int speedCurve, double PR, double& rMFP, double& rUCs, double& rEta)
{
	int tempRow = 0;
	
	// Find the two points on the steady curve that PR lies between and interpolate
	while(steadyPRMFP[speedCurve][tempRow][lblPR] < PR && tempRow < numSteadyPoints[speedCurve] - 1) ++tempRow;
	// PR lies between tempRow and tempRow-1
/*
	// --------------------
	// Linear interpolation
	if(tempRow==0) tempRow=1; // If PR out of range, interpolate using first two points
	
	// Interpolate these points
	return
		steadyPRMFP[MFPlabel][tempRow-1]
		+ 
		((PR - steadyPRMFP[PRlabel][tempRow-1])/(steadyPRMFP[PRlabel][tempRow] - steadyPRMFP[PRlabel][tempRow-1])
		*(steadyPRMFP[MFPlabel][tempRow] - steadyPRMFP[MFPlabel][tempRow-1]));
	// --------------------
*/
	// --------------------
	// Quadratic interpolation
	// Select three most appropriate points and interpolate quadratically; y = Ax^2 + Bx + C
	int point1, point2, point3;
	if(tempRow==0)
	{
		if(numSteadyPoints[speedCurve] >= 3){point1 = tempRow; point2 = tempRow + 1; point3 = tempRow + 2;}
		else 
		{
			if(numSteadyPoints[speedCurve] >= 2){point1 = tempRow; point2 = tempRow + 1; point3 = point2;}
			else {point1 = tempRow; point2 = point1; point3 = point1;}
		}
	}
	else
	{
		if(tempRow == numSteadyPoints[speedCurve] - 1)
		{
			if(numSteadyPoints[speedCurve] >= 3){point3 = tempRow; point2 = tempRow - 1; point1 = tempRow - 2;}
			else 
			{
				if(numSteadyPoints[speedCurve] >= 2){point2 = tempRow; point1 = tempRow - 1; point3 = point1;}
				else {point1 = tempRow; point2 = point1; point3 = point1;}
			}
		}
		else {point2 = tempRow; point1 = tempRow - 1; point3 = tempRow + 1;}
	}

	double x1, x2, x3, y1, y2, y3, z1, z2, z3, zz1, zz2, zz3;
	double Ay, By, Cy, Az, Bz, Cz, Azz, Bzz, Czz;
	x1 = steadyPRMFP[speedCurve][point1][lblPR]; y1 = steadyPRMFP[speedCurve][point1][lblMFP]; z1 = steadyPRMFP[speedCurve][point1][lblUCs]; zz1 = steadyPRMFP[speedCurve][point1][lblEta];
	x2 = steadyPRMFP[speedCurve][point2][lblPR]; y2 = steadyPRMFP[speedCurve][point2][lblMFP]; z2 = steadyPRMFP[speedCurve][point2][lblUCs]; zz2 = steadyPRMFP[speedCurve][point2][lblEta];
	x3 = steadyPRMFP[speedCurve][point3][lblPR]; y3 = steadyPRMFP[speedCurve][point3][lblMFP]; z3 = steadyPRMFP[speedCurve][point3][lblUCs]; zz3 = steadyPRMFP[speedCurve][point3][lblEta];
	Ay = ((y3-y2)/((x3-x2)*(x3-x1))) - ((y1-y2)/((x1-x2)*(x3-x1)));
	By = ((y1-y2) + Ay*(pow(x2,2)-pow(x1,2)))/(x1-x2);//	B = ((y3-y2) + A*(pow(x2,2)-pow(x3,2)))/(x2-x3);
	Cy = y2 - Ay*pow(x2,2) - By*x2;//	C = y1 - A*pow(x1,2) - B*x1;//	C = y3 - A*pow(x3,2) - B*x3;
	Az = ((z3-z2)/((x3-x2)*(x3-x1))) - ((z1-z2)/((x1-x2)*(x3-x1)));
	Bz = ((z1-z2) + Az*(pow(x2,2)-pow(x1,2)))/(x1-x2);
	Cz = z2 - Az*pow(x2,2) - Bz*x2;
	Azz = ((zz3-zz2)/((x3-x2)*(x3-x1))) - ((zz1-zz2)/((x1-x2)*(x3-x1)));
	Bzz = ((zz1-zz2) + Azz*(pow(x2,2)-pow(x1,2)))/(x1-x2);
	Czz = zz2 - Azz*pow(x2,2) - Bzz*x2;
	rMFP = Ay*pow(PR,2) + By*PR + Cy;
	rUCs = Az*pow(PR,2) + Bz*PR + Cz;
	rEta = Azz*pow(PR,2) + Bzz*PR + Czz;

	cout << "x1 = " << x1 << "\tx2 = " << x2 << "\tx3 = " << x3 << "\tPR = " << PR << endl;
	cout << "y1 = " << y1 << "\ty2 = " << y2 << "\ty3 = " << y3 << "\trMFP = " << rMFP << endl;
	cout << "z1 = " << z1 << "\tz2 = " << z2 << "\tz3 = " << z3 << "\trUCs = " << rUCs << endl;
	cout << "zz1 = " << zz1 << "\tzz2 = " << zz2 << "\tzz3 = " << zz3 << "\trEta = " << rEta << endl;
	cout << endl;
	return; 
}
void CAPLDev::LinearlyInterpolateSteadyCharacteristics(int speedCurve, double PR, double& rMFP, double& rUCs, double& rEta)
{
	// Find the two points on the steady curve that PR lies between and interpolate
	int tempRow = 1;
	while(steadyPRMFP[speedCurve][tempRow][lblPR] < PR && tempRow < numSteadyPoints[speedCurve] - 1) ++tempRow;	// PR lies between tempRow and tempRow-1

	Interpolate(steadyPRMFP[speedCurve][tempRow-1][lblPR], steadyPRMFP[speedCurve][tempRow-1][lblMFP], 
				steadyPRMFP[speedCurve][tempRow][lblPR], steadyPRMFP[speedCurve][tempRow][lblMFP], PR, rMFP);

	Interpolate(steadyPRMFP[speedCurve][tempRow-1][lblPR], steadyPRMFP[speedCurve][tempRow-1][lblUCs], 
				steadyPRMFP[speedCurve][tempRow][lblPR], steadyPRMFP[speedCurve][tempRow][lblUCs], PR, rUCs);
/*
	cout << "steadyPRMFP[speedCurve][tempRow-1][lblPR] = " << steadyPRMFP[speedCurve][tempRow-1][lblPR] << endl;
	cout << "steadyPRMFP[speedCurve][tempRow-1][lblUCs] = " << steadyPRMFP[speedCurve][tempRow-1][lblUCs] << endl;
	cout << "steadyPRMFP[speedCurve][tempRow][lblPR] = " << steadyPRMFP[speedCurve][tempRow][lblPR] << endl;
	cout << "steadyPRMFP[speedCurve][tempRow][lblUCs] = " << steadyPRMFP[speedCurve][tempRow][lblUCs] << endl;
	cout << "PR = " << PR << endl;
	cout << "rUCs = " << rUCs << endl;
	cout << endl;
//*/
	Interpolate(steadyPRMFP[speedCurve][tempRow-1][lblPR], steadyPRMFP[speedCurve][tempRow-1][lblEta], 
				steadyPRMFP[speedCurve][tempRow][lblPR], steadyPRMFP[speedCurve][tempRow][lblEta], PR, rEta);

	return; 
}
void CAPLDev::InterpolateLossCoefficients(CProperties* pPpt, double time, double curveValue, double xValue, double &rYValue, double &rZValue, bool &rEXTRAPOLATE, int lblCurve, int lblXtemp, int lblYtemp, int lblZtemp, bool SHOW_WARNINGS, int admissionMapNo)
//InterpolateLossCoefficients(pPpt, curveValue==PARAM_S, xValue==PARAM_R, rYValue==K, rZValue==eta, lblSpeed, lblM2File, lblKFile, lblKEtaFile);// Interpolate value of K, eta
// Linearly interpolate three-dimensional array [curveValue][xValue][yValue][zValue] given curveValue and xValue to estimate yValue, zValue
{
	if(numLossSpeeds[admissionMapNo] < 2){pPpt->Out("Need at least two curves to linearly interpolate (numLossSpeeds[admissionMapNo="); pPpt->Out(admissionMapNo); pPpt->Out("]="); pPpt->Out(numLossSpeeds[admissionMapNo]); pPpt->Out(")\n");}
	int tempCurve, tempRow;
	int tempCurveBelow, tempCurveAbove, tempRowBelow, tempRowAbove;
	double curveValueBelow, curveValueAbove;
	double xCurveBelowPointBelow, yCurveBelowPointBelow, zCurveBelowPointBelow;
	double xCurveBelowPointAbove, yCurveBelowPointAbove, zCurveBelowPointAbove;
	double xCurveAbovePointBelow, yCurveAbovePointBelow, zCurveAbovePointBelow;
	double xCurveAbovePointAbove, yCurveAbovePointAbove, zCurveAbovePointAbove;
	double yBelow, yAbove, zBelow, zAbove;

	// Find the two curves that curveValue lies between
	tempCurve = 1; 
	tempRow = 0;
	while(K_loss_file[admissionMapNo][tempCurve][tempRow][lblCurve] < curveValue && tempCurve < numLossSpeeds[admissionMapNo] - 1) ++tempCurve;
	// curveValue lies between [tempCurve] and [tempCurve-1], or lies outside range
	tempCurveBelow = tempCurve - 1;
	tempCurveAbove = tempCurve;
	curveValueBelow = K_loss_file[admissionMapNo][tempCurveBelow][tempRow][lblCurve];
	curveValueAbove = K_loss_file[admissionMapNo][tempCurveAbove][tempRow][lblCurve];
	
	// Find pair of points on tempCurveBelow that xValue lies between
	tempRowBelow = 1;
	while(K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow][lblXtemp] < xValue && tempRowBelow < numLossPoints[admissionMapNo][tempCurveBelow] - 1) ++tempRowBelow;
	// xValue lies between [tempRowBelow] and [tempRowBelow-1], or lies outside range
	xCurveBelowPointBelow = K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow - 1][lblXtemp];
	yCurveBelowPointBelow = K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow - 1][lblYtemp];
	zCurveBelowPointBelow = K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow - 1][lblZtemp];
	xCurveBelowPointAbove = K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow][lblXtemp];
	yCurveBelowPointAbove = K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow][lblYtemp];
	zCurveBelowPointAbove = K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow][lblZtemp];

	// Find pair of points on tempCurveAbove that xValue lies between
	tempRowAbove = 1;
	while(K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove][lblXtemp] < xValue && tempRowAbove < numLossPoints[admissionMapNo][tempCurveAbove] - 1) ++tempRowAbove;
	// xValue lies between [tempRowAbove] and [tempRowAbove-1], or lies outside range
	xCurveAbovePointBelow = K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove - 1][lblXtemp];
	yCurveAbovePointBelow = K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove - 1][lblYtemp];
	zCurveAbovePointBelow = K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove - 1][lblZtemp];
	xCurveAbovePointAbove = K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove][lblXtemp];
	yCurveAbovePointAbove = K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove][lblYtemp];
	zCurveAbovePointAbove = K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove][lblZtemp];

	bool PERMIT_EXTRAPOLATION;
	PERMIT_EXTRAPOLATION = true;
	//PERMIT_EXTRAPOLATION = false;
	rEXTRAPOLATE = false;

	// ----------------------------------------------------------------------------------------------------

	{
		double xValueExtrapBelow;
		bool belowEXTRAPOLATE; belowEXTRAPOLATE = false;

		// The following prevents extrapolation and instead uses the closest value on the curve
		if(xCurveBelowPointBelow < xCurveBelowPointAbove)	// xValue==PARAM_R==M2
		{
			if(xValue > xCurveBelowPointAbove)
			{
				if(!PERMIT_EXTRAPOLATION) xValue = xCurveBelowPointAbove;
				xValueExtrapBelow = xCurveBelowPointAbove;
				rEXTRAPOLATE = true;
				belowEXTRAPOLATE = true;
			}
			else
			{
				if(xValue < xCurveBelowPointBelow)
				{
					if(!PERMIT_EXTRAPOLATION) xValue = xCurveBelowPointBelow;
					xValueExtrapBelow = xCurveBelowPointBelow;
					rEXTRAPOLATE = true;
					belowEXTRAPOLATE = true;
				}
			}
		}
		else //xCurveBelowPointBelow >= xCurveBelowPointAbove
		{
			if(xValue < xCurveBelowPointAbove)
			{
				if(!PERMIT_EXTRAPOLATION) xValue = xCurveBelowPointAbove;
				xValueExtrapBelow = xCurveBelowPointAbove;
				rEXTRAPOLATE = true;
				belowEXTRAPOLATE = true;
			}
			else
			{
				if(xValue > xCurveBelowPointBelow)
				{
					if(!PERMIT_EXTRAPOLATION) xValue = xCurveBelowPointBelow;
					xValueExtrapBelow = xCurveBelowPointBelow;
					rEXTRAPOLATE = true;
					belowEXTRAPOLATE = true;
				}
			}
		}
		
		// Linearly interpolate for yValue on CurveBelow; YValue==K
		yBelow	= yCurveBelowPointBelow + 
					(((yCurveBelowPointAbove - yCurveBelowPointBelow)/(xCurveBelowPointAbove - xCurveBelowPointBelow))	// Gradient
					*(xValue - xCurveBelowPointBelow));																	// Rise

		// Linearly interpolate for zValue on CurveBelow; ZValue==eta
		zBelow	= zCurveBelowPointBelow + 
					(((zCurveBelowPointAbove - zCurveBelowPointBelow)/(xCurveBelowPointAbove - xCurveBelowPointBelow))	// Gradient
					*(xValue - xCurveBelowPointBelow));																	// Rise
		
		if(belowEXTRAPOLATE) // Check extrapolated value is reasonable
		{
			if(yBelow<=0)
			{
				//Re-interpolate within known data
				yBelow	= yCurveBelowPointBelow + 
					(((yCurveBelowPointAbove - yCurveBelowPointBelow)/(xCurveBelowPointAbove - xCurveBelowPointBelow))	// Gradient
					*(xValueExtrapBelow - xCurveBelowPointBelow));														// Rise
			}
			if(zBelow<=0 || zBelow>1)
			{
				//Re-interpolate within known data
				zBelow	= zCurveBelowPointBelow + 
					(((zCurveBelowPointAbove - zCurveBelowPointBelow)/(xCurveBelowPointAbove - xCurveBelowPointBelow))	// Gradient
					*(xValueExtrapBelow - xCurveBelowPointBelow));														// Rise
			}
		}	
	}

	// ----------------------------------------------------------------------------------------------------

	{
		double xValueExtrapAbove;
		bool aboveEXTRAPOLATE; aboveEXTRAPOLATE = false;

		// The following prevents extrapolation and instead uses the closest value on the curve
		if(xCurveAbovePointBelow < xCurveAbovePointAbove)
		{
			if(xValue > xCurveAbovePointAbove)
			{
				if(!PERMIT_EXTRAPOLATION) xValue = xCurveAbovePointAbove;
				xValueExtrapAbove = xCurveAbovePointAbove;
				rEXTRAPOLATE = true;
				aboveEXTRAPOLATE = true;
			}
			else
			{
				if(xValue < xCurveAbovePointBelow)
				{
					if(!PERMIT_EXTRAPOLATION) xValue = xCurveAbovePointBelow;
					xValueExtrapAbove = xCurveAbovePointBelow;
					rEXTRAPOLATE = true;
					aboveEXTRAPOLATE = true;
				}
			}
		}
		else //xCurveAbovePointBelow >= xCurveAbovePointAbove
		{
			if(xValue < xCurveAbovePointAbove)
			{
				if(!PERMIT_EXTRAPOLATION) xValue = xCurveAbovePointAbove;
				xValueExtrapAbove = xCurveAbovePointAbove;
				rEXTRAPOLATE = true;
				aboveEXTRAPOLATE = true;
			}
			else
			{
				if(xValue > xCurveAbovePointBelow)
				{
					if(!PERMIT_EXTRAPOLATION) xValue = xCurveAbovePointBelow;
					xValueExtrapAbove = xCurveAbovePointBelow;
					rEXTRAPOLATE = true;
					aboveEXTRAPOLATE = true;
				}
			}
		}

		// Linearly interpolate for yValue on CurveAbove; YValue==K
		yAbove	= yCurveAbovePointBelow + 
					(((yCurveAbovePointAbove - yCurveAbovePointBelow)/(xCurveAbovePointAbove - xCurveAbovePointBelow))	// Gradient
					*(xValue - xCurveAbovePointBelow));																	// Rise

		// Linearly interpolate for zValue on CurveAbove; ZValue==eta
		zAbove	= zCurveAbovePointBelow + 
					(((zCurveAbovePointAbove - zCurveAbovePointBelow)/(xCurveAbovePointAbove - xCurveAbovePointBelow))	// Gradient
					*(xValue - xCurveAbovePointBelow));																	// Rise

		if(aboveEXTRAPOLATE) // Check extrapolated value is reasonable
		{
			if(yAbove<=0)
			{
				//Re-interpolate within known data
				yAbove	= yCurveAbovePointBelow + 
					(((yCurveAbovePointAbove - yCurveAbovePointBelow)/(xCurveAbovePointAbove - xCurveAbovePointBelow))	// Gradient
					*(xValueExtrapAbove - xCurveAbovePointBelow));														// Rise									
			}
			if(zAbove<=0 || zAbove>1)
			{
				//Re-interpolate within known data
				zAbove	= zCurveAbovePointBelow + 
					(((zCurveAbovePointAbove - zCurveAbovePointBelow)/(xCurveAbovePointAbove - xCurveAbovePointBelow))	// Gradient
					*(xValueExtrapAbove - xCurveAbovePointBelow));														// Rise
			}
		}	
	}

	// ----------------------------------------------------------------------------------------------------

	// Linearly interpolate for yValue between curves based on curveValue
	rYValue	= yBelow + (((yAbove - yBelow)/(curveValueAbove - curveValueBelow))										// Gradient
				*(curveValue - curveValueBelow));																	// Rise
	
	// Linearly interpolate for zValue between curves based on curveValue
	rZValue	= zBelow + (((zAbove - zBelow)/(curveValueAbove - curveValueBelow))										// Gradient
				*(curveValue - curveValueBelow));																	// Rise

/*
if(time>0.028)
{
cout << endl;
cout << "xCurveBelowPointBelow = " << xCurveBelowPointBelow << endl;
cout << "yCurveBelowPointBelow = " << yCurveBelowPointBelow << endl;
cout << "zCurveBelowPointBelow = " << zCurveBelowPointBelow << endl;
cout << "xCurveBelowPointAbove = " << xCurveBelowPointAbove << endl;
cout << "yCurveBelowPointAbove = " << yCurveBelowPointAbove << endl;
cout << "zCurveBelowPointAbove = " << zCurveBelowPointAbove << endl;
cout << endl;
cout << "xCurveAbovePointBelow = " << xCurveAbovePointBelow << endl;
cout << "yCurveAbovePointBelow = " << yCurveAbovePointBelow << endl;
cout << "zCurveAbovePointBelow = " << zCurveAbovePointBelow << endl;
cout << "xCurveAbovePointAbove = " << xCurveAbovePointAbove << endl;
cout << "yCurveAbovePointAbove = " << yCurveAbovePointAbove << endl;
cout << "zCurveAbovePointAbove = " << zCurveAbovePointAbove << endl;
cout << endl;
cout << "xValue = " << xValue << endl;
cout << "rYValue = " << rYValue << endl;
cout << "rZValue = " << rZValue << endl;
}
*/

/*
	if(true//SHOW_WARNINGS 
		//&& xValue>0.5
		)
	{
		if(xValue > K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow][lblXtemp] && xValue > K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow-1][lblXtemp])
		{
			pPpt->Out("CAPLDev::InterpolateLossCoefficients: Extrapolation warning:\n");
			pPpt->Out("\txValue ("); pPpt->Out(xValue); pPpt->Out(") > K_loss_file[admissionMapNo][tempCurveBelow="); pPpt->Out(tempCurveBelow); pPpt->Out("][tempRowBelow="); pPpt->Out(tempRowBelow); pPpt->Out("][lblXtemp]="); pPpt->Out(K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow][lblXtemp]); pPpt->Out("\n");
			pPpt->Out("\txValue ("); pPpt->Out(xValue); pPpt->Out(") > K_loss_file[admissionMapNo][tempCurveBelow="); pPpt->Out(tempCurveBelow); pPpt->Out("][tempRowBelow-1="); pPpt->Out(tempRowBelow-1); pPpt->Out("][lblXtemp]="); pPpt->Out(K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow-1][lblXtemp]); pPpt->Out("\n");
			pPpt->Out("\n");
		}

		if(xValue < K_loss_file[admissionMapNo][tempCurveBelow][tempRow][lblXtemp] && xValue < K_loss_file[admissionMapNo][tempCurveBelow][tempRow-1][lblXtemp])
		{
			pPpt->Out("CAPLDev::InterpolateLossCoefficients: Extrapolation warning:\n");
			pPpt->Out("\txValue ("); pPpt->Out(xValue); pPpt->Out(") < K_loss_file[admissionMapNo][tempCurveBelow="); pPpt->Out(tempCurveBelow); pPpt->Out("][tempRowBelow-1="); pPpt->Out(tempRowBelow-1); pPpt->Out("][lblXtemp]="); pPpt->Out(K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow-1][lblXtemp]); pPpt->Out("\n");
			pPpt->Out("\txValue ("); pPpt->Out(xValue); pPpt->Out(") < K_loss_file[admissionMapNo][tempCurveBelow="); pPpt->Out(tempCurveBelow); pPpt->Out("][tempRowBelow="); pPpt->Out(tempRowBelow); pPpt->Out("][lblXtemp]="); pPpt->Out(K_loss_file[admissionMapNo][tempCurveBelow][tempRowBelow][lblXtemp]); pPpt->Out("\n");
			pPpt->Out("\n");
		}

		if(xValue > K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove][lblXtemp] && xValue > K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove-1][lblXtemp])
		{
			pPpt->Out("CAPLDev::InterpolateLossCoefficients: Extrapolation warning:\n");
			pPpt->Out("\txValue ("); pPpt->Out(xValue); pPpt->Out(") > K_loss_file[admissionMapNo][tempCurveAbove="); pPpt->Out(tempCurveAbove); pPpt->Out("][tempRowAbove="); pPpt->Out(tempRowAbove); pPpt->Out("][lblXtemp]="); pPpt->Out(K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove][lblXtemp]); pPpt->Out("\n");
			pPpt->Out("\txValue ("); pPpt->Out(xValue); pPpt->Out(") > K_loss_file[admissionMapNo][tempCurveAbove="); pPpt->Out(tempCurveAbove); pPpt->Out("][tempRowAbove-1="); pPpt->Out(tempRowAbove-1); pPpt->Out("][lblXtemp]="); pPpt->Out(K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove-1][lblXtemp]); pPpt->Out("\n");
			pPpt->Out("\n");
		}

		if(xValue < K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove][lblXtemp] && xValue < K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove-1][lblXtemp])
		{
			pPpt->Out("CAPLDev::InterpolateLossCoefficients: Extrapolation warning:\n");
			pPpt->Out("\txValue ("); pPpt->Out(xValue); pPpt->Out(") < K_loss_file[admissionMapNo][tempCurveAbove="); pPpt->Out(tempCurveAbove); pPpt->Out("][tempRowAbove-1="); pPpt->Out(tempRowAbove-1); pPpt->Out("][lblXtemp]="); pPpt->Out(K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove-1][lblXtemp]); pPpt->Out("\n");
			pPpt->Out("\txValue ("); pPpt->Out(xValue); pPpt->Out(") < K_loss_file[admissionMapNo][tempCurveAbove="); pPpt->Out(tempCurveAbove); pPpt->Out("][tempRowAbove="); pPpt->Out(tempRowAbove); pPpt->Out("][lblXtemp]="); pPpt->Out(K_loss_file[admissionMapNo][tempCurveAbove][tempRowAbove][lblXtemp]); pPpt->Out("\n");
			pPpt->Out("\n");
		}

		pPpt->Out("For speed = "); pPpt->Out(curveValue); pPpt->Out("\tM2 = "); pPpt->Out(xValue); pPpt->Out("\tinterp. K = "); pPpt->Out(rYValue); pPpt->Out("\tinterp. eta = "); pPpt->Out(rZValue); pPpt->Out("\n");
		pPpt->Out("\n");
    //exit(1);
	}
//*/
	return;
}
void CAPLDev::RecordTurbineOperatingPoint(CProperties* pPpt)
// ====================================================================================================
// Calculate turbine flow and efficiency parameters
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".RecordTurbineOperatingPoint\n");}

	// Reset values
	m_stage[LIMB_INLET] = 0; m_stage[LIMB_EXIT] = 0; m_stage[ROTOR] = 0;
	m_abs_sum[LIMB_INLET] = 0; m_abs_sum[LIMB_EXIT] = 0; m_abs_sum[ROTOR] = 0;
	p_sum[LIMB_INLET] = 0; p_sum[LIMB_EXIT] = 0; p_sum[ROTOR] = 0;
	PR_sum = 0;
	UCs_sum = 0;
	p0_sum = 0;
	mT0_sum = 0;
	PR_stage_mean = 0;
	MFP_stage = 0;
	UCs_stage_mean = 0;
	Ws_stage_inlet = 0;
	UCs_stage_shaft = 0;
	Ws_stage_shaft = 0;
	Wact_stage_shaft = 0;
	dh0_turb = 0;
	

	p_exit = pExitEndEnv->pBN[ONE_SIDE]->p_dash*pPpt->PREF;	// Record stage exit static pressure (bar)

	
	CPipe* pPipe_temp;				// Declare temporary pipe pointer to be used for either transmissive or end envirionment entry boundaries
	CNode** pBN_temp;				// Declare temporary boundary node pointer to be used for either transmissive or end envirionment entry boundaries
	pBN_temp = new CNode* [nLocs];	// Two locations in each limb, plus one just upstream of the rotor
	
	int limbID, loc;
	for(limbID=0; limbID<nEntries; ++limbID)					// For each entry
	{
		// Set pointers to appropriate locations
		if(transmEnt)
		{
			pPipe_temp = pEntryTransm[limbID]->pPipe[ONE_SIDE];
			pBN_temp[LIMB_INLET/*0*/] = pEntryTransm[limbID]->pBN[ONE_SIDE];
			pBN_temp[LIMB_EXIT/*1*/] = &(pPipe_temp->Node[pPipe_temp->N-1]);	// Appropriate only for the single entry or simple twin entry models		
		}
		else
		{
			pPipe_temp = pEntryEndEnv[limbID]->pPipe[ONE_SIDE];
			pBN_temp[LIMB_INLET/*0*/] = pEntryEndEnv[limbID]->pBN[ONE_SIDE];
			pBN_temp[LIMB_EXIT/*1*/] = &(pPipe_temp->Node[pPipe_temp->N-1]);	// Appropriate only for the single entry or simple twin entry models
		}
		pBN_temp[ROTOR/*2*/] = pBN[0];	// Irrespective of limb. The nominally upstream boundary node for this boundary condition
		//pBN_temp[ROTOR/*2*/] = &(pPipe[0]->Node[0]);	// Irrespective of limb. The nominally upstream boundary node for this boundary condition
		//pBN_temp[ROTOR/*2*/] = pRotorInletJunction->pBN[pRotorInletJunction->Get_MAX_FLOW_BRANCH()];

		// Location measurements (LIMB_INLET==0, LIMB_EXIT==1, ROTOR==2)
		for(loc=0; loc<nLocs; ++loc)
		{
			p[limbID][loc]/*bar*/ = pBN_temp[loc]->p_dash*pPpt->PREF;/*bar*/
			T[limbID][loc] = pBN_temp[loc]->T;
			m[limbID][loc] = pBN_temp[loc]->mdot;
			u[limbID][loc] = pBN_temp[loc]->U*pPipe_temp->AREF;
			p0[limbID][loc]/*bar*/ = TotalPressureBar(pPpt, p[limbID][loc]/*bar*/, T[limbID][loc], u[limbID][loc]);
			T0[limbID][loc] = TotalTemperature(pPpt, T[limbID][loc], u[limbID][loc]);
		}
		
		// Stage performance parameters, i.e., [LIMB_INLET]
		PR[limbID] = p0[limbID][LIMB_INLET]/p_exit; // (bar/bar)
		MFP[limbID] = m[limbID][LIMB_INLET]*sqrt(T0[limbID][LIMB_INLET])/p0[limbID][LIMB_INLET]/*bar*/;
		tipVel[limbID] = PI*(PARAM_S*sqrt(T0[limbID][LIMB_INLET]))*(tipDiameter*1e-3); // v=omega*r=(2*pi*Nrpm/60)*(tipDia/2)==(pi*Nrps)*tipDia=pi*(SP*sqrt(T01))*tipDia
		
		// Prevent division by zero
		if(fabs(p0[limbID][LIMB_INLET]/*bar*/ - p_exit/*bar*/)>pPpt->ZERO_TOL)
		{
			UCs[limbID] = tipVel[limbID]/sqrt(2*pPpt->cpAir(T[limbID][LIMB_INLET])*T0[limbID][LIMB_INLET]*(1-pow(p_exit/p0[limbID][LIMB_INLET],(pPpt->gammaAir(T[limbID][LIMB_INLET])-1)/pPpt->gammaAir(T[limbID][LIMB_INLET]))));

			Ws[limbID] = m[limbID][LIMB_INLET]*pow(sqrt(2*pPpt->cpAir(T[limbID][LIMB_INLET])*T0[limbID][LIMB_INLET]*(1-pow(p_exit/p0[limbID][LIMB_INLET],(pPpt->gammaAir(T[limbID][LIMB_INLET])-1)/pPpt->gammaAir(T[limbID][LIMB_INLET])))),2)/2;
		}
		else{UCs[limbID] = 0; Ws[limbID] = 0;}

		m_stage[LIMB_INLET] += m[limbID][LIMB_INLET];				// Sum mass flow rate at limb inlet
		m_stage[LIMB_EXIT] += m[limbID][LIMB_EXIT];					// Sum mass flow rate at limb exit
		m_stage[ROTOR] += m[limbID][ROTOR];							// Sum mass flow rate at rotor
		m_abs_sum[LIMB_INLET] += fabs(m[limbID][LIMB_INLET]);
		m_abs_sum[LIMB_EXIT] += fabs(m[limbID][LIMB_EXIT]);
		m_abs_sum[ROTOR] += fabs(m[limbID][ROTOR]);
		p_sum[LIMB_INLET]/*bar*/ += p[limbID][LIMB_INLET]/*bar*/;	// Sum pressure
		p_sum[LIMB_EXIT]/*bar*/ += p[limbID][LIMB_EXIT]/*bar*/;		// Sum pressure
		p_sum[ROTOR]/*bar*/ += p[limbID][ROTOR]/*bar*/;				// Sum pressure
		p0_sum += p0[limbID][LIMB_INLET];							// Sum p0, only needed at loc = LIMB_INLET thus far
		PR_sum += PR[limbID];										// Sum PR
		MFP_stage += MFP[limbID];									// Sum MFP
		UCs_sum += UCs[limbID];										// Sum UCs
		Ws_stage_inlet += Ws[limbID];								// Sum Ws
	} // End of FOR, for each entry

//	// lambda based on mass flow rates (but this falls down when mass flow rate in one limb is negative)
//	if(m_stage[LIMB_INLET]==0) lambda[LIMB_INLET] = 0.5;		// Both zero, but therefore equal admission
//	else lambda[LIMB_INLET] = m[0/*==FIRST LIMB*/][LIMB_INLET]/m_stage[LIMB_INLET];	// Instantaneous mass admission ratio near rotor
//	if(m_stage[LIMB_EXIT]==0) lambda[LIMB_EXIT] = 0.5;			// Both zero, but therefore equal, admission
//	else lambda[LIMB_EXIT] = m[0/*==FIRST LIMB*/][LIMB_EXIT]/m_stage[LIMB_EXIT];	// Instantaneous mass admission ratio near rotor


	// lambda based on absolute mass flow rates
	if(m_abs_sum[LIMB_INLET]==0) lambda[LIMB_INLET] = 0.5;		// Total mass flow is zero, therefore equal admission
	else lambda[LIMB_INLET] = fabs(m[0][LIMB_INLET])/m_abs_sum[LIMB_INLET];	// Instantaneous mass admission ratio at entry

/*
	if(m_abs_sum[LIMB_EXIT]==0) lambda[LIMB_EXIT] = 0.5;		// Total mass flow is zero, therefore equal admission
	else lambda[LIMB_EXIT] = fabs(m[0][LIMB_EXIT])/m_abs_sum[LIMB_EXIT];	// Instantaneous mass admission ratio at rotor
*/

	// lambda based on static pressures (no problem with -ve pressures, but constant pressure junction means this is 0.5 ALWAYS!)
//	lambda[LIMB_INLET] = p[0/*==FIRST LIMB*/][LIMB_INLET]/p_sum[LIMB_INLET];// Instantaneous admission ratio near rotor
//	lambda[LIMB_EXIT] = p[0/*==FIRST LIMB*/][LIMB_EXIT]/p_sum[LIMB_EXIT];	// Instantaneous admission ratio near rotor

	if(averagePR) PR_stage_mean = PR_sum/nEntries;				// Calculate mean pressure ratio across all entries
	else PR_stage_mean = PR[entryPR];							// Otherwise choose which pressure ratio

	if(massAverageMFP && nEntries>1)
	{
		for(limbID=0; limbID<nEntries; ++limbID)	mT0_sum += (m[limbID][LIMB_INLET]/m_stage[LIMB_INLET])*T0[limbID][LIMB_INLET];	// Sum mass averaged T0, can only perform after m_stage[LIMB_INLET] completed 
		MFP_stage = m_stage[LIMB_INLET]* sqrt(mT0_sum)/(p0_sum/nEntries);		// Overwrite MFP_stage with mass averaged MFP for multiple entries turbine, i.e., Eqn. 4.20 in Romagnoli (2010th)
	}

	
	UCs_stage_mean = UCs_sum/nEntries;							// Calculate mean velocity ratio across all entries

	// Close to the shaft, i.e., across this APL rotor boundary condition
	//double p2_shaft;
	//p2_shaft = pBN[1/*DOWNSTREAM*/]->p_dash*pPpt->PREF; // bar
	//p2_shaft = pBN[1/*DOWNSTREAM*/]->p_dash*pPpt->PREF*1e5; // Pa
	
	
	
	//mRotorInletJunction = pRotorInletJunction->Get_m_dot(pRotorInletJunction->Get_MAX_FLOW_BRANCH());
	//mRotorInletJunction = pRotorInletJunction->pBN[pRotorInletJunction->Get_MAX_FLOW_BRANCH()]->mdot;
	limbID = 0;							// Limb values are identical for the rotor
	mRotorInletJunction = m[limbID][ROTOR];

/*	double m_dot_inner, T_inner, T0_inner, p_inner, p0_inner, u_inner;
	double m_dot_outer, T_outer, T0_outer, p_outer, p0_outer, u_outer;
	double m_dot_inner_following, m_dot_outer_following;
	double p2_shaft;
	
/*
cout << pInnerInlet->ID << endl;
cout << pInnerFollowing->ID << endl;
cout << pOuterInlet->ID << endl;
cout << pOuterFollowing->ID << endl;
exit(1);
//*/

	int whichNodeInner, whichNodeOuter;
	whichNodeInner = pInnerInlet->N-1/*EVEN END*/;
	whichNodeOuter = pInnerInlet->N-1/*EVEN END*/;
	//whichNodeInner = 0/*ODD END*/;
	//whichNodeOuter = 0/*ODD END*/;

	m_dot_inner = pInnerInlet->Node[whichNodeInner].mdot;
	u_inner = pInnerInlet->Node[whichNodeInner].U*pInnerInlet->AREF;
	T_inner = pInnerInlet->Node[whichNodeInner].T;
	T0_inner = TotalTemperature(pPpt, T_inner, u_inner);
	p_inner/*bar*/ = pInnerInlet->Node[whichNodeInner].p_dash*pPpt->PREF;
	p0_inner/*bar*/ = TotalPressureBar(pPpt, p_inner/*bar*/, T_inner, u_inner);
	m_dot_inner_following = pInnerFollowing->Node[0/*ODD END*/].mdot;

	m_dot_outer = pOuterInlet->Node[whichNodeOuter].mdot;
	u_outer = pOuterInlet->Node[whichNodeOuter].U*pOuterInlet->AREF;
	T_outer = pOuterInlet->Node[whichNodeOuter].T;
	T0_outer = TotalTemperature(pPpt, T_outer, u_outer);
	p_outer/*bar*/ = pOuterInlet->Node[whichNodeOuter].p_dash*pPpt->PREF;
	p0_outer/*bar*/ = TotalPressureBar(pPpt, p_outer/*bar*/, T_outer, u_outer);
	m_dot_outer_following = pOuterFollowing->Node[0/*ODD END*/].mdot;

	if(ABS_LAMBDA)	// Use adsolute mass flow for local admission ratio calc.
	{
		if( (fabs(m_dot_inner) + fabs(m_dot_outer))<pPpt->ZERO_TOL ) lambda[LIMB_EXIT] = 0.5; // Total mass flow is zero, therefore equal admission
		//if( (fabs(m_dot_inner) + fabs(m_dot_outer))<pPpt->ZERO_TOL ) lambda[LIMB_EXIT] = 0.5; // Total mass flow is zero, therefore equal admission
		else
		{
			lambda[LIMB_EXIT] = fabs(m_dot_inner)/(fabs(m_dot_inner) + fabs(m_dot_outer));	// Instantaneous mass admission ratio at rotor
			//lambda[LIMB_EXIT] = fabs(m_dot_inner - m_dot_inner_following)/(fabs(m_dot_inner - m_dot_inner_following) + fabs(m_dot_outer - m_dot_outer_following));	// Instantaneous mass admission ratio at rotor

			/*		
			cout << "m_dot_inner = " << m_dot_inner << endl;
			cout << "m_dot_inner_following = " << m_dot_inner_following << endl;
			cout << "m_dot_outer = " << m_dot_outer << endl;
			cout << "m_dot_outer_following = " << m_dot_outer_following << endl;
			cout << "lambda[LIMB_EXIT] = " << lambda[LIMB_EXIT] << endl;
			//*/
		}
	}
	else	//	Use actual mass flow for local admission ratio calc.
	{
		// =========================================================
		// Revised definition for instantaneous mass admission ratio
		// =========================================================
		lambda[LIMB_EXIT] = m_dot_inner/(m_dot_inner + m_dot_outer);	// Using actual m_dot values instead of the absolute
		if(LAMBDA_THRESHOLD)											// Apply admission ratio threshold to limit its value within 0--1, since negative mass flow may presence
		{
			if(lambda[LIMB_EXIT] < 0) lambda[LIMB_EXIT] = 0;
			if(lambda[LIMB_EXIT] > 1) lambda[LIMB_EXIT] = 1;
		}
		//cout<<lambda[LIMB_EXIT]<<endl;
	}


//lambda[LIMB_EXIT] = lambda[LIMB_INLET]; 
// Though this will be phase shifted, this is the only place where mass flow will differ 
// (because of the constant pressure junctions)

//p2_shaft = pBN[1]->p_dash*pPpt->PREF*1e5; // Pa
p2_shaft = pBN[1]->p_dash*pPpt->PREF; // bar

///*
	if(fabs(p0[limbID][ROTOR] - p2_shaft)>pPpt->ZERO_TOL)
	{	
//		Pa:
//		Ws_stage_shaft = 
//			m_dot_inner*pow(sqrt(2*pPpt->cpAir(T_inner)*T0_inner*(1-pow(p2_shaft/p0_inner,(pPpt->gammaAir(T_inner)-1)/pPpt->gammaAir(T_inner)))),2)/2
//		+
//			m_dot_outer*pow(sqrt(2*pPpt->cpAir(T_outer)*T0_outer*(1-pow(p2_shaft/p0_outer,(pPpt->gammaAir(T_outer)-1)/pPpt->gammaAir(T_outer)))),2)/2;

//		bar:
//		Ws_stage_shaft = m[limbID][ROTOR]
//			*(pPpt->cpAir(T[limbID][ROTOR])*T0[limbID][ROTOR]
//				*(1-pow(p2_shaft/p0[limbID][ROTOR],(pPpt->gammaAir(T[limbID][ROTOR])-1)/pPpt->gammaAir(T[limbID][ROTOR]))));



		UCs_stage_shaft = tipVel[limbID]/sqrt(2*pPpt->cpAir(T[limbID][ROTOR])*T0[limbID][ROTOR]*(1-pow(p2_shaft/p0[limbID][ROTOR],(pPpt->gammaAir(T[limbID][ROTOR])-1)/pPpt->gammaAir(T[limbID][ROTOR]))));

		//Ws_stage_shaft = mRotorInletJunction*pow(sqrt(2*pPpt->cpAir(T[limbID][ROTOR])*T0[limbID][ROTOR]*(1-pow(p2_shaft/p0[limbID][ROTOR],(pPpt->gammaAir(T[limbID][ROTOR])-1)/pPpt->gammaAir(T[limbID][ROTOR])))),2)/2;
		Ws_stage_shaft = mRotorInletJunction*pPpt->cpAir(T[limbID][ROTOR])*T0[limbID][ROTOR]*(1-pow(p2_shaft/p0[limbID][ROTOR],(pPpt->gammaAir(T[limbID][ROTOR])-1)/pPpt->gammaAir(T[limbID][ROTOR])));

/*
cout << "mRotorInletJunction = " << mRotorInletJunction << endl;
cout << "pPpt->cpAir(T[limbID][ROTOR]) = " << pPpt->cpAir(T[limbID][ROTOR]) << endl;
cout << "T0[limbID][ROTOR] = " << T0[limbID][ROTOR] << endl;
cout << "p2_shaft = " << p2_shaft << endl;
cout << "p0[limbID][ROTOR] = " << p0[limbID][ROTOR] << endl;
cout << "pPpt->gammaAir(T[limbID][ROTOR]) = " << pPpt->gammaAir(T[limbID][ROTOR]) << endl;
cout << "UCs_stage_shaft = " << UCs_stage_shaft << endl;
cout << "Ws_stage_shaft = " << Ws_stage_shaft << endl;
cout << endl;
*/
	}
	else Ws_stage_shaft = 0;

//cout << "Ws_stage_shaft = " << Ws_stage_shaft << endl << endl;
//*/

	//////////////////////////////////////////////
/*
	// Sum isentropic powers in all rotor inlet junction branches with positive flow (flow towards the junction)
	
	p2_shaft = pBN[1]->p_dash*pPpt->PREF; // bar
	Ws_stage_shaft = 0; // Reset
	int j;
	for(j=0; j<pRotorInletJunction->NPIPES; ++j)
	{
		tempNode = pRotorInletJunction->pBN[j];
		double tempMdot, tempP0;
		tempMdot = pRotorInletJunction->Get_m_dot(j);

		//if(tempMdot > 0) // If flow is toward rotor inlet junction
		{
//cout << "tempMdot[branch " << j << "] = " << tempMdot << endl;			
			tempP0 = TotalPressurePa(pPpt, tempNode->p_dash*pPpt->PREF*1e5, 
				tempNode->T, tempNode->U*pRotorInletJunction->pPipe[j]->AREF)/1e5; // bar
			if(fabs(tempP0 - p2_shaft) > pPpt->ZERO_TOL) // Prevents division by zero
			{

//				Ws_stage_shaft += 
//					tempMdot
//					*pow(sqrt(2*pPpt->cpAir(tempNode->T)
//					*TotalTemperature(pPpt, tempNode->T, tempNode->U*pRotorInletJunction->pPipe[j]->AREF)
//					*(1-pow(p2_shaft/tempP0,(pPpt->gammaAir(tempNode->T)-1)/pPpt->gammaAir(tempNode->T)))),2)/2;


				Ws_stage_shaft += 
					fabs(tempMdot)
					//tempMdot
					*pPpt->cpAir(tempNode->T)
					*TotalTemperature(pPpt, tempNode->T, tempNode->U*pRotorInletJunction->pPipe[j]->AREF)
					*(1-pow(p2_shaft/tempP0,(pPpt->gammaAir(tempNode->T)-1)/pPpt->gammaAir(tempNode->T)));

//cout << "Ws_stage_shaft[branch " << j << "] = " << Ws_stage_shaft << endl;
			}
		}
	}
//cout << "Ws_stage_shaft total = " << Ws_stage_shaft << endl << endl;
//*/
	//////////////////////////////////////////////

	Wact_stage_shaft = eta_OVERALL*Ws_stage_shaft;		// Calulate rotor power using values recorded near the rotor
	
	dh0_turb = eta_OVERALL*Ws_stage_inlet/mRotorInletJunction;	// Calulate specific rotor stagnation entahlpy change using (eta_OVERALL*inlet stage Ws) & mass flow rate recorded near the rotor
}

double CAPLDev::MeanlineCalc(CProperties* pPpt, double time, int timestep)
{
	mstarttime = clock();	cout<<"mStart = "<<mstarttime<<"\t";	// This is the meanline loop start time.
	DONE = false;
	// This loop will only be run for LEFT_TO_RIGHT flow, therefore need not to consider the reverse flow scenario.
	// These loss equations aren't valid for reverse flow.

	// Obtain rotor upstream flow state, i.e., p_us, T_us, u_us, rho_us, a_us & k_us.
	// ----------------------------------------------------------------------------------------------------
	double T_us = pBN[UPSTREAM]->T;										// Upstream static temperature (K)
	double k_us = pPpt->gammaAir(T_us);									// Upstream ratio of specific heat
	double p_us = pow(AA_1, (2*k_us)/(k_us-1)) * pPpt->PREF * 1e5;		// Upstream static pressure (Pa)
	double a_us = pBN[UPSTREAM]->A* (EX ? pPpt->AREFi : pPpt->AREFe);	// Upstream speed of sound (m/s)
	double u_us = pBN[UPSTREAM]->U* (EX ? pPpt->AREFi : pPpt->AREFe);	// Upstream flow velocity (m/s)
	double rho_us = p_us/(pPpt->R_air*T_us);							// Upstream static density (kg/m^3)
	double F_us = pBN[UPSTREAM]->f;										// Upstream pipe cross-sectional area (m^2)
	double mdot_us = rho_us*u_us*F_us;									// Upstream mass flow rate (kg/s)
	//double M_us = pBN[UPSTREAM]->M;
	
	// Routine check:
	//pPpt->Out("T_us = "); pPpt->Out(T_us); pPpt->Out("\t p_us = "); pPpt->Out(p_us); pPpt->Out("\n");
	//double a_us = sqrt(k_us*pPpt->R_air*T_us);
	//double u_us = (lambda_in[UPSTREAM]-a_us)/((k_us-1)/2);
	//pPpt->Out("a_us = "); pPpt->Out(a_us); pPpt->Out(", "); pPpt->Out(a_temp); pPpt->Out("\t k_us = "); pPpt->Out(k_us); pPpt->Out("\n");
	//pPpt->Out("lambda_in[UPSTREAM] = "); pPpt->Out(lambda_in[UPSTREAM]); pPpt->Out("\t lambda_in[DOWNSTREAM] = "); pPpt->Out(lambda_in[DOWNSTREAM]); pPpt->Out("\n");
	//pPpt->Out("u_us = "); pPpt->Out(u_us); pPpt->Out("\t rho_us = "); pPpt->Out(rho_us); pPpt->Out("\n\n");
	

	// Obtain rotor downstream parameter currently known, i.e., U_turb_ds & F_ds.
	// ----------------------------------------------------------------------------------------------------
	double U_turb_ds = ((55.9/1000)/2)*(2*PI*791.71);//((exDiameter)/1000/2)*(2*PI*N_turbine);		// Rotor exducer tangential velocity (m/s)
	double F_ds = pBN[DOWNSTREAM]->f;								// Downstream pipe cross-sectional area (m^2)

	// Upstream velocity triangle components
	// ----------------------------------------------------------------------------------------------------
	double Cm_us = u_us;											// Meridional component of upstream absolute velocity (m/s)
	double Cteta_us = ((-0.0146)*(pow(Cm_us,2) )) + 3.8067*Cm_us;	// Tangential component of upstream absolute velocity (m/s)
	double C_us = sqrt(pow(Cm_us,2) + pow(Cteta_us,2));				// Upstream absolute velocity (m/s)
	double alpha_us = atan(Cteta_us/Cm_us);							// Upstream absolute flow angle (rad)
	
	double U_turb_us = ((83.58/1000)/2)*(2*PI*791.71);//((tipDiameter/1000)/2)*(2*PI*N_turbine);		// Rotor inducer tangential velocity (m/s)
	double Wteta_us = U_turb_us - fabs(Cteta_us);					// Tangential component of upstream relative velocity (m/s)
	double W_us = sqrt(pow(Wteta_us,2) + pow(Cm_us,2));				// Upstream relative velocity (m/s) 
	double beta_us;
	if(Wteta_us < pPpt->ZERO_TOL) {beta_us = fabs(acos(Cm_us/ W_us));}	// Upstream relative flow angle (rad)
	else {beta_us = -1*fabs(acos(Cm_us/ W_us));}
	
	// Upstream relative flow parameter
	// ----------------------------------------------------------------------------------------------------
	double Mr_us = fabs(W_us/a_us);									// Upstream relative Mach no.
	double T0r_us = T_us*(1 + ((k_us-1)/2)*pow(Mr_us,2));						// Upstream relative stagnation temperature (K)
	double p0r_us = p_us*pow((1 + ((k_us-1)/2)*pow(Mr_us,2)), k_us/(k_us-1));	// Upstream relative stagnation temperature (K)
	double Wcr_us = sqrt((2*k_us*pPpt->R_air*T0r_us/*T0_us*/)/ (k_us+1));		// Upstream critical relative velocity (m/s) -- Eq. B70 in Meitner and Glassman (1983)

	double K_inc = 1.4;
	double K_p = 0.3;
	double in_betab = 0.3491;
	double ex_betab = -0.6981;

	// Components of Incidence loss from upstream velocity triangle
	// ----------------------------------------------------------------------------------------------------
	//double i_opt = atan((-1.98*tan(alpha_us))/ (Z*(1 - 1.98*Z)));
	double i_opt = atan((-1.98*tan(alpha_us))/ (12*(1 - 1.98*12)));
	double beta_diff = beta_us-in_betab-i_opt;
	// Original incidence loss equation from Mizumachi et al. (1979):
	//	if(beta_diff < (PI/4)) {L_inc' = 0.5*K_inc*pow(W_ds*sin(beta_diff),2);}
	//	else {L_inc' = 0.5*K_inc*pow(W_ds,2)*(0.5+fabs(beta_diff)-(PI/4) );}
	// The W_ds^2 term is moved outside incidence loss eq to ease the calculation, hence:
	double L_inc;
	if(beta_diff <= (PI/4)) {L_inc = 0.5*K_inc*pow(sin(beta_diff),2);}
	else {L_inc = 0.5*K_inc*(0.5+fabs(beta_diff)-(PI/4));}
	
	// Components of Passage loss across rotor (as a function of upstream & downstream relative velocities)
	// ----------------------------------------------------------------------------------------------------
	// Original passage loss equation from Futral and Wasserbauer (1965):
	//	L_p' = (K_p/2)*(pow(W_us,2)*( (cos(fabs(beta_us-i_opt)))^2) + pow(W_ds,2) )
	// This equation is broken up to ease the calculation:
	double L_pa = (K_p/2);
	double L_pb = L_pa*pow(cos(fabs(beta_us-i_opt)),2);
	
	// The changes in relative stagnation temperature, critical velocity and ideal relative
	//  stagnation pressure due to radius change
	// ----------------------------------------------------------------------------------------------------
	
	double T0r_ds = T0r_us+(((k_us-1)/ (2*k_us*pPpt->R_air))*(pow(U_turb_ds,2) - pow(U_turb_us,2)) );	// Downstream relative stagnation temperature (K) -- Eq. B83
	double Wcr_ds = Wcr_us* sqrt(T0r_ds/T0r_us);					// Downstream critical relative velocity (m/s) -- Eq. B84
	double p0r_ds_id = p0r_us*((T0r_ds/T0r_us),(k_us/(k_us-1)) );	// Ideal downstream relative stagnation relative pressure (Pa) -- Eq. B85
	
	
	/*if(time>0.1)
	{
		// Routine check:
		pPpt->Out("p_us = "); pPpt->Out(p_us); pPpt->Out("\t T_us = "); pPpt->Out(T_us); pPpt->Out("\n");
		pPpt->Out("W_us = "); pPpt->Out(W_us); pPpt->Out("\t U_turb_us = "); pPpt->Out(U_turb_us); pPpt->Out("\n");
		pPpt->Out("p0r_us = "); pPpt->Out(p0r_us); pPpt->Out("\t T0r_us = "); pPpt->Out(T0r_us); pPpt->Out("\t Wcr_us = "); pPpt->Out(Wcr_us); pPpt->Out("\n");
		
		pPpt->Out("p0r_ds_id = "); pPpt->Out(p0r_ds_id); pPpt->Out("\t T0r_ds = "); pPpt->Out(T0r_ds); pPpt->Out("\t Wcr_ds = "); pPpt->Out(Wcr_ds); pPpt->Out("\n");
	}*/

	
	// Iterate for p_ds such that mass flow rate across rotor wheel is conserved:
	// ----------------------------------------------------------------------------------------------------
	double p_ds;					// Downstream static pressure (Pa)
	p_ds = p_us;					// Initial assumption. Note p_ds must be lower than p_us or else the device becomes compressor!
	double mdot_localError = 1e6;	// Assume initial mass flow rate error
	double del_p_ds = 100;			// Assume initial del_p_ds
	int pcounter = 0;

	do
	{
		++pcounter;
		double W_ds_id = sqrt(pow(Wcr_ds,2)*((k_us+1)/ (k_us-1))*(1 - pow((p_ds/p0r_ds_id),((k_us-1)/k_us)) ));	// Ideal downstream relative velocity (m/s) -- Eq. B90
		// (L_inc' + L_p') = (pow(W_ds_id,2) - pow(W_ds,2)), hence:
		double W_ds = sqrt((pow(W_ds_id,2) - (L_pb + L_inc)*pow(W_us,2))/ (1+L_pa));	// Downstream relative velocity (m/s)
		double Ca_ds = W_ds*cos(ex_betab);												// Axial component of downstream absolute velocity (m/s)
		double Wteta_ds = W_ds*fabs(sin(ex_betab));										// Tangential component of downstream relative velocity (m/s)
		double Cteta_ds, alpha_ds/*, C_ds*/; 
		if(fabs(U_turb_ds)-fabs(Wteta_ds) <= pPpt->ZERO_TOL)
		{
			Cteta_ds = 0;
			alpha_ds = 0;	// Zero exit swirl
			//C_ds = Ca_ds; 
		}	
		else
		{
			if(fabs(U_turb_ds)-fabs(Wteta_ds) > pPpt->ZERO_TOL)
			{
				Cteta_ds = U_turb_ds - Wteta_ds;	// Positive Cteta_ds
				alpha_ds = atan(Cteta_ds/Ca_ds);	// Positive alpha_ds
				//C_ds = sqrt(pow(Cteta_ds,2) + pow(Ca_ds,2));
			}
			else
			{
				Cteta_ds = U_turb_ds - Wteta_ds;	// Negative Cteta_ds
				alpha_ds = atan(Cteta_ds/Ca_ds);	// Negative alpha_ds
				//C_ds = sqrt(pow(Cteta_ds,2) + pow(Ca_ds,2));
			}
		}
		double T_ds = T0r_ds - ( ((k_us-1)/(2*k_us*pPpt->R_air))*pow(W_ds,2) );		// Downstream static temperature (K)
		double k_ds = pPpt->gammaAir(T_ds);											// Downstream ratio of specific heat
		double Mr_ds = W_ds/ sqrt(k_ds*pPpt->R_air*T_ds);							// Downstream relative Mach no.
		double rho_ds = p_ds/(pPpt->R_air*T_ds);									// Downstream flow density (kg/m^3)
		double mdot_ds_d = rho_ds*W_ds*F_ds;										// Derived downstream mass flow rate (kg/s)
		
		mdot_localError = mdot_ds_d - mdot_us;
		if(mdot_localError > pPpt->ZERO_TOL)
		{
			/*pPpt->Out("Warning! mdot_ds_d > mdot_us. mdot_localError = "); pPpt->Out(mdot_localError);
			pPpt->Out("Continue? (y/n) \n");
			char cont;
			cin >> cont;
			pPpt->Out("\n");
			if(cont=='n' || cont=='N')
			{
				pPpt->Out("User answers NO. Exiting... \n");
				exit(1);
			}*/
			p_ds = p_ds - del_p_ds;
		}
		else 
		{
			/*pPpt->Out("Warning! mdot_ds_d > mdot_us. mdot_localError = "); pPpt->Out(mdot_localError);
			pPpt->Out("Continue? (y/n) \n");
			char cont;
			cin >> cont;
			pPpt->Out("\n");
			if(cont=='n' || cont=='N')
			{
				pPpt->Out("User answers NO. Exiting... \n");
				exit(1);
			}*/
			p_ds = p_ds + del_p_ds;
		}

		if(pcounter>1000) {DONE = true;}
		/*if(pcounter>500)
		{
			pPpt->Out("\n pcounter is > 500, = "); pPpt->Out(pcounter); pPpt->Out("\n");
			pPpt->Out("p_ds is = "); pPpt->Out(p_ds); pPpt->Out(" mdot_localError is = "); pPpt->Out(mdot_localError); pPpt->Out("\n"); 
			pPpt->Out("Continue? (y/n) \n");
			char cont;
			cin >> cont;
			pPpt->Out("\n");
			if(cont=='n' || cont=='N')
			{
				pPpt->Out("User answers NO. Exiting... \n");
				exit(1);
			}
		}*/
	}while(mdot_localError >= pPpt->ZERO_TOL && !DONE);
	
	// Routine check:
	//pPpt->Out("p_ds = "); pPpt->Out(p_ds); pPpt->Out("\t T_ds = "); pPpt->Out(T_ds); pPpt->Out("\n");
	//pPpt->Out("W_ds = "); pPpt->Out(W_ds); pPpt->Out("\t U_turb_ds = "); pPpt->Out(U_turb_ds); pPpt->Out("\n");
		
	double dh0_main = (U_turb_us*Cteta_us) - (U_turb_ds*Cteta_ds);		// Main stream stagnation enthalpy change
	
	// Disk friction loss
	// ----------------------------------------------------------------------------------------------------
	//double L_df = (0.02125* pow(U_turb_us,2)* pow(rho_us,2)* pow( ((tipDiameter/1000)/2),2) )/ 
	//				(mdot_us* pow( ((rho_us*U_turb_us*( (tipDiameter/1000)/2))/ pPpt->ViscosityAir(T_us)),0.2) );		// Eq. B129
	double L_df = (0.02125* pow(U_turb_us,2)* pow(rho_us,2)* pow( ((83.58/1000)/2),2) )/ 
					(mdot_us* pow( ((rho_us*U_turb_us*( (83.58/1000)/2))/ pPpt->ViscosityAir(T_us)),0.2) );

	// Clearance flow loss
	// ----------------------------------------------------------------------------------------------------
	//double L_cl = (2* dh0_main*ex_clrHeight/ ex_tipDiameter)/ (1 - (ex_hubDiameter/ex_tipDiameter));		// Eq. B130
	double L_cl = (2* dh0_main*0.43/ 78.65)/ (1 - (31.1/78.65));

	double dh0_turb = dh0_main-L_df-L_cl;		// Stagnation enthalpy change across rotor wheel
	double Work_turb = mdot_us*dh0_turb;		// Actual turbine output power (Watt)
	double a_ds = sqrt(pPpt->gammaAir(T_ds)*pPpt->R_air*T_ds);		// Downstream flow sonic velocity (m/s) 




	//if(time>0.1)
	//{
		// Routine check:
		//pPpt->Out("Cm_us = "); pPpt->Out(Cm_us); pPpt->Out("\t Cteta_us = "); pPpt->Out(Cteta_us); pPpt->Out("\n");
		//pPpt->Out("alpha_us = "); pPpt->Out(alpha_us); pPpt->Out("\t beta_us = "); pPpt->Out(beta_us); pPpt->Out("\n");
		//pPpt->Out("i_opt = "); pPpt->Out(i_opt); pPpt->Out("\n\n");
		//pPpt->Out("p_us = "); pPpt->Out(p_us); pPpt->Out("\t T_us = "); pPpt->Out(T_us); pPpt->Out("\n");
		//pPpt->Out("W_us = "); pPpt->Out(W_us); pPpt->Out("\t U_turb_us = "); pPpt->Out(U_turb_us); pPpt->Out("\n");
		//pPpt->Out("p0r_us = "); pPpt->Out(p0r_us); pPpt->Out("\t T0r_us = "); pPpt->Out(T0r_us); pPpt->Out("\t Wcr_us = "); pPpt->Out(Wcr_us); pPpt->Out("\n");
		
		//pPpt->Out("p0r_ds_id = "); pPpt->Out(p0r_ds_id); pPpt->Out("\t T0r_ds = "); pPpt->Out(T0r_ds); pPpt->Out("\t Wcr_ds = "); pPpt->Out(Wcr_ds); pPpt->Out("\n");
		
		//pPpt->Out("p_ds = "); pPpt->Out(p_ds); pPpt->Out("\t T_ds = "); pPpt->Out(T_ds); pPpt->Out("\n");
		//pPpt->Out("W_ds = "); pPpt->Out(W_ds); pPpt->Out("\t U_turb_ds = "); pPpt->Out(U_turb_ds); pPpt->Out("\n");
	//}
	
	// Just to complete the function for now
	double K_meanline;	
	K_meanline = ((p_us-p_ds)/p_us*pow(pBN[UPSTREAM]->M,2)); 
	
	mendtime = clock();	cout<<"mEnd = "<<mendtime<<"\t";	// This is the meanline loop end time.
	melapsedtime = melapsedtime+(mendtime-mstarttime);	pPpt->Out("mElapsed time = "); pPpt->Out(melapsedtime); pPpt->Out("\n");

	//double d; d = 0;
	return K_meanline;
	//return d;
}