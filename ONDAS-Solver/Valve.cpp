// Valve.cpp: implementation of the CValve class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Tools.h"
#include "Valve.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CValve::CValve()
{

}

CValve::~CValve()
{
	delete [] vt_raw;
	delete [] vt;
	delete [] Cd;
}

// Copy constructor
CValve::CValve(const CValve& inValve)
{
	// Configuration
	// =============
	ID = inValve.ID;
	EX = inValve.EX;
	open = inValve.open;

	end = inValve.end;
	pPipe = inValve.pPipe;
	pBN = inValve.pBN;
	pCLIN = inValve.pCLIN;
	pCLOUT = inValve.pCLOUT;
	pend_flow = inValve.pend_flow;
	pPathLine = inValve.pPathLine;

	FP = inValve.FP;
	
	// Algorithm control
	// =================
	num_loops_current = inValve.num_loops_current;
	num_loops_max = inValve.num_loops_max;

	// Valve timing related variables
	// ==============================
	ref_dia = inValve.ref_dia;
	lash = inValve.lash;
	cam_timing_ang = inValve.cam_timing_ang;
	cam_timing_loc = inValve.cam_timing_loc;

	//lift_multiplier = inValve.lift_multiplier;
	//area_multiplier = inValve.area_multiplier;
	opening_multiplier = inValve.opening_multiplier;
	ang_multiplier = inValve.ang_multiplier;

	datapoints = inValve.datapoints;
//	vt_raw = inValve.vt_raw;
//	vt = inValve.vt;
//	datapoints_flow = inValve.datapoints_flow;
//	Cd = inValve.Cd;
	
	THETA = inValve.THETA;
	LIFT = inValve.LIFT;
	AREA = inValve.AREA;

	CAD = inValve.CAD;
	TARGET = inValve.TARGET;
	RECORDED = inValve.RECORDED;
	ERROR = inValve.ERROR;
	ERROR_GRAD = inValve.ERROR_GRAD;
	SIGN_CHANGE = inValve.SIGN_CHANGE;
	CORR_FACTOR = inValve.CORR_FACTOR;
	PSI_UNCORR = inValve.PSI_UNCORR;
	PSI_CORR = inValve.PSI_CORR;
	PSI_CORR_PREV = inValve.PSI_CORR_PREV;

	REF = inValve.REF;
	FORWARD_CD = inValve.FORWARD_CD;
	REVERSE_CD = inValve.REVERSE_CD;

	raw_lift = inValve.raw_lift;
	final_lift = inValve.final_lift;
	valve_Cf_or_Cd = inValve.valve_Cf_or_Cd;
	eff_area = inValve.eff_area;
	dmdt = inValve.dmdt;
}

CValve& CValve::operator=(const CValve& inValve)
{
	if(this != &inValve)
	{
		// Configuration
		// =============
		ID = inValve.ID;
		EX = inValve.EX;
		open = inValve.open;

		end = inValve.end;
		pPipe = inValve.pPipe;
		pBN = inValve.pBN;
		pCLIN = inValve.pCLIN;
		pCLOUT = inValve.pCLOUT;
		pend_flow = inValve.pend_flow;
		pPathLine = inValve.pPathLine;

		FP = inValve.FP;
	
		// Algorithm control
		// =================
		num_loops_current = inValve.num_loops_current;
		num_loops_max = inValve.num_loops_max;

		// Valve timing related variables
		// ==============================
		ref_dia = inValve.ref_dia;
		lash = inValve.lash;
		cam_timing_ang = inValve.cam_timing_ang;
		cam_timing_loc = inValve.cam_timing_loc;

		//lift_multiplier = inValve.lift_multiplier;
		//area_multiplier = inValve.area_multiplier;
		opening_multiplier = inValve.opening_multiplier;
		ang_multiplier = inValve.ang_multiplier;

		datapoints = inValve.datapoints;
//		vt_raw = inValve.vt_raw;
//		vt = inValve.vt;
//		datapoints_flow = inValve.datapoints_flow;
//		Cd = inValve.Cd;
	
		THETA = inValve.THETA;
		LIFT = inValve.LIFT;
		AREA = inValve.AREA;

		CAD = inValve.CAD;
		TARGET = inValve.TARGET;
		RECORDED = inValve.RECORDED;
		ERROR = inValve.ERROR;
		ERROR_GRAD = inValve.ERROR_GRAD;
		SIGN_CHANGE = inValve.SIGN_CHANGE;
		CORR_FACTOR = inValve.CORR_FACTOR;
		PSI_UNCORR = inValve.PSI_UNCORR;
		PSI_CORR = inValve.PSI_CORR;
		PSI_CORR_PREV = inValve.PSI_CORR_PREV;

		REF = inValve.REF;
		FORWARD_CD = inValve.FORWARD_CD;
		REVERSE_CD = inValve.REVERSE_CD;

		raw_lift = inValve.raw_lift;
		final_lift = inValve.final_lift;
		valve_Cf_or_Cd = inValve.valve_Cf_or_Cd;
		eff_area = inValve.eff_area;
		dmdt = inValve.dmdt;
	}
	return *this;
}

void CValve::Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rVPIPES, int** &rVPIPES_ENDS, double* &rENDCORR,
						int id, bool ex, int npipes, CEngine* EngPtr, std::string param_dir/*, char* vt_dir*/, int assyid, int cyl_id, string parent_assy_res_dir, 
						int nPipesInAssy)
//--------------------------------------------------//
// Poppet valve	initialisation						//
// ---------------------------						//
// Loads values from properties and connects the	//
// valve to the relevant pipe.						//
//--------------------------------------------------//
{
	// Label labels
	// ------------
	THETA = 0; LIFT = 1; AREA = 1; // Lift and area are the same on purpose
	
	REF = 0; FORWARD_CD = 1; REVERSE_CD = 2;

	CAD = 0; TARGET = 1; RECORDED = 2; ERROR = 3; ERROR_GRAD = 4; SIGN_CHANGE = 5;
	CORR_FACTOR = 6; PSI_UNCORR = 7; PSI_CORR = 8; PSI_CORR_PREV = 9;// Active valve arrays


	// Member variables to be initialised
	// ----------------------------------
	CYL_ID = cyl_id;	// ID of cylinder on which the valve is located
	pEng = EngPtr;
	num_loops_max = 0;
	num_loops_current = 0;
	dmdt = 0;
	final_lift = 0;
	valve_Cf_or_Cd = 0;
	eff_area = 0;

//std::cout << ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->VT_DIR/*vt_dir*/), pEng->EV_FILE) << std::endl;
//exit(1);

	// Valve timing related variables
	// ------------------------------
	if(ex) {
		EX = true;
		ref_dia = pEng->REF_DIA_EXH;
		L_OR_A = pEng->EV_L_OR_A;
		lash = pEng->LASH_EXH;
		cam_timing_ang = pEng->CAM_TNG_ANG_EXH;
		cam_timing_loc = pEng->CAM_TNG_LOC_EXH;
		onum = pEng->ONUM_EXH;
		const_ref_area = pEng->CONST_REF_AREA_EXH;
		//lift_multiplier = pEng->LIFT_MULT_EXH;
		//area_multiplier = pEng->AREA_MULT_EXH;
		opening_multiplier = pEng->EV_OPENING_MULT;
		ang_multiplier = pEng->ANG_MULT_EXH;
	}
	else {
		EX = false;
		ref_dia = pEng->REF_DIA_INT;
		L_OR_A = pEng->IV_L_OR_A;
		lash = pEng->LASH_INT;
		cam_timing_ang = pEng->CAM_TNG_ANG_INT;
		cam_timing_loc = pEng->CAM_TNG_LOC_INT;
		onum = pEng->ONUM_INT;
		const_ref_area = pEng->CONST_REF_AREA_INT;
		//lift_multiplier = pEng->LIFT_MULT_INT;
		//area_multiplier = pEng->AREA_MULT_INT;
		opening_multiplier = pEng->IV_OPENING_MULT;
		ang_multiplier = pEng->ANG_MULT_INT;
	}	

	// Configuration
	if(EX || (!EX && pEng->IPIPE))
	{
		InitialiseGen(pPpt, pPipes, rPipe, rVPIPES, rVPIPES_ENDS, rENDCORR, id, ex, npipes, assyid, "CValve", parent_assy_res_dir);
		
		// Boundary name
		NAME = new int [NPIPES];
		if (EX) NAME[ONE_SIDE] = VALVE; //EXHVALVE;
		else NAME[ONE_SIDE] = VALVE; //INTVALVE;
	}
	else // Perform initialisation for valve that exists but doesn't connect to pipe
	{
		// Identification
		ID = id;
		EX = ex;
		AssyID = assyid;
		
		// Construct the correct assembly results directory
/*		
		std::string dir_str;
		dir_str = "res_assembly";
		
		if(AssyID>=10) dir_str += int(AssyID/10) + 48;
		dir_str += (AssyID - int(AssyID/10)*10) + 48;
		dir_str += "\\";
		//RES_DIR = StringToChar(pPpt->case_dir + dir_str); // Set the results directory for this object
		RES_DIR = StringToChar(pPpt->case_res_dir + dir_str); // Set the results directory for this object
		//	std::cout << dir_str << std::endl;
		//	cout << RES_DIR << endl;
*/
		RES_DIR = StringToChar(parent_assy_res_dir);
	}

	// Load valve timing
	// -----------------
	if(EX)
	{
		LoadValveTiming(ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->VT_DIR), pEng->EV_FILE));
		if(id==0) pEng->EVO = VO; // VO is read in LoadLiftArrays; if this is the first valve, set engine EVO
		if(L_OR_A) LoadFlowArrays(ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->VT_DIR), pEng->EVFA_FILE));
	}
	else 
	{
		LoadValveTiming(ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->VT_DIR), pEng->IV_FILE));
		if(id==0)
		{
			pEng->IVO = VO; // VO is read in LoadLiftArrays; if this is the first valve, set engine IVO
			pEng->IVC = VC; // VC is read in LoadLiftArrays; if this is the first valve, set engine IVC
		}
		if(L_OR_A) LoadFlowArrays(ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->VT_DIR), pEng->IVFA_FILE));
	}

	//if (!CONSTANT) period = (1 / f); // Otherwise already read in directly from file

/*
	// Labels
	TIME = 0;
	TIME_IN_PERIOD = 1;
	TARGET_P = 2;
	APPLIED_P = 3;
	RECORDED_P = 4;
	NEXT_P = 5;
	THIS_P = 6;
	TARGET_T = 7;
	TARGET_T0 = 8;
	APPLIED_T = 9;
	RECORDED_T = 10;
	NEXT_T = 11;
	THIS_T = 12;
	TARGET_V = 13;
	APPLIED_V = 14;
	RECORDED_V = 15;
	NEXT_V = 16;
	THIS_V = 17;
	TARGET_MFR = 18;
	APPLIED_MFR = 19;
	RECORDED_MFR = 20;
	NEXT_MFR = 21;
	THIS_MFR = 22;
	VEL_GRAD = 23;
	P_VEL = 24;
	P_VEL_GRAD = 25;
	UNDER_EST = 26;
	UNDER_EST_PREV = 27;
	UNDER_EST_FACT = 28;
	K_LOSS = 29;

	//WAVE_TIME = 0;
	//WAVE_PARAMETER = 1;
	//WAVE_MFR = 1;

	// AVT
	AVTError = 0;
	AVTErrorPrev = 0;
	AVTErrorGrad = 0;

	target_applied_recorded = new double* [pEng->nSamples];
	for (int i = 0; i < pEng->nSamples; ++i){
		target_applied_recorded[i] = new double[30];
		// time, time_in_period, target, applied, recorded, next

		// Set appropriate times and unit values in these arrays
		target_applied_recorded[i][TIME] = 0; // +i * (period / nSamples);
		target_applied_recorded[i][TIME_IN_PERIOD] = 0; // +i * (period / nSamples);

		target_applied_recorded[i][TARGET_P] = 1;
		target_applied_recorded[i][APPLIED_P] = 1;
		target_applied_recorded[i][RECORDED_P] = 1;
		target_applied_recorded[i][NEXT_P] = 1;
		target_applied_recorded[i][THIS_P] = 1;

		target_applied_recorded[i][TARGET_T] = 0;
		target_applied_recorded[i][TARGET_T0] = 0;
		target_applied_recorded[i][APPLIED_T] = 0;
		target_applied_recorded[i][RECORDED_T] = 0;
		target_applied_recorded[i][NEXT_T] = 0;
		target_applied_recorded[i][THIS_T] = 0;

		target_applied_recorded[i][TARGET_V] = 0;
		target_applied_recorded[i][APPLIED_V] = 0;
		target_applied_recorded[i][RECORDED_V] = 0;
		target_applied_recorded[i][NEXT_V] = 0;
		target_applied_recorded[i][THIS_V] = 0;

		target_applied_recorded[i][TARGET_MFR] = 0;
		target_applied_recorded[i][APPLIED_MFR] = 0;
		target_applied_recorded[i][RECORDED_MFR] = 0;
		target_applied_recorded[i][NEXT_MFR] = 0;
		target_applied_recorded[i][THIS_MFR] = 0;

		target_applied_recorded[i][VEL_GRAD] = 0;
		target_applied_recorded[i][P_VEL] = 0;
		target_applied_recorded[i][P_VEL_GRAD] = 0;
		target_applied_recorded[i][UNDER_EST] = 0;
		target_applied_recorded[i][UNDER_EST_PREV] = 0;
		target_applied_recorded[i][UNDER_EST_FACT] = 1;
		target_applied_recorded[i][K_LOSS] = 1;
	}
*/

	SONIC = false;

	// Set up files
	PRINT_MOVIE_FILE = false; // Default

	std::string bcname_str = "res_cyl";
	bcname_str += IntToString(CYL_ID);
	if(EX) bcname_str += "_exh"; else bcname_str += "_int";
	bcname_str += "valve";
	SetupFiles(pPpt, bcname_str);
}

void CValve::InitialiseVolumeValve(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rVPIPES, int** &rVPIPES_ENDS, double* &rENDCORR,
						int id, bool ex, int npipes, std::string param_dir, int assyid, int vol_id, string parent_assy_res_dir)
//--------------------------------------------------//
// Poppet valve	initialisation						//
// ---------------------------						//
// Loads values from properties and connects the	//
// valve to the relevant pipe.						//
//--------------------------------------------------//
{
/*
	// Label labels
	// ------------
	THETA = 0; LIFT = 1;
	REF = 0; FORWARD_CD = 1; REVERSE_CD = 2;
*/
	// Member variables to be initialised
	// ----------------------------------
	VOL_ID = vol_id;	// ID of volume on which the valve is located
	//pEng = EngPtr;
/*
	dmdt = 0;
	final_lift = 0;
	valve_Cf_or_Cd = 0;
	eff_area = 0;
*/
//std::cout << ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->VT_DIR/*vt_dir*/), pEng->EV_FILE) << std::endl;
//exit(1);

	// Valve timing related variables
	// ------------------------------
	if(ex)
	{
		EX = true;
/*
		ref_dia = pEng->REF_DIA_EXH;
		lash = pEng->LASH_EXH;
		cam_timing_ang = pEng->CAM_TNG_ANG_EXH;
		onum = pEng->ONUM_EXH;
		lift_multiplier = pEng->LIFT_MULT_EXH;
		ang_multiplier = pEng->ANG_MULT_EXH;
		L_OR_A = pEng->EV_L_OR_A;
		const_ref_area = pEng->CONST_REF_AREA_EXH;
		LoadLiftArrays(ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->VT_DIR), pEng->EV_FILE));
		if(id==0) pEng->EVO = VO; // VO is read in LoadLiftArrays; if this is the first valve, set engine EVO
		LoadFlowArrays(ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->VT_DIR/), pEng->EVFA_FILE));
*/
	}
	else 
	{
		EX = false;
/*
		ref_dia = pEng->REF_DIA_INT;
		lash = pEng->LASH_INT;
		cam_timing_ang = pEng->CAM_TNG_ANG_INT;
		//cam_timing_ang = pEng->IVO;
		onum = pEng->ONUM_INT;
		lift_multiplier = pEng->LIFT_MULT_INT;
		ang_multiplier = pEng->ANG_MULT_INT;
		L_OR_A = pEng->IV_L_OR_A;
		const_ref_area = pEng->CONST_REF_AREA_INT;
		LoadLiftArrays(ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->VT_DIR), pEng->IV_FILE));
		if(id==0)
		{
			pEng->IVO = VO; // VO is read in LoadLiftArrays; if this is the first valve, set engine IVO
			pEng->IVC = VC; // VC is read in LoadLiftArrays; if this is the first valve, set engine IVC
		}
		LoadFlowArrays(ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->VT_DIR), pEng->IVFA_FILE));
*/
	}	

	// Configuration
//	if(EX || (!EX && pEng->IPIPE))
//	{
		InitialiseGen(pPpt, pPipes, rPipe, rVPIPES, rVPIPES_ENDS, rENDCORR, id, ex, npipes, assyid, "CValve", parent_assy_res_dir);
		
		// Boundary name
		NAME = new int [NPIPES];
		NAME[ONE_SIDE] = VOLUMEVALVE;
//	}
/*
	else // Perform initialisation for valve that exists but doen't connect to pipe
	{
		// Identification
		ID = id;
		EX = ex;
		AssyID = assyid;
		
		// Construct the correct assembly results directory
		std::string dir_str;
		dir_str = "res_assembly";
		
		if(AssyID>=10) dir_str += int(AssyID/10) + 48;
		dir_str += (AssyID - int(AssyID/10)*10) + 48;
		dir_str += "\\";
		//RES_DIR = StringToChar(pPpt->case_dir + dir_str); // Set the results directory for this object
		RES_DIR = StringToChar(pPpt->case_res_dir + dir_str); // Set the results directory for this object
		//	std::cout << dir_str << std::endl;
		//	cout << RES_DIR << endl;
	}
*/
/*
	SONIC = false;

	// Set up files
	PRINT_MOVIE_FILE = false; // Default

	std::string bcname_str = "res_cyl";
	bcname_str += IntToString(CYL_ID);
	if(EX) bcname_str += "_exh"; else bcname_str += "_int";
	bcname_str += "valve";
	SetupFiles(pPpt, bcname_str);
*/
}

void CValve::ConfigureExtra(CProperties* pPpt, CPipe** pExhaustSystem, CPipe** pIntakeSystem, int nExPipes, int nInPipes)
{
	// Record area of pipe to which valve is attached
	FP = pBN[ONE_SIDE]->f_dash * pPpt->fref;


	// AVT
	// ----------------------------------------------------------------------------------------------------

	if (this->EX) ACTIVE = pEng->EV_ACTIVE; else ACTIVE = pEng->IV_ACTIVE;		

	if (ACTIVE) { // Only set up if this valve is used for AVT

		// Is this valve's target location in the exhaust or intake system
		if (this->EX) {
			TARGET_PIPE_EX = pEng->EV_TARGET_PIPE_EX;
			target_pipe_num = pEng->EV_target_pipe_num;
			target_pipe_loc = pEng->EV_target_pipe_loc;
		}
		else {
			TARGET_PIPE_EX = pEng->IV_TARGET_PIPE_EX;
			target_pipe_num = pEng->IV_target_pipe_num;
			target_pipe_loc = pEng->IV_target_pipe_loc;
		}

		if (TARGET_PIPE_EX) target_pipe_len = pExhaustSystem[target_pipe_num]->length;
		else target_pipe_len = pIntakeSystem[target_pipe_num]->length;

		if (this->EX) {
			if (pEng->EV_ACTIVE) {

				// Load AVT target waveform
				LoadActiveValveTarget(ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->VT_DIR), pEng->EV_ACTIVE_FILE_P));

				// Initialise target tapping

				// Check is chosen target pipe is actually available
				pEng->nExPipesInAssy = nExPipes;
				pEng->nInPipesInAssy = nInPipes;
				if ((TARGET_PIPE_EX && target_pipe_num >= pEng->nExPipesInAssy) || (!TARGET_PIPE_EX && target_pipe_num >= pEng->nInPipesInAssy)) {	// If not, exit safely and print reasonto screen

					pPpt->Out("In the input file for "); pPpt->Out(pEng->Identify()); pPpt->Out(" target_pipe_num is set to "); pPpt->Out(target_pipe_num);
					pPpt->Out(" but this pipe does not exist in the assembly. There are only ");
					if (TARGET_PIPE_EX) pPpt->Out(pEng->nExPipesInAssy); else pPpt->Out(pEng->nInPipesInAssy);
					pPpt->Out(" pipe(s) in the assembly. Exiting.\n");
					exit(1);
				}

				// Locate the two nodes either side of the location of this tapping
				int S = 0;
				bool located;

				if (TARGET_PIPE_EX) pTargetPipeSys = pExhaustSystem; // Target location is in exhaust system
				else pTargetPipeSys = pIntakeSystem;				
//				cout << "pTargetPipeSys[0]->N = " << pTargetPipeSys[0]->N << endl; exit(1);

				if (pTargetPipeSys[target_pipe_num]->N > 1) {

					SINGLE_NODED_PIPE = false;
					located = false;
					while (!located) {
						if (pTargetPipeSys[target_pipe_num]->Node[S + 1].x > (target_pipe_loc * pTargetPipeSys[target_pipe_num]->length) ||
							fabs(pTargetPipeSys[target_pipe_num]->Node[S + 1].x - (target_pipe_loc * pTargetPipeSys[target_pipe_num]->length)) < 1e-6)
							located = true;
						else ++S;
					};

					// Then the location is between nodes S and S+1	
					pAVTNodeLeft = &(pTargetPipeSys[target_pipe_num]->Node[S]);
					pAVTNodeRight = &(pTargetPipeSys[target_pipe_num]->Node[S + 1]);
				}
				else
				{
					SINGLE_NODED_PIPE = true;
					// For single node pipes make these the same:
					pAVTNodeLeft = &pTargetPipeSys[target_pipe_num]->Node[S];
					pAVTNodeRight = &pTargetPipeSys[target_pipe_num]->Node[S];
				}
				//cout << "AVTMeasureNode.p_dash * pPpt->PREF = " << AVTMeasureNode.p_dash * pPpt->PREF << " bar\n"; exit(1);

				// Initialize counters and flags
				AVTCumulativeVal = 0; 
				AVTInCycleStepCount = 0;
				recorded_cyc_av = 1; // 0;
				cyc_av_error = 0;
				FIRST_CYCLE = true;
				cyc_av_factor = 5; // 2; // 3;// 5;//  2.5; //10
				//AVTFactorSF = 0; // 0.2;//  0.25;
				cyc_av_factor_orig = cyc_av_factor;
				//AVTFactorSFOrig = AVTFactorSF; 
				cyc_av_error = 0;
				cyc_av_error_prev = 0;
				cyc_av_error_grad = 0;
				CYC_AV_MATCHED = false;
				CYC_AV_ERR_SIGN_CHANGE = false;

				PHASE_MATCHING = false;
				PHASE_MATCHED = false;
				phase_shift = 0;

				plenum_press = pEng->pres;
				plenum_press_prev = plenum_press;
			}
		}
		else {
			if (pEng->IV_ACTIVE) {

				// FILL THIS IN ONCE EXHAUST AVT IS COMPLETE
	/*
				LoadActiveValveTarget(ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->VT_DIR), pEng->IV_ACTIVE_FILE_P));
				pEng->nInPipesInAssy = nPipesInAssy; // Initialise tapping
				if (pEng->target_in_pipe >= pEng->nInPipesInAssy) {	// Check is chosen target pipe is actually available
					// If not, exit safely and print reasonto screen
					pPpt->Out("In the input file for ");
					pPpt->Out(pEng->Identify());
					pPpt->Out(" target_in_pipe is set to ");
					pPpt->Out(pEng->target_in_pipe);
					pPpt->Out(" but this pipe does not exist in the assembly. There are only ");
					pPpt->Out(pEng->nInPipesInAssy);
					pPpt->Out(" pipe(s) in the assembly. Exiting.\n");
					exit(1);
				}
	*/
			}
		}
	}
}

void CValve::ListProperties(CProperties* pPpt, char* prefix)
{
	pPpt->Out(Underline(Identify(), "=", "\t\t")); pPpt->Out("\n");
	pPpt->Out(prefix); pPpt->Out("\tReference diameter (ref_dia)\t\t=\t"); pPpt->Out(ref_dia); pPpt->Out(" mm\n");
	pPpt->Out(prefix); pPpt->Out("\tInitial effective area (eff_area)\t=\t"); pPpt->Out(eff_area); pPpt->Out(" m^2\n");
	pPpt->Out(prefix); pPpt->Out("\tAdjoining pipe area (FP)\t\t=\t"); pPpt->Out(FP); pPpt->Out(" m^2\n");
	pPpt->Out(prefix); pPpt->Out("\tArea ratio (psi_val)\t\t\t=\t"); pPpt->Out(psi_val); pPpt->Out("\n");
	pPpt->Out("\n");
	if (EX) {
		if (pEng->EV_ACTIVE) {
			pPpt->Out(Underline("Exhaust AVT setup", "-", "\t\t"));
			pPpt->Out(prefix); pPpt->Out("\tTarget average (target_cyc_av)\t=\t"); pPpt->Out(target_cyc_av); pPpt->Out(" bar\n");
			pPpt->Out("\n");
		}
	}
	else {
		if (pEng->IV_ACTIVE) {
			pPpt->Out(Underline("Intake AVT setup", "-", "\t\t"));
			pPpt->Out(prefix); pPpt->Out("\tTarget average (target_cyc_av)\t=\t"); pPpt->Out(target_cyc_av); pPpt->Out(" bar\n");
			pPpt->Out("\n");
		}
	}
	pPpt->Out("\n");
}

void CValve::LoadValveTiming(char* InputFile)
//--------------------------------------------------//
// Poppet valve	vt arrays							//
// ------------------------							//
// Loads lift or area arrays from file, i.e.  		//
// lift (mm) or area (m^2) as a function			//
// of theta (degrees).								//
//--------------------------------------------------//
{
	int c = 0;
	//float temp;
	double temp;
	int num_columns = 2;
	int col, row;
	datapoints = 0;
	
	FILE* stream;
	stream = fopen(InputFile, "r");

	if(stream == NULL)
	{
		printf("Error opening valve timing file\n");
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
			fscanf(stream, "%lf", &temp);		// Runs over the theta value
			fscanf(stream, "%lf", &temp);		// Runs over the vt value
//cin >> pause;
		}while(fscanf(stream, "%l")!=EOF);
//cout << "Lift arrays datapoints = " << datapoints << endl;
//cin >> pause;
		// Can now dimension timing arrays
		vt_raw = new double* [num_columns]; 
		for(col=0; col<num_columns; ++col) vt_raw[col] = new double [datapoints];

		fseek(stream, 0L, SEEK_SET);	// Reset pointer to beginning of file

		for (row=0; row<datapoints; row++)
		{
			fscanf(stream, "%d", &c);	// Scans past the line number
			for (col=0; col<num_columns; col++)
			{
				// Finds the angle and the lift or area values
				fscanf(stream, "%lf", &vt_raw[col][row]);	
//cout << "vt_raw[col][row] = " << vt_raw[col][row] << endl;
			}
		}
	}
	fclose(stream);

	// Identify raw angle of start of valve opening
	row = 0;
	while(vt_raw[LIFT/*or AREA*/][row]==0 && row<datapoints) ++row;
	double start_opening_angle, start_opening_value;
	if(row>0)
	{
		start_opening_angle = vt_raw[THETA][row-1];
		start_opening_value = vt_raw[LIFT/*or AREA*/][row-1];
	}
	else
	{
		start_opening_angle = vt_raw[THETA][row];
		start_opening_value = vt_raw[LIFT/*or AREA*/][row];
	}

	// Identify raw angle of end of valve closing - start from current row
	while(vt_raw[LIFT/*or AREA*/][row]>0 && row<datapoints) ++row;
	double end_closing_angle, end_closing_value;
	end_closing_angle = vt_raw[THETA][row];
	end_closing_value = vt_raw[LIFT/*or AREA*/][row];

	// Identify raw angle of raw maximum opening
	double max_opening_angle = 0;
	double max_opening_value = 0;
	for(row=0; row<datapoints; ++row)
	{
		if(vt_raw[LIFT/*or AREA*/][row]>max_opening_value)
		{
			max_opening_angle = vt_raw[THETA][row];
			max_opening_value = vt_raw[LIFT/*or AREA*/][row];
		}
	}

	// Now adjust raw values of theta and lift or area
	vt = new double* [num_columns];
	for(col=0; col<num_columns; ++col) vt[col] = new double [datapoints];

	for(row=0; row<datapoints; ++row)
	{
		// Shift raw angle data according to cam_tng_ang and cam_tng_loc
		switch(cam_timing_loc)
		{
		case 0:// cam_tng_ang gives start of valve opening
			vt[THETA][row] = vt_raw[THETA][row] + (cam_timing_ang - start_opening_angle);
			break;
		case 1:// cam_tng_ang gives maximum opening
			vt[THETA][row] = vt_raw[THETA][row] + (cam_timing_ang - max_opening_angle);
			break;
		case 2:// cam_tng_ang gives end of valve closing
			vt[THETA][row] = vt_raw[THETA][row] + (cam_timing_ang - end_closing_angle);
			break;
		default:// cam_tng_ang gives maximum opening
			cam_timing_loc = 1;
			vt[THETA][row] = vt_raw[THETA][row] + (cam_timing_ang - max_opening_angle);
			break;
		}
		// Now adjust the spread of the opening period with ang_multiplier
		vt[THETA][row] = ((vt[THETA][row] - cam_timing_ang)*ang_multiplier) + cam_timing_ang;

		if(L_OR_A) {
			// lift = lift array * lift multiplier - valve lash (all in mm);
			//if(vt_raw[LIFT][row]*lift_multiplier - lash < 0) vt[LIFT][row] = 0;
			//else vt[LIFT][row] = vt_raw[LIFT][row]*lift_multiplier - lash; 
			if(vt_raw[LIFT][row]*opening_multiplier - lash < 0) vt[LIFT][row] = 0;
			else vt[LIFT][row] = vt_raw[LIFT][row]* opening_multiplier - lash;
		}
		else // File contains area values (m^2)
			//vt[AREA][row] = vt_raw[AREA][row]*area_multiplier;
			vt[AREA][row] = vt_raw[AREA][row] * opening_multiplier;
	}

	// Calculate adjusted valve opening and closing
	switch(cam_timing_loc)
	{
	case 0:// cam_tng_ang gives start of valve opening
		VO = start_opening_angle + (cam_timing_ang - start_opening_angle);
		VC = end_closing_angle + (cam_timing_ang - start_opening_angle);
		break;
	case 1:// cam_tng_ang gives maximum opening
		VO = start_opening_angle + (cam_timing_ang - max_opening_angle);
		VC = end_closing_angle + (cam_timing_ang - max_opening_angle);
		break;
	case 2:// cam_tng_ang gives end of valve closing
		VO = start_opening_angle + (cam_timing_ang - end_closing_angle);
		VC = end_closing_angle + (cam_timing_ang - end_closing_angle);
		break;
	default:// cam_tng_ang gives maximum opening
		VO = start_opening_angle + (cam_timing_ang - max_opening_angle);
		VC = end_closing_angle + (cam_timing_ang - max_opening_angle);
		break;
	}
	// Now adjust for any spread
	VO = ((VO - cam_timing_ang)*ang_multiplier) + cam_timing_ang;
	VC = ((VC - cam_timing_ang)*ang_multiplier) + cam_timing_ang;
			
	if(pEng->PRE_vt) {
		cout << Underline("Raw valve timing", "-", "\t");
		cout << "\tRaw angle of start of valve opening, start_opening_angle\t=\t" << start_opening_angle << Deg() << endl;
		cout << "\tStart of opening value, start_opening_value\t\t\t=\t" << start_opening_value;
		if(L_OR_A) cout << " mm"; else cout << " m^2"; cout << endl;
		cout << endl;
		cout << "\tRaw angle of maximum opening, max_opening_angle\t\t\t=\t" << max_opening_angle << Deg() << endl;
		cout << "\tMaximum opening value, max_opening_value\t\t\t=\t" << max_opening_value;
		if(L_OR_A) cout << " mm"; else cout << " m^2"; cout << endl;
		cout << endl;
		cout << "\tRaw angle of end of valve closing, end_closing_angle\t\t=\t" << end_closing_angle << Deg() << endl;
		cout << "\tEnd of closing value, end_closing_value\t\t\t\t=\t" << end_closing_value;
		if(L_OR_A) cout << " mm"; else cout << " m^2"; cout << endl;
		cout << endl;

		cout << Underline("Adjusted valve timing", "-", "\t");
		cout << "\tInput file = " << InputFile << endl << endl;
		cout << "\tCrank angle (" << Deg() << "CA)\t";
		if(L_OR_A) cout << "Lift (mm)"; else cout << "Area (m^2)" << endl;
		for(row=0; row<datapoints; ++row) cout << "\t" << vt[THETA][row] << "\t\t\t" << vt[LIFT/*or AREA*/][row] << endl;
		cout << endl;
		
		cout << "\tValve opening (VO) = " << VO << endl;
		cout << "\tValve closing (VC) = " << VC << endl;
		cout << endl;
	}
}

void CValve::LoadFlowArrays(char* InputFile)
//--------------------------------------------------//
// Poppet valve	flow arrays							//
// ------------------------							//
// Loads flow arrays from file, i.e. forward and 	//
// reverse coefficients of drag for a range of		//
// reference lift/diameter values.					//
//--------------------------------------------------//
{
	int c = 0;
	double temp;
	int num_columns = 3;
	int col, row;
	datapoints_flow = 0;
	
	FILE* stream;
	stream = fopen(InputFile, "r");

	if(stream == NULL)
	{
		printf("Error opening valve flow arrays file\n");
		exit(1);
	}
	else
	{
		fseek(stream, 0L, SEEK_SET);	// Set pointer to beginning of file
		do
		{
			fscanf(stream, "%d", &c);	// Scans for the line number (ints), puts it in c
//cout << "c = " << c << endl;
			if(c>datapoints_flow) datapoints_flow = c;
			fscanf(stream, "%lf", &temp);		// Runs over the ref value
//cout << "temp = " << temp << endl;
			fscanf(stream, "%lf", &temp);		// Runs over the forward Cd value
//cout << "temp = " << temp << endl;
			fscanf(stream, "%lf", &temp);		// Runs over the reverse Cd value
//cout << "temp = " << temp << endl;
//cout << endl;
//cin >> pause;
		}while(fscanf(stream, "%l")!=EOF);
		//cout << "Flow arrays datapoints = " << datapoints_flow << endl;
		//cin >> pause;
		
		// Can now dimension flow arrays
		Cd = new double* [num_columns];
		for(col=0; col<num_columns; ++col) Cd[col] = new double [datapoints_flow];

		fseek(stream, 0L, SEEK_SET);	// Reset pointer to beginning of file

bool PRINT_DATA = false;//false;//true;
if(PRINT_DATA)
{
	cout << endl;
	cout << InputFile << ":" << endl;
	cout << "REF(LIFT/DIA)\tFORWARD CD\tREVERSE CD\n";
}
		for(row=0; row<datapoints_flow; row++)
		{
			fscanf(stream, "%d", &c);	// Scans past the line number
			// Finds the reference value, the forward Cd value, and the reverse Cd value
			for (col=0; col<num_columns; col++)
			{
				fscanf(stream, "%lf", &Cd[col][row]);
			}
if(PRINT_DATA) cout << Cd[REF][row] << "\t\t" << Cd[FORWARD_CD][row] << "\t\t" << Cd[REVERSE_CD][row] << endl;
		}
	}
	fclose(stream);
}

void CValve::LoadActiveValveTarget(char* InputFile)
//--------------------------------------------------//
// Poppet valve	vt arrays							//
// ------------------------							//
// Loads taget pressure from file, i.e.  			//
// pressure (bar) as a function						//
// of theta (CA degrees).							//
//--------------------------------------------------//
{
	int c = 0;
	double temp;
	int col, row;

	int num_columns = 10; // CAD, TARGET, RECORDED, ERROR, ERROR_GRAD, SIGN_CHANGE, CORR_FACTOR, PSI_UNCORR, PSI_CORR, PSI_CORR_PREV
	int num_columns_in_file = 2; // CAD, TARGET
	numAVTDataPoints = 0;

//int pause;

	FILE* stream;
	stream = fopen(InputFile, "r");

	if (stream == NULL) {
		printf("Error opening active valve target file\n");
		exit(1);
	}
	else{
		fseek(stream, 0L, SEEK_SET);						// Set pointer to beginning of file
		do{
			fscanf(stream, "%d", &c);						// Scans for the line number (ints), puts it in c
			if (c > numAVTDataPoints) numAVTDataPoints = c; // Number of AVT datapoints
			fscanf(stream, "%lf", &temp);					// Runs over the CAD value
			fscanf(stream, "%lf", &temp);					// Runs over the target value
		}while (fscanf(stream, "%l") != EOF);
		
		// Can now dimension arrays
		avt_raw = new double* [numAVTDataPoints];  //avt_raw = new double* [num_columns];
		for (row = 0; row < numAVTDataPoints; ++row) { //for (col = 0; col < num_columns; ++col) {
			avt_raw[row] = new double[num_columns]; //avt_raw[col] = new double[numAVTDataPoints];
		}

		fseek(stream, 0L, SEEK_SET);						// Reset pointer to beginning of file

		for(row = 0; row < numAVTDataPoints; row++){
			fscanf(stream, "%d", &c);	// Scans past the line number
			
			// Read the CAD and TARGET values from the file
			for (col = 0; col < num_columns_in_file; col++){
				fscanf(stream, "%lf", &avt_raw[row][col]);  //fscanf(stream, "%lf", &avt_raw[col][row]);
//cout << "avt_raw[col=" << col << "][row=" << row << " ] = " << avt_raw[col][row] << endl;
//cout << "avt_raw[row=" << row << "][col=" << col << "] = " << avt_raw[row][col] << "\t";
			}

			// Initialize the other columns that are not read from file
			avt_raw[row][RECORDED] = 0;
			avt_raw[row][ERROR] = 0;
			avt_raw[row][ERROR_GRAD] = 0;
			avt_raw[row][SIGN_CHANGE] = double(false);
			avt_raw[row][CORR_FACTOR] = 0.1; // 1; // Unity = no correction
			avt_raw[row][PSI_UNCORR] = 0;
			avt_raw[row][PSI_CORR] = 0;
			avt_raw[row][PSI_CORR_PREV] = 0;

/*
for (col = 0; col < num_columns; col++) {
	cout << "avt_raw[row=" << row << "][col=" << col << "] = " << avt_raw[row][col] << "\t";
}
cout << endl;
//exit(1);
//*/
		}
	}
	fclose(stream);
//exit(1);

	// Identify raw angle of start of pressure target
	//row = 0;
	//double start_opening_angle, start_opening_value;
	//start_opening_angle = avt_raw[CAD][row];
	//start_opening_value = avt_raw[TARGET][row];

	// Calculate target cycle average, peak, peak location etc.
	
	double cumulative_TARGET = 0;
	target_max_row = 0;
	target_max_ang = 0;
	target_max_val = 0;
	for (row = 0; row < numAVTDataPoints; row++) {

		cumulative_TARGET += avt_raw[row][TARGET]; //cumulative_TARGET += avt_raw[TARGET][row];
		
		if (avt_raw[row][TARGET] > target_max_val) {
			target_max_row = row;
			target_max_ang = avt_raw[row][CAD];
			target_max_val = avt_raw[row][TARGET];
		}
												   
//cout << "cumulative_TARGET = " << cumulative_TARGET << endl;
//cout << endl;
	}
	target_cyc_av = cumulative_TARGET/numAVTDataPoints;
/*
cout << "target_cyc_av = " << target_cyc_av << endl;
cout << "target_max_val = " << target_max_val << " bar" << endl; 
cout << "target_max_loc lies at avt_raw[row=" << target_max_loc << "]" << endl;
exit(1);
//cin >> pause;
*/






/*
	// Now populate p_des_app_seen[][PARAMETER_DESIRED]
	for (int i_fill = 0; i_fill < pEng->nSamples; ++i_fill)
	{
		row = 0;
		do ++row;
		while (parameter_wave[WAVE_TIME][row] < p_des_app_seen[i_fill][TIME_IN_PERIOD] && row < datapoints_parameter - 1);

		if (target_applied_recorded[i_fill][TIME_IN_PERIOD] <= phi * period)
		{
			// Interpolate between [row-1] and [row]
			target_applied_recorded[i_fill][PARAMETER_DESIRED] = parameter_wave[WAVE_PARAMETER][row - 1] +
				((parameter_wave[WAVE_PARAMETER][row] - parameter_wave[WAVE_PARAMETER][row - 1]) /
					(parameter_wave[WAVE_TIME][row] - parameter_wave[WAVE_TIME][row - 1]))
				* (target_applied_recorded[i_fill][TIME_IN_PERIOD] - parameter_wave[WAVE_TIME][row - 1]);
		}
		else
		{
			//cout << "Outside of phi*period" << endl;
			target_applied_recorded[i_fill][PARAMETER_DESIRED] = raw_low;
		}

		// USEFUL for debugging:
		//if(ID==0)
		//cout << p_des_app_seen[i_fill][TIME_IN_PERIOD] << "\t" << p_des_app_seen[i_fill][PARAMETER_DESIRED] << endl;
*/


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
//	}
}

void CValve::EffectiveArea(CProperties* pPpt, double ca)
//--------------------------------------------------//
// Poppet valve	effective area						//
// ---------------------------						//
// Interpolates for lift then calculates			//
// effective area for the given crank angle.		//
//--------------------------------------------------//
{
	int row;
	double lookup;
	double lookup_ref;

	// angle in theta array IS ca:
	lookup = ca;

	row = 0;
	//if(lookup<vt_raw[THETA][0] || lookup>vt_raw[THETA][datapoints-1])
	if(lookup<vt[THETA][0] || lookup>vt[THETA][datapoints-1])
	{
//cout << "Valve timing interpolation: lookup out of range. Returning zero area." << endl;
		raw_lift = 0;
		final_lift = 0;
/*
		if(lookup>126 && EX && pPipe[ONE_SIDE]->ID==0)
		{
			cout << "lookup=" << lookup << endl;
			cout << "vt[THETA][0]=" << vt[THETA][0] << endl;
			cout << "vt[THETA][datapoints-1]=" << vt[THETA][datapoints-1] << endl;
			cout << "- out of range" << endl;
			cout << endl;
		}
*/
	}
	else
	{
//if(this->EX) cout << "Exhaust Valve lookup = " << lookup << endl;
		// Find the two datapoints either side of lookup
		//while(!(lookup>=vt[THETA][row] && lookup<vt[THETA][row+1]))
		row = 0;
		while((lookup>=vt[THETA][row] && row < datapoints-1) || row==0)
		{
//			if(this->EX) cout << "EV vt[THETA=" << THETA << "][row=" << row << "] = " << vt[THETA][row] << endl;
			++row;
		}
		// Interpolate between row and row-1
		raw_lift = vt_raw[LIFT][row-1] + ((lookup - vt_raw[THETA][row-1])/(vt_raw[THETA][row] - vt_raw[THETA][row-1])
											*(vt_raw[LIFT][row] - vt_raw[LIFT][row-1]));

		final_lift = vt[LIFT][row-1] + ((lookup - vt[THETA][row-1])/(vt[THETA][row] - vt[THETA][row-1])
										  *(vt[LIFT][row] - vt[LIFT][row-1]));
/*
		if(lookup>126 && EX && pPipe[ONE_SIDE]->ID==0)
		{
			cout << "vt[THETA][row] = " << vt[THETA][row] << endl;
			cout << "lookup = " << lookup << endl;
			cout << "vt[THETA][row-1] = " << vt[THETA][row-1] << endl;
			cout << "vt[LIFT][row] = " << vt[LIFT][row] << endl;
			cout << "vt[LIFT][row-1] = " << vt[LIFT][row-1] << endl;
			cout << "final_lift = " << final_lift << endl;
			cout << endl;
		}
*/
/*
		// Interpolate for lift. vt[THETA][row] is the lookup, vt[LIFT][row] is the lift(mm).
		raw_lift = vt_raw[LIFT][row] + (lookup - vt_raw[THETA][row])/(vt_raw[THETA][row+1] - vt_raw[THETA][row])
								*(vt_raw[LIFT][row+1] - vt_raw[LIFT][row]);

		final_lift = lift[LIFT][row] + (lookup - vt[THETA][row])/(vt[THETA][row+1] - vt[THETA][row])
								*(vt[LIFT][row+1] - vt[LIFT][row]);
*/
	}

	if(!L_OR_A) // If data is area values rather than lift values
	{
//cout << "raw_lift = " << raw_lift << endl;
//cout << "final_lift = " << final_lift << endl;
		eff_area = final_lift;
		valve_Cf_or_Cd = 0;
	}
	else // Lift data is lift values
	{
//if(this->EX) cout << "raw_lift = " << raw_lift << endl;
//if(this->EX) cout << "Exhaust Valve final_lift = " << final_lift << endl;
		// Now interpolate for the Cd value for this reference value (lift/dia)
		lookup_ref = final_lift/ref_dia;
//		lookup_ref = raw_lift/ref_dia;

		row = 0;
		if(lookup_ref<Cd[REF][0] || lookup_ref>Cd[REF][datapoints_flow-1])
		{
//if(this->EX) cout << "Exhaust Valve Cd interpolation: lookup_ref out of range. Returning zero area." << endl;
			valve_Cf_or_Cd = 0;
			eff_area = 0;
		}
		else
		{
			while(!(lookup_ref>=Cd[REF][row] && lookup_ref<Cd[REF][row+1])) ++row;

			int CD_DIRECTION = FORWARD_CD;
			// Interpolate for Cd. Cd[REF][row] is the lookup_ref, Cd[CD_DIRECTION][row] is the Cd
			valve_Cf_or_Cd = Cd[CD_DIRECTION][row] + (lookup_ref - Cd[REF][row])/(Cd[REF][row+1] - Cd[REF][row])
									*(Cd[CD_DIRECTION][row+1] - Cd[CD_DIRECTION][row]);
	
			// Can now calculate effective area
			if(const_ref_area)
				eff_area = onum * valve_Cf_or_Cd * (PI/4)*pow(ref_dia*1e-3,2);
			else // Curtain area:
				eff_area = onum * valve_Cf_or_Cd * ( PI*(ref_dia*1e-3)*(final_lift*1e-3) );
		}
	} 

//if(this->EX) cout << "Exhaust Valve eff_area = " << eff_area << endl;
}

void CValve::EffectiveArea(CProperties* pPpt)
//--------------------------------------------------//
// Poppet valve	effective area						//
// ---------------------------						//
// Sets effective area.								//
//--------------------------------------------------//
{
	onum = 1; // Flow area multiplier (default=1)
	valve_Cf_or_Cd = 1; // Disharge coefficient
	eff_area = onum * valve_Cf_or_Cd * (PI/4)*pow(ref_dia*1e-3,2);

//cout << "Volume valve eff_area = " << eff_area << endl;
//cout << endl;
}

void CValve::Poppet_H(CProperties* pPpt, int timestep, double time, double &rDMDT, double psi, double AREF, double PREF, double PC, double TC, double Fp)
//--------------------------------------------------//
// Homentropic poppet valve							//
// ---------------------------------------			//
// Inflow algorithm based on Benson					//
// 													//
// Outflow will revert to a homentropic nozzle		//
//--------------------------------------------------//
{
	double lambda_in = (*(pCLIN[ONE_SIDE]))[R+1];

	// Test if valve open (psi > 0)
	if(psi<=0.0) // Valve closed - treat as closed end
	{
		open = false;
//		rlambda_out = lambda_in;
		(*(pCLOUT[ONE_SIDE]))[R+1] = lambda_in;
		rDMDT = 0.0;
		dmdt = rDMDT;
		return;
	}
	else
	{
		open = true;
		if(psi>1.0) psi=1.0;// In case valve area is greater than pipe area
	}

	double K1, K2, K3, S, pi_alpha, C, D, lambda_out_temp, T;	// Declare constants
	// Test flow direction
	T = TC;
	K1 = (pPpt->gammaAir(T)-1.0)/(2.0*pPpt->gammaAir(T));
	pi_alpha = lambda_in*pow((PREF/PC), K1);
	if(pi_alpha<1.0)	// Outflow from cylinder through valve into pipe
	{
		// Use cylinder gas conditions
		T = TC;
		K1 = (pPpt->gammaAir(T)-1.0)/(2.0*pPpt->gammaAir(T));
		K2 = 2.0/(pPpt->gammaAir(T)-1.0);
		K3 = 1.0/K1;
		S = sqrt(2.0/(pPpt->gammaAir(T)+1.0));
		pi_alpha = lambda_in*pow((PREF/PC), K1);
		C = 1.0/pi_alpha;
		D = pi_alpha;

//		rpipe_flow = INFLOW;
		*(pend_flow[ONE_SIDE]) = INFLOW;
	}
	else
	{
		if(pi_alpha>1.0)
		{
			// Use pipe gas conditions
			T = pBN[ONE_SIDE]->T;
			K1 = (pPpt->gammaAir(T)-1.0)/(2.0*pPpt->gammaAir(T));
			K2 = 2.0/(pPpt->gammaAir(T)-1.0);
			K3 = 1.0/K1;
			S = sqrt(2.0/(pPpt->gammaAir(T)+1.0));
			pi_alpha = lambda_in*pow((PREF/PC), K1);
			C = 1.0/pi_alpha;
			D = pi_alpha;

//			rpipe_flow = OUTFLOW;
			*(pend_flow[ONE_SIDE]) = OUTFLOW;
		// Inflow into the cylinder, OUTFLOW from the pipe, use nozzle subroutine		
		}
		else 
//			rpipe_flow = NOFLOW;
			*(pend_flow[ONE_SIDE]) = NOFLOW;
	}

//	if(rpipe_flow==NOFLOW)
	if(*(pend_flow[ONE_SIDE])==NOFLOW)
	{
//		rlambda_out = lambda_in;
		(*(pCLOUT[ONE_SIDE]))[R+1] = lambda_in;
		rDMDT = 0.0;
		dmdt = rDMDT;
		return;
	}
	
//	if(rpipe_flow==INFLOW)
	if(*(pend_flow[ONE_SIDE])==INFLOW)
	{
		// Line 2
		// Test for sonic/subsonic flow
		double U = S*(sqrt(1.0 + (pow(pPpt->gammaAir(T), 2) - 1.0)*pow(psi, 2)) - 1.0)/(psi*(pPpt->gammaAir(T)-1.0));
		C = sqrt(1.0 - pow(U, 2)/K2)/S;
		double pi_cr = S - U/(K2*C);	// Calculate PI critical		
		C = 1.0/pi_alpha;				// Reset entropy
		double FNPI, B, XPI, DXPI;
		if(pi_alpha<=pi_cr) // Sonic flow
		{
//			cout << "Sonic flow" << endl;
			XPI = 0.5*(pi_alpha + S);
			DXPI = 0.25*(S - pi_alpha);
			double E = psi*pow(2.0/(pPpt->gammaAir(T)+1.0), (pPpt->gammaAir(T)+1.0)/(2.0*(pPpt->gammaAir(T)-1.0)));
			
			do
			{
				B = (XPI - pi_alpha)*C;
				// Equation 6.51, f2(PI). Converged when this = 0.
				FNPI = pow(XPI, K3) - E*((1.0 - K2*pow(B, 2))/(K2*B)); 
				C = sqrt(1.0 - K2*pow(B, 2))/XPI; // Entropy change
				D = 1.0/C;
				if(FNPI>0.0) XPI = XPI - DXPI;
				else XPI = XPI + DXPI;
				DXPI = 0.5*DXPI;
			}while(fabs(FNPI)>=pPpt->POPPET_H_TOL1/*0.00001*/ && DXPI>=pPpt->POPPET_H_TOL2/*0.000001*/);
		}
		else // Subsonic flow
		{
//			cout << "Subsonic flow" << endl;
			XPI = 0.5*(1.0 + S);
			DXPI = 0.25*(1.0 - S);
			
			do
			{
				B = (XPI - pi_alpha)*C;
				// Equation 6.47, f1(PI). Converged when this = 0.
				FNPI = (psi/XPI)*sqrt(K2*(1.0/pow(XPI, 2) - 1.0)) - (B*K2/(1.0 - pow(B, 2)*K2)); 
				C = sqrt(1.0 - K2*pow(B, 2))/XPI; // Entropy change
				D = 1.0/C;
				if(FNPI>0.0) XPI = XPI + DXPI;
				else XPI = XPI - DXPI;
				DXPI = 0.5*DXPI;
			}while(fabs(FNPI)>=pPpt->POPPET_H_TOL1/*0.00001*/ && DXPI>=pPpt->POPPET_H_TOL2/*0.000001*/);
		}			
		lambda_out_temp = 2.0*XPI*pow(PC/PREF, K1) - lambda_in;
//		rlambda_out = lambda_out_temp;
		(*(pCLOUT[ONE_SIDE]))[R+1] = lambda_out_temp;
	}
	else // Outflow from pipe into cylinder through valve, use nozzle subroutine
	{
		// Send psi as the value for phi to the nozzle subroutine
//old		lambda_out_temp = common_HN(pi_alpha, psi, PREF, PC, *(pend_flow[ONE_SIDE]), true, pPpt);
		lambda_out_temp = common_HN_code(pPpt, timestep, time, lambda_in, psi, PC, pBN[ONE_SIDE]->CHOKED, T); // Run nozzle procedure

		// Apply conversion
		lambda_out_temp = lambda_out_temp*pow(PC/PREF, K1);
		lambda_in = pi_alpha*pow(PC/PREF, K1);  // Why convert it if is doesn't get used
//		rlambda_out = lambda_out_temp;
		(*(pCLOUT[ONE_SIDE]))[R+1] = lambda_out_temp;
	}
	
	// Calculate following pipe (p) conditions, in order to calculate mass flow rate 
	double Ap = 0.5*(lambda_in + lambda_out_temp);		// A = a/AREF
	double Up = (lambda_in - lambda_out_temp)/(pPpt->gammaAir(T)-1.0);	// U = u/AREF
	double Pp = pow(Ap, 2.0*pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1.0));		// P = p/PREF
	// For PREF in bar, Fp in m^2, then mass flow rate is:
	rDMDT = pPpt->gammaAir(T)*(PREF/AREF)*Fp*(Up/pow(Ap, 2))*Pp*1.0E5;
	dmdt = rDMDT;
	return;
}

double CValve::Poppet_NH(CProperties* pPpt, double psi, double P0, double T0, double AREF, int cyl_id, bool UPDATE_PIPE, int timestep, double time)
//--------------------------------------------------//
// Non-homentropic end environment					//
// -------------------------------					//
// For flow to/from stagnation conditions P0, T0	//
//--------------------------------------------------//
{
	double lambda_in_n, lambda_out_n, AA_n, PIp, T;
	double *lambda_in_c, *lambda_out_c, *AA_c;
	int *pipe_flow;
	lambda_in_c = new double [NPIPES]; lambda_out_c = new double [NPIPES]; AA_c = new double [NPIPES];	pipe_flow = new int [NPIPES];
	bool* CHOKED; CHOKED = new bool [NPIPES]; 
	CHOKED[ONE_SIDE] = false; // Reset here for all cases

	lambda_in_n = (*(pCLIN[ONE_SIDE]))[R+1];
	lambda_out_n = (*(pCLOUT[ONE_SIDE]))[R+1];
	AA_n = pBN[ONE_SIDE]->AA[R+1];

	// Iterative flow direction test
	T = T0; // Initial guess - assume inflow so use T0 for temperature
	bool FINISHED = false;
	do
	{
		PIp = pow(P0/pPpt->PREF, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T)));
		if(fabs(lambda_in_n/AA_n - PIp) < pPpt->ZERO_TOL || psi==0) // ZERO FLOW
		{
			T = T0; // Arbitrary value for zero flow
			FINISHED = true;
		}
		else
		{
			if(lambda_in_n/AA_n < PIp) // INFLOW - T should be T0	
			{
				if(fabs(T - T0) < pPpt->ZERO_TOL) FINISHED = true;
				else T = T0;
			}
			else // OUTFLOW - T should be pBN[ONE_SIDE]->T
			{
				if(fabs(T - pBN[ONE_SIDE]->T) < pPpt->ZERO_TOL) FINISHED = true;
				else T = pBN[ONE_SIDE]->T;
			}
		}
	}while(!FINISHED);

	pipe_flow_old[ONE_SIDE] = *(pend_flow[ONE_SIDE]);
	
	if(psi <= 0)
	{
		open = false;
	}
	else
	{
		if(psi > 1.0) psi = 1.0;
		open = true;
	}

	//if(fabs(lambda_in_n/AA_n - PIp) < pPpt->ZERO_TOL || !open)
	if( fabs( (lambda_in_n/AA_n)*pow(1/(P0/pPpt->PREF),(pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T))) - 1 ) < pPpt->ZERO_TOL || !open)
	{
		// NOFLOW
		pipe_flow[ONE_SIDE] = NOFLOW;

		lambda_in_c[ONE_SIDE] = lambda_in_n;
		lambda_out_c[ONE_SIDE] = lambda_in_n; // Treat as closed
		AA_c[ONE_SIDE] = AA_n;

		dmdt = 0.0;
	}
	else
	{
		//if(lambda_in_n/AA_n < PIp)	
		if( (lambda_in_n/AA_n)*pow(1/(P0/pPpt->PREF),(pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T))) < 1.0)	
		{
			// INFLOW (into pipe)
			pipe_flow[ONE_SIDE] = INFLOW;

//			cout << GetBoundaryName(NAME[ONE_SIDE]) << "[" << ID << "]" << ": calling CBoundary::common_NHI_code";// << endl;
//			cout << "cyl_id = " << cyl_id << endl;
//			cout << "psi = " << psi << endl;
//			cout << endl;

			common_NHI_code(pPpt, lambda_in_n, lambda_out_n, AA_n, 
						lambda_in_c[ONE_SIDE], lambda_out_c[ONE_SIDE], AA_c[ONE_SIDE], psi, P0, T0, 
						CHOKED[ONE_SIDE], SONIC, T, timestep, time, true/*true=constant pressure (valve) model or false=pressure loss (port) model*/, 1e-6/*inflow main loop tolerance*/);

//			cout << " - DONE" << endl;
		}
		else 
		{

			// OUTFLOW (into cylinder)
			pipe_flow[ONE_SIDE] = OUTFLOW;

//			cout << GetBoundaryName(NAME[ONE_SIDE]) << "[" << ID << "]" << ": calling CBoundary::common_NHN_code";// << endl;
			//cout << GetBoundaryName(NAME[ONE_SIDE]) << "[" << ID << "]" << ": calling CBoundary::NHNozzle";// << endl;
			
			// Need new values for lambda_out_c only
			lambda_in_c[ONE_SIDE] = lambda_in_n;
			AA_c[ONE_SIDE] = AA_n;

			// Run nozzle procedure; for open ends set phi = 1
			lambda_out_c[ONE_SIDE] = common_NHN_code(pPpt, lambda_in_n, AA_n, psi, P0, CHOKED[ONE_SIDE], T, time); 
			// Use newer version (better, N-R convergence):
			//lambda_out_c[ONE_SIDE] = NHNozzle(pPpt, lambda_in_n, lambda_out_n, AA_n, psi, P0, CHOKED[ONE_SIDE]);

//			cout << " - DONE" << endl;

//			cout << "lambda_in_c[ONE_SIDE] = " << lambda_in_c[ONE_SIDE] << endl;
//			cout << "lambda_out_c[ONE_SIDE] = " << lambda_out_c[ONE_SIDE] << endl;
//			cout << endl;
		}

		// For PREF in bar, FP in m^2, then mass flow rate is (Eq. 7.169 of Benson):
		dmdt = (4*pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1))*(pPpt->PREF*(pBN[ONE_SIDE]->f_dash*pPpt->fref)/AREF)*(((*(pCLIN[ONE_SIDE]))[R+1] - (*(pCLOUT[ONE_SIDE]))[R+1])/pow((*(pCLIN[ONE_SIDE]))[R+1] + (*(pCLOUT[ONE_SIDE]))[R+1],2))
						*pow(((*(pCLIN[ONE_SIDE]))[R+1] + (*(pCLOUT[ONE_SIDE]))[R+1])/(2*pBN[ONE_SIDE]->AA[R+1]),(2*pPpt->gammaAir(T))/(pPpt->gammaAir(T)-1))*1.0E5;
	}
	if(UPDATE_PIPE) common_UPDATE_NH(pPpt, lambda_in_c, lambda_out_c, AA_c, pipe_flow, CHOKED);
	delete [] lambda_in_c;
	delete [] lambda_out_c;
	delete [] AA_c;
	delete [] pipe_flow;
	delete [] CHOKED;
	return dmdt;
}

void CValve::Equations
(CProperties* pPpt, double &lambda_in_c, double &lambda_out_c, double &AA_c, double &AA_calc, 
 double lambda_in_n, double AA_n, double Ac, double &U, double &C, double &PpbyPc, 
 double psi, double rc, bool &rCHOKED, bool &ERROR, double T)
//--------------------------------------------------//
// Equations for non-homentropic poppet valve		//
// ------------------------------------------		//
// Set of equations called multiple time from		//
// non-homentropic poper valve function.			//
//--------------------------------------------------//
{	
	lambda_in_c = 2*lambda_in_n*(AA_c/(AA_c + AA_n)) + lambda_out_c*((AA_c - AA_n)/(AA_c + AA_n));
//	lambda_in_c = ((lambda_in_c + lambda_out_c)/2)*(1 - (AA_n/AA_c)) + lambda_in_n;
	if(lambda_in_c < 0)
	{
		ERROR = true;
		return;
	}

	// (6) ======================================================
	double temp = (pow(pPpt->gammaAir(T),2)-1)*pow(Ac,2) + 2*(1-pPpt->gammaAir(T))*pow(lambda_in_c,2);
	if(temp<0)
	{
		ERROR = true;//cout << "temp < 0!" << endl;
//		return;
	}
	else ERROR = false;
	lambda_out_c = ((3-pPpt->gammaAir(T))/(pPpt->gammaAir(T)+1))*lambda_in_c + (2/(pPpt->gammaAir(T)+1))*sqrt(temp);

	// (7) ======================================================
//cout << "Ac = " << Ac << endl;
	U = (lambda_out_c - lambda_in_c)/((pPpt->gammaAir(T)-1)*Ac);
//	U = (lambda_out_c - lambda_in_c)/((pPpt->gammaAir(T)-1)*1);

	// (8) ======================================================
	C = (((pPpt->gammaAir(T)-1)/2)*pow(U,2))/pow((1 - ((pPpt->gammaAir(T)-1)/2)*pow(U,2)),2);

	// (9, 10) ======================================================
	if((4*C)/(pow(pPpt->gammaAir(T),2)-1) - pow(psi,2)>=0)	// Sonic flow in the VALVE, VALVE is choked
	{
		SONIC = true;
	
		// If there is choked flow in the pipe u_p = a_p, or U_p = sqrt(2/(pPpt->gammaAir(T) + 1));
		// This makes PpbyPc simplify to:
		// PpbyPc = psi * pow(2/(pPpt->gammaAir(T)+1), pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1)); // eq. 7.148
		// But can set U_pipe = sqrt(2/(pPpt->gammaAir(T) + 1)), and the above eqn will suffice.

		if(U >= sqrt(2/(pPpt->gammaAir(T) + 1))) // Choked flow in PIPE
		{
			rCHOKED = true;
			U = sqrt(2/(pPpt->gammaAir(T) + 1)); // The limit on U in the pipe now that PIPE is choked
		}
		
		PpbyPc = psi * pow(2/(pPpt->gammaAir(T)+1), (pPpt->gammaAir(T)+1)/(2*(pPpt->gammaAir(T)-1))) * (1 - ((pPpt->gammaAir(T)-1)/2)*pow(U,2)) * (1/U); // eq. 7.142
	}
	else // Subsonic flow in the valve
	{

		if(U >= sqrt(2/(pPpt->gammaAir(T) + 1))) // Choked flow in PIPE
		{
			rCHOKED = true;
			U = sqrt(2/(pPpt->gammaAir(T) + 1)); // The limit on U in the pipe now that PIPE is choked
		}
		C = (((pPpt->gammaAir(T)-1)/2)*pow(U,2))/pow((1 - ((pPpt->gammaAir(T)-1)/2)*pow(U,2)),2);

		SONIC = false;
		PpbyPc = pow((1/(2*C))*( psi*sqrt(pow(psi,2) + 4*C) - pow(psi,2) ), pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1));
	}



/*	
//	Not needed - intrinsic to method. Have checked with graph comparisons - makes no difference
	// Test for choked flow
//	double Choked = ((pPpt->gammaAir(T)+1)/(3-pPpt->gammaAir(T)))*lambda_in_c;
	double Choked = ((pPpt->gammaAir(T)+1)/(3-pPpt->gammaAir(T)))*lambda_in_n;
//	cout << "Choked = " << Choked << endl;
	if(lambda_out_c>Choked)
//	PpbyPc = psi * pow(2/(pPpt->gammaAir(T)+1), pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1)); // eq. 7.148
	{lambda_out_c = Choked; cout << "Invoking choked flow limitation on valve[" << ID << "]\n";}
*/




	// (11) ======================================================
	AA_calc = ((lambda_in_c + lambda_out_c)/2)*pow(1/PpbyPc, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T)))*pow(1/rc, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T)));
	if(AA_calc < 0)
	{
		ERROR = true;
		return;
	}
}
			

void CValve::NHP(CProperties* pPpt, double &rlambda_in, double &rlambda_out, double lambda_out_old, 
 double &rAA, CPathLine &rPathLine, double &rDMDT, double psi, double K, double AREF, 
 double PREF, double PC, double TC, double AC, double Fp, double ca, 
 int &rpipe_flow, int end, double XPIPE, int timestep, double time)
 //--------------------------------------------------//
// Non-homentropic port								//
// --------------------								//
// Inflow algorithm based on Benson's constant		//
// pressure port model.								//
// Outflow will revert to a non-homentropic nozzle	//
//--------------------------------------------------//
{
	double lambda_in_n, lambda_out_n, AA_n, PIp, T;
	double *lambda_in_c, *lambda_out_c, *AA_c;
	int *pipe_flow;
	lambda_in_c = new double [NPIPES]; lambda_out_c = new double [NPIPES]; AA_c = new double [NPIPES];	pipe_flow = new int [NPIPES];
	bool* CHOKED; CHOKED = new bool [NPIPES]; 
	CHOKED[ONE_SIDE] = false; // Reset here for all cases

	lambda_in_n = (*(pCLIN[ONE_SIDE]))[R+1];
	lambda_out_n = (*(pCLOUT[ONE_SIDE]))[R+1];
	AA_n = pBN[ONE_SIDE]->AA[R+1];

//	double lambda_in_n = rlambda_in;		// Store the uncorrected value
//	double AA_n = rAA;						// Store the uncorrected value

	// Iterative flow direction test
	T = TC; // Initial guess - assume inflow so use TC for temperature
	bool FINISHED = false;
	do
	{
		PIp = pow(PC/pPpt->PREF, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T)));
		if(fabs(lambda_in_n/AA_n - PIp) < pPpt->ZERO_TOL || psi==0) // ZERO FLOW
		{
			T = TC; // Arbitrary value for zero flow
			FINISHED = true;
		}
		else
		{
			if(lambda_in_n/AA_n < PIp) // INFLOW - T should be TC	
			{
				if(fabs(T - TC) < pPpt->ZERO_TOL) FINISHED = true;
				else T = TC;
			}
			else // OUTFLOW - T should be pBN[ONE_SIDE]->T
			{
				if(fabs(T - pBN[ONE_SIDE]->T) < pPpt->ZERO_TOL) FINISHED = true;
				else T = pBN[ONE_SIDE]->T;
			}
		}
	}while(!FINISHED);

	// Test if valve open (psi > 0)
	if(psi<=0.0) // Valve closed - treat as closed end
	{
		open = false;
		rlambda_out = lambda_in_n; // No update needed to rlambda_in or rAA
		rDMDT = 0.0;
		return;
	}
	else open = true;

	// If valve area is greater than pipe area, set psi=1.0
	if(psi>1.0) psi=1.0;

	// Set up constants
	double K1 = (pPpt->gammaAir(T)-1.0)/(2.0*pPpt->gammaAir(T));
	double K3 = 1.0/K1;

	// Establish flow direction due to cylinder
	double lambda_in_star = (lambda_in_n/AA_n)*pow(PREF/PC, K1);

	int rpipe_flow_old = rpipe_flow; // Store old value

	if(lambda_in_star<1.0) rpipe_flow = INFLOW; 
	// Outflow from the cylinder, INFLOW into the pipe; entropy adjustment required
	// Change the flow direction indicator at the appropriate end of the relevant pipe
	// by reference
	else
	{ 
		if(lambda_in_star>1.0) 	rpipe_flow = OUTFLOW; 
		// Inflow into the cylinder, OUTFLOW from the pipe; no entropy adjustment		
		else rpipe_flow = NOFLOW;
	}
/*
	// Check whether direction agrees with pipe pathlines
	if(rpipe_flow_old==rpipe_flow) ;
		//cout << "Directions match" << endl;
	else 
	{
		cout << "Flow direction mismatch!" << endl;
		cout << "rpipe_flow before getting altered = " << rpipe_flow_old << endl;
		cout << "rpipe_flow = " << rpipe_flow << endl;
		cout << "lambda_in_n = " << lambda_in_n << endl;
		cout << "AA_n = " << AA_n << endl;
		cout << "pow(PREF/PC, K1) = " << pow(PREF/PC, K1) << endl;
		cout << "lambda_in_star = " << lambda_in_star << endl;
		char pause;
//		cin >> pause;//exit(1);
	}
*/	
	if(rpipe_flow==NOFLOW)
	{
		rlambda_out = lambda_in_n;
		// No update needed to rlambda_in or rAA
		rDMDT = 0.0;
		return;
	}
	
	if(rpipe_flow==INFLOW)	
	// INFLOW (from cylinder, through valve, to pipe)
	// If rpipe_flow has just been changed above, this still updates the first/last
	// pathlines, whether specially created or not
	{
		// Benson Sec. 7.5. Outflow from a cylinder to a pipe through a valve: constant
		// pressure model.
		// Known variables: lambda_in_n, AA_n, PC, AC, psi
		// Required variables: lambda_in_c, lambda_out_c, AA_c, Pp (pipe pressure)
		// NB AA_c here is AA_est in book
		bool converged, converged2, less_than, switched;
		int counter=0;
		int counter2;
		double U, RR, PpbyPc, AA_calc;

		// (2) -------------------------------------------------------
		double Ac = AC/AREF;// AC is not dimensionless!
		double AA_is = Ac*pow(PREF/PC, K1);	
		cout << "Ac = " << Ac << endl;
		cout << "AA_is = " << AA_is << endl;
		double lambda_in_c = lambda_in_n;			// Use uncorrected as initial estimate
		double AA_c;	// Corrected entropy level in the pipe immediately after the
						// gas enters the pipe from the valve
						// AA_n (passed in) is the entropy level in the pipe immediately
						// before the gas enters the pipe from the valve
		double del_AA;	// Entropy level adjustment to enable conergence
		double lambda_out_c = lambda_out_old;	// Use previous result fot initial estimate

		// (3) -------------------------------------------------------
		if((PC/PREF)<1.0) // Replace later with rc??
		{AA_c = AA_is; del_AA = 0.5;}
		else
		{AA_c = AA_is; 
		del_AA = (Ac - AA_is)/2;}
	
		// Since the initial value of AA_c is isentropic, AA_c will always be
		// greater than AA_calc to start with, so...
		less_than = false;
		switched = false;
//if(ca>224) cout << "here" << endl;
		converged = false;
		do
		{
			++counter;
			counter2=0;
			if(counter>150) cout << "counter = " << counter << endl;
		// (4) -------------------------------------------------------
			converged2 = false;
			do
			{
				++counter2;
				lambda_in_c = lambda_in_n + ((lambda_in_c + lambda_out_c)/2)*(1 - AA_n/AA_c);
		// (5) -------------------------------------------------------
				if((Ac - lambda_in_c)<0)
				{
					AA_c = AA_is;
					del_AA = del_AA/2;
					lambda_in_c = lambda_in_n + ((lambda_in_c + lambda_out_c)/2)*(1 - AA_n/AA_c); 
					// Must repeat calculation
					converged2 = true;	// do loop once only
				}
				else converged2 = true;
			}while(!converged2);
//if(ca > 136) cout << "here6aa" << endl;
		// (6) -------------------------------------------------------
			lambda_out_c = ((3-pPpt->gammaAir(T))/(pPpt->gammaAir(T)+1))*lambda_in_c 
						   + (2/(pPpt->gammaAir(T)+1))*sqrt( 
								(pow(pPpt->gammaAir(T),2)-1)*pow(Ac,2) + 2*(1-pPpt->gammaAir(T))*pow(lambda_in_c,2) );
		// (7) -------------------------------------------------------
			U = (lambda_out_c - lambda_in_c)/((pPpt->gammaAir(T)-1)*Ac);
		// (8) -------------------------------------------------------
			RR = pow(psi,2)*pow((1 - ((pPpt->gammaAir(T)-1)/2)*pow(U,2)), (pPpt->gammaAir(T)-3)/(pPpt->gammaAir(T)-1))/((pPpt->gammaAir(T)-1)*pow(U,2));
		// (9, 10) ---------------------------------------------------
			// Subsonic/sonic test
			if( (2/(pPpt->gammaAir(T)+1))*pow((1 - ((pPpt->gammaAir(T)-1)/2)*pow(U,2)), 2/(pPpt->gammaAir(T)-1))*pow(U,2) - pow(psi,2) >0 )	
			// Sonic flow in the ports
			{
				if(U<=sqrt(2.0/(pPpt->gammaAir(T)+1.0)))	// Choked flow in the pipe
				{
//					cout << "Choked flow in pipe" << endl;
		//			U = sqrt(2.0/(pPpt->gammaAir(T)+1.0));	// Do this and 7.142 simplifies correctly
		//			PpbyPc = (sqrt((4*C)/(pow(pPpt->gammaAir(T),2)-1)))*pow(2/(pPpt->gammaAir(T)-1), (pPpt->gammaAir(T)+1)/(2*(pPpt->gammaAir(T)-1)))*(1 - ((pPpt->gammaAir(T)-1)/2)*pow(U,2))*(1/U); // eq. 7.142		
		//			PpbyPc = psi*pow(2/(pPpt->gammaAir(T)-1), (pPpt->gammaAir(T)+1)/(2*(pPpt->gammaAir(T)-1)))*(1 - ((pPpt->gammaAir(T)-1)/2)*pow(U,2))*(1/U); // eq. 7.142
					PpbyPc = psi*pow(2/(pPpt->gammaAir(T)+1), pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1));
//					choked_pipe=true;
				}
				else
//*/
				{
//					cout << "Choked flow in valve" << endl;
		//			U = sqrt(2/(pPpt->gammaAir(T)+1));
		//			psi = sqrt((4*C)/(pow(pPpt->gammaAir(T),2)-1));
		//			PpbyPc = (sqrt((4*C)/(pow(pPpt->gammaAir(T),2)-1)))/*psi*/*pow(2/(pPpt->gammaAir(T)-1), (pPpt->gammaAir(T)+1)/(2*(pPpt->gammaAir(T)-1)))*(1 - ((pPpt->gammaAir(T)-1)/2)*pow(U,2))*(1/U); // eq. 7.142
					PpbyPc = psi*pow(2/(pPpt->gammaAir(T)-1), (pPpt->gammaAir(T)+1)/(2*(pPpt->gammaAir(T)-1)))*(1 - ((pPpt->gammaAir(T)-1)/2)*pow(U,2))*(1/U); // eq. 7.142
		//			PpbyPc = (sqrt((4*C)/(pow(pPpt->gammaAir(T),2)-1)))*pow(2/(pPpt->gammaAir(T)+1), pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1));
		//			PpbyPc = pow(2/(pPpt->gammaAir(T)+1), pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1));
				}
//				choked_valve=true;
			}
			else	// Subsonic
			{
				PpbyPc = pow((R*( sqrt(1 + (2/RR)*(1 - ((pPpt->gammaAir(T)-1)/2)*pow(U,2))) - 1 )), pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1));
//				choked_valve=false;
//				choked_pipe=false;
			}
		// (11) ------------------------------------------------------
			AA_calc = ((lambda_in_c + lambda_out_c)/2)
						*pow(1/PpbyPc, K1)
						*pow(PREF/PC, K1);		// Replace later with 1/rc??
		// (12) ------------------------------------------------------
			if(fabs(AA_calc - AA_c)<0.0001)
				converged = true;
				// lambda_in = lambda_in_c, lambda_out = lambda_out_c, AA = AA_c
			else
			{
				if(AA_c<AA_calc)								//if((AA_calc - AA_c)>0)
				{	
					// My convergence method will not affect the converged value
					// it only helps to speed up convergence
					if(!less_than) switched=true; else switched=false;
					less_than = true;

					if(switched)
					{AA_c = AA_c + del_AA; del_AA = del_AA/2;}				
					// Halve the step size when we shoot too far and pass AA_calc
					else AA_c = AA_c + del_AA;	
					// If AA_c < AA_calc last time as well, continue with the same step size	
				}
				else // AA_c > AA_calc
				{
					if(less_than) switched=true; else switched=false;
					less_than = false;

					if(switched)
					{AA_c = AA_c - del_AA; del_AA = del_AA/2;}				
					// Halve the step size when we shoot too far and pass AA_calc				
					else AA_c = AA_c - del_AA;	
					// If AA_c > AA_calc last time as well, continue with the same step size
				}
			}
/*			
			cout << "AA_calc = " << AA_calc << endl;
			cout << "AA_c = " << AA_c << endl;
			cout << "fabs(AA_calc - AA_c) = " << fabs(AA_calc - AA_c) << endl;
			cout << "del_AA = " << del_AA << endl;
*/
///*
				if(ca>136)
				{
					char pause;
					cout << "counter = " << counter << endl;
					cout << "(pow(pPpt->gammaAir(T),2)-1)*pow(Ac,2) + 2*(1-pPpt->gammaAir(T))*pow(lambda_in_c,2) = " <<
						(pow(pPpt->gammaAir(T),2)-1)*pow(Ac,2) + 2*(1-pPpt->gammaAir(T))*pow(lambda_in_c,2) << endl;
					cout << "Ac = " << Ac << endl;
					cout << "U = " << U << endl;
					cout << "RR = " << RR << endl;
				//	cout << "(4*C)/(pow(pPpt->gammaAir(T),2)-1) - pow(psi,2) = " << (4*C)/(pow(pPpt->gammaAir(T),2)-1) - pow(psi,2) << endl;
//					if(choked_valve) cout << "Choked valve" << endl;
//					else cout << "Subsonic flow" << endl;
				//	cout << "K = " << K << endl;
				//	cout << "K1 = " << K1 << endl;
					cout << "psi = " << psi << endl;
//					cout << "psi_cr = " << sqrt((4*C)/(pow(pPpt->gammaAir(T),2)-1)) << endl;
				//	cout << "lambda_in_n = " << lambda_in_n << endl;
				//	cout << "AA_n = " << AA_n << endl;
				//	cout << "AA_is = " << AA_is << endl;
					cout << "AA_c = " << AA_c << endl;
					cout << "AA_calc = " << AA_calc << endl;
					cout << "lambda_in_c = " << lambda_in_c << endl;
					cout << "lambda_out_c = " << lambda_out_c << endl;
					cout << "PpbyPc = " << PpbyPc << endl;
					cout << "Pp = " << PpbyPc * PC << endl;
					cout << "del_AA = " << del_AA << endl;
					cout << "(Ac - lambda_in_c) = " << (Ac - lambda_in_c) << endl;
					cin >> pause;
				}
//*/
		}while(!converged && counter<10000);
//		char pause;
//		cout << "counter = " << counter << endl;
//		cout << "paused..." << endl; 
//		cin >> pause;
//if(ca>224) cout << "here fin" << endl;
		// Must also save the corrected lambda_in, and entropy level in the pipe
		// as well as returning the usual, now corrected, lambda_out
		rAA = AA_c;						// Correct pipe referenced AA; NEEDED????
		
		if((pPipe[ONE_SIDE]->METHOD==pPipe[ONE_SIDE]->MMOC && !pPpt->HOMENTROPIC) || 
		   (pPipe[ONE_SIDE]->METHOD==pPipe[ONE_SIDE]->W_ALPHA_BETA && pPpt->COMBINED_WAB_MOC))
			rPathLine.AA = AA_c;			// Update path line referenced AA
		
		if(rpipe_flow_old != rpipe_flow) 
		// If old flow direction was not INFLOW then there was no new path line created
		// so adjust existing one by setting its XK value to the appropriate end
		{
			if((pPipe[ONE_SIDE]->METHOD==pPipe[ONE_SIDE]->MMOC && !pPpt->HOMENTROPIC) || 
			   (pPipe[ONE_SIDE]->METHOD==pPipe[ONE_SIDE]->W_ALPHA_BETA && pPpt->COMBINED_WAB_MOC))
			{
				if(end==ODD){rPathLine.XK = 0; rPathLine.XK_old = rPathLine.XK;}
				else{rPathLine.XK = XPIPE; rPathLine.XK_old = rPathLine.XK;}
			}
		}
		else 
		// Check that a new path line HAS been created correctly
		// This is INFLOW only anyway, no need to allow for NOFLOW
		{
			if((pPipe[ONE_SIDE]->METHOD==pPipe[ONE_SIDE]->MMOC && !pPpt->HOMENTROPIC) || 
			   (pPipe[ONE_SIDE]->METHOD==pPipe[ONE_SIDE]->W_ALPHA_BETA && pPpt->COMBINED_WAB_MOC))
			{
				if((end==ODD && rPathLine.XK!=0) || (end==EVEN && rPathLine.XK!=XPIPE))
					cout << "Valve: should've created new path line but its XK is not right!!" << endl;
			}
		}
		
		rlambda_in = lambda_in_c;//*pow(PC/PREF, K1);		// Correct pipe referenced lambda_in
		rlambda_out = lambda_out_c;//*pow(PREF/PC, K1);		// Return lambda_out by reference
	}
	else 
	// OUTFLOW to cylinder through valve from pipe, use nozzle subroutine
	// Should therefore only use OUTFLOW section of nozzle subroutine
	// NOFLOW taken care of above
	{
		// Same boundary conditions are used for this case as for homentropic
		// flow i.e., isentropic flow is assumed from the pipe pressure to the
		// cylinder pressure. Then the partially open end boundary (nozzle) can
		// applied with the back pressure set equal to the cylinder pressure
		// and psi = phi;
/*		
		double PIA = lambda_in_n*pow((PREF/PC), K1);	// I dunno why we do this
		
		// Send PIA as lambda_in, and psi as the value for PHI to the NOZZLE subroutine, PC as the value of Pb
		double lambda_out_temp = 
			common_NHN(pPpt, PIA, rAA, rPathLine, psi, PREF, PC, rpipe_flow, end, XPIPE, true, pPipe[ONE_SIDE]->EX, ID, false, timestep, time);
		
		lambda_out_temp = lambda_out_temp*pow(PC/PREF, K1); // Apply conversion
		rlambda_in = PIA*pow(PC/PREF, K1);
		rlambda_out = lambda_out_temp;
		// No update needed to rAA, or rAAK, these remain as they were
*/
		// OUTFLOW
		pipe_flow[ONE_SIDE] = OUTFLOW;
			
		// Need new values for lambda_out_c only
		lambda_in_c[ONE_SIDE] = lambda_in_n;
		AA_c[ONE_SIDE] = AA_n;

		// Run nozzle procedure; phi = psi
		lambda_out_c[ONE_SIDE] = common_NHN_code(pPpt, lambda_in_n, AA_n, psi, PC, CHOKED[ONE_SIDE], T, time);

		common_UPDATE_NH(pPpt, lambda_in_c, lambda_out_c, AA_c, pipe_flow, CHOKED);
	}


	// For PREF in bar, Fp in m^2, then mass flow rate is:
	// Eq. 7.169 of Benson
	rDMDT = (4*pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1))*(PREF*Fp/AREF)*((rlambda_in - rlambda_out)/pow(rlambda_in + rlambda_out,2))
						*pow((rlambda_in + rlambda_out)/(2*rAA),K3)*1.0E5; // *1.0E5 correct??
	return;
}

void CValve::PrintToScreen(CProperties* pPpt)
{
	if(EX) cout << Underline("Exhaust Valve", "-", ID, "\t\t");
	else cout << Underline("Intake Valve", "-", ID, "\t\t");

	if(open) {
		cout << "\t\tValve open:\n";
//		cout << "\t\tLift (final_lift)\t\t=\t" << final_lift << " mm\n";
//		cout << "\t\tEff. Area (eff_area)\t\t=\t" << eff_area << " m^2\n";
//		if(!pPpt->HOMENTROPIC && SONIC) cout << "\t\tCHOKED flow in valve\n";
	}
	else cout << "\t\tValve closed\n";

	if(L_OR_A) {
		cout << "\t\tLift, final_lift\t\t=\t" << final_lift << " mm" << endl;
		cout << "\t\tDischarge coefficient, valve_Cf_or_Cd\t=\t" << valve_Cf_or_Cd << endl;
		cout << "\t\tEffective area, eff_area\t=\t" << eff_area << " m^2" << endl;
		cout << "\t\tArea ratio, eff_area/FP\t\t=\t" << eff_area/FP << endl;
	}
	else {
		cout << "\t\tEffective area, eff_area\t=\t" << eff_area << " m^2" << endl;
		cout << "\t\tArea ratio, eff_area/FP\t\t=\t" << eff_area/FP << endl;
	}

	if(!pPpt->HOMENTROPIC && SONIC) cout << "\t\tCHOKED flow in valve\n";
}

void CValve::PrintToScreenVolume(CProperties* pPpt)
{
	cout << Underline("Valve", "-", ID, "\t");

	if(open)
	{
		cout << "\tValve open:\n";
	}
	else cout << "\tValve closed\n";

	cout << "\tEffective area (eff_area)\t\t\t=\t" << eff_area << " m^2\n";
cout << "\tArea ratio (psi_val)\t\t\t\t=\t" << psi_val << "\n";

cout << endl;
}

void CValve::PrintToFile(CProperties* pPpt, int timestep, double time, double ca_elapsed, double ca)
// ============================================================ //
// Prints instantaneous pipe object data to file.				//
// This function is called from the main function.				//
// ============================================================ //
{
	//	if(ntappings>0)
	//	{
	int f;
	f = 0;

	if (timestep == 0)
	{
		std::string temp_str = "RESULTS FILE FOR ";
		if (EX) temp_str += "EXHAUST"; else temp_str += "INTAKE";
		temp_str += " VALVE [";
		temp_str += IntToString(ID);
		temp_str += "], CYLINDER [";
		temp_str += IntToString(CYL_ID);
		temp_str += "]";
		//std::cout << temp_str << std::endl;

		fprintf(OUTPUT_FILE, "%s\n", Underline(StringToChar(temp_str), "-"));

		fprintf(OUTPUT_FILE, "%s", "Time(s)");
		if (!pPpt->CONTINUOUS) // Periodic
			fprintf(OUTPUT_FILE, "\t%s%c%s\t%s%c%s", "Elapsed (", Deg(), "CA)", "Crank angle (", Deg(), "CA)");

		fprintf(OUTPUT_FILE, "\t%s\t%s\t%s\t%s", "Lift (mm)", "Discharge Coeff. ()", "Effective Area (mm^2)", "Mass Flow Rate (kg/s)");

		if (ACTIVE) {
			fprintf(OUTPUT_FILE, "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
				"Target value (units)", "Recorded value at target location (units)", "Error (units)", "Gradient (units/s)", 
				"Target cycle average (units)", "Recorded cycle average ss(units)", "Cycle average error (units)", "Cycle average error gradient (units/s)");
		}

		fprintf(OUTPUT_FILE, "\n");
	}
	/*
			if(time>=print_from_time)// Only start recording data as specified
			{
	*/
	if (timestep % pPpt->freq == 0) // Print data at the specified sampling frequency
	{
		fprintf(OUTPUT_FILE, "%f", time);

		if (!pPpt->CONTINUOUS) // Periodic
			fprintf(OUTPUT_FILE, "\t%f\t%f", ca_elapsed, ca);

		fprintf(OUTPUT_FILE, "\t%f\t%f\t%f\t%f", final_lift, valve_Cf_or_Cd, eff_area * 1e6/*convert to mm^2*/, dmdt);

		if (ACTIVE) {
			//fprintf(OUTPUT_FILE, "\t%s\t%s\t%s\t%s", "Target value (units)", "Recorded value at target location (units)", "Error (units)", "Gradient (units/s)");
			fprintf(OUTPUT_FILE, "\t%f\t%f\t%f\t%f", AVTargetValue, AVTCurrentValue, AVTError, AVTErrorGrad);
			fprintf(OUTPUT_FILE, "\t%f\t%f\t%f\t%f", target_cyc_av, recorded_cyc_av, cyc_av_error, cyc_av_error_grad);
		}

		fprintf(OUTPUT_FILE, "\n");

	}
	/*
			}
	*/
	//	}
}

void CValve::RecordAVT(CProperties* pPpt, int timestep, double time_or_CADelapsed, double time_step_length) //, CPipe* Pipe)
{
	int pause;
	if (pEng->NEW_CYCLE){

		if (!FIRST_CYCLE) {

			// Find cycle average, maximum value, and location of maximum value in the recorded waveform			
			recorded_cyc_av = AVTCumulativeVal / AVTInCycleStepCount;
			recorded_max_row = 0;
			recorded_max_ang = 0;
			recorded_max_val = 0;
			for (int row = 0; row < numAVTDataPoints; row++) {
				if (avt_raw[row][RECORDED] > recorded_max_val) {
					recorded_max_row = row;
					recorded_max_ang = avt_raw[row][CAD];
					recorded_max_val = avt_raw[row][RECORDED];
				}
			}

			phase_lag_between_peaks = recorded_max_ang - target_max_ang; // Positive values indicate recorded peak location is lagging target peak location
			
			cyc_av_error_prev = cyc_av_error; // Record previous cycle value
			cyc_av_error = recorded_cyc_av - target_cyc_av; // [bar]
			
			if (fabs(plenum_press - plenum_press_prev) > 1e-6) {
				cyc_av_error_grad = (cyc_av_error - cyc_av_error_prev) / (plenum_press - plenum_press_prev); // [(value - target) bar / plenum bar] // / time_step_length;
			}
			else cyc_av_error_grad = 1;// 0;

			if (cyc_av_error / cyc_av_error_prev < 0) CYC_AV_ERR_SIGN_CHANGE = true;
			 
			if (fabs((cyc_av_error / target_cyc_av) * 100) < pEng->tol || PHASE_MATCHING) CYC_AV_MATCHED = true;
			else CYC_AV_MATCHED = false;

cout << endl;
cout << endl;
cout << "START OF NEW CYCLE\n";
cout << endl;
cout << "Target max lies at avt_raw[row=" << target_max_row << "][CAD=" << avt_raw[target_max_row][CAD] << "] deg with value = " << avt_raw[target_max_row][TARGET] << " bar" << endl;
cout << "target_max_ang\t\t=\t" << target_max_ang << " deg" << endl;
cout << "target_max_val\t\t=\t" << target_max_val << " bar" << endl;
cout << endl;
cout << "Recorded max lies at avt_raw[row=" << recorded_max_row << "][CAD=" << avt_raw[recorded_max_row][CAD] << "] deg with value = " << avt_raw[recorded_max_row][RECORDED] << " bar" << endl;
cout << "recorded_max_ang\t=\t" << recorded_max_ang << " deg" << endl;
cout << "recorded_max_val\t=\t" << recorded_max_val << " bar" << endl;
cout << endl;
cout << "phase_lag_between_peaks\t=\t" << phase_lag_between_peaks << " deg" << endl;
cout << endl;
cout << endl;
cout << "target_cyc_av\t\t=\t" << target_cyc_av << " bar" << endl;
cout << "recorded_cyc_av\t\t=\t" << recorded_cyc_av << " bar" << endl;
cout << "cyc_av_error\t\t=\t" << cyc_av_error << " bar" << endl;
cout << "cyc_av_error_grad\t=\t" << cyc_av_error_grad << endl;
			
			if (CYC_AV_ERR_SIGN_CHANGE) cout << "CYC_AV_ERR_SIGN_CHANGE\t=\t" << TrueOrFalse(CYC_AV_ERR_SIGN_CHANGE) << endl;

			if (CYC_AV_ERR_SIGN_CHANGE && !CYC_AV_MATCHED) {
				cyc_av_factor /= 2;
				//AVTFactorSF /= 2;
				CYC_AV_ERR_SIGN_CHANGE = false;
			}

			if (CYC_AV_MATCHED || PHASE_MATCHING) {

				if (CYC_AV_MATCHED) {

cout << "CYCLE AVERAGE MATCH ACHIEVED " << endl;
cout << "at a plenum pressure\t=\t" << pEng->pres << " bar" << endl;
cout << endl;
				// Reset factors
//				cyc_av_factor = cyc_av_factor_orig; // Don't reset!?
				//AVTFactorSF = AVTFactorSFOrig;
				}

				// Next do phase matching
				PHASE_MATCHING = true;
				phase_shift -= 0.5*phase_lag_between_peaks;
				
				if (fabs(phase_lag_between_peaks) < 1.0) { // If peaks align within 1 degree CA

					PHASE_MATCHED = true;
					PHASE_MATCHING = false; 
					CYC_AV_MATCHED = false; // Go back to cycle average matching

					cout << "PHASE MATCHED!" << endl;
					cout << endl;
				}
				else {
					cout << "PHASE_MATCHING = true" << endl;
					cout << "phase_lag_between_peaks\t=\t" << phase_lag_between_peaks << " deg" << endl;
					cout << "phase_shift\t\t=\t" << phase_shift << " deg" << endl;
				}	
			}
			else{

cout << "Old plenum pressure\t=\t" << pEng->pres << " bar" << endl;
//cout << "Old pEng->scale_factor_ev = \t" << pEng->scale_factor_ev << endl;
				
				// Adjust plenum pressure in order to move value_cyc_av towards target_cyc_av
				//pEng->pres += -1 * cyc_av_factor * cyc_av_error;
				plenum_press_prev = plenum_press; // Store previous value
				if (fabs(cyc_av_error_grad) > 1e-6)
				{
					// Don't allow the new plenum pressure to go negative!
					double test_val;
					do {
						test_val = plenum_press + -1 * cyc_av_factor * (cyc_av_error * (1 / cyc_av_error_grad));
						if (test_val < 0) cyc_av_factor *= 0.5; // Half the movement	
					} 					
					while(test_val <= 0);
					plenum_press = test_val;
				}
				else plenum_press += 0;

				pEng->pres = plenum_press;

cout << "cyc_av_factor\t\t=\t" << cyc_av_factor << endl;
cout << "Plenum pressure rise\t=\t" << plenum_press - plenum_press_prev << " bar" << endl;
cout << "New plenum pressure\t=\t" << pEng->pres << " bar" << endl;

				//pEng->scale_factor_ev += -1 * AVTFactorSF * cyc_av_error;
				//cout << "New pEng->scale_factor_ev = " << pEng->scale_factor_ev << endl;
			}


			// Reset counters for the next cycle
			AVTCumulativeVal = 0;
			AVTInCycleStepCount = 0;
		}
		else {
			FIRST_CYCLE = false; // Can turn this off now
		}
	}
	
	// Need to update the interpolatation between the two nodes either side of the target tapping upon every timestep
	if (!SINGLE_NODED_PIPE) {

		AVTMeasureNode = *pAVTNodeLeft + (*pAVTNodeRight - *pAVTNodeLeft) *
			(((target_pipe_loc * target_pipe_len) - pAVTNodeLeft->x) / (pAVTNodeRight->x - pAVTNodeLeft->x));
	}
	else AVTMeasureNode = *pAVTNodeLeft; // For single node pipes

	AVTCurrentValue = AVTMeasureNode.p_dash * pPpt->PREF;
	AVTCumulativeVal += AVTCurrentValue;

	++AVTInCycleStepCount;





	// What fraction through pulse are we?
	time_in_period = fmod(time_or_CADelapsed / (360 * (this->pEng->cycle / 2)) * this->pEng->period, this->pEng->period);
	//CAD_in_period = fmod(time_or_CADelapsed, 360 * (this->pEng->cycle/2));

	fraction_cycle = time_in_period / this->pEng->period;
	/*
		cout << "time_or_CADelapsed = " << time_or_CADelapsed << "\ttime_in_period = " << time_in_period //<< "s\tCAD_in_period = " << CAD_in_period
			<< " deg\tthis->pEng->ca = " << this->pEng->ca << " deg\tfraction_of_cycle = " << fraction_cycle << endl;
	*/

	// Find the relevant data point in avt_raw for the current crank angle
	int row = 0;
	do { ++row; } // Hence first row evaluated is [1], on purpose
	while (avt_raw[row][CAD] < pEng->ca && row < numAVTDataPoints - 1);


	// Interpolate for target value at time_in_period between [row-1] and [row]
	//cout << "pEng->ca = " << pEng->ca << "\t lies between avt_raw[row-1=" << row-1 << "][CAD=" << CAD << "] = " << avt_raw[row-1][CAD] << " and avt_raw[row=" << row << "][CAD=" << CAD << "] = " << avt_raw[row][CAD] << endl;
	AVTargetValue = avt_raw[row - 1][TARGET] 
					+ (pEng->ca - avt_raw[row - 1][CAD]) 
					* ((avt_raw[row][TARGET] - avt_raw[row - 1][TARGET]) / (avt_raw[row][CAD] - avt_raw[row - 1][CAD]));

	//cout << "AVTargetValue = " << AVTargetValue << endl;
	
	AVTErrorPrev = avt_raw[row][ERROR];
	AVTError = AVTCurrentValue - AVTargetValue;
	
	if (AVTError / AVTErrorPrev < 0) avt_raw[row][SIGN_CHANGE] = double(true);
	else avt_raw[row][SIGN_CHANGE] = double(false);

	if (fabs(avt_raw[row][PSI_CORR] - avt_raw[row][PSI_CORR_PREV]) > 1e-6) { // Prevent div by 0

		AVTErrorGrad = (AVTError - AVTErrorPrev) / (avt_raw[row][PSI_CORR] - avt_raw[row][PSI_CORR_PREV]);
	}
	else AVTErrorGrad = 0;
	
	avt_raw[row][RECORDED] = AVTCurrentValue;
	avt_raw[row][ERROR] = AVTError;
	avt_raw[row][ERROR_GRAD] = AVTErrorGrad;

/*
	cout << "avt_raw[row=" << row << "][CAD=" << CAD << "] = " << avt_raw[row][CAD] << "\t";
	cout << "avt_raw[row=" << row << "][TARGET=" << TARGET << "] = " << avt_raw[row][TARGET] << "\t";
	cout << "avt_raw[row=" << row << "][RECORDED=" << RECORDED << "] = " << avt_raw[row][RECORDED] << "\t";
	cout << "avt_raw[row=" << row << "][ERROR=" << ERROR << "] = " << avt_raw[row][ERROR] << "\t";
	cout << "avt_raw[row=" << row << "][ERROR_GRAD=" << ERROR_GRAD << "] = " << avt_raw[row][ERROR_GRAD] << "\t";
	cout << "avt_raw[row=" << row << "][SIGN_CHANGE=" << SIGN_CHANGE << "] = " << TrueOrFalse(DoubleToBool(avt_raw[row][SIGN_CHANGE])) << "\t";
	
	cout << "avt_raw[row=" << row << "][CORR_FACTOR=" << CORR_FACTOR << "] = " << avt_raw[row][CORR_FACTOR] << "\t";
	cout << "avt_raw[row=" << row << "][PSI_UNCORR=" << PSI_UNCORR << "] = " << avt_raw[row][PSI_UNCORR] << "\t";
	cout << "avt_raw[row=" << row << "][PSI_CORR=" << PSI_CORR << "] = " << avt_raw[row][PSI_CORR] << "\t";
	cout << "avt_raw[row=" << row << "][PSI_CORR_PREV=" << PSI_CORR_PREV << "] = " << avt_raw[row][PSI_CORR_PREV] << "\t";
	cout << endl;
	//exit(1);
//*/

}
