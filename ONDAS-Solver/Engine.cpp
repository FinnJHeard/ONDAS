// Engine.cpp: implementation of the CEngine class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Engine.h"
#include "Tools.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CEngine::CEngine()
{

}

CEngine::~CEngine()
{

}

// Copy constructor
CEngine::CEngine(const CEngine& inEng)
{
	// Configuration
	// =============
	ID = inEng.ID;
	
	cycle = inEng.cycle;
	reveng = inEng.reveng;

	dcyl = inEng.dcyl;
	stroke = inEng.stroke;
	conrod = inEng.conrod;
	cr = inEng.cr;

	pcr = inEng.pcr;
	Tcr = inEng.Tcr;
	EVO = inEng.EVO;
	IVO = inEng.IVO;

	IVOL = inEng.IVOL;
	IAIR = inEng.IAIR;
	IPIPE = inEng.IPIPE;

	// Simulation variables
	// ====================
	rotation = inEng.rotation;
	ca = inEng.ca;
	ca_old = inEng.ca_old;
	ca_elapsed = inEng.ca_elapsed;
	degrees_elapsed = inEng.degrees_elapsed;
	ca_start = inEng.ca_start;
	rev = inEng.rev;
	NEW_CYCLE = inEng.NEW_CYCLE;
}

CEngine& CEngine::operator=(const CEngine& inEng)
{
	if(this != &inEng)
	{
		// Configuration
		// =============
		ID = inEng.ID;
	
		cycle = inEng.cycle;
		reveng = inEng.reveng;

		dcyl = inEng.dcyl;
		stroke = inEng.stroke;
		conrod = inEng.conrod;
		cr = inEng.cr;

		pcr = inEng.pcr;
		Tcr = inEng.Tcr;
		EVO = inEng.EVO;
		IVO = inEng.IVO;

		IVOL = inEng.IVOL;
		IAIR = inEng.IAIR;
		IPIPE = inEng.IPIPE;

		// Simulation variables
		// ====================
		rotation = inEng.rotation;
		ca = inEng.ca;
		ca_old = inEng.ca_old;
		ca_elapsed = inEng.ca_elapsed;
		degrees_elapsed = inEng.degrees_elapsed;
		ca_start = inEng.ca_start;
		rev = inEng.rev;
		NEW_CYCLE = inEng.NEW_CYCLE;
	}
	return *this;
}
/*
void CEngine::ReadProperties(CProperties* pPpt, int id)
{
	ReadInput(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->ENG_FILE));
}
*/
void CEngine::Initialise(CProperties* pPpt, int id, std::string param_dir, char* eng_file, int ncyls, int ninpipes, int assyid, string parent_assy_res_dir, string calling_object_str)
{
	if (pPpt->SHOW_calls) { pPpt->Out("CEngine.Initialise ("); pPpt->Out(calling_object_str); pPpt->Out(")\n"); }

	this->ENG_FILE = eng_file;
	this->NCYLS = ncyls;
	ReadInput(pPpt, ConstructString(pPpt, param_dir, this->ENG_FILE/*eng_file*/), ninpipes);
	
	// Apply engine properties
	// =======================
	ID = id;
	AssyID = assyid;

	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Initialise\n");}

	// Initialise conditions
	// ================
	if(ca_start>=720)
	{
		cout << "ca_start was " << ca_start << endl;
		ca_start = ca_start - 720*( (ca_start/720) - (fmod(ca_start, 720)/720) ); 
		// Calculates how many full 720 degrees in ca_start and subtracts to get a value 0-720
		cout << "ca_start should be 0-720 degrees!\n";
		cout << "new ca_start is now " << ca_start << endl;
	}

	rev = 1;
	if(ca_start>=360)
	{
		rotation = ca_start - 360;
		++rev;
	}
	else rotation = ca_start;

	del_theta = 0;
	ca = 0 + ca_start;
	ca_elapsed = 0 + ca_start;
	degrees_elapsed = 0;

	// Work and power
	// ==============
	del_WK_total = 0;		// Zero the work done for the present iteration
	P_inst = 0;				// Zero the instantaneous power
	CYCLE_WK = 0;			// Zero the cumulative cycle work
	PREV_CYCLE_WK = 0;		// Zero the previous cycle work
	PREV_CYCLE_POWER = 0;	// Zero the power calulated over most recently complete cycle
	CYCLE_START_TIME = 0;	// Zero current cycle start time

	PREV_CYCLE_FUEL = 0;
	CYCLE_FUEL = 0;
	eta = 1;
	PREV_CYCLE_SFC = 0;


	// AVT setup
	period = 1 / ( reveng / (this->cycle / 2));
	//cout << "reveng = " << reveng << " s^-1" << endl;
	//cout << "period = " << period << " s" << endl;
	//exit(1);
}

void CEngine::RunBoundary(CProperties* pPpt, CTime &rMyTime, double DELZe, double DELZi, double time)
{
	Update(rMyTime, DELZe, DELZi, pPpt, time);
}

void CEngine::Update(CTime &rMyTime, double DELZe, double DELZi, CProperties* pPpt, double time)
{
	// Set up time step
	double del_t;
	double del_t_e = DELZe*pPpt->xref/pPpt->AREFe;
	double del_t_i = DELZi*pPpt->xref/pPpt->AREFi;

	// Check that del_t_e and del_t_i are the same
	if(fabs(del_t_e-del_t_i)<1e-12) 
	{
		del_t = del_t_e; 
		//cout << "del_t_e = del_t_i. OK." << endl;
	}
	else
	{
		del_t = del_t_e; 
//		cout << "Cylinder: del_t_e != del_t_i. Using del_t_e." << endl;
	}

	double del_theta_rad = 2.0*PI*reveng*del_t;
	del_theta = del_theta_rad*180.0/PI;

	rotation += del_theta;
	ca_old = ca;
	ca += del_theta;
	ca_elapsed += del_theta;
	degrees_elapsed += del_theta; 
	
	if(rotation>360.0)	// 360 degrees per rev. whether 4- or 2- stroke
	{
		rotation = rotation - 360.0; // Adjust rotation accordingly
		++rev; // Increment rev. counter
	}

	if(ca>180.0*cycle)	
	{
		ca = ca - 180.0*cycle; // Adjust ca accordingly
		NEW_CYCLE = true;	// Flag a new cycle
//		++cycle; // Increment cycle counter
		
		// Calculate cycle power for cycle just completed
		PREV_CYCLE_POWER = CYCLE_WK/(time - CYCLE_START_TIME);
		
		// Reset cycle variables
		CYCLE_START_TIME = time;		// Set new cycle start time
		
		PREV_CYCLE_WK = CYCLE_WK;		// Save previous cycle work
		CYCLE_WK = 0;
		
		PREV_CYCLE_FUEL = CYCLE_FUEL;	// Save previous cycle fuel	
		CYCLE_FUEL = 0;
		
		eta = PREV_CYCLE_WK/(PREV_CYCLE_FUEL*lhv); // Thermal efficiency
		PREV_CYCLE_SFC = (PREV_CYCLE_FUEL*1000)/(PREV_CYCLE_WK/(60*60*1000)); // (Indicated) specific fuel consumption (g/kWh)
	}
	else {
		NEW_CYCLE = false;
	}
//*/

	// Work and power
	// ==============
	del_WK_total = 0;			// Zero the work done for the present iteration
	P_inst = 0;					// Zero the instantaneous power
}

void CEngine::ReadInput(CProperties* pPpt, char *InputFile, int ninpipes)
{
	int last_entry;
	CommonReadInput(InputFile, pPpt->NUM_PARAMS, this->labels, this->values, this->strings, last_entry);
	for(int r=0; r<last_entry+1; ++r)
	{
		// ====================================================================================================
		// Parameter file for Engine []
		// ====================================================================================================

		// Optional object description (max. 500 characters - use underscores for spaces)
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "strDesc") == 0) strDesc = strings[r];

		// Operation
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "ca_start") == 0) ca_start = values[r];
		if(strcmp(labels[r], "cycle") == 0) cycle = values[r];
		if(strcmp(labels[r], "reveng") == 0) reveng = values[r]/60;	// Conversion to s^-1 from min^-1

		// Cylinder geometry
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "dcyl") == 0) dcyl = values[r]/1000; // Conversion from mm to m
		if(strcmp(labels[r], "stroke") == 0) stroke = values[r]/1000;
		if(strcmp(labels[r], "conrod") == 0) conrod = values[r]/1000;
		if(strcmp(labels[r], "cr") == 0) cr = values[r];

		// Special cases
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "IVOL") == 0) IVOL = DoubleToBool(values[r]);
		if(strcmp(labels[r], "IAIR") == 0) IAIR = DoubleToBool(values[r]);
		if(strcmp(labels[r], "IPIPE") == 0) IPIPE = DoubleToBool(values[r]);
		
			// Air receiver (IPIPE=false)
			// ----------------------------------------------------------------------------------------------------
			if(strcmp(labels[r], "PA") == 0) PA = values[r];
			if(strcmp(labels[r], "TA") == 0) TA = values[r];

		// Cylinder model
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "MODEL") == 0)
		{
			switch(int(values[r]))
			{
			case 0:
				MODEL = CONST_P;
				break;
			case 1:
				MODEL = CONST_VOL;
				break;
			case 2:
				MODEL = READ_P_T;
				break;
			case 3:
				MODEL = RESET_P_T;
				break;
			case 4:
				MODEL = WATSON;
				break;
			default:
				cout << "CEngine::ReadInput [" << ID << "]: unknown cylinder model (" << int(values[r]) << "). Selecting default RESET_P_T\n";
				MODEL = RESET_P_T;
				break;
			}
		}
		
		// (0) Constant conditions (cylinder acts as a reservoir)
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "pres") == 0) pres = values[r];
		if(strcmp(labels[r], "Tres") == 0) Tres = values[r];
		if(strcmp(labels[r], "Vres") == 0) Vres = values[r]/1000;		// Convert to m^3

		// (1) Constant volume, variable conditions (cylinder acts as a volume)
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "vconst") == 0) vconst = values[r]/1000;	// Convert to m^3
		if(strcmp(labels[r], "pinit") == 0) pinit = values[r];
		if(strcmp(labels[r], "Tinit") == 0) Tinit = values[r];

		// (2) Read cylinder pressure & temperature from file
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "cyl_file") == 0) cyl_file = strings[r];

		// (3) Cylinder pressure & temperature reset at EVO
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "RESET_EVO") == 0) RESET_EVO = DoubleToBool(values[r]);
		if(strcmp(labels[r], "pcr") == 0) pcr = values[r];
		if(strcmp(labels[r], "Tcr") == 0) Tcr = values[r];

		// (4) Single-zone heat release model
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "x") == 0) x = int(values[r]);
		if(strcmp(labels[r], "y") == 0) y = int(values[r]);
		if(strcmp(labels[r], "lhv") == 0) lhv = values[r];
		if(strcmp(labels[r], "ign_ca") == 0) ign_ca = values[r];
		if(strcmp(labels[r], "comb_stop") == 0) comb_stop = values[r];
		if(strcmp(labels[r], "inj_start") == 0) inj_start = values[r];
		if(strcmp(labels[r], "inj_stop") == 0) inj_stop = values[r];
		if(strcmp(labels[r], "m_f0") == 0) m_f0 = values[r];
		if(strcmp(labels[r], "comb_a") == 0) comb_a = values[r];
		if(strcmp(labels[r], "comb_b") == 0) comb_b = values[r];
		if(strcmp(labels[r], "comb_c") == 0) comb_c = values[r];
		if(strcmp(labels[r], "ADD_MFB") == 0) ADD_MFB = DoubleToBool(values[r]);
		
		// Cylinder heat transfer model
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "HT_MODEL") == 0)
		{
			switch(int(values[r]))
			{
			case 0:
				HT_MODEL = NONE;
				break;
			case 1:
				HT_MODEL = ANNAND;
				break;
			case 2:
				HT_MODEL = WOSCHNI;
				break;
			default:
				cout << "CEngine::ReadInput [" << ID << "]: unknown cylinder heat transfer model (" << int(values[r]) << "). Selecting default ANNAND\n";
				HT_MODEL = ANNAND;
				break;
			}
		}
		
		// (1) Annand convective heat transfer
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "a") == 0) a = values[r];
		if(strcmp(labels[r], "b") == 0) b = values[r];
		if(strcmp(labels[r], "Tw") == 0) Tw = values[r];

		// (2) Woschni correlation
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "T_c") == 0) T_c = values[r];
		if(strcmp(labels[r], "h_cc") == 0) h_cc = values[r];
		if(strcmp(labels[r], "k") == 0) k = values[r];
		if(strcmp(labels[r], "tw") == 0) tw = values[r];
		if(strcmp(labels[r], "epsilon") == 0) epsilon = values[r];
		if(strcmp(labels[r], "WOS_TOL") == 0) WOS_TOL = values[r];
	
		// Valve and general AVT parameters
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "PRE_vt") == 0) PRE_vt = DoubleToBool(values[r]);
		if (strcmp(labels[r], "nSamples") == 0) nSamples = int(values[r]);
		if (strcmp(labels[r], "tol") == 0) tol = double(values[r]);
		if (strcmp(labels[r], "SUSPEND") == 0) SUSPEND = DoubleToBool(values[r]);

		// Exhaust valve parameters
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NEXVALVES") == 0) NEXVALVES = int(values[r]);
		if (strcmp(labels[r], "REF_DIA_EXH") == 0) REF_DIA_EXH = values[r];
		if(strcmp(labels[r], "EV_WAVEFORM") == 0) EV_WAVEFORM = int(values[r]);
//		if(strcmp(labels[r], "EV_SWITCH") == 0) EV_SWITCH = DoubleToBool(values[r]);
//		if (strcmp(labels[r], "EV_L_OR_A") == 0) EV_L_OR_A = DoubleToBool(values[r]);
		if (strcmp(labels[r], "EV_L_OR_A") == 0) EV_L_OR_A = bool(values[r]);
		if (strcmp(labels[r], "EV_MIN_OPENING") == 0) EV_MIN_OPENING = values[r];
		if (strcmp(labels[r], "EV_MAX_OPENING") == 0) EV_MAX_OPENING = values[r];
		
		// If exhaust valve acts as a switch (open/close), i.e., EV_WAVEFORM == 0
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "EV_OPEN") == 0) EV_OPEN = DoubleToBool(values[r]);
//		if(strcmp(labels[r], "EV_CALC") == 0) EV_CALC = DoubleToBool(values[r]);
//		if(strcmp(labels[r], "EV_SQUARE") == 0) EV_SQUARE = DoubleToBool(values[r]);
		
		// If exhaust valve lift or area profile is triangular or square or sinusoidal or a Fourier series, i.e., EV_WAVEFORM == 1 or 2 or 3 or 4
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "EV_BEGIN_OPENING") == 0) EV_BEGIN_OPENING= values[r];
		if (strcmp(labels[r], "EV_END_or_phi") == 0) EV_END_or_phi = bool(values[r]);
		if (strcmp(labels[r], "EV_END_CLOSING") == 0) EV_END_CLOSING = values[r];
		if (strcmp(labels[r], "EV_phi") == 0) EV_phi = values[r];

		// If exhaust valve lift or area profile is read from file, i.e., EV_WAVEFORM == 5
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "EV_FILE") == 0) EV_FILE = strings[r];
		if (strcmp(labels[r], "CAM_TNG_ANG_EXH") == 0) CAM_TNG_ANG_EXH = values[r];
		if (strcmp(labels[r], "CAM_TNG_LOC_EXH") == 0) CAM_TNG_LOC_EXH = int(values[r]);

		// If exhaust valve data are lift values, i.e., EV_L_OR_A == 0
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "LASH_EXH") == 0) LASH_EXH = values[r];
		if (strcmp(labels[r], "EVFA_FILE") == 0) EVFA_FILE = strings[r];
		if (strcmp(labels[r], "CONST_REF_AREA_EXH") == 0) CONST_REF_AREA_EXH = DoubleToBool(values[r]);
		if (strcmp(labels[r], "ONUM_EXH") == 0) ONUM_EXH = values[r];

//		// Else exhaust valve data are area values, i.e., EV_L_OR_A == 1
//		// ----------------------------------------------------------------------------------------------------
		
//		// Exhaust valve flow or discharge coefficients (can be used with either lift or area values)
//		// ----------------------------------------------------------------------------------------------------
//		if (strcmp(labels[r], "CONST_REF_AREA_EXH") == 0) CONST_REF_AREA_EXH = DoubleToBool(values[r]);
//		if(strcmp(labels[r], "ONUM_EXH") == 0) ONUM_EXH = values[r];

		// Variable valve timing and lift
		// ----------------------------------------------------------------------------------------------------
//		if(strcmp(labels[r], "LIFT_MULT_EXH") == 0) LIFT_MULT_EXH = values[r];
//		if(strcmp(labels[r], "AREA_MULT_EXH") == 0) AREA_MULT_EXH = values[r];
		if (strcmp(labels[r], "EV_OPENING_MULT") == 0) EV_OPENING_MULT = values[r];
		if (strcmp(labels[r], "ANG_MULT_EXH") == 0) ANG_MULT_EXH = values[r];
		if (strcmp(labels[r], "scale_factor_ev") == 0) scale_factor_ev = values[r];
		
		// Active Valve Train
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "EV_ACTIVE") == 0) EV_ACTIVE = DoubleToBool(values[r]);
		if (strcmp(labels[r], "EV_TARGET_PIPE_EX") == 0) EV_TARGET_PIPE_EX = bool(values[r]);
		if (strcmp(labels[r], "EV_target_pipe_num") == 0) EV_target_pipe_num = int(values[r]);
		if (strcmp(labels[r], "EV_target_pipe_loc") == 0) EV_target_pipe_loc = values[r];
		if (strcmp(labels[r], "EV_ACTIVE_FILE_P") == 0) EV_ACTIVE_FILE_P = strings[r];
		
		// Intake valve parameters
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "NINVALVES") == 0) NINVALVES = int(values[r]);
		if (strcmp(labels[r], "REF_DIA_INT") == 0) REF_DIA_INT = values[r];
		if (strcmp(labels[r], "IV_WAVEFORM") == 0) IV_WAVEFORM = int(values[r]);
//		if (strcmp(labels[r], "IV_SWITCH") == 0) IV_SWITCH = DoubleToBool(values[r]);
//		if (strcmp(labels[r], "IV_L_OR_A") == 0) IV_L_OR_A = DoubleToBool(values[r]);
		if (strcmp(labels[r], "IV_L_OR_A") == 0) IV_L_OR_A = bool(values[r]);
		if (strcmp(labels[r], "IV_MIN_OPENING") == 0) IV_MIN_OPENING = values[r];
		if (strcmp(labels[r], "IV_MAX_OPENING") == 0) IV_MAX_OPENING = values[r];

		// If intake valve acts as a switch (open/close), i.e., IV_WAVEFORM == 0
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "IV_OPEN") == 0) IV_OPEN = DoubleToBool(values[r]);
//		if (strcmp(labels[r], "IV_CALC") == 0) IV_CALC = DoubleToBool(values[r]);
//		if (strcmp(labels[r], "IV_SQUARE") == 0) IV_SQUARE = DoubleToBool(values[r]);

		// If intake valve lift or area profile is triangular or square or sinusoidal or a Fourier series, i.e., IV_WAVEFORM == 1 or 2 or 3 or 4
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "IV_BEGIN_OPENING") == 0) IV_BEGIN_OPENING = values[r];
		if (strcmp(labels[r], "IV_END_or_phi") == 0) IV_END_or_phi = bool(values[r]);
		if (strcmp(labels[r], "IV_END_CLOSING") == 0) IV_END_CLOSING = values[r];
		if (strcmp(labels[r], "IV_phi") == 0) IV_phi = values[r];

		// If intake valve lift or area profile is read from file, i.e., IV_WAVEFORM == 5
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "IV_FILE") == 0) IV_FILE = strings[r];
		if (strcmp(labels[r], "CAM_TNG_ANG_INT") == 0) CAM_TNG_ANG_INT = values[r];
		if (strcmp(labels[r], "CAM_TNG_LOC_INT") == 0) CAM_TNG_LOC_INT = int(values[r]);

		// If intake valve data are lift values, i.e., EV_L_OR_A == 1 == true
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "LASH_INT") == 0) LASH_INT = values[r];
		if (strcmp(labels[r], "IVFA_FILE") == 0) IVFA_FILE = strings[r];
		if (strcmp(labels[r], "CONST_REF_AREA_INT") == 0) CONST_REF_AREA_INT = DoubleToBool(values[r]);
		if (strcmp(labels[r], "ONUM_INT") == 0) ONUM_INT = values[r];

//		// Else intake valve data are area values, i.e., IV_L_OR_A == 0 == false
//		// ----------------------------------------------------------------------------------------------------

//		// Intake valve flow or discharge coefficients (can be used with either lift or area values)
//		// ----------------------------------------------------------------------------------------------------
//		if (strcmp(labels[r], "CONST_REF_AREA_INT") == 0) CONST_REF_AREA_INT = DoubleToBool(values[r]);
//		if (strcmp(labels[r], "ONUM_INT") == 0) ONUM_INT = values[r];

		// Variable valve timing and lift
		// ----------------------------------------------------------------------------------------------------
//		if(strcmp(labels[r], "LIFT_MULT_INT") == 0) LIFT_MULT_INT = values[r];
//		if(strcmp(labels[r], "AREA_MULT_INT") == 0) AREA_MULT_INT = values[r];
		if (strcmp(labels[r], "IV_OPENING_MULT") == 0) IV_OPENING_MULT = values[r];
		if (strcmp(labels[r], "ANG_MULT_INT") == 0) ANG_MULT_INT = values[r];
		if (strcmp(labels[r], "scale_factor_iv") == 0) scale_factor_iv = values[r];

		// Active Valve Train
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "IV_ACTIVE") == 0) IV_ACTIVE = DoubleToBool(values[r]);
		if (strcmp(labels[r], "IV_TARGET_PIPE_EX") == 0) IV_TARGET_PIPE_EX = bool(values[r]);
		if (strcmp(labels[r], "IV_target_pipe_num") == 0) IV_target_pipe_num = int(values[r]);
		if (strcmp(labels[r], "IV_target_pipe_loc") == 0) IV_target_pipe_loc = values[r];
		if (strcmp(labels[r], "IV_ACTIVE_FILE_P") == 0) IV_ACTIVE_FILE_P = strings[r];
	}

	if(pPpt->HOMENTROPIC && HT_MODEL != NONE)
	{
		HT_MODEL = NONE;
		cout << "Automatically selecting no cylinder heat transfer model (NONE) for this homentropic simulation\n";
	}
	if(MODEL==CONST_P && HT_MODEL != NONE)
	{
		HT_MODEL = NONE;
		cout << "Automatically selecting no cylinder heat transfer model (NONE) for this cylinder (CONST_P)\n";
	}
	if(MODEL==CONST_VOL && HT_MODEL != NONE)
	{
		HT_MODEL = NONE;
		cout << "Automatically selecting no cylinder heat transfer model (NONE) for this cylinder (CONST_VOL)\n";
	}
	if(MODEL==READ_P_T && HT_MODEL != NONE)
	{
		HT_MODEL = NONE;
		cout << "Automatically selecting no cylinder heat transfer model (NONE) for this cylinder (READ_P_T)\n";
	}

	// Make sure enough intake pipes exist to allow intake side of cylinders to be intake pipes
	if(NCYLS>0 && ninpipes==0 && IPIPE) {
		cout << "At least one intake pipe is required to have an intake pipe as intake to cylinder, i.e. for IPIPE==" << TrueOrFalse(IPIPE) << endl;
		cout << "Adjust " << InputFile << " to satisfy these requirements, and re-run" << endl;
		exit(1);
	}

	// Initialise coeffs
	// ----------------------------------------------------------------------------------------------------
	coeffs = k*h_cc/(k + h_cc*tw);		
			
	// Initialise combustion parameters before ListProperties
	// ----------------------------------------------------------------------------------------------------
	double z = y/x;										// Fuel composition ratio of Hydrogen to Carbon.
	phiFAstoich = (12.011 + 1.008*z)/(34.56*(4 + z));	// Calculate stoichiometric fuel/air ratio (need only do once)

	// Initialise overall engine work and power variables before ListProperties
	// ----------------------------------------------------------------------------------------------------
	del_WK_total = 0;
	P_inst = 0;
	CYCLE_START_TIME = 0;
	CYCLE_WK = 0;
	PREV_CYCLE_WK = 0;
	PREV_CYCLE_POWER = 0;
	CYCLE_FUEL = 0;
	PREV_CYCLE_FUEL = 0;
	eta = 0;
	PREV_CYCLE_SFC = 0;
}

void CEngine::ListProperties(CProperties* pPpt)
{
	if (pPpt->SHOW_calls) { pPpt->Out(Identify()); pPpt->Out(".ListProperties\n"); }

	// ====================================================================================================
	// Parameter file for Engine []
	// ====================================================================================================

	pPpt->Out(Underline(Identify(), "=", "\t", strDesc));
	pPpt->Out("\n");

	cout << Underline("Operation", "-", "\t");
	// ----------------------------------------------------------------------------------------------------
	cout << "\tEngine parameter file, ENG_FILE\t\t\t=\t" << /*pPpt->*/ENG_FILE << "\n";
	cout << "\tNo. of cylinders, NCYLS\t\t\t\t=\t" << /*pPpt->*/NCYLS << "\n";
	cout << "\tCycle type, cycle\t\t\t\t=\t" << cycle << "-stroke\n";
	cout << "\tEngine speed, reveng\t\t\t\t=\t" << reveng * 60 << " rev/min\n";
	cout << "\tCA at sim. start, ca_start\t\t\t=\t" << ca_start << Deg() << "CA" << endl;
	cout << endl;

	cout << Underline("Cylinder geometry", "-", "\t");
	// ----------------------------------------------------------------------------------------------------
	cout << "\tBore, dcyl\t\t\t\t\t=\t" << dcyl * 1000 << " mm\n";
	cout << "\tStroke, stroke\t\t\t\t\t=\t" << stroke * 1000 << " mm\n";
	cout << "\tCylinder displacement\t\t\t\t=\t" << (0.25 * PI * pow(dcyl, 2) * stroke) * 1000 << " l\n";
	cout << "\tTotal displacement\t\t\t\t=\t" << /*pPpt->*/NCYLS * (0.25 * PI * pow(dcyl, 2) * stroke) * 1000 << " l\n";
	cout << "\tConrod length, conrod\t\t\t\t=\t" << conrod * 1000 << " mm\n";
	cout << "\tCompression ratio, cr\t\t\t\t=\t" << cr << "\n";
	//cout << "\tExhaust valve open, EVO\t\t\t=\t" << EVO << Deg() << "CA" << endl;
	//cout << "\tIntake valve open, IVO\t\t\t\t=\t" << IVO << Deg() << "CA" << endl;
	//cout << "\tIntake valve close, IVC\t\t\t=\t" << IVC << Deg() << "CA" << endl;
	cout << endl;

	cout << Underline("Special cases", "-", "\t");
	// ----------------------------------------------------------------------------------------------------
	if (IVOL)	cout << "\tVariable volume cylinder, IVOL\t\t\t=\t" << TrueOrFalse(IVOL) << "\n";
	else		cout << "\tConstant volume cylinder, IVOL\t\t\t=\t" << TrueOrFalse(IVOL) << "\n";

	if (IAIR) {
		cout << "\tIntake valve present, IAIR\t\t\t=\t" << TrueOrFalse(IAIR) << "\n";
		if (IPIPE) cout << "\t- connects to pipe, IPIPE\t\t\t=\t" << TrueOrFalse(IPIPE) << "\n";
		else {
			cout << "\t- connects to receiver, IPIPE\t\t\t=\t" << TrueOrFalse(IPIPE) << "\n";
			cout << "\t  with:" << "\n";
			cout << "\t  - stagnation pressure, PA\t\t\t=\t" << PA << " bar" << endl;
			cout << "\t  - stagnation temperature, TA\t\t\t=\t" << TA << " K" << endl;
		}
	}
	else cout << "\tNo intake to cylinder, IAIR\t\t\t=\t" << TrueOrFalse(IAIR) << "\n";
	cout << endl;

	cout << Underline("Cylinder model", "-", "\t");
	// ----------------------------------------------------------------------------------------------------
	switch (MODEL) {
	case CONST_P:
		cout << "\tConstant pressure cylinder" << endl;
		cout << "\t- pressure, pres\t\t\t\t=\t" << pres << " bar" << endl;
		cout << "\t- temperature, Tres\t\t\t\t=\t" << Tres << " K" << endl;
		cout << "\t- volume, Vres\t\t\t\t\t=\t" << Vres * 1000 << " litres" << endl;
	case CONST_VOL:
		cout << "\tConstant volume cylinder, variable conditions" << endl;
		cout << "\t- constant cylinder volume, vconst\t\t=\t" << vconst * 1000 << " litres" << endl;
		cout << "\t- initial pressure, pinit\t\t\t=\t" << pinit << " bar" << endl;
		cout << "\t- initial temperature, Tinit\t\t\t=\t" << Tinit << " K" << endl;
		break;
	case READ_P_T:
		cout << "\tRead cylinder pressure and temperature from file" << endl;
		cout << "\tP & T trace file (cyl_file)\t\t\t=\t" << cyl_file << endl;
		break;
	case RESET_P_T:
		if (RESET_EVO) cout << "\tResetting cylinder conditions at EVO" << endl;
		else cout << "\tResetting cylinder conditions at TDC" << endl;
		cout << "\t- reset pressure (pcr)\t\t\t\t=\t" << pcr << " bar\n";
		cout << "\t- reset temperature (Tcr)\t\t\t=\t" << Tcr << " K\n";
		break;
	case WATSON:
		cout << "\tWatson single-zone heat release model" << endl;
		cout << "\tFuel composition CxHy, where: \n";
		cout << "\t- C proportion, x\t\t\t\t=\t" << x << " \n";
		cout << "\t- H proportion, y\t\t\t\t=\t" << y << " \n";
		cout << "\t- Stoichiometric fuel/air ratio (phiFAstoich)\t=\t" << phiFAstoich << " \n";
		cout << "\t- lower heating value, lhv\t\t\t=\t" << lhv << " J/kg\n";
		cout << "\tIgnition, ign_ca\t\t\t\t=\t" << ign_ca << Deg() << "CA\n";
		cout << "\tCombustion assumed complete at, comb_stop\t=\t" << comb_stop << Deg() << "CA\n";
		cout << "\tStart of fuel injection, inj_start\t\t=\t" << inj_start << Deg() << "CA\n";
		cout << "\tEnd of fuel injection, inj_stop\t\t=\t" << inj_stop << Deg() << "CA\n";
		cout << "\tTotal fuel mass inj. per cycle per cyl., m_f0\t=\t" << m_f0 << " kg\n";
		cout << "\tCombustion parameter a, comb_a\t\t\t=\t" << comb_a << " \n";
		cout << "\tCombustion parameter b, comb_b\t\t\t=\t" << comb_b << " \n";
		cout << "\tCombustion parameter c, comb_c\t\t\t=\t" << comb_c << " \n";
		if (ADD_MFB) cout << "\tAdding burned fuel mass to cyl. mass, ADD_MFB\t=\t" << TrueOrFalse(ADD_MFB) << endl;
		else cout << "\tNot adding burned fuel to cyl. mass, ADD_MFB\t=\t" << TrueOrFalse(ADD_MFB) << endl;
		break;
	default:
		MODEL = RESET_P_T;
		cout << "\tCalculating P & T normally\n";
		cout << "\tNo combustion; ";
		if (RESET_EVO) cout << "resetting cylinder at EVO\n";
		else cout << "resetting cylinder at TDC\n";
		cout << "\t- reset pressure (pcr)\t\t\t\t=\t" << pcr << " bar\n";
		cout << "\t- reset temperature (Tcr)\t\t\t=\t" << Tcr << " K\n";
		break;
	}
	cout << endl;

	cout << Underline("Cylinder heat transfer model", "-", "\t");
	// ----------------------------------------------------------------------------------------------------
	switch (HT_MODEL)
	{
	case NONE:
		cout << "\tNo cylinder heat transfer";
		if (pPpt->HOMENTROPIC) cout << " - homentropic simulation";
		cout << "\n";
		break;
	case ANNAND:
		cout << "\tAnnand convective heat transfer" << endl;
		cout << "\t- charge motion/engine design coefficient (a)\t=\t" << a << " \n";
		cout << "\t- combustion coefficient (normal = 0.7) (b)\t=\t" << b << " \n";
		cout << "\tCylinder wall temperature (Tw)\t\t\t=\t" << Tw << " K\n";
		break;
	case WOSCHNI:
		cout << "\tWoschni correlation" << endl;
		cout << "\tCoolant temperature (T_c)\t\t\t=\t" << T_c << " K\n";
		cout << "\tCoolant heat transfer coefficient (h_cc)\t=\t" << h_cc << " W/(m^2.K)\n";
		cout << "\tLiner wall thermal conductivity (k)\t\t=\t" << k << " W/(m.K)\n";
		cout << "\tLiner wall thickness (tw)\t\t\t=\t" << tw << " m\n";
		cout << "\tHeat transfer combined coefficients (coeffs)\t=\t" << coeffs << " \n";
		cout << "\tEmissivity (peak conditions) (epsilon)\t\t=\t" << epsilon << " \n";
		cout << "\tWoschni radiative loop tolerance (WOS_TOL)\t=\t" << WOS_TOL << " \n";
		break;
	default:
		HT_MODEL = ANNAND;
		break;
	}
	cout << endl;

	pPpt->Out(Underline("Valve parameters", "-", "\t"));
	// ----------------------------------------------------------------------------------------------------
	if (PRE_vt) { pPpt->Out("\tShowing valve timing in preamble, PRE_vt\t=\t"); pPpt->Out(TrueOrFalse(PRE_vt)); pPpt->Out("\n"); }
	else { pPpt->Out("\tNot showing valve timing in preamble, PRE_vt\t=\t"); pPpt->Out(TrueOrFalse(PRE_vt)); pPpt->Out("\n"); }
	pPpt->Out("\tNo. exhaust valves per cylinder, NEXVALVES\t=\t"); pPpt->Out(NEXVALVES); pPpt->Out("\n");
	pPpt->Out("\tNo. intake valves per cylinder, NINVALVES\t=\t"); pPpt->Out(NINVALVES); pPpt->Out("\n");

	if (EV_ACTIVE || IV_ACTIVE) {

		pPpt->Out("\n");
		pPpt->Out(Underline("General AVT parameters", "-", "\t"));
		// ----------------------------------------------------------------------------------------------------
		pPpt->Out("\t\tNumber of samples per pulse, nSamples\t=\t"); pPpt->Out(nSamples); pPpt->Out("\n");
		pPpt->Out("\t\tPulse match tolerance, tol\t\t=\t"); pPpt->Out(tol); pPpt->Out("%\n");
		if (SUSPEND) pPpt->Out("\t\tSuspend simulation upon match, SUSPEND\t=\t");
		else pPpt->Out("\t\tDo not suspend simulation upon match, SUSPEND\t=\t");
		pPpt->Out(TrueOrFalse(SUSPEND)); pPpt->Out("\n");
		//pPpt->Out("\tCycles between turning on/off, nMatchCycles\t=\t"); pPpt->Out(nMatchCycles); pPpt->Out("\n");
		//pPpt->Out("\n");
	}
	pPpt->Out("\n");
	
	if (NEXVALVES > 0) {
		pPpt->Out(Underline("Exhaust valve parameters", "-", "\t\t"));
		pPpt->Out("\t\tReference diameter, REF_DIA_EXH\t\t=\t"); pPpt->Out(REF_DIA_EXH); pPpt->Out(" mm\n");
		switch (EV_WAVEFORM) {
		case 0:
			pPpt->Out("\t\tValve acts as switch, EV_WAVEFORM\t=\t"); pPpt->Out(EV_WAVEFORM); pPpt->Out("\n");
			if (EV_OPEN) pPpt->Out("\t\t- fully open, EV_OPEN\t\t\t=\t");
			else pPpt->Out("\t\t- fully closed, EV_OPEN\t\t\t=\t");
			pPpt->Out(TrueOrFalse(EV_OPEN)); pPpt->Out("\n");
			break;
		case 1:
			pPpt->Out("\t\tSquare wave opening, EV_WAVEFORM\t=\t"); pPpt->Out(EV_WAVEFORM); pPpt->Out("\n");
			break;
		case 2:
			pPpt->Out("\t\tTriangular wave opening, EV_WAVEFORM\t=\t"); pPpt->Out(EV_WAVEFORM); pPpt->Out("\n");
			break;
		case 3:
			pPpt->Out("\t\tSinusoidal wave opening, EV_WAVEFORM\t=\t"); pPpt->Out(EV_WAVEFORM); pPpt->Out("\n");
			break;
		case 4:
			pPpt->Out("\t\tFourier series wave opening, EV_WAVEFORM\t=\t"); pPpt->Out(EV_WAVEFORM); pPpt->Out("\n");
			break;
		case 5:
			pPpt->Out("\t\tFile wave opening, EV_WAVEFORM\t=\t"); pPpt->Out(EV_WAVEFORM); pPpt->Out("\n");
			break;
		default:
			pPpt->Out("\t\tUnknown wave opening, EV_WAVEFORM\t=\t"); pPpt->Out(EV_WAVEFORM); pPpt->Out(". Exiting.\n");
			exit(1);
			break;
		}

		// Switch on exhaust valve data units
		// ----------------------------------------------------------------------------------------------------
		switch (EV_L_OR_A) {
		case 0:
			pPpt->Out("\t\tValve data are lift (mm), EV_L_OR_A\t=\t"); pPpt->Out(EV_L_OR_A); pPpt->Out("\n");
			pPpt->Out("\t\tValve lash, LASH_EXH\t\t\t=\t"); pPpt->Out(LASH_EXH); pPpt->Out(" mm\n");
			pPpt->Out("\t\tFlow/discharge coeff vs l:d, EVFA_FILE\t=\t"); pPpt->Out(EVFA_FILE); pPpt->Out("\n");

			// Exhaust valve flow or discharge coefficients (can only be used with lift data since l/d is the lookup)
			// ----------------------------------------------------------------------------------------------------
			if (CONST_REF_AREA_EXH) pPpt->Out("\t\tConst ref area, CONST_REF_AREA_EXH\t=\t");
			else pPpt->Out("\t\tCurtain ref area, CONST_REF_AREA_EXH\t=\t");
			pPpt->Out(TrueOrFalse(CONST_REF_AREA_EXH)); pPpt->Out("\n");
			pPpt->Out("\t\tFlow/discharge coeff multplr, ONUM_EXH\t=\t"); pPpt->Out(ONUM_EXH);
			if (ONUM_EXH == 1.0) pPpt->Out(" (def.)\n");

			break;
		case 1:
			pPpt->Out("\t\tProfile is flow area (m^2), EV_L_OR_A\t=\t"); pPpt->Out(EV_L_OR_A); pPpt->Out("\n");
			break;
		case 2:
			pPpt->Out("\t\tProfile is psi (valve:pipe), EV_L_OR_A\t=\t"); pPpt->Out(EV_L_OR_A); pPpt->Out("\n");
			break;
		default:
			//
			break;
		}

		if (EV_WAVEFORM != 5)
		{
			pPpt->Out("\t\tMin lift/area/psi value, EV_MIN_OPENING\t=\t"); pPpt->Out(EV_MIN_OPENING); pPpt->Out("\n");
			pPpt->Out("\t\tMax lift/area/psi value, EV_MAX_OPENING\t=\t"); pPpt->Out(EV_MAX_OPENING); pPpt->Out("\n");

			// If exhaust valve lift or area profile is square or triangular or sinusoidal or a Fourier series, i.e., EV_WAVEFORM == 1 or 2 or 3 or 4
			// ----------------------------------------------------------------------------------------------------
			if (EV_WAVEFORM == 1 || EV_WAVEFORM == 2 || EV_WAVEFORM == 3 || EV_WAVEFORM == 4) {
				pPpt->Out("\t\tBegin opening angle, EV_BEGIN_OPENING\t=\t"); pPpt->Out(EV_BEGIN_OPENING); pPpt->Out("\n");
				
				if (EV_END_or_phi == 0) {
					pPpt->Out("\t\tEV_END_CLOSING decides closing, EV_END_or_phi\t=\t"); pPpt->Out(EV_END_or_phi); pPpt->Out("\n");
					pPpt->Out("\t\tEnd closing angle, EV_END_CLOSING\t=\t"); pPpt->Out(EV_END_CLOSING); pPpt->Out("\n");
				}
				else {
					if (EV_END_or_phi == 1) {
						pPpt->Out("\t\tEV_phi decides closing, EV_END_or_phi\t=\t"); pPpt->Out(EV_END_or_phi); pPpt->Out("\n");
						pPpt->Out("\t\tOpen fraction of wavelength (0-1), EV_phi\t=\t"); pPpt->Out(EV_phi); pPpt->Out("\n");
						//	double phi_orig;			// phi(0-1) based on original excitation
					}
				}
				pPpt->Out("\n");
			}
			else
			{
				if (EV_WAVEFORM == 0)
				{ 
					//
				}
				else
				{
					pPpt->Out("\t\tUnknown wave opening, EV_WAVEFORM\t=\t"); pPpt->Out(EV_WAVEFORM); pPpt->Out(". Exiting.\n");
					exit(1);
				}
			}
		}
		else
		{
			// Exhaust valve lift or area profile is read from file, i.e., EV_WAVEFORM == 5
			// ----------------------------------------------------------------------------------------------------
			pPpt->Out("\t\tValve timing file, EV_FILE\t\t=\t"); pPpt->Out(EV_FILE); pPpt->Out("`n");
			pPpt->Out("\t\tCam timing angle, CAM_TNG_ANG_EXH\t=\t"); pPpt->Out(CAM_TNG_ANG_EXH); pPpt->Out(Deg()); pPpt->Out("\n");
			pPpt->Out("\t\t- ");
			switch (CAM_TNG_LOC_EXH) {
			case 0:
				pPpt->Out("start of opening");
				break;
			case 1:
				pPpt->Out("maximum opening");
				break;
			case 2:
				pPpt->Out("end of closing");
				break;
			default:
				CAM_TNG_LOC_EXH = 1;
				pPpt->Out("maximum opening");
				break;
			}
			pPpt->Out(", CAM_TNG_LOC_EXH\t=\t"); pPpt->Out(CAM_TNG_LOC_EXH); pPpt->Out("\n");
		}
		
/*
		if (EV_L_OR_A) {
			pPpt->Out("\t\tValve data are lift (mm), EV_L_OR_A\t=\t"); pPpt->Out(TrueOrFalse(EV_L_OR_A)); pPpt->Out("\n");
//			cout << "\t\tLift m'plier, LIFT_MULT_EXH\t\t=\t" << LIFT_MULT_EXH << endl;
			pPpt->Out("\t\tValve lash, LASH_EXH\t\t\t=\t"); pPpt->Out(LASH_EXH); pPpt->Out(" mm\n");
			pPpt->Out("\t\tFlow/discharge coeff vs l:d, EVFA_FILE\t=\t"); pPpt->Out(EVFA_FILE); pPpt->Out("\n");
		}
		else {
			pPpt->Out("\t\tValve data are area (m^2), EV_L_OR_A\t=\t"); pPpt->Out(TrueOrFalse(EV_L_OR_A)); pPpt->Out("\n");
//			cout << "\t\tArea m'plier, AREA_MULT_EXH\t\t=\t" << AREA_MULT_EXH << endl;
		}
*/

		// Variable valve timing and lift
		// ----------------------------------------------------------------------------------------------------
		pPpt->Out("\t\tOpening multiplier, EV_OPENING_MULT\t=\t"); pPpt->Out(EV_OPENING_MULT); pPpt->Out("\n");
		pPpt->Out("\t\tOpen duration multiplier, ANG_MULT_EXH\t=\t"); pPpt->Out(ANG_MULT_EXH); pPpt->Out("\n");
		pPpt->Out("\t\tOpen duration compress, scale_factor_ev\t=\t"); pPpt->Out(scale_factor_ev); pPpt->Out("\n");
		
/*
		if(EV_SWITCH) {
			cout << "\t\tExhaust valve acts as switch, EV_SWITCH\t=\t" << TrueOrFalse(EV_SWITCH) << endl;	
			if(EV_OPEN) cout << "\t\t- fully open, EV_OPEN\t\t\t=\t" << TrueOrFalse(EV_OPEN) << endl;
			else cout << "\t\t- fully closed, EV_OPEN\t\t\t=\t" << TrueOrFalse(EV_OPEN) << endl;
		}
		else {
			if(EV_CALC) {
				cout << "\t\tCalculate valve timing, EV_CALC\t\t=\t" << TrueOrFalse(EV_CALC) << endl;
				cout << "\t\tSquare wave opening, EV_SQUARE\t\t=\t" << TrueOrFalse(EV_SQUARE) << endl;
				cout << "\t\tBegin valve opening, EV_BEGIN\t\t=\t" << EV_BEGIN << endl;
				cout << "\t\tEnd valve closing, EV_END\t\t=\t" << EV_END << endl;
				cout << "\t\tMaximum psi value, EV_PSI_MAX\t\t=\t" << EV_PSI_MAX << endl;
				cout << "\t\tMinimum psi value, EV_PSI_MIN\t\t=\t" << EV_PSI_MIN << endl;
			}
			else {
				cout << "\t\tInterpolate valve timing, EV_CALC\t=\t" << TrueOrFalse(EV_CALC) << endl;
				cout << "\t\tValve timing file, EV_FILE\t\t=\t" << EV_FILE << endl;
				cout << "\t\tCam timing angle, CAM_TNG_ANG_EXH\t=\t" << CAM_TNG_ANG_EXH << Deg() << endl;
				cout << "\t\t- ";
				switch(CAM_TNG_LOC_EXH) {
				case 0:
					cout << "start of opening";
					break;
				case 1:
					cout << "maximum opening";
					break;
				case 2:
					cout << "end of closing";
					break;
				default:
					CAM_TNG_LOC_EXH = 1;
					cout << "maximum opening";
					break;
				}
				cout << ", CAM_TNG_LOC_EXH\t=\t" << CAM_TNG_LOC_EXH << endl;
				cout << "\t\tAngle multiplier, ANG_MULT_EXH\t\t=\t" << ANG_MULT_EXH << endl;
				if(EV_L_OR_A) {
					cout << "\t\tData are lift values (mm), EV_L_OR_A\t=\t" << TrueOrFalse(EV_L_OR_A) << endl;
					cout << "\t\tLift m'plier, LIFT_MULT_EXH\t\t=\t" << LIFT_MULT_EXH << endl;
					cout << "\t\tLash, LASH_EXH\t\t\t\t=\t" << LASH_EXH << " mm" << endl;
					cout << "\t\tFlow array file, EVFA_FILE\t\t=\t" << EVFA_FILE << endl;
					if(CONST_REF_AREA_EXH)	cout << "\t\tConstant ref. area, CONST_REF_AREA_EXH\t=\t" << TrueOrFalse(CONST_REF_AREA_EXH) << endl;
					else					cout << "\t\tCurtain ref. area, CONST_REF_AREA_EXH\t=\t" << TrueOrFalse(CONST_REF_AREA_EXH) << endl;
					cout << "\t\tFlow area multiplier, ONUM_EXH\t\t=\t" << ONUM_EXH;
					if(ONUM_EXH==1.0) cout << " (def.)"; cout << endl;
					cout << "\t\tRef. dia., REF_DIA_EXH\t\t\t=\t" << REF_DIA_EXH << " mm" << endl;
				}
				else {
					cout << "\t\tData are area values (m^2), EV_L_OR_A\t=\t" << TrueOrFalse(EV_L_OR_A) << endl;
					cout << "\t\tArea m'plier, AREA_MULT_EXH\t\t=\t" << AREA_MULT_EXH << endl;
				}
			}
		}
*/

		if (EV_ACTIVE) {
			pPpt->Out("\n");
			pPpt->Out(Underline("Exhaust AVT setup", "-", "\t\t"));
			pPpt->Out("\t\tActive exhaust valve, EV_ACTIVE\t\t=\t"); pPpt->Out(TrueOrFalse(EV_ACTIVE)); pPpt->Out("\n");
			pPpt->Out("\t\tTarget tapping in pipe, EV_target_pipe_num\t=\t");
			if (EV_TARGET_PIPE_EX) pPpt->Out("Exhaust["); else pPpt->Out("Intake[");
			pPpt->Out(EV_target_pipe_num); pPpt->Out("]\n");
			pPpt->Out("\t\t- at a fraction of, EV_target_pipe_loc\t=\t"); pPpt->Out(EV_target_pipe_loc); pPpt->Out("\n");
			pPpt->Out("\t\tTarget pressure file, EV_ACTIVE_FILE_P\t=\t"); pPpt->Out(EV_ACTIVE_FILE_P); pPpt->Out("\n");
			//pPpt->Out("\t\tTarget cycle average, target_cycle_average\t=\t"); pPpt->Out(this->target_cycle_average); pPpt->Out("\n");
		}
		pPpt->Out("\n");
	}	

	if(NINVALVES>0) {
		pPpt->Out(Underline("Intake valve parameters", "-", "\t\t"));
		pPpt->Out("\t\tReference diameter, REF_DIA_INT\t\t=\t"); pPpt->Out(REF_DIA_INT); pPpt->Out(" mm\n");
		switch (IV_WAVEFORM) {
		case 0:
			pPpt->Out("\t\tValve acts as switch, IV_WAVEFORM\t=\t"); pPpt->Out(IV_WAVEFORM); pPpt->Out("\n");
			if (IV_OPEN) pPpt->Out("\t\t- fully open, IV_OPEN\t\t\t=\t");
			else pPpt->Out("\t\t- fully closed, IV_OPEN\t\t\t=\t");
			pPpt->Out(TrueOrFalse(IV_OPEN)); pPpt->Out("\n");
			break;
		case 1:
			pPpt->Out("\t\tSquare wave opening, IV_WAVEFORM\t=\t"); pPpt->Out(IV_WAVEFORM); pPpt->Out("\n");
			break;
		case 2:
			pPpt->Out("\t\tTriangular wave opening, IV_WAVEFORM\t=\t"); pPpt->Out(IV_WAVEFORM); pPpt->Out("\n");
			break;
		case 3:
			pPpt->Out("\t\tSinusoidal wave opening, IV_WAVEFORM\t=\t"); pPpt->Out(IV_WAVEFORM); pPpt->Out("\n");
			break;
		case 4:
			pPpt->Out("\t\tFourier series wave opening, IV_WAVEFORM\t=\t"); pPpt->Out(IV_WAVEFORM); pPpt->Out("\n");
			break;
		case 5:
			pPpt->Out("\t\tFile wave opening, IV_WAVEFORM\t=\t"); pPpt->Out(IV_WAVEFORM); pPpt->Out("\n");
			break;
		default:
			pPpt->Out("\t\tUnknown wave opening, IV_WAVEFORM\t=\t"); pPpt->Out(IV_WAVEFORM); pPpt->Out(". Exiting.\n");
			exit(1);
			break;
		}

		// Switch on intake valve data units
		// ----------------------------------------------------------------------------------------------------
		switch (IV_L_OR_A) {
		case 0:
			pPpt->Out("\t\tValve data are lift (mm), IV_L_OR_A\t=\t"); pPpt->Out(IV_L_OR_A); pPpt->Out("\n");
			pPpt->Out("\t\tValve lash, LASH_INT\t\t\t=\t"); pPpt->Out(LASH_INT); pPpt->Out(" mm\n");
			pPpt->Out("\t\tFlow/discharge coeff vs l:d, IVFA_FILE\t=\t"); pPpt->Out(IVFA_FILE); pPpt->Out("\n");

			// Intake valve flow or discharge coefficients (can only be used with lift data since l/d is the lookup)
			// ----------------------------------------------------------------------------------------------------
			if (CONST_REF_AREA_INT) pPpt->Out("\t\tConst ref area, CONST_REF_AREA_INT\t=\t");
			else pPpt->Out("\t\tCurtain ref area, CONST_REF_AREA_INT\t=\t");
			pPpt->Out(TrueOrFalse(CONST_REF_AREA_INT)); pPpt->Out("\n");
			pPpt->Out("\t\tFlow/discharge coeff multplr, ONUM_INT\t=\t"); pPpt->Out(ONUM_INT);
			if (ONUM_INT == 1.0) pPpt->Out(" (def.)\n");

			break;
		case 1:
			pPpt->Out("\t\tProfile is flow area (m^2), IV_L_OR_A\t=\t"); pPpt->Out(IV_L_OR_A); pPpt->Out("\n");
			break;
		case 2:
			pPpt->Out("\t\tProfile is psi (valve:pipe), IV_L_OR_A\t=\t"); pPpt->Out(IV_L_OR_A); pPpt->Out("\n");
			break;
		default:
			//
			break;
		}

		if (IV_WAVEFORM != 5)
		{
			pPpt->Out("\t\tMin lift/area/psi value, IV_MIN_OPENING\t=\t"); pPpt->Out(IV_MIN_OPENING); pPpt->Out("\n");
			pPpt->Out("\t\tMax lift/area/psi value, IV_MAX_OPENING\t=\t"); pPpt->Out(IV_MAX_OPENING); pPpt->Out("\n");

			// If intake valve lift or area profile is square or triangular or sinusoidal or a Fourier series, i.e., IV_WAVEFORM == 1 or 2 or 3 or 4
			// ----------------------------------------------------------------------------------------------------
			if (IV_WAVEFORM == 1 || IV_WAVEFORM == 2 || IV_WAVEFORM == 3 || IV_WAVEFORM == 4) {
				pPpt->Out("\t\tBegin opening angle, IV_BEGIN_OPENING\t=\t"); pPpt->Out(IV_BEGIN_OPENING); pPpt->Out("\n");
		
				if (IV_END_or_phi == 0) {
					pPpt->Out("\t\tIV_END_CLOSING decides closing, IV_END_or_phi\t=\t"); pPpt->Out(IV_END_or_phi);
					pPpt->Out("\t\tEnd closing angle, IV_END_CLOSING\t=\t"); pPpt->Out(IV_END_CLOSING); pPpt->Out("\n");
				}
				else {
					if (IV_END_or_phi == 1) {
						pPpt->Out("\t\tIV_phi decides closing, IV_END_or_phi\t=\t"); pPpt->Out(IV_END_or_phi); pPpt->Out("\n");
						pPpt->Out("\t\tOpen fraction of wavelength (0-1), IV_phi\t=\t"); pPpt->Out(IV_phi); pPpt->Out("\n"); pPpt->Out("\n");
						//	double phi_orig;			// phi(0-1) based on original excitation
					}
				}
				pPpt->Out("\n");

				pPpt->Out("\t\tEnd closing angle, IV_END_OPENING\t=\t"); pPpt->Out(IV_END_CLOSING); pPpt->Out("\n");
				pPpt->Out("\t\tOpen fraction of wavelength (0-1), IV_phi\t=\t"); pPpt->Out(IV_phi); pPpt->Out("\n");
			}
			else
			{
				if (IV_WAVEFORM == 0)
				{ 
					// 
				}
				else
				{
					pPpt->Out("\t\tUnknown wave opening, IV_WAVEFORM\t=\t"); pPpt->Out(IV_WAVEFORM); pPpt->Out(". Exiting.\n");
					exit(1);
				}
			}
		}
		else 
		{
			// If intake valve lift or area profile is read from file, i.e., IV_WAVEFORM == 5
			// ----------------------------------------------------------------------------------------------------
			pPpt->Out("\t\tValve timing file, IV_FILE\t\t=\t"); pPpt->Out(IV_FILE); pPpt->Out("`n");
			pPpt->Out("\t\tCam timing angle, CAM_TNG_ANG_INT\t=\t"); pPpt->Out(CAM_TNG_ANG_INT); pPpt->Out(Deg()); pPpt->Out("\n");
			pPpt->Out("\t\t- ");
			switch (CAM_TNG_LOC_INT) {
			case 0:
				pPpt->Out("start of opening");
				break;
			case 1:
				pPpt->Out("maximum opening");
				break;
			case 2:
				pPpt->Out("end of closing");
				break;
			default:
				CAM_TNG_LOC_INT = 1;
				pPpt->Out("maximum opening");
				break;
			}
			pPpt->Out(", CAM_TNG_LOC_INT\t=\t"); pPpt->Out(CAM_TNG_LOC_INT); pPpt->Out("\n");
		}

		/*
		// If intake valve data are lift values, i.e., IV_L_OR_A == 1 == true
		// ----------------------------------------------------------------------------------------------------
		if (IV_L_OR_A) {
			pPpt->Out("\t\tValve data are lift (mm), IV_L_OR_A\t=\t"); pPpt->Out(TrueOrFalse(IV_L_OR_A)); pPpt->Out("\n");
			//			cout << "\t\tLift m'plier, LIFT_MULT_INT\t\t=\t" << LIFT_MULT_INT << endl;
			pPpt->Out("\t\tValve lash, LASH_INT\t\t\t=\t"); pPpt->Out(LASH_INT); pPpt->Out(" mm\n");
			pPpt->Out("\t\tFlow/discharge coeff vs l:d, IVFA_FILE\t=\t"); pPpt->Out(IVFA_FILE); pPpt->Out("\n");
		}
		else {
			pPpt->Out("\t\tValve data are area (m^2), IV_L_OR_A\t=\t"); pPpt->Out(TrueOrFalse(IV_L_OR_A)); pPpt->Out("\n");
			//			cout << "\t\tArea m'plier, AREA_MULT_INT\t\t=\t" << AREA_MULT_INT << endl;
		}

		// Intake valve flow or discharge coefficients (can be used with either lift or area values)
		// ----------------------------------------------------------------------------------------------------
		if (CONST_REF_AREA_INT) pPpt->Out("\t\tConst ref area, CONST_REF_AREA_INT\t=\t");
		else pPpt->Out("\t\tCurtain ref area, CONST_REF_AREA_INT\t=\t");
		pPpt->Out(TrueOrFalse(CONST_REF_AREA_INT)); pPpt->Out("\n");
		pPpt->Out("\t\tFlow/discharge coeff multplr, ONUM_INT\t=\t"); pPpt->Out(ONUM_INT);
		if (ONUM_INT == 1.0) pPpt->Out(" (def.)\n");
*/
		// Variable valve timing and lift
		// ----------------------------------------------------------------------------------------------------
		pPpt->Out("\t\tOpening multiplier, IV_OPENING_MULT\t=\t"); pPpt->Out(IV_OPENING_MULT); pPpt->Out("\n");
		pPpt->Out("\t\tOpen duration multiplier, ANG_MULT_INT\t=\t"); pPpt->Out(ANG_MULT_INT); pPpt->Out("\n");
		pPpt->Out("\t\tOpen duration compress, scale_factor_iv\t=\t"); pPpt->Out(scale_factor_iv); pPpt->Out("\n");
/*
		if(IV_SWITCH) {
			cout << "\t\tIntake valve acts as switch, IV_SWITCH\t=\t" << TrueOrFalse(IV_SWITCH) << endl;
			if(IV_OPEN) cout << "\t\t- fully open, IV_OPEN\t\t\t=\t" << TrueOrFalse(IV_OPEN) << endl;
			else cout << "\t\t- fully closed, IV_OPEN\t\t\t=\t" << TrueOrFalse(IV_OPEN) << endl;
		}
		else {
			if(IV_CALC) {
				cout << "\t\tCalculate valve timing, IV_CALC\t\t=\t" << TrueOrFalse(IV_CALC) << endl;
				cout << "\t\tSquare wave opening, IV_SQUARE\t\t=\t" << TrueOrFalse(IV_SQUARE) << endl;
				cout << "\t\tBegin valve opening, IV_BEGIN\t\t=\t" << IV_BEGIN << endl;
				cout << "\t\tEnd valve closing, IV_END\t\t=\t" << IV_END << endl;
				cout << "\t\tMaximum psi value, IV_PSI_MAX\t\t=\t" << IV_PSI_MAX << endl;
				cout << "\t\tMinimum psi value, IV_PSI_MIN\t\t=\t" << IV_PSI_MIN << endl;
			}
			else {
				cout << "\t\tInterpolate valve timing, IV_CALC\t=\t" << TrueOrFalse(IV_CALC) << endl;
				cout << "\t\tValve timing file, IV_FILE\t\t=\t" << IV_FILE << endl;
				cout << "\t\tCam timing angle, CAM_TNG_ANG_INT\t=\t" << CAM_TNG_ANG_INT << Deg() << endl;
				cout << "\t\t- ";
				switch(CAM_TNG_LOC_INT) {
				case 0:
					cout << "start of opening";
					break;
				case 1:
					cout << "maximum opening";
					break;
				case 2:
					cout << "end of closing";
					break;
				default:
					CAM_TNG_LOC_INT = 1;
					cout << "maximum opening";
					break;
				}
				cout << ", CAM_TNG_LOC_INT\t=\t" << CAM_TNG_LOC_INT << endl;
				cout << "\t\tAngle multiplier, ANG_MULT_INT\t\t=\t" << ANG_MULT_INT << endl;
				if(IV_L_OR_A) {
					cout << "\t\tData are lift values (mm), IV_L_OR_A\t=\t" << TrueOrFalse(IV_L_OR_A) << endl;
					cout << "\t\tLift m'plier, LIFT_MULT_INT\t\t=\t" << LIFT_MULT_INT << endl;
					cout << "\t\tLash, LASH_INT\t\t\t\t=\t" << LASH_INT << " mm" << endl;
					if(CONST_REF_AREA_INT)	cout << "\t\tConstant ref. area, CONST_REF_AREA_INT\t=\t" << TrueOrFalse(CONST_REF_AREA_INT) << endl;
					else					cout << "\t\tCurtain ref. area, CONST_REF_AREA_INT\t=\t" << TrueOrFalse(CONST_REF_AREA_INT) << endl;
					cout << "\t\tFlow area multiplier, ONUM_INT\t\t=\t" << ONUM_INT;
					if(ONUM_INT==1.0) cout << " (def.)"; cout << endl;
					cout << "\t\tRef. dia., REF_DIA_INT\t\t\t=\t" << REF_DIA_INT << " mm" << endl;
				}
				else {
					cout << "\t\tData are area values (m^2), IV_L_OR_A\t=\t" << TrueOrFalse(IV_L_OR_A) << endl;
					cout << "\t\tArea m'plier, AREA_MULT_INT\t\t=\t" << AREA_MULT_INT << endl;
				}
			}
		}	
		pPpt->Out("\n");
*/

		if (IV_ACTIVE) {
			pPpt->Out(Underline("Intake AVT setup", "-", "\t\t"));
			pPpt->Out("\t\tActive intake valve, IV_ACTIVE\t\t=\t"); pPpt->Out(TrueOrFalse(IV_ACTIVE)); pPpt->Out("\n");
			pPpt->Out("\t\tTarget tapping in pipe, IV_target_pipe_num\t=\t");
			if (IV_TARGET_PIPE_EX) pPpt->Out("Exhaust["); else pPpt->Out("Intake[");
			pPpt->Out(IV_target_pipe_num); pPpt->Out("]\n");
			pPpt->Out("\t\t- at a fraction of, IV_target_pipe_loc\t=\t"); pPpt->Out(IV_target_pipe_loc); pPpt->Out("\n");
			pPpt->Out("\t\tTarget pressure file, IV_ACTIVE_FILE_P\t=\t"); pPpt->Out(IV_ACTIVE_FILE_P); pPpt->Out("\n");
			//pPpt->Out("\t\tTarget cycle average, target_cycle_average\t=\t"); pPpt->Out(this->target_cycle_average); pPpt->Out("\n");
		}
		pPpt->Out("\n");
	}
	pPpt->Out("\n");
}

void CEngine::PrintToScreen()
{
	cout << Underline(Identify(), "=");
	cout << "Rev no.\t\t\t\t\t\t=\t" << rev << "\n";
	cout << "Rotation\t\t\t\t\t=\t" << rotation << Deg() << "CA\n";
	cout << "Crank angle\t\t\t\t\t=\t" << ca << Deg() << "CA\n";
	cout << "Elapsed ca\t\t\t\t\t=\t" << ca_elapsed << Deg() << "CA\n";
	cout << "Degrees since start\t\t\t\t=\t" << degrees_elapsed << Deg() << "CA" << endl;;
	cout << endl;
	cout << Underline("Overall engine work and power", "-");
	cout << "Work done most recent iteration (del_WK_total)\t=\t" << del_WK_total << " J\n";
	cout << "Instantaneous power (P_inst)\t\t\t=\t" << P_inst << " W\n";
	cout << "Current cycle work (CYCLE_WK)\t\t\t=\t" << CYCLE_WK << " J\n";
	cout << "Previous cycle work (PREV_CYCLE_WK)\t\t=\t" << PREV_CYCLE_WK << " J\n";
	cout << "Previous cycle power (PREV_CYCLE_POWER)\t\t=\t" << PREV_CYCLE_POWER << " W\n";
	cout << endl;
	if(MODEL==WATSON)
	{
		cout << Underline("Fuel Economy", "-");
		cout << "Current cycle fuel inj. (CYCLE_FUEL)\t\t=\t" << CYCLE_FUEL << " kg\n";
		cout << "Previous cycle fuel inj. (PREV_CYCLE_FUEL)\t=\t" << PREV_CYCLE_FUEL << " kg\n";
		cout << "Previous cycle thermal efficiency (eta)\t\t=\t" << eta << " \n";
		cout << "Previous cycle s.f.c. (PREV_CYCLE_SFC)\t\t=\t" << PREV_CYCLE_SFC << " g/kWh\n";
		cout << endl;
	}
}

char* CEngine::Identify()
// ============================================================ //
// Returns identification of the current object					//
// ============================================================ //
{
	std::string sz;
	sz = "Assembly [";
	sz += IntToString(AssyID);
	sz += "], ";
	sz += "Engine [";
	sz += IntToString(ID);
	sz += "]";

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str());
	return szz;
}
