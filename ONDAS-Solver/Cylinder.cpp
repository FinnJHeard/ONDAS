/// Cylinder.cpp: implementation of the CCylinder class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Cylinder.h"
#include "Engine.h"
#include "Tools.h"
#include "Valve.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCylinder::CCylinder()
{
	
}

CCylinder::~CCylinder()
{

}

// Copy constructor
CCylinder::CCylinder(const CCylinder& inCyl)
{
	// Configuration
	// =============
	ID = inCyl.ID;
//	AREF = inCyl.AREF;
//	npipes = inCyl.npipes;

//	for(int p=0; p<2; ++p)
//	{
//		end[p] = inCyl.end[p];
//		pPipe[p] = inCyl.pPipe[p];
//		pBN[p] = inCyl.pBN[p]; // Valves now
//		pCLIN[p] = inCyl.pCLIN[p];
//		pCLOUT[p] = inCyl.pCLOUT[p];
//		pend_flow[p] = inCyl.pend_flow[p];
//		pPathLine[p] = inCyl.pPathLine[p];
//	}

	ExhaustValve = inCyl.ExhaustValve;
	IntakeValve = inCyl.IntakeValve;
	
	// Read from file
	// =========
//	char** labels;	// Parameter list of labels in no order
//	double* values;	// Parameter list of values in same order as labels

	offset = inCyl.offset;
	
	// Other
	// =====
	TC = inCyl.TC;
	CONRAT = inCyl.CONRAT;
	FCYL = inCyl.FCYL;
	THETA = inCyl.THETA;

	
	FNN = inCyl.FNN;
	X = inCyl.X;
	VC = inCyl.VC;
	PC = inCyl.PC;
	AC = inCyl.AC;
	MC = inCyl.MC;
	DXDT = inCyl.DXDT;
	DTHDT = inCyl.DTHDT;
	DXDTH = inCyl.DXDTH;
	DVCDT = inCyl.DVCDT;
	DPCDT = inCyl.DPCDT;

//	NCYCLES = inCyl.NCYCLES;
	
//	DMEDT = inCyl.DMEDT;
//	DMIDT = inCyl.DMIDT;

//	DT = inCyl.DT;

//	a = inCyl.a;
//	b = inCyl.b;
//	Tw = inCyl.Tw;

	WAIT = inCyl.WAIT;

	datapoints_cyl = inCyl.datapoints_cyl;
//	cyl_data = inCyl.cyl_data
}

CCylinder& CCylinder::operator=(const CCylinder& inCyl)
{
	if(this != &inCyl)
	{
		// Configuration
		// =============
		ID = inCyl.ID;
//		AREF = inCyl.AREF;
//		npipes = inCyl.npipes;

//		for(int p=0; p<2; ++p)
//		{
//			end[p] = inCyl.end[p];
//			pPipe[p] = inCyl.pPipe[p];
//			pBN[p] = inCyl.pBN[p]; // Valves now	
//			pCLIN[p] = inCyl.pCLIN[p];
//			pCLOUT[p] = inCyl.pCLOUT[p];
//			pend_flow[p] = inCyl.pend_flow[p];
//			pPathLine[p] = inCyl.pPathLine[p];
//		}
	
		ExhaustValve = inCyl.ExhaustValve;
		IntakeValve = inCyl.IntakeValve;

		offset = inCyl.offset;
	
		// Other
		// =====
		TC = inCyl.TC;
		CONRAT = inCyl.CONRAT;
		FCYL = inCyl.FCYL;
		THETA = inCyl.THETA;
	
		FNN = inCyl.FNN;
		X = inCyl.X;
		VC = inCyl.VC;
		PC = inCyl.PC;
		AC = inCyl.AC;
		MC = inCyl.MC;
		DXDT = inCyl.DXDT;
		DTHDT = inCyl.DTHDT;
		DXDTH = inCyl.DXDTH;
		DVCDT = inCyl.DVCDT;
		DPCDT = inCyl.DPCDT;

//		NCYCLES = inCyl.NCYCLES;
	
//		DMEDT = inCyl.DMEDT;
//		DMIDT = inCyl.DMIDT;

//		DT = inCyl.DT;

//		a = inCyl.a;
//		b = inCyl.b;
//		Tw = inCyl.Tw;

		WAIT = inCyl.WAIT;

		datapoints_cyl = inCyl.datapoints_cyl;
//		cyl_data = inCyl.cyl_data;
	}
	return *this;
}

void CCylinder::Initialise(CProperties* pPpt, 
						   CPipe** pPipes,
						   CPipe* &rExhaustPipe, /*CPipe* &rIntakePipe,*/ int** &rCYLPIPES, int** &rCYLPIPES_ENDS, double* &rENDCORR,
							int id, bool ex, int npipes, CEngine* EngPtr, std::string param_dir/*, char* cyl_dir*//*, char* vt_dir*/, 
							int assyid, string parent_assy_res_dir, string calling_object_str, int nExPipesInAssy, int nInPipesInAssy)
{
	if (pPpt->SHOW_calls) { pPpt->Out("CCylinder.Initialise ("); pPpt->Out(calling_object_str); pPpt->Out(")\n"); }

	// General initialization of exhaust side
	InitialiseGen(pPpt, pPipes, rExhaustPipe, rCYLPIPES, rCYLPIPES_ENDS, rENDCORR, id, ex, npipes, assyid, "CCylinder", parent_assy_res_dir);

	// General initialization of intake side
	// NOT PRESENT YET

	LoadCylinderData(ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->CYL_DIR), EngPtr->cyl_file));
	
	// Apply cylinder properties
	EnginePtr = EngPtr;
	ID = id;
	AssyID = assyid;				// ID of assembly on which boundary belongs
	
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Initialise\n");}

	std::string bcname_str = "CYLINDER";
	ReadInput(pPpt, ConstructString(pPpt, param_dir, bcname_str, ID));

	// Boundary name
	//int p;
	//NPIPES = npipes;
	//NAME = new int[NPIPES];
	//for (p = 0; p < NPIPES; ++p) NAME[p] = CYLINDER;

	bool intake_fitted;

	if(EnginePtr->IPIPE) intake_fitted = true; // Intake pipe fitted
	else intake_fitted = false;

	THETA = this->EnginePtr->ca + offset;
	double THETAR = THETA*PI/180.0;

	PC = EnginePtr->pcr;
	TC = EnginePtr->Tcr;
	EVO_SET = false;

	NCYCLES = 0;
	W_cig = 0;
	W_cin = 0;

	if(!EnginePtr->IVOL){	// Constant volume cylinder, read volume, PC, TC
		if(EnginePtr->MODEL==CONST_P){
			PC = EnginePtr->pres;
			TC = EnginePtr->Tres;
			VC = EnginePtr->Vres;
		}
		else{
			if(EnginePtr->MODEL==CONST_VOL){
				VC = EnginePtr->vconst;
				PC = EnginePtr->pinit;
				TC = EnginePtr->Tinit;
			}
			else{
				PC = EnginePtr->pcr;
				TC = EnginePtr->Tcr;
			}
		}
		AC = sqrt(pPpt->gammaAir(TC)*287.0*TC);
		MC = 1.0E5*PC*VC/(287.0*TC);
		DVCDT = 0;
		DPCDT = -pPpt->gammaAir(TC)*DVCDT*PC/VC;
	}
	else{
		if(EnginePtr->MODEL==CONST_P){
			PC = EnginePtr->pres; 
			// Allow TC to vary according to variable volume
		}
		// Variable volume cylinder
		// Initialise values
		CONRAT = 2.0*EnginePtr->conrod/EnginePtr->stroke;
		FNN = sqrt(pow(CONRAT, 2) - pow(sin(THETAR), 2));
		X = 0.5*EnginePtr->stroke*(1.0 + CONRAT - FNN - cos(THETAR));
		FCYL = 0.25*PI*pow(EnginePtr->dcyl, 2);
		VC = FCYL*(X + EnginePtr->stroke/(EnginePtr->cr - 1.0));
		AC = sqrt(pPpt->gammaAir(TC)*287.0*TC);
		MC = 1.0E5*PC*VC/(287.0*TC);
		DXDT = 0.5*EnginePtr->stroke*sin(THETAR)*(1.0 + cos(THETAR)/FNN);
		DTHDT = 2.0*PI*EnginePtr->reveng;
		DVCDT = FCYL*DXDT*DTHDT;
		DPCDT = -pPpt->gammaAir(TC)*DVCDT*PC/VC;
	}

	if(THETA>EnginePtr->EVO) WAIT = true;
	else WAIT = false;
/*
	// Mean entropy levels (used in homentropic simulations only)
	sum_mc = 0;

//	aA_1_mean = 0;
	sum_mc_aA_1 = 0;
	aA_1 = 0;

//	aA_2_mean = 0;
	sum_mc_aA_2 = 0;
	aA_2 = 0;
	P1 = 0;

	calculate_entropy = false;
*/
	
	int i;

	// Setup exhaust valve(s)
	// ----------------------
	CPipe** pPipesExhaust;
	int** EVPIPES; int** EVPIPES_ENDS; double* EVENDCORR;
	pPipesExhaust = new CPipe* [EnginePtr->NEXVALVES];
	EVPIPES = new int* [EnginePtr->NEXVALVES]; EVPIPES_ENDS = new int* [EnginePtr->NEXVALVES]; EVENDCORR = new double [EnginePtr->NEXVALVES];	
	
	for(i=0; i<EnginePtr->NEXVALVES; ++i) {EVPIPES[i] = new int [1]; EVPIPES_ENDS[i] = new int [1];}
	
	for(i=0; i<EnginePtr->NEXVALVES; ++i){
		pPipesExhaust[i] = pPipes[i];
		EVPIPES[i][ONE_SIDE] = rCYLPIPES[ID][i];
		EVPIPES_ENDS[i][ONE_SIDE] = rCYLPIPES_ENDS[ID][i];
		EVENDCORR[i] = rENDCORR[ID];
	}
	ExhaustValve = new CValve [EnginePtr->NEXVALVES];
	for(i=0; i<EnginePtr->NEXVALVES; ++i)
		ExhaustValve[i].Initialise(pPpt, pPipesExhaust, rExhaustPipe, EVPIPES, EVPIPES_ENDS, EVENDCORR, i, true, 1, EngPtr, param_dir/*, vt_dir*/, 
			AssyID, ID, parent_assy_res_dir, nExPipesInAssy);

	DMEDT = new double [EnginePtr->NEXVALVES];
	for(i=0; i<EnginePtr->NEXVALVES; ++i) DMEDT[i] = 0.0;
/*
	// Setup intake valve(s)
	// ---------------------
	CPipe** pPipesIntake;
	int** IVPIPES; int** IVPIPES_ENDS; double* IVENDCORR; 
	pPipesIntake = new CPipe* [EnginePtr->NINVALVES];
	IVPIPES = new int* [EnginePtr->NINVALVES]; IVPIPES_ENDS = new int* [EnginePtr->NINVALVES]; IVENDCORR = new double [EnginePtr->NINVALVES];
	for(i=0; i<EnginePtr->NINVALVES; ++i)
	{IVPIPES[i] = new int [1]; IVPIPES_ENDS[i] = new int [1];}
	for(i=0; i<EnginePtr->NINVALVES; ++i)
	{
		pPipesIntake[i] = pPipes[i+EnginePtr->NEXVALVES];
		IVPIPES[i][ONE_SIDE] = rCYLPIPES[ID][i+EnginePtr->NEXVALVES];
		IVPIPES_ENDS[i][ONE_SIDE] = rCYLPIPES_ENDS[ID][i+EnginePtr->NEXVALVES];
		IVENDCORR[i] = rENDCORR[ID];
	}
	
	IntakeValve = new CValve [EnginePtr->NINVALVES];
	for(i=0; i<EnginePtr->NINVALVES; ++i)
		IntakeValve[i].Initialise(pPpt, pPipesIntake,
								  rIntakePipe, IVPIPES, IVPIPES_ENDS, IVENDCORR, i, false, 1, EngPtr, param_dir, AssyID, ID, parent_assy_res_dir);
	
	DMIDT = new double [EnginePtr->NINVALVES];
	for(i=0; i<EnginePtr->NINVALVES; ++i) DMIDT[i] = 0.0;
*/

/*	Now done as valve data is read
	// Calculate IVC
	// -------------
	EnginePtr->IVC = 360; // Default value
	if(EnginePtr->IAIR)
	{	for(int i=0; i<EnginePtr->NINVALVES; ++i)
		{	if(EnginePtr->IVO + IntakeValve[i].Get_VC() > EnginePtr->IVC) EnginePtr->IVC = EnginePtr->IVO + IntakeValve[i].Get_VC();}}
*/

	// Motored values
	// ==============
	if(EnginePtr->HT_MODEL == WOSCHNI) 
	{
		PC_mot = EnginePtr->pcr;
		TC_mot = EnginePtr->Tcr;

		if(!EnginePtr->IVOL)	// Constant volume cylinder, read volume, PC, TC
		{
			PC_mot = EnginePtr->pcr;
			TC_mot = EnginePtr->Tcr;
			AC_mot = sqrt(pPpt->gammaAir(TC)*287.0*TC);
		}
		else
		{	// Variable volume cylinder
			AC_mot = sqrt(pPpt->gammaAir(TC)*287.0*TC_mot);
			MC_mot = 1.0E5*PC_mot*VC/(287.0*TC_mot);
			DPC_motDT = -pPpt->gammaAir(TC)*DVCDT*PC_mot/VC;	
		}

		DMEDT_mot = new double [EnginePtr->NEXVALVES];
		for(i=0; i<EnginePtr->NEXVALVES; ++i) DMEDT_mot[i] = 0.0;

		DMIDT_mot = new double [EnginePtr->NINVALVES];
		for(i=0; i<EnginePtr->NINVALVES; ++i) DMIDT_mot[i] = 0.0;

		MOT_SET = false;
		PC_mot_ref = PC_mot; TC_mot_ref = TC; VC_mot_ref = VC;
	}

	// Combustion
	// ==========
	COMB_RESET = false;
	INJ_PREV_ITER = false;
	m_finj = 0;		// Fuel mass injected so far this cycle (kg)
	m_fub = 0;		// Mass of unburnt fuel in cylinder (kg)
	m_fb = 0;		// Mass of fuel burned so far this cycle (kg)
	del_m_fb = 0;	// Mass of fuel burned this iteration (kg)
	QLHVrate = 0;	// Combustion heat release rate (J/s)
	phiFAoverall = 1;// Overall fuel/air equivalence ratio

	// Zero cycle variables
	W_cig = 0;
	W_cig_prev = 0;
	W_cin = 0;
	W_cin_prev = 0;

	W_cig_mot = 0;
	W_cig_mot_prev = 0;
	W_cin_mot = 0;
	W_cin_mot_prev = 0;

	SetupFiles(pPpt, parent_assy_res_dir);
}

void CCylinder::Configure(CProperties* pPpt, CPipe** pExhaustSystem, CPipe** pIntakeSystem, int nExPipes, int nInPipes)
{
	int i;
	// Configure valves
	for(i=0; i<EnginePtr->NEXVALVES; ++i)
		ExhaustValve[i].Configure(pPpt);

	if(EnginePtr->IAIR && EnginePtr->IPIPE)
		for(i=0; i<EnginePtr->NINVALVES; ++i)
			IntakeValve[i].Configure(pPpt);

	for(i=0; i<EnginePtr->NEXVALVES; ++i)
		ExhaustValve[i].ConfigureExtra(pPpt, pExhaustSystem, pIntakeSystem, nExPipes, nInPipes);

	if(EnginePtr->IAIR && EnginePtr->IPIPE)
		for(i=0; i<EnginePtr->NINVALVES; ++i)
			IntakeValve[i].ConfigureExtra(pPpt, pExhaustSystem, pIntakeSystem, nExPipes, nInPipes);
}

void CCylinder::RunBoundary(CProperties* pPpt, CTime &rMyTime, double DELZe, double DELZi, double ca, double del_ca, int timestep, double time)
{
	Update(rMyTime, DELZe, DELZi, pPpt, /*ca, del_ca,*/ timestep, time);
}

void CCylinder::Update(CTime &rMyTime, double DELZe, double DELZi, CProperties* pPpt, /*double ca, double del_ca,*/ int timestep, double time)
{
	int i;
	// Check and set up time step
	double DT;
	if(EnginePtr->IPIPE)
	{
		if( (DELZe*pPpt->xref/pPpt->AREFe) != (DELZi*pPpt->xref/pPpt->AREFi) ) {cout << "Cylinder[" << ID << "] CCylinder::Update: exhaust and intake derived time steps do not match" << endl; exit(1);}
		else DT = DELZe*pPpt->xref/pPpt->AREFe; // Could use either though
	}
	else DT = DELZe*pPpt->xref/pPpt->AREFe; // No intake so use exhaust values

	THETA += EnginePtr->del_theta;		// Update member THETA from engine

	// Test for start of new cycle
	// ===========================
	if(THETA>180.0*EnginePtr->cycle)	// Test for start of new cycle 
	{
		++NCYCLES;
		THETA = THETA - 180.0*EnginePtr->cycle;	// Reset THETA

		// Print cycle info to file
		fprintf(FILE_CYCLE,"%f\t%f\t%f\t%f\n", time, EnginePtr->ca_elapsed, THETA, W_cin);
	
		// Record and reset cycle variables
		W_cig_prev = W_cig;
		W_cig = 0;
		W_cin_prev = W_cin;
		W_cin = 0;

		W_cig_mot_prev = W_cig_mot;
		W_cig_mot = 0;
		W_cin_mot_prev = W_cin_mot;
		W_cin_mot = 0;
		
		WAIT = false;					// No need to WAIT as we have started a new cycle from the beginning

		switch(EnginePtr->MODEL)
		{
		case CONST_P:
			break;
		case CONST_VOL:
			break;
		case READ_P_T:
			break;
		case RESET_P_T:
			if(!EnginePtr->RESET_EVO)	// If we reset conditions at TDC and not EVO
			{
				// Reset release conditions for TDC
				double TDC_RAD = 0;		// TDC in radians
				double IVC_RAD = (EnginePtr->IVO + IntakeValve[0].Get_VC())*PI/180.0; // VC should be same so use [0]
				FNN = sqrt(pow(CONRAT, 2) - pow(sin(TDC_RAD), 2)); X = 0.5*EnginePtr->stroke*(1.0 + CONRAT - FNN - cos(TDC_RAD));
				double VCC_TDC = FCYL*(X + EnginePtr->stroke/(EnginePtr->cr - 1.0));
				FNN = sqrt(pow(CONRAT, 2) - pow(sin(IVC_RAD), 2)); X = 0.5*EnginePtr->stroke*(1.0 + CONRAT - FNN - cos(IVC_RAD));
				double VCC_IVC = FCYL*(X + EnginePtr->stroke/(EnginePtr->cr - 1.0));			
				TC = (EnginePtr->pcr*VCC_TDC)*TC/(PC*VCC_IVC); // Conservation of mass in cylinder during closed valves period
				PC = EnginePtr->pcr;
			}
			else EVO_SET = false;		// Flag to reset conditions once at EVO
			break;
		case WATSON:
			COMB_RESET = false;			// Since new cycle, THETA not > IVC so can now flag the reset
			break;
		default:
			EnginePtr->MODEL = RESET_P_T;
			cout << "CCylinder::Update [" << ID << "]: Unknown cylinder model; selecting RESET_P_T\n";
			break;
		}

		if(EnginePtr->HT_MODEL==WOSCHNI) MOT_SET = false;				// Since new cycle, THETA not > IVC so can now flag the reset
	}

	// Calculate cylinder pressure
	// ===========================
	double DQDT = 0; double DQDT_mot = 0;
	switch(EnginePtr->MODEL)
	{
	case CONST_P:
		PC = EnginePtr->pres; 
		if(!EnginePtr->IVOL) TC = EnginePtr->Tres;	// TC must be allowed to vary if variable volume cylinder
		break;
	case CONST_VOL:
		PC = PC + DPCDT*DT;							// Calculate new pressure due to mass flows (no volume change)
		break;
	case READ_P_T:
		Interpolate(pPpt, THETA);	// Sets PC, TC
		break;
	case RESET_P_T:
		PC = PC + DPCDT*DT; // For both RESET_P and WATSON cylinder models
		if(EnginePtr->RESET_EVO && THETA>=EnginePtr->EVO && !EVO_SET)// If we reset conditions at EVO and not TDC
		{
//cout << "EVO\n\n\n";
			// Reset release conditions
			double EVO_RAD = EnginePtr->EVO*PI/180.0;		// EVO in radians
			double IVC_RAD = (EnginePtr->IVO + IntakeValve[0].Get_VC())*PI/180.0; // VC should be same so use [0]
			FNN = sqrt(pow(CONRAT, 2) - pow(sin(EVO_RAD), 2)); X = 0.5*EnginePtr->stroke*(1.0 + CONRAT - FNN - cos(EVO_RAD));
			double VCC_EVO = FCYL*(X + EnginePtr->stroke/(EnginePtr->cr - 1.0));
			FNN = sqrt(pow(CONRAT, 2) - pow(sin(IVC_RAD), 2)); X = 0.5*EnginePtr->stroke*(1.0 + CONRAT - FNN - cos(IVC_RAD));
			double VCC_IVC = FCYL*(X + EnginePtr->stroke/(EnginePtr->cr - 1.0));			
			TC = (EnginePtr->pcr*VCC_EVO)*TC/(PC*VCC_IVC); // Conservation of mass in cylinder during closed valves period
			PC = EnginePtr->pcr;
// Fix
if(NCYCLES<=1)
{
PC = EnginePtr->pcr;
TC = EnginePtr->Tcr;
MC = ((PC*1e5)/(pPpt->R_air*TC))*VCC_EVO;
}
//
			EVO_SET = true;
		}	
		break;
	case WATSON:
		{
			PC = PC + DPCDT*DT; // For both RESET_P and WATSON cylinder models
			// Reset combustion variables prior to combustion			
			if(THETA>EnginePtr->IVC && AllClosed() && !COMB_RESET) // Double check all valves are closed
			{
				double phiFAactual = EnginePtr->m_f0/MC;			// Assumes cylinder mass at IVC is only air
				phiFAoverall = phiFAactual/EnginePtr->phiFAstoich;	// Calculate overall fuel/air equivalence ratio
				INJ_PREV_ITER = false;	m_finj = 0;	m_fub = 0; m_fb = 0; del_m_fb = 0; QLHVrate = 0;							
				COMB_RESET = true;
			}
			DQDT += CombustionDQDT(pPpt, DT); // Calculate combustion heat release
		}
		break;
	default:
		EnginePtr->MODEL = RESET_P_T;
		cout << "CCylinder::Update [" << ID << "]: Unknown cylinder model; selecting RESET_P_T\n";
		break;
	}

	// Cylinder heat transfer
	// ======================
	switch(EnginePtr->HT_MODEL)
	{
	case NONE:
		DQDT += 0;
		break;
	case ANNAND:
		{
			double SPBAR = 2*EnginePtr->stroke*EnginePtr->reveng;	// Mean piston speed
			double RC_temp =  MC/VC;
			double Re_temp = RC_temp*EnginePtr->dcyl*SPBAR/pPpt->ViscosityAir(TC);
			double kq = pPpt->cpAir(TC)*pPpt->ViscosityAir(TC)/0.7; // Conductivity (Pr = 0.7)                                                  
			
			// Calculate perimetered/walled area of the gas
			double ROOF = FCYL;	// Let chamber roof area = piston area
			double FW = ((X + EnginePtr->stroke/(EnginePtr->cr - 1.0))*PI*EnginePtr->dcyl) + FCYL + ROOF;		
			// Total area = cylindrical walls + piston surface area + cylinder head area
			
			// Annand expression for heat transfer rate		
			DQDT += (EnginePtr->a*pow(Re_temp, EnginePtr->b)/EnginePtr->dcyl)*kq*FW*(EnginePtr->Tw - TC);
		}
		break;
	case WOSCHNI:
		{
			PC_mot = PC_mot + DPC_motDT*DT; // For both RESET_P and WATSON cylinder models

			double C1, C2;
			if(EnginePtr->ign_ca>360) // Ignition before TDC (arbitrary dividing point)
			{
				if(THETA>=EnginePtr->IVC && THETA<EnginePtr->ign_ca)
				{C1 = 2.28; C2 = 0;} // Compression period
				else
				{
					if(THETA>=EnginePtr->EVO && THETA<EnginePtr->IVC)
					{C1 = 6.18; C2 = 0;} // Gas exchange period
					else {C1 = 2.28; C2 = 3.24e-3;} // Combustion & expansion period
				}
			}
			else // Ignition after TDC
			{
				if(THETA>=EnginePtr->IVC || THETA<EnginePtr->ign_ca)
				{C1 = 2.28; C2 = 0;} // Compression period
				else
				{
					if(THETA>=EnginePtr->EVO && THETA<EnginePtr->IVC)
					{C1 = 6.18; C2 = 0;} // Gas exchange period
					else {C1 = 2.28; C2 = 3.24e-3;} // Combustion & expansion period
				}
			}
		
			double SPBAR = 2*EnginePtr->stroke*EnginePtr->reveng;	// Mean piston speed
			double T_g = TC; // Mean gas temperature
			double T_g_mot = TC_mot; // Motored mean gas temperature
			double Vd = EnginePtr->stroke*FCYL;						// Swept (displaced) volume
			double w = C1*SPBAR + C2*Vd*TC_mot_ref/(PC_mot_ref*VC_mot_ref)*(PC - PC_mot);	
			// Local average gas velocity in cylinder (pressure units don't matter - cancel out)
			double w_mot = C1*SPBAR; // Motored average gas velocity  - putting PC = PC_mot gives this

			double h_cg = 3.26*pow(EnginePtr->dcyl, -0.2)*pow(PC*100, 0.8)*pow(T_g, -0.55)*pow(w, 0.8);	// Woschni's correlation; put pressure into kPa
			double h_cg_mot = 3.26*pow(EnginePtr->dcyl, -0.2)*pow(PC_mot*100, 0.8)*pow(T_g_mot, -0.55)*pow(w_mot, 0.8);	// Motored version
			
			double T_wg = (h_cg*T_g + EnginePtr->coeffs*EnginePtr->T_c)/(h_cg + EnginePtr->coeffs); // Use as initial estimate for iteration
			double T_wg_mot = (h_cg_mot*T_g_mot + EnginePtr->coeffs*EnginePtr->T_c)/(h_cg_mot + EnginePtr->coeffs);
		
			double T_wg_old, T_wg_mot_old;
			do // Iterate due to radiative terms
			{				
				T_wg_old = T_wg;
				T_wg = (h_cg*T_g + pPpt->sigma*EnginePtr->epsilon*pow(T_g, 4) + EnginePtr->coeffs*EnginePtr->T_c)
					  /(h_cg + EnginePtr->coeffs + pPpt->sigma*EnginePtr->epsilon*pow(T_wg, 3));
				T_wg_mot_old = T_wg_mot;
				T_wg_mot = (h_cg_mot*T_g_mot + pPpt->sigma*EnginePtr->epsilon*pow(T_g_mot, 4) + EnginePtr->coeffs*EnginePtr->T_c)
						  /(h_cg_mot + EnginePtr->coeffs + pPpt->sigma*EnginePtr->epsilon*pow(T_wg_mot, 3));
			}
			while(fabs(T_wg - T_wg_old)/T_wg_old > EnginePtr->WOS_TOL);

			double THETAR = THETA*PI/180;
			double FCH = FCYL;	// Cylinder head surface area, m_Fch (m^2).
			double FCC = FCH + FCYL + (PI*EnginePtr->dcyl*EnginePtr->stroke/2)*(CONRAT + 1 - cos(THETAR) - sqrt(pow(CONRAT, 2) - pow(sin(THETAR),2)));
			// Instantaneous combustion chamber surface area
			
			double Qcoolrate = h_cg*FCC*(T_g - T_wg);
			DQDT -= Qcoolrate;	// Substract from overall heat transfer to cylinder

			double Qcoolrate_mot = h_cg_mot*FCC*(T_g_mot - T_wg_mot);
			DQDT_mot -= Qcoolrate_mot;	// Substract from overall motored heat transfer to cylinder
		}
		break;
	default:
		EnginePtr->HT_MODEL = ANNAND;
		cout << "CCylinder::Update [" << ID << "]: Unknown cylinder heat transfer model; selecting ANNAND\n";
		break;
	}

WAIT = false;


	// Valve mass fluxes
	// =================
	bool UPDATE_PIPE = true; // Update pipe values with valve calls (as opposed to motored valve calls)
	
	
	// Valve opening profile
	// ----------------------------------------------------------------------------------------------------
	int WAVEFORM;
	bool VALVE_OPEN;
	int END_or_phi;
	double end_closing;
	double psi_temp;
	double cycle = EnginePtr->cycle;
	double cycle_length = (cycle / 2) * 360;  // Cycle length is usually 720 or 360
	double begin;
	double phase; // = 180; // Default value (hardcoded on purpose)
	double min; 
	double max;
	double amp;
	double mean;
	double phi;
	double sf;
	double end;
	bool CYCLE_EXCEEDED = false;


	if(NCYCLES>0) // Don't run valves on first cycle (else may start simulation in middle of valve open)
	{
		// Set new exhaust valve area

		for(i=0; i<EnginePtr->NEXVALVES; ++i) {

			WAVEFORM = EnginePtr->EV_WAVEFORM;
			VALVE_OPEN = EnginePtr->EV_OPEN;
			END_or_phi = EnginePtr->EV_END_or_phi;
			end_closing = EnginePtr->EV_END_CLOSING;
			min = EnginePtr->EV_MIN_OPENING;
			max = EnginePtr->EV_MAX_OPENING;
			begin = EnginePtr->EV_BEGIN_OPENING;
			phase = 180; // Default value (hardcoded on purpose)
			phi = EnginePtr->EV_phi;
			sf = EnginePtr->scale_factor_ev;
			amp = (max - min) / 2;
			mean = (max + min) / 2;

			// Put next code as separate function?
			// -----------------------------------

			if (WAVEFORM != 5) // If valve opening is not read from file, i.e., calculated
			{
				if (WAVEFORM == 0) { // Valve is either open or closed
				
					if (VALVE_OPEN) psi_temp = max; // Set to maximum opening
					else psi_temp = min; // Set to minimum opening
				}
				else {

					if (true) //if (NCYCLES == 1 || !(ExhaustValve[i].ACTIVE)) // If ACTIVE, only calculate psi waveform on 1st open cycle
					{ 

						// If valve opening profile is square or triangular or sinusoidal or a Fourier series, i.e., WAVEFORM == 1 or 2 or 3 or 4
						if (WAVEFORM == 1 || WAVEFORM == 2 || WAVEFORM == 3 || WAVEFORM == 4) {

							if (END_or_phi == 0) {		// End angle is specified (but still subject to scale factor)

								end = end_closing;
							}
							else {
								if (END_or_phi == 1) {	// Use phi fraction to determine closing point

									end = begin + phi * cycle_length;
								}
								else {
									// UNDEFINED AS YET
								}
							}

							// Apply scale factor
							end = begin + (end - begin) / sf;
							double end_corr = end; // Initiate end_corr to be same as end for now

							if (end > cycle_length) {

								CYCLE_EXCEEDED = true; // Flag the overflow
								end_corr -= cycle_length; // Now the end of the opening is before the beginning (on a cycle basis)
							}
							else CYCLE_EXCEEDED = false;


							if ((!CYCLE_EXCEEDED && (THETA < begin || THETA > end)) || (CYCLE_EXCEEDED && (THETA > end_corr && THETA < begin))) {
								// If THETA lies outside the range where the opening profile is varying, set it to the minimum opening for all cases
								psi_temp = min;
							}
							else {
								// Square wave opening
								if (WAVEFORM == 1) psi_temp = max;

								else {

									if (WAVEFORM == 2) { // Triangular wave opening

										double middle = (begin + end) / 2;
										double triangle_grad_mag = fabs((max - min) / (middle - begin));

										if (!CYCLE_EXCEEDED) { // End of valve closing does not exceed cycle length

											if (begin <= end) { // Conventional labelling

												if (THETA <= middle) // The angle of the peak value (i.e., the middle)

													psi_temp = min + triangle_grad_mag * (THETA - begin);

												else
													psi_temp = min + triangle_grad_mag * (end - THETA);
											}
											else {
												// UNCONVENTIONAL begin AND end LABELLING
											}
										}
										else { // End of valve closing would exceed cycle length

											// Now two sub cases depending on where the middle point lies
											double middle_corr = middle;
											bool MIDDLE_CORRECTED = false;

											if (middle >= cycle_length) {
												MIDDLE_CORRECTED = true;
												middle_corr -= cycle_length; // Corrected middle point
											}

											if (!MIDDLE_CORRECTED) { // Middle point not corrected

												if (THETA <= middle && THETA >= begin) {

													psi_temp = min + triangle_grad_mag * (THETA - begin);
												}
												else {

													if (THETA >= middle)

														psi_temp = max - triangle_grad_mag * (THETA - middle);

													else

														psi_temp = min + triangle_grad_mag * (end_corr - THETA);
												}
											}
											else { // Corrected middle point

												if (THETA <= middle && THETA >= begin) {

													psi_temp = min + triangle_grad_mag * (THETA - begin);
												}
												else {

													if (THETA >= middle_corr)

														psi_temp = max - triangle_grad_mag * (THETA - middle_corr);

													else

														psi_temp = max - triangle_grad_mag * (middle_corr - THETA);
												}
											}
										}
									}
									else {

										if (EnginePtr->EV_WAVEFORM == 3) { // Sinusoidal opening

											// Sine profile: value = amp* SIN(  ( (THETA - (begin + phase * phi / sf) ) / (cycle/2) / phi * sf) * (PI() / 180)  ) + mean
											//psi_temp = amp * sin(((THETA - (begin + phase * phi / sf)) / (cycle / 2) / phi * sf) * (PI / 180)) + mean;
											psi_temp = amp * sin(((THETA - (begin + (phase + ExhaustValve[i].phase_shift) * phi / sf)) / (cycle / 2) / phi * sf) * (PI / 180)) + mean;
											// THETA and PHASE in degrees
										}
										else {
											// Fourier series opening
											//
										}
									}
								}
							}
						}
						else
						{
							pPpt->Out(this->Identify()); pPpt->Out(" - unknown WAVEFORM setting. Exiting.");
							exit(1);
						}

						if (ExhaustValve[i].ACTIVE) // If ACTIVE, now record the calculated psi to the array for next cycle
						{
							// Find the relevant data point in avt_raw for the current crank angle
							int row = 0;
							do { ++row; } // Hence first row evaluated is [1], on purpose
							while (ExhaustValve[i].avt_raw[row][ExhaustValve[i].CAD] < this->EnginePtr->ca && row < ExhaustValve[i].numAVTDataPoints - 1);

							ExhaustValve[i].avt_raw[row][ExhaustValve[i].PSI_UNCORR] = psi_temp;
							ExhaustValve[i].avt_raw[row][ExhaustValve[i].PSI_CORR] = psi_temp;;
						}

					}
					else { // ACTIVE and NCYCLES >= 2

						// Find the relevant data point in avt_raw for the current crank angle
						int row = 0;
						do { ++row; } // Hence first row evaluated is [1], on purpose
						while (ExhaustValve[i].avt_raw[row][ExhaustValve[i].CAD] < this->EnginePtr->ca && row < ExhaustValve[i].numAVTDataPoints - 1);

						//double error = ExhaustValve[i].avt_raw[row][ExhaustValve[i].ERROR];
						//double error_grad = ExhaustValve[i].avt_raw[row][ExhaustValve[i].ERROR_GRAD]; // Error gradient with lift
						//double psi_current = ExhaustValve[i].avt_raw[row][ExhaustValve[i].PSI_CORR];
						//double psi_prev = ExhaustValve[i].avt_raw[row][ExhaustValve[i].PSI_CORR_PREV];

						psi_temp = ExhaustValve[i].avt_raw[row][ExhaustValve[i].PSI_CORR]; // Set psi_temp as value from last cycle (since not calculated any more)
						
						ExhaustValve[i].avt_raw[row][ExhaustValve[i].PSI_UNCORR] = psi_temp; // The unmodified psi (lift or area)

						if (DoubleToBool(ExhaustValve[i].avt_raw[row][ExhaustValve[i].SIGN_CHANGE]))
							ExhaustValve[i].avt_raw[row][ExhaustValve[i].CORR_FACTOR] *= 0.5; // Half the correction factor if the error has changed sign

						// AVT: Adjust instantaneous opening according to error once cycle average pressure in matched well enough
						if (ExhaustValve[i].ACTIVE && ExhaustValve[i].CYC_AV_MATCHED && ExhaustValve[i].PHASE_MATCHED) { // if (ExhaustValve[i].ACTIVE && ExhaustValve[i].CYC_AV_MATCHED) { // && PHASE_MATCHED

///*
							// Row-by-row psi adjustment
							if (fabs(ExhaustValve[i].avt_raw[row][ExhaustValve[i].ERROR_GRAD]) > 1e-6) { // Prevent div by 0

								double test_val;
								int counter = 0;
								bool VIABLE = true;
								do {
									++counter;
									// Must prevent negative opening (lift or area)!
									test_val = psi_temp - //+ //- //+ 

										ExhaustValve[i].avt_raw[row][ExhaustValve[i].CORR_FACTOR]
										* (ExhaustValve[i].avt_raw[row][ExhaustValve[i].ERROR] - 0)
										/ ExhaustValve[i].avt_raw[row][ExhaustValve[i].ERROR_GRAD];

									//if (test_val < 0) ExhaustValve[i].avt_raw[row][ExhaustValve[i].CORR_FACTOR] *= 0.5;
									if (test_val < 0 ||
										(test_val > psi_temp && test_val / psi_temp > 2) ||
										(psi_temp > test_val && psi_temp / test_val > 2)) // Limit maximum adjustment
									{
										ExhaustValve[i].avt_raw[row][ExhaustValve[i].CORR_FACTOR] *= 0.5;
										VIABLE = false;
									}
									else VIABLE = true;

									if (counter > 20) cout << "counter = " << counter << endl;
								} while (!VIABLE && counter < 100);

								if (VIABLE) psi_temp = test_val; //else psi_temp = psi_temp;
							}
							else psi_temp += 0;
//*/
						}
						
						ExhaustValve[i].avt_raw[row][ExhaustValve[i].PSI_CORR_PREV] = ExhaustValve[i].avt_raw[row][ExhaustValve[i].PSI_CORR];
						ExhaustValve[i].avt_raw[row][ExhaustValve[i].PSI_CORR] = psi_temp;
/*
							cout << endl;
							//cout << "[" << row << "][CAD]=" << ExhaustValve[i].avt_raw[row][ExhaustValve[i].CAD] << "\t";
							//cout << "[" << row << "][TARGET]=" << ExhaustValve[i].avt_raw[row][ExhaustValve[i].TARGET] << "\t";
							//cout << "[" << row << "][RECORDED]=" << ExhaustValve[i].avt_raw[row][ExhaustValve[i].RECORDED] << "\t";
							//cout << "[" << row << "][ERROR]=" << ExhaustValve[i].avt_raw[row][ExhaustValve[i].ERROR] << "\t";
							//cout << "[" << row << "][ERROR_GRAD]=" << ExhaustValve[i].avt_raw[row][ExhaustValve[i].ERROR_GRAD] << "\t";
							//cout << "[" << row << "][SIGN_CHANGE]=" << TrueOrFalse(DoubleToBool(ExhaustValve[i].avt_raw[row][ExhaustValve[i].SIGN_CHANGE])) << "\t";
							cout << "[" << row << "][CORR_FACTOR]=" << ExhaustValve[i].avt_raw[row][ExhaustValve[i].CORR_FACTOR] << "\t";
							cout << "[" << row << "][PSI_UNCORR]=" << ExhaustValve[i].avt_raw[row][ExhaustValve[i].PSI_UNCORR] << "\t";
							cout << "[" << row << "][PSI_CORR]=" << ExhaustValve[i].avt_raw[row][ExhaustValve[i].PSI_CORR] << "\t";
							//cout << "[" << row << "][PSI_CORR_PREV]=" << ExhaustValve[i].avt_raw[row][ExhaustValve[i].PSI_CORR_PREV] << "\t";
							cout << endl;
							//exit(1);
//*/
					}
				}




				// Now set the effective flow area in the valve
				if (EnginePtr->EV_L_OR_A == 0) { // If exhaust valve data are lift values, i.e., EV_L_OR_A == 0, and psi_temp = lift
					
					ExhaustValve[i].Set_final_lift(psi_temp); // [mm]

					if (EnginePtr->CONST_REF_AREA_EXH) {	// Use constant reference area (i.e., values are flow coefficient) (1) 
						// Applying flow coefficient and a reference area based on the (constant) valve face diameter
						double flow_coeff = 0;
						if (psi_temp > 0) {
							flow_coeff = 1; // THIS NEEDS TO BE INTERPOLATED FROM FILE *****
						}
						else {
							flow_coeff = 0;
						}

						ExhaustValve[i].Set_valve_Cf_or_Cd(flow_coeff);
						ExhaustValve[i].Set_eff_area(flow_coeff * ((PI / 4) * pow(this->EnginePtr->REF_DIA_EXH, 2)) / 1e6); // [m^2]
					}
					else { // Or curtain reference area (i.e., values are discharge coefficient) (0) to interpret file values
						// Applying discharge coefficient and a reference area based on the (varying) valve curtain area
						
						double discharge_coeff = 1; // THIS NEEDS TO BE INTERPOLATED FROM FILE *****

						ExhaustValve[i].Set_valve_Cf_or_Cd(discharge_coeff);
						ExhaustValve[i].Set_eff_area(discharge_coeff * (psi_temp * PI * this->EnginePtr->REF_DIA_EXH) / 1e6); // [m^2]	
					}
				}
				else
				{
					if (EnginePtr->EV_L_OR_A == 1) {
						// Else exhaust valve data are effective flow area values, i.e., EV_L_OR_A == 1, and psi_temp = effective flow area (m^2)
						ExhaustValve[i].Set_eff_area(psi_temp);
					}
					else {
						// Else exhaust valve data are psi values (valve to pipe area ratio), i.e., EV_L_OR_A == 2 and psi_temp = psi
						ExhaustValve[i].Set_eff_area(psi_temp * ExhaustValve[i].Get_FP()); // FP is the area of the adjoining pipe
					}
				}
			}
			else ExhaustValve[i].EffectiveArea(pPpt, THETA * EnginePtr->scale_factor_ev); // Valve opening interpolated from file. Interpolates an up-to-date area.
			
			// In any case:
			if(ExhaustValve[i].Get_eff_area()<=0.0) {DMEDT[i]=0.0; ExhaustValve[i].Set_open(false);}
			else ExhaustValve[i].Set_open(true);
		}

		// Get intake valve area
		for (i = 0; i < EnginePtr->NINVALVES; ++i)
		{
/*
			double psi_temp;
			double iv_psi_max = EnginePtr->IV_MAX_OPENING;
			double iv_psi_min = EnginePtr->IV_MIN_OPENING;

			//if(EnginePtr->IV_SWITCH) // Valve is either fully open or fully closed

			if (EnginePtr->IV_WAVEFORM != 5) // If valve opening is not read from file, i.e., calculated
			{
				if (EnginePtr->IV_WAVEFORM == 0) // Valve is either open or closed
				{
					//if(EnginePtr->IV_OPEN) IntakeValve[i].Set_eff_area(IntakeValve[i].Get_FP()); // Fully open; set a psi of 1
					//else IntakeValve[i].Set_eff_area(0); // Fully closed; set a psi of 0
					//if (EnginePtr->IV_OPEN) IntakeValve[i].Set_eff_area(IntakeValve[i].Get_FP() * iv_psi_max); // Set to maximum opening
					//else IntakeValve[i].Set_eff_area(IntakeValve[i].Get_FP() * iv_psi_min); // Set to minimum opening
					if (EnginePtr->IV_OPEN) psi_temp = iv_psi_max; // Set to maximum opening
					else psi_temp = iv_psi_min; // Set to minimum opening
				}
				else {
					//if(EnginePtr->IV_CALC) // Calculate valve timing
					// If intake valve lift or area profile is triangular or square or sinusoidal or a Fourier series, i.e., IV_WAVEFORM == 1 or 2 or 3 or 4
					if (EnginePtr->IV_WAVEFORM == 1 || EnginePtr->IV_WAVEFORM == 2 || EnginePtr->IV_WAVEFORM == 3 || EnginePtr->IV_WAVEFORM == 4) {

						// Calculate valve timing
						double iv_begin = EnginePtr->IV_BEGIN_OPENING;
						double iv_end; 
						
						if (EnginePtr->IV_END_or_phi == 0) { // Use IV_END_CLOSING to determine closing point
							
							iv_end = EnginePtr->IV_END_CLOSING;
						}
						else {
							if (EnginePtr->IV_END_or_phi == 1) { // Use IV_phi fraction to determine closing point

								iv_end = iv_begin + EnginePtr->IV_phi * (360 * (EnginePtr->cycle / 2));
							}
							else {
							// UNDEFINED AS YET
							}
						}

						// Apply scale factor
						iv_end = iv_begin + (iv_end - iv_begin) / EnginePtr->scale_factor_iv;

						if ((iv_begin <= iv_end && (THETA <= iv_begin || THETA >= iv_end)) || (iv_begin > iv_end && (THETA <= iv_begin && THETA >= iv_end)))
							psi_temp = iv_psi_min; // Outside period of variation
						else {
							//if(EnginePtr->IV_SQUARE) 
							if (EnginePtr->IV_WAVEFORM == 1) psi_temp = iv_psi_max; // Square wave opening
							else {
								// Triangular wave opening
								if (EnginePtr->IV_WAVEFORM == 2) {
									if (iv_begin <= iv_end) {
										if (THETA <= (iv_begin + iv_end) / 2)
											psi_temp = iv_psi_min + (iv_psi_max - iv_psi_min) / ((iv_end - iv_begin) / 2) * (THETA - iv_begin);
										else
											psi_temp = iv_psi_min + (iv_psi_max - iv_psi_min) / ((iv_end - iv_begin) / 2) * (iv_end - THETA);
									}
									else {
										double end_temp = iv_end + (EnginePtr->cycle / 2) * 360;
										double middle = (iv_begin + end_temp) / 2;
										if (middle > (EnginePtr->cycle / 2) * 360) middle -= (EnginePtr->cycle / 2) * 360;

										if (middle > iv_begin) {
											if (THETA <= middle && THETA >= iv_begin) {
												psi_temp = iv_psi_min + ((iv_psi_max - iv_psi_min) / (middle - iv_begin)) * (THETA - iv_begin);
											}
											else {
												if (THETA >= middle)
													psi_temp = iv_psi_max - ((iv_psi_max - iv_psi_min) / (middle - iv_begin)) * (THETA - middle);
												else
													psi_temp = iv_psi_max - ((iv_psi_max - iv_psi_min) / (middle - iv_begin)) * ((THETA + (EnginePtr->cycle / 2) * 360) - middle);
											}
										}
										else {}
										if (THETA >= middle && THETA <= iv_end) {
											psi_temp = iv_psi_max - ((iv_psi_max - iv_psi_min) / (iv_end - middle)) * (THETA - middle);
										}
										else {
											if (THETA <= middle)
												psi_temp = iv_psi_min + ((iv_psi_max - iv_psi_min) / (iv_end - middle)) * ((THETA + (EnginePtr->cycle / 2) * 360) - iv_begin);
											else
												psi_temp = iv_psi_min + ((iv_psi_max - iv_psi_min) / (iv_end - middle)) * (THETA - iv_begin);
										}
									}
								}
								else {
									if (EnginePtr->IV_WAVEFORM == 3) { // Sinusoidal opening

										// TO DO
									}
									else { // Fourier series opening
										//
									}
								}
							}
						}
					}
					else
					{
						pPpt->Out(this->Identify()); pPpt->Out(" - unknown IV_WAVFORM setting. Exiting.");
						exit(1);
					}
				}

				// Now set the effective flow area in the valve
				if (EnginePtr->IV_L_OR_A == 0) { // If intake valve data are lift values, i.e., IV_L_OR_A == 0, and psi_temp = lift
					
					IntakeValve[i].Set_final_lift(psi_temp); // [mm]

					if (EnginePtr->CONST_REF_AREA_INT) {	// Use constant reference area (i.e., values are flow coefficient) (1) 
						// Applying flow coefficient and a reference area based on the (constant) valve face diameter
						double flow_coeff = 0;
						if (psi_temp > 0) {
							flow_coeff = 1; // THIS NEEDS TO BE INTERPOLATED FROM FILE
						}
						else {
							flow_coeff = 0;
						}
						IntakeValve[i].Set_eff_area(flow_coeff * ((PI / 4) * pow(this->EnginePtr->REF_DIA_INT, 2)) / 1e6); // [m^2]
					}
					else { // Or curtain reference area (i.e., values are discharge coefficient) (0) to interpret file values
						// Applying discharge coefficient and a reference area based on the (varying) valve curtain area
						double discharge_coeff = 1; // Allow this to be interpolated from a file
						IntakeValve[i].Set_eff_area(discharge_coeff * (psi_temp * PI * this->EnginePtr->REF_DIA_INT) / 1e6); // [m^2]	
					}
				}
				else
				{
					if (EnginePtr->IV_L_OR_A == 1) {
						// Else intake valve data are effective flow area values, i.e., IV_L_OR_A == 1, and psi_temp = effective flow area
						IntakeValve[i].Set_eff_area(psi_temp);
					}
					else {
						// Else intake valve data are psi values (valve to pipe area ratio), i.e., IV_L_OR_A == 2 and psi_temp = psi
						IntakeValve[i].Set_eff_area(psi_temp * IntakeValve[i].Get_FP()); // FP is the area of the adjoining pipe
					}
				}
			}
			else IntakeValve[i].EffectiveArea(pPpt, THETA * EnginePtr->scale_factor_iv); // Valve opening interpolated from file. Interpolates an up-to-date area. 

			// In any case:
			if (IntakeValve[i].Get_eff_area() <= 0.0) { DMIDT[i] = 0.0; IntakeValve[i].Set_open(false); }
			else IntakeValve[i].Set_open(true);
*/
		}

/*
		for(i=0; i<EnginePtr->NINVALVES; ++i)
		{
			if(EnginePtr->IAIR) // If there is an intake valve
			{
				if(EnginePtr->IV_SWITCH) // Valve is either fully or open or fully closed
				{
					if(EnginePtr->IV_OPEN) IntakeValve[i].Set_eff_area(IntakeValve[i].Get_FP()); // Fully open; set a psi of 1
					else IntakeValve[i].Set_eff_area(0); // Fully closed; set a psi of 0
				}
				else // Use valve timing
				{
					if(EnginePtr->IV_CALC) // Calculate valve timing
					{		
						double psi_temp;
						double iv_begin = EnginePtr->IV_BEGIN;
						double iv_end = EnginePtr->IV_END;
						double iv_psi_max = EnginePtr->IV_PSI_MAX;
						double iv_psi_min = EnginePtr->IV_PSI_MIN;
	
						if( (iv_begin <= iv_end && (THETA <= iv_begin || THETA >= iv_end))
							||
							(iv_begin > iv_end && (THETA <= iv_begin && THETA >= iv_end))
							)
							psi_temp = iv_psi_min; // Outside period of variation
						else
						{
							if(EnginePtr->IV_SQUARE) // Square wave opening
							{
								psi_temp = iv_psi_max;
							}
							else // Triangular wave opening
							{
								if(iv_begin <= iv_end) 
								{
									if(THETA <= (iv_begin + iv_end)/2)
										psi_temp = iv_psi_min + (iv_psi_max - iv_psi_min)/((iv_end - iv_begin)/2)*(THETA - iv_begin);
									else
										psi_temp = iv_psi_min + (iv_psi_max - iv_psi_min)/((iv_end - iv_begin)/2)*(iv_end - THETA);
								}
								else
								{					
									double end_temp = iv_end + (EnginePtr->cycle/2)*360;
									double middle = (iv_begin + end_temp)/2;
									if(middle > (EnginePtr->cycle/2)*360) middle -= (EnginePtr->cycle/2)*360;
		
									if(middle > iv_begin)
									{
										if(THETA <= middle && THETA >= iv_begin)
										{
											psi_temp = iv_psi_min + ((iv_psi_max - iv_psi_min)/(middle - iv_begin))*(THETA - iv_begin);
										}
										else
										{
											if(THETA >= middle)
												psi_temp = iv_psi_max - ((iv_psi_max - iv_psi_min)/(middle - iv_begin))*(THETA - middle);
											else
												psi_temp = iv_psi_max - ((iv_psi_max - iv_psi_min)/(middle - iv_begin))*((THETA + (EnginePtr->cycle/2)*360) - middle);
										}
									}
									else
									{
										if(THETA >= middle && THETA <= iv_end)
										{
											psi_temp = iv_psi_max - ((iv_psi_max - iv_psi_min)/(iv_end - middle))*(THETA - middle);
										}
										else
										{
											if(THETA <= middle)
												psi_temp = iv_psi_min + ((iv_psi_max - iv_psi_min)/(iv_end - middle))*((THETA + (EnginePtr->cycle/2)*360) - iv_begin);
											else
												psi_temp = iv_psi_min + ((iv_psi_max - iv_psi_min)/(iv_end - middle))*(THETA - iv_begin);
										}
									}
								}
							}
						}
						IntakeValve[i].Set_eff_area(psi_temp*IntakeValve[i].Get_FP()); // Apply psi value
					}
					else IntakeValve[i].EffectiveArea(pPpt, THETA); // Interpolates an up-to-date area
				}
	
				if(IntakeValve[i].Get_eff_area()<=0.0) {DMIDT[i]=0.0; IntakeValve[i].Set_open(false);}
				else IntakeValve[i].Set_open(true);
			}
			else
			{ 
				IntakeValve[i].Set_open(false);
				DMIDT[i]=0.0;
			}
		}
*/
	}

	// Enter exhaust valve boundary
	double* AEXH; AEXH = new double [EnginePtr->NEXVALVES];
	for(i=0; i<EnginePtr->NEXVALVES; ++i)
	{
		double psi_exh = ExhaustValve[i].Get_eff_area()/ExhaustValve[i].Get_FP();
		if(psi_exh>1) psi_exh=1; // Though eff_area can be greater than FP, the ratio cannot be greater than 1
		if(WAIT) psi_exh = 0; // Waiting until next cycle - simulate closed valve by setting psi = 0

		if(!pPpt->HOMENTROPIC) DMEDT[i] = ExhaustValve[i].Poppet_NH(pPpt, psi_exh, PC, TC, pPpt->AREFe, ID, UPDATE_PIPE, timestep, time);
		else ExhaustValve[i].Poppet_H(pPpt, timestep, time, DMEDT[i], psi_exh, pPpt->AREFe, pPpt->PREF, PC, TC, ExhaustValve[i].Get_FP());

		AEXH[i] = ((ExhaustValve[i].Get_CLOUT() + ExhaustValve[i].Get_CLIN())/2)*pPpt->AREFe;
	}

	// Enter intake valve boundary
	double* AAIR; AAIR = new double [EnginePtr->NINVALVES];
	for(i=0; i<EnginePtr->NINVALVES; ++i)
	{
		if(EnginePtr->IAIR)// && !WAIT)
		{
			if(EnginePtr->IPIPE)
			{
				double psi_air = IntakeValve[i].Get_eff_area()/IntakeValve[i].Get_FP();
				if(psi_air>1) psi_air=1; // Though eff_area can be greater than FP, the ratio cannot be greater than 1
				if(WAIT) psi_air = 0; // Waiting until next cycle - simulate closed valve by setting psi = 0

				// Intake system is a pipe with valve, not a receiver
				if(!pPpt->HOMENTROPIC) DMIDT[i] = IntakeValve[i].Poppet_NH(pPpt, psi_air, PC, TC, pPpt->AREFi, ID, UPDATE_PIPE, timestep, time);
				else IntakeValve[i].Poppet_H(pPpt, timestep, time, DMIDT[i], psi_air, pPpt->AREFi, pPpt->PREF, PC, TC, IntakeValve[i].Get_FP());
///				
				AAIR[i] = ((IntakeValve[i].Get_CLOUT() + IntakeValve[i].Get_CLIN())/2)*pPpt->AREFi;
			}
			else
			{
				double AREA_AIR = IntakeValve[i].Get_eff_area();
				if(WAIT) AREA_AIR = 0; // Waiting until next cycle - simulate closed valve by setting psi = 0

				// Intake system is the receiver
				AAIR[i] = sqrt(pPpt->gammaAir(EnginePtr->TA)*287.0*EnginePtr->TA);

				// If we have an open intake receiver valve
				if(IntakeValve[i].Get_eff_area()>0)	AirReceiver(pPpt, DMIDT[i], EnginePtr->PA, EnginePtr->TA, PC, TC, AAIR[i], AC, AREA_AIR);
				else DMIDT[i]=0.0;
			}
		}
	}

	// Sum mass fluxes for each set of valves
	double DMEDT_TOTAL = 0;	for(i=0; i<EnginePtr->NEXVALVES; ++i) DMEDT_TOTAL += DMEDT[i];
	double DMIDT_TOTAL = 0;	for(i=0; i<EnginePtr->NINVALVES; ++i) DMIDT_TOTAL += DMIDT[i];
	
	MC = MC + (DMIDT_TOTAL + DMEDT_TOTAL)*DT;	
	// New mass of air in cylinder - here both DMEDT_TOTAL and DMIDT_TOTAL are positive when entering cylinder

	// Calculate cylinder volume and DVC/DT
	if(!EnginePtr->IVOL) DVCDT = 0.0; // Constant volume cylinder
	else
	{
		double THETAR = THETA*PI/180;
		FNN = sqrt(pow(CONRAT, 2) - pow(sin(THETAR), 2));
		X = 0.5*EnginePtr->stroke*(1.0 + CONRAT - FNN - cos(THETAR));
		VC = FCYL*(X + EnginePtr->stroke/(EnginePtr->cr - 1.0));
		DXDTH = 0.5*EnginePtr->stroke*sin(THETAR)*(1.0 + cos(THETAR)/FNN);
		DTHDT = 2.0*PI*EnginePtr->reveng;
		DVCDT = FCYL*DXDTH*DTHDT;
	}
	if( (EnginePtr->MODEL != READ_P_T) && (EnginePtr->MODEL != CONST_P) ) TC = 1.0E5*PC*VC/(287.0*MC); // Do not calculate TC if read from file
	AC = sqrt(pPpt->gammaAir(TC)*287.0*TC);

	// Check direction of flow for exhaust and select appropriate speed of sound
	double* AE; AE = new double [EnginePtr->NEXVALVES];
	for(i=0; i<EnginePtr->NEXVALVES; ++i)
	{	
		if(DMEDT[i]<0.0) AE[i] = AC;
		if(DMEDT[i]>0.0) AE[i] = AEXH[i];
		if(DMEDT[i]==0.0) AE[i] = 0.0;
	}

	// Check direction of flow for intake and select appropriate speed of sound
	double* AI; AI = new double [EnginePtr->NINVALVES];
	for(i=0; i<EnginePtr->NINVALVES; ++i)
	{
		if(DMIDT[i]<0.0) AI[i] = AC;
		if(DMIDT[i]>0.0) AI[i] = AAIR[i];
		if(DMIDT[i]==0.0) AI[i] = 0.0;
	}
	
	// Calculate pressure rates of change due to individual mass fluxes
	double DPCDT_DMEDT = 0; for(i=0; i<EnginePtr->NEXVALVES; ++i) DPCDT_DMEDT += DMEDT[i]*pow(AE[i], 2)/(VC*1.0E5);
	double DPCDT_DMIDT = 0; for(i=0; i<EnginePtr->NINVALVES; ++i) DPCDT_DMIDT += DMIDT[i]*pow(AI[i], 2)/(VC*1.0E5);

	// Calculate overall pressure rate of change
	DPCDT = -pPpt->gammaAir(TC)*DVCDT*PC/VC					// Volume change component
				+ DPCDT_DMEDT								// Exhaust mass flux component
				+ DPCDT_DMIDT								// Intake mass flux component
				+ ((pPpt->gammaAir(TC)-1)/(VC*1.0E5))*DQDT		// Heat transfer component (includes heat release due to combustion) 
				;

	// Now calculate increment in work done (PC*VC)
	double del_WK = (DVCDT/fabs(DVCDT))		// Work is positive when cyclinder volume expands (e.g. power stroke)
					*(PC*1e5)*VC;			// (Pa * m^3) = J or Nm

	double del_WK_mot = (DVCDT/fabs(DVCDT))	
						*(PC_mot*1e5)*VC;		

	if(THETA <= 180.0 || THETA >= 540.0)
	{
		// Increment gross indicated work per cycle (compression and expansion strokes only)
		W_cig += del_WK;	
		W_cig_mot += del_WK_mot;
	}
	// Increment net indicated work per cycle (all strokes)
	W_cin += del_WK;
	W_cin_mot += del_WK_mot;
	
	EnginePtr->del_WK_total += del_WK;		// Add to overall engine del_WK
	EnginePtr->CYCLE_WK += del_WK;			// Add to cumulative engine CYCLE_WK
	EnginePtr->P_inst = EnginePtr->del_WK_total/DT; // Recalulate engine instantaneous power
	
	// Print to files
	if(NCYCLES >= print_from_cycle - 1)
	{
		PrintToFile(pPpt, timestep, time, EnginePtr->ca_elapsed);
		if(PRINT_TO_FILE_MOV && timestep%mov_freq==0)
			PrintToFileMovie(pPpt, timestep, time, EnginePtr->ca_elapsed);
	}

	// Motored calculations	(Woschni heat transfer only)
	// =================================================
	if(EnginePtr->HT_MODEL == WOSCHNI) // Need to calculate motored values for Woschni heat transfer correlation
	{
		// Reset combustion variables prior to combustion			
		if(THETA>EnginePtr->IVC && AllClosed() && !MOT_SET) // Double check all valves are closed
		{
			PC_mot_ref = PC; TC_mot_ref = TC; VC_mot_ref = VC;	
			// At IVC get reference values for mean local gas velocity calculation
			MOT_SET = true;
		}

		// Valve mass fluxes
		// =================
		UPDATE_PIPE = false; // Do not update pipe values with motored valve calls
		
		for(i=0; i<EnginePtr->NEXVALVES; ++i){if(ExhaustValve[i].Get_eff_area()<=0.0) DMEDT_mot[i]=0.0;} // Already found exhaust valve area
		
		for(i=0; i<EnginePtr->NINVALVES; ++i)
		{
			if(EnginePtr->IAIR) // If there is an intake valve
				{if(IntakeValve[i].Get_eff_area()<=0.0) DMIDT_mot[i]=0.0;} // Already found intake valve area
			else DMIDT_mot[i]=0.0;
		}
		
		// Enter exhaust valve boundary - motored
		for(i=0; i<EnginePtr->NEXVALVES; ++i)
		{
			double psi_exh = ExhaustValve[i].Get_eff_area()/ExhaustValve[i].Get_FP();
			if(WAIT) psi_exh = 0; // Waiting until next cycle - simulate closed valve by setting psi = 0
			
			DMEDT_mot[i] = ExhaustValve[i].Poppet_NH(pPpt, psi_exh, PC_mot, TC_mot, pPpt->AREFe, ID, UPDATE_PIPE, timestep, time); 
			// Woschni requires non-homentropic method
		}
		
		// Enter intake valve boundary - motored
		for(i=0; i<EnginePtr->NINVALVES; ++i)
		{
			if(EnginePtr->IAIR)// && !WAIT)
			{
				if(EnginePtr->IPIPE)
				{
					double psi_air = IntakeValve[i].Get_eff_area()/IntakeValve[i].Get_FP();
					if(WAIT) psi_air = 0; // Waiting until next cycle - simulate closed valve by setting psi = 0
					
					// Intake system is a pipe with valve, not a receiver
					DMIDT_mot[i] = IntakeValve[i].Poppet_NH(pPpt, psi_air, PC_mot, TC_mot, pPpt->AREFi, ID, UPDATE_PIPE, timestep, time);
					// Woschni requires non-homentropic method
				}
				else
				{
					double AREA_AIR = IntakeValve[i].Get_eff_area();
					if(WAIT) AREA_AIR = 0; // Waiting until next cycle - simulate closed valve by setting psi = 0
					
					// If we have an open intake receiver valve
					if(IntakeValve[i].Get_eff_area()>0)	AirReceiver(pPpt, DMIDT_mot[i], EnginePtr->PA, EnginePtr->TA, PC_mot, TC_mot, AAIR[i], AC_mot, AREA_AIR);
					else DMIDT_mot[i]=0.0;
				}
			}
		}

		// Sum motored mass fluxes for each set of valves
		double DMEDT_mot_TOTAL = 0;	for(i=0; i<EnginePtr->NEXVALVES; ++i) DMEDT_mot_TOTAL += DMEDT_mot[i];
		double DMIDT_mot_TOTAL = 0;	for(i=0; i<EnginePtr->NINVALVES; ++i) DMIDT_mot_TOTAL += DMIDT_mot[i];
		
		MC_mot = MC_mot + (DMIDT_mot_TOTAL + DMEDT_mot_TOTAL)*DT;	
		// New motored mass of air in cylinder - here both DMEDT_mot_TOTAL and DMIDT_mot_TOTAL are positive when entering cylinder
		
		TC_mot = 1.0E5*PC_mot*VC/(287.0*MC_mot);
		AC_mot = sqrt(pPpt->gammaAir(TC_mot)*287.0*TC_mot);

		// Check direction of flow for exhaust and select appropriate motored speed of sound
		double* AE_mot; AE_mot = new double [EnginePtr->NEXVALVES];
		for(i=0; i<EnginePtr->NEXVALVES; ++i)
		{	
			if(DMEDT_mot[i]<0.0) AE_mot[i] = AC_mot;
			if(DMEDT_mot[i]>0.0) AE_mot[i] = AEXH[i];
			if(DMEDT_mot[i]==0.0) AE_mot[i] = 0.0;
		}
		
		// Check direction of flow for intake and select appropriate motored speed of sound
		double* AI_mot; AI_mot = new double [EnginePtr->NINVALVES];
		for(i=0; i<EnginePtr->NINVALVES; ++i)
		{
			if(DMIDT_mot[i]<0.0) AI_mot[i] = AC_mot;
			if(DMIDT_mot[i]>0.0) AI_mot[i] = AAIR[i];
			if(DMIDT_mot[i]==0.0) AI_mot[i] = 0.0;
		}

		// Calculate motored pressure rates of change due to individual mass fluxes
		double DPC_motDT_DMEDT = 0; for(i=0; i<EnginePtr->NEXVALVES; ++i) DPC_motDT_DMEDT += DMEDT_mot[i]*pow(AE_mot[i], 2)/(VC*1.0E5);
		double DPC_motDT_DMIDT = 0; for(i=0; i<EnginePtr->NINVALVES; ++i) DPC_motDT_DMIDT += DMIDT_mot[i]*pow(AI_mot[i], 2)/(VC*1.0E5);
		
		// Calculate overall motored pressure rate of change
		DPC_motDT = -pPpt->gammaAir(TC_mot)*DVCDT*PC_mot/VC	// Volume change component
			+ DPC_motDT_DMEDT							// Exhaust mass flux component
			+ DPC_motDT_DMIDT							// Intake mass flux component
			+ ((pPpt->gammaAir(TC_mot)-1)/(VC*1.0E5))*DQDT_mot // Heat transfer component (no combustion heat release - motored) 
			;
/*
		// Now calculate increment in work done (PC*VC)
		double del_WK_mot = (DVCDT/fabs(DVCDT))	
							*(PC_mot*1e5)*VC;		

		if(THETA <= 180.0 || THETA >= 540.0)
		{
			// Increment gross indicated work per cycle (compression and expansion strokes only)
			W_cig_mot += del_WK_mot;
		}
		// Increment net indicated work per cycle (all strokes)
		W_cin_mot += del_WK_mot;
*/
	}
	return;
}

double CCylinder::CombustionDQDT(CProperties* pPpt, double DT)
{	
	if(AllClosed()) // If all valves are closed
	{
		double m_finj_old, m_finj_this_iter;

		// Obtain amount of fuel injected so far this cycle
		if(THETA >= EnginePtr->inj_start || THETA < EnginePtr->inj_stop)	// If inside the injection period
		{
			// Calculate length of injection period in terms of crank angle
			double THETA_ELAPSED_INJ, INJ_LENGTH;

			if(EnginePtr->inj_start >= 360) // Arbitrary value to decide that inj_start is before TDC
			{
				if(THETA >= 360) THETA_ELAPSED_INJ = THETA - EnginePtr->inj_start; // THETA is before TDC
				else THETA_ELAPSED_INJ = THETA + (720 - EnginePtr->inj_start); // THETA is after TDC

				INJ_LENGTH = (720 - EnginePtr->inj_start) + EnginePtr->inj_stop;
			}
			else
			{
				THETA_ELAPSED_INJ = THETA - EnginePtr->inj_start;	// If inj_start is after TDC.
				INJ_LENGTH = EnginePtr->inj_stop - EnginePtr->inj_start;
			}

			m_finj_old = m_finj;
			m_finj = (THETA_ELAPSED_INJ/INJ_LENGTH)*EnginePtr->m_f0;	// Fuel mass injected so far this cycle
			m_finj_this_iter = m_finj - m_finj_old;						// Fuel mass injected this iteration

			EnginePtr->CYCLE_FUEL += m_finj_this_iter;					// Add to cumulative engine CYCLE_FUEL

			INJ_PREV_ITER = true;										// Injection occured this iteration
		}
		else // Valves closed but outside the injection period
		{
			if(INJ_PREV_ITER) // Have just finished the injection period
			{
				m_finj_old = m_finj;
				m_finj = EnginePtr->m_f0;	// Inject any remaining fuel the first iteraion following the end of injection period
				m_finj_this_iter = m_finj - m_finj_old;
				
				INJ_PREV_ITER = false;
			}
			else m_finj_this_iter = 0;
		}

		double m_fb_old = m_fb;
		m_fb = SingleZoneWatson(pPpt);	// Returns new mass of fuel to be burned up to current point
		del_m_fb = m_fb - m_fb_old;		// Mass of fuel burned this iteration (kg)
		if(del_m_fb<0) del_m_fb = 0;

		if(del_m_fb > m_fub)			// If the mass of fuel predicted to be burned is greater than that available...
		{
			del_m_fb = m_fub;			// ...burn all that is available
			m_fb = m_fb_old + del_m_fb;	// Recalculate the mass of fuel burned
		}
		m_fub = m_finj - m_fb;	// Resolve the mass of unburnt fuel still in the cylinder	
		
//		m_air_used = (EnginePtr->x + EnginePtr->y/4)*del_m_fb;	

		if(EnginePtr->ADD_MFB) MC += del_m_fb;	// Add mass of fuel burned to MC if desired
		// Subtract the mass of air required for stoichiometric combustion (assumes MC is mass of air)
		QLHVrate = del_m_fb*EnginePtr->lhv/DT;	// Calculate combustin heat release rate
	}
	else QLHVrate = 0;
	return QLHVrate;
}

double CCylinder::SingleZoneWatson(CProperties* pPpt)
{
	double m_fb_temp;
	if(THETA >= EnginePtr->ign_ca || THETA < EnginePtr->comb_stop) // If inside the combustion period
	{
		// Calculate length of combustion period in terms of crank angle
		double THETA_ELAPSED_COMB, COMB_LENGTH;
		
		if(EnginePtr->ign_ca >= 360) // Arbitrary value to decide that ign_ca is before TDC
		{
			if(THETA >= 360) THETA_ELAPSED_COMB = THETA - EnginePtr->ign_ca; // THETA is before TDC
			else THETA_ELAPSED_COMB = THETA + (720 - EnginePtr->ign_ca); // THETA is after TDC

			COMB_LENGTH = (720 - EnginePtr->ign_ca) + EnginePtr->comb_stop;
		}
		else
		{
			THETA_ELAPSED_COMB = THETA - EnginePtr->ign_ca;	// If ign_ca is after TDC.
			COMB_LENGTH = EnginePtr->comb_stop - EnginePtr->ign_ca;
		}
	
		double t_prime = THETA_ELAPSED_COMB/COMB_LENGTH;	// Fraction through combustion period
		double t_id = 1000*(((EnginePtr->ign_ca - EnginePtr->inj_start)/360)/EnginePtr->reveng); // Ignition delay in milliseconds
		double beta = 1 - EnginePtr->comb_a*pow(phiFAoverall, EnginePtr->comb_b)/pow(t_id, EnginePtr->comb_c);
		double K1 = 2 + 1.25e-8*pow((t_id*(EnginePtr->reveng*60)), 2.4); // Engine speed converted to rpm for K1.
		double K2 = 5000;
		double K3 = 14.2/pow(phiFAoverall, 0.644);
		double K4 = 0.79*pow(K3, 0.25);
		double f1 = 1 - pow((1 - pow(t_prime, K1)), K2);
		double f2 = 1 - exp(-K3*pow(t_prime, K4));
//cout << "exp(1) = " << exp(1) << endl;
//cout << "t_prime = " << t_prime << endl;
//cout << "f2 = " << f2 << endl;	
		m_fb_temp = ((beta*f1 + (1 - beta)*f2)*EnginePtr->m_f0);///t_prime;
	}
	else m_fb_temp = 0;	// Valves closed but outside the combustion period
	return m_fb_temp;
}

bool CCylinder::AllClosed(void)
{
	int i;
	bool ALL_CLOSED = true;
	for(i=0; i<EnginePtr->NEXVALVES; ++i) if(ExhaustValve[i].Get_open()) ALL_CLOSED = false;
	for(i=0; i<EnginePtr->NINVALVES; ++i) if(IntakeValve[i].Get_open()) ALL_CLOSED = false;
	return ALL_CLOSED;
}

/*
void CCylinder::CylinderUpdate(double dt,int cylnum,CPipe* InPipe,CPipe* ExPipe,
								CProperties* ppt, bool stopped)
{
	ncyl=ppt->ncyl;
	if (cylnum<(ncyl-1))	// For all cylinders except for the last.
	{pipenum=2*cylnum+1;}
	if (cylnum==(ncyl-1))	// For the last cylinder, or the only cylinder.
	{pipenum=2*cylnum;}

	in_flow_status=3;
	ex_flow_status=3;

	// Combustion chamber area:
	m_Fcc = m_Fch + Fc + (PI*ppt->B*ppt->S/2)*(m_R + 1 - cos(thetar) - sqrtawc(pow(m_R,2) - pow(sin(thetar),2))); 

	if( ((fmod(thetad,720) < ppt->EVO || fmod(thetad,720) >= EVC) && (fmod(thetad,720) < ppt->IVO || fmod(thetad,720) >= IVC)) )
	// If both valves are closed. If thetad is < 0 then fortunately we also get into this section.
	{
		if(openlastrun==true)CylinderOncePerCycle(ppt, stopped);
		IVOpen=false;
		InF1=0;
		EVOpen=false;
		OutF1=0;

		s = ppt->S/2*cos(thetar) + sqrtawc(pow(ppt->L,2) - pow(ppt->S/2,2)*pow(sin(thetar),2)); // Current crank centre to piston pin displacement.
		ss = ssmin + (s_TDC - s); // Current displacement of piston face from roof of cylinder downwards.
	
		// Instantaneous piston speed:
		dsdt = -1*ppt->S*ppt->N*PI*sin(thetar) * (1 + cos(thetar)/sqrtawc(pow(2*ppt->L/ppt->S, 2) - pow(sin(thetar),2)));
		dssdt = -dsdt;
		Sp_bar=2*ppt->S*ppt->N;		// Mean piston speed.
		V=ss*Fc;					// Current combustion chamber volume.
		dVdt=dssdt*Fc;

		// Flag the point at which TDC is passed:
		if((fmod(thetad,720) - fmod(thetad_old,720))<0){crossed=true;}

		// Obtain amount of fuel injected so far this cycle:
		if( ((fmod(thetad,720)>=ppt->inj_start || fmod(thetad,720)<ppt->inj_stop))
			&& thetad>ppt->EVO )		// If inside the injection period and past EVO the first time:
		{
			if(ppt->inj_start>=360) // Arbitrary value to decide that ignition is before TDC.
			{
				if(crossed==false){thetad_elapsed_inj=fmod(thetad,720) - ppt->inj_start;} // While before TDC.
				else{thetad_elapsed_inj=fmod(thetad,720) + (720 - ppt->inj_start);}
			}
			else{thetad_elapsed_inj=fmod(thetad,720) - ppt->inj_start;}	// If injection start is after TDC.
	
			m_mfinj_total_old = m_mfinj_total;
			m_mfinj_total = (thetad_elapsed_inj/dthetad_inj)*mf0;
			m_mfinj=m_mfinj_total-m_mfinj_total_old;
			injecting_last_run=true;
		}
		else // Valves closed but outside the injection period:
		{
			if(injecting_last_run==true)
			{
				m_mfinj_total_old = m_mfinj_total;
				m_mfinj_total=mf0;  // Inject any remaining fuel the first iteraion following the end of injection period.
				m_mfinj=m_mfinj_total-m_mfinj_total_old;
				injecting_last_run=false;
			}
			else{m_mfinj=0;}
		}

		m_mfub = m_mfinj_total - m_mfb;							// Resolve the mass of unburnt fuel still in the cylinder:	

		if(thetad>ppt->EVO)CylinderCombustion(dt, ppt);	
		// Sets new mass of fuel burned up to now, and the increment to be burned this iteration.

		if(m_dmfb>m_mfub)				// If the mass of fuel predicted to be burned is greater than that available...
		{
			m_dmfb = m_mfub;			// ...burn all that is available.
			m_mfb = m_mfb_old + m_dmfb; // Need to recalculate the mass of fuel burned.
			if(m_dmfb==0)				// If there is no fuel available:
			{fuel_status=1;}			// NO FUEL TO BURN
			else {fuel_status=2;}		// else BURNING ALL FUEL AVAILABLE
		}
		else							// If there is the fuel predicted available...
		{
			if(fmod(thetad,720)<ppt->ignition && fmod(thetad,720)>IVC)	// If we are between IVC and ignition:
			{fuel_status=5;}			// VALVES CLOSED BUT COMBUSTION NOT YET STARTED
			else
			{
				if(m_dmfb==0)			
				{fuel_status=3;}		// COMBUSTION STOPPED		// This isn't exactly right!!!!!
				else {fuel_status=4;}	// BURNING SOME FUEL
			}
		}
		
		// Corresponding air used to combine with the increment in fuel burnt (stoichiometric calc.):
		m_air_stoich = (ppt->x + ppt->y/4)*(m_dmfb);	// Mass of air required for stoich combustion.
		m_air -= m_air_stoich;							// Subtract the mass of air required for combustion.

		// Cylinder cooling:
		if( (fmod(thetad,720)>=IVC && fmod(thetad,720)<ppt->ignition) ){C2=0;}else{C2=3.24E-3;} 
		// Within compression period (up to ignition) set C2=0. In combustion and expansion periods C2=3.24E-3.

		m_w = 2.28*Sp_bar + C2*Vd*m_Tr/(m_pr*m_Vr)*(m_pc - m_pcm); 
	//	m_w = 0.028*Sp_bar;// + C2*Vd*m_Tr/(m_pr*m_Vr)*(pc - pcm); 
		// Local average gas velocity in cylinder. C1=2.28 outside gas exchange process, C2=0 in gas exchange process.
		m_wm = 2.28*Sp_bar; // Put pc = pcm in above expression!

		m_hcg = 3.26*pow(ppt->B,-0.2)*pow(m_pc/1000,0.8)*pow(m_Tg,-0.55)*pow(m_w,0.8);	// pc must be in kPa! Use OLD pressure.
		m_hcgm = 3.26*pow(ppt->B,-0.2)*pow(m_pcm/1000,0.8)*pow(m_Tgm,-0.55)*pow(m_wm,0.8);	// Motored version.
		//m_hcg = 0;
		//m_hcgm = 0;

		m_Twg = (m_hcg*m_Tg + m_coeffs*ppt->m_Tc)/(m_hcg + m_coeffs); // Use as initial estimate for iteration
		//m_Twgm = (m_hcgm*m_Tgm + m_coeffs*ppt->m_Tc)/(m_hcgm + m_coeffs);
		//cout << "Twg initial = " << m_Twg << "\n";
		//cout << "hcg initial = " << m_hcg << "\n";
		//cout << "Tg initial = " << m_Tg << "\n";
		//cout << "Tc initial = " << ppt->m_Tc << "\n";
		//cout << "coeffs initial = " << m_coeffs << "\n";
		
		// Now with radiative term, so have to iterate
		int counter=0;
		do
		{
			counter++;

			m_Twg_old = m_Twg;
			m_Twg = (m_hcg*m_Tg + ppt->m_sigma*ppt->m_epsilon*pow(m_Tg,4) + m_coeffs*ppt->m_Tc)
					/(m_hcg + m_coeffs + ppt->m_sigma*ppt->m_epsilon*pow(m_Twg,3));
			
			m_Twgm_old = m_Twgm;
			m_Twgm = (m_hcgm*m_Tgm + ppt->m_sigma*ppt->m_epsilon*pow(m_Tgm,4) + m_coeffs*ppt->m_Tc)
					/(m_hcgm + m_coeffs + ppt->m_sigma*ppt->m_epsilon*pow(m_Twgm,3));
		}
		while(100*fabs(m_Twg - m_Twg_old)/m_Twg_old > 0.0001);

		//cout << "counter = " << counter << "\n";

		if( (fmod(thetad,720)>=ppt->ignition || fmod(thetad,720)<ppt->combustion_stop) )
		{
			m_Qht+=m_hcg*m_Fcc*(m_Tg - m_Twg)*dt;
			m_Qch+=m_dmfb*ppt->LHV;
		}
		Qcoolrate=m_hcg*m_Fcc*(m_Tg - m_Twg);
		QLHVrate=m_dmfb*ppt->LHV/dt;

		// Calculate new pressure rate change:
		dpcdt = (m_dmfb*ppt->LHV/dt 
					- m_hcg*m_Fcc*(m_Tg - m_Twg)
					- ppt->mGammaCom/(ppt->mGammaCom-1)*m_pc*dVdt) * (ppt->mGammaCom-1)/V;

		// Motored pressure rate calculation:
		dpcmdt = (0 
					- m_hcgm*m_Fcc*(m_Tgm - m_Twgm)
					- ppt->mGammaCom/(ppt->mGammaCom-1)*m_pcm*dVdt) * (ppt->mGammaCom-1)/V; // Is this exactly right?

		// Calculate and check total mass of 'stuff' in cylinder (air + unburnt fuel + burnt fuel).
		m_total += m_mfinj;		// The only mass flow.

		// Primitive variables calculation:
		m_pc+=dpcdt*dt;
		m_pcm+=dpcmdt*dt;						// Motored pressure.
		m_Tg=m_pc*V/(m_total*ppt->mR);		// Note this R is for air, not air+fuel! Consider this point later.
		m_Tgm=m_pcm*V/(m_totalm*ppt->mR);		// Motored temperature.
			
		
		// Calculate incoming Riemann variables:
		InPipeClin=InPipe[pipenum].newReiman(dt,2,ppt);
		ExPipeClin=ExPipe[pipenum].newReiman(dt,0,ppt);
		// Calculate outgoing Riemann variables:  
		InPipeClout=InPipeClin;		// Since closed b.c.
		ExPipeClout=ExPipeClin;		// Since closed b.c.
		// Update boundary conditions:
		InPipe[pipenum].ReimanUpdate(InPipeClin,InPipeClout,InPipe[pipenum].GridPts-1,ppt);
		ExPipe[pipenum].ReimanUpdate(ExPipeClin,ExPipeClout,0,ppt);
	}
	else
	// One or both valves are open or we have not passed EVO for the first time. Marched forward using Euler integration (note, this
	// is an inaccurate method, if neccessary change to something like a 4th order Runga Kutta Method).
	{	
		// Calculate incoming Riemann variables from both intake and exhaust sides:
		ExPipeClin=ExPipe[pipenum].newReiman(dt,0,ppt);
		InPipeClin=InPipe[pipenum].newReiman(dt,2,ppt);

		//if ( ((CylinderGetAngleFromTDC()<ppt->IVO) || (CylinderGetAngleFromTDC()>=IVC)) )
		if ( ((fmod(thetad,720)<ppt->IVO) || (fmod(thetad,720)>=IVC)) )
		{
			IVOpen=false;
			InF1=0;
			InPipeClout=InPipeClin;
			in_flow_status=0;
		}
		else
		{
			IVOpen=true;
			InF1=MyInValve->ValveTellMeArea(thetad-ppt->IVO);
			InF2=InPipe[pipenum].Area[InPipe[pipenum].GridPts-1];
//			MyInValve->ValveTellMeClout(InPipeClin,pc,ppt->mAtmosPressure,
//				ppt->aref,InF2,ppt->mGamma,InF1/InF2,dmidt,InPipeClout);
			MyInValve->BensonValve(InPipeClout, dmidt, InPipeClin, InF1/InF2, ppt, m_pc, InF2);
			in_flow_status=MyInValve->ValveGetflow_status();

			// Motored:
			MyInValve->BensonValve(InPipeCloutm, dmimdt, InPipeClin, InF1/InF2, ppt, m_pcm, InF2);
			// Only make use of dmimdt in the above expression
		}
		if ( ((fmod(thetad,720)<ppt->EVO) || (fmod(thetad,720)>=EVC)) )
		{
			EVOpen=false;
			OutF1=0;
			ExPipeClout=ExPipeClin;
			ex_flow_status=0;		// NO FLOW
		}
		else
		{
			EVOpen=true;;
			OutF1=MyExValve->ValveTellMeArea(thetad-ppt->EVO);
			OutF2=ExPipe[pipenum].Area[0];
//			MyExValve->ValveTellMeClout(ExPipeClin,pc,ppt->mAtmosPressure,
//				ppt->aref,OutF2,ppt->mGamma,OutF1/OutF2,dmedt,ExPipeClout);
			MyExValve->BensonValve(ExPipeClout, dmedt, ExPipeClin, OutF1/OutF2, ppt, m_pc, OutF2);
			ex_flow_status=MyExValve->ValveGetflow_status();
			
			// Motored:
			MyExValve->BensonValve(ExPipeCloutm, dmemdt, ExPipeClin, OutF1/OutF2, ppt, m_pcm, OutF2);
			// Only make use of dmimdt in the above expression
		}
		
		InPipe[pipenum].ReimanUpdate(InPipeClin,InPipeClout,InPipe[pipenum].GridPts-1,ppt);
		ExPipe[pipenum].ReimanUpdate(ExPipeClin,ExPipeClout,0,ppt);
			
		if (dmidt>0.0){ai=SpeedOfSound(m_Tg);}else{ai=SpeedOfSound(InPipe[pipenum].T[InPipe[pipenum].GridPts-1]);}
		if (dmedt>0.0){ae=SpeedOfSound(m_Tg);}else{ae=SpeedOfSound(ExPipe[pipenum].T[0]);}

		// Motored values:
		if (dmimdt>0.0){aim=SpeedOfSound(m_Tgm);}else{aim=SpeedOfSound(InPipe[pipenum].T[InPipe[pipenum].GridPts-1]);}
		if (dmemdt>0.0){aem=SpeedOfSound(m_Tgm);}else{aem=SpeedOfSound(ExPipe[pipenum].T[0]);}

		s = ppt->S/2*cos(thetar) + sqrtawc(pow(ppt->L,2) - pow(ppt->S/2,2)*pow(sin(thetar),2)); // Current crank centre to piston pin displacement.
		ss = ssmin + (s_TDC - s); // Current displacement of piston face from roof of cylinder downwards.
		// Instantaneous piston speed:
		dsdt = -1*ppt->S*ppt->N*PI*sin(thetar) * (1 + cos(thetar)/sqrtawc(pow(2*ppt->L/ppt->S, 2) - pow(sin(thetar),2)));
		//Differentiating the equation for 'ss' gives:
		dssdt = -dsdt;
		Sp_bar=2*ppt->S*ppt->N;						// Mean piston speed.
		V=ss*Fc; // Current combustion chamber volume.
		dVdt=dssdt*Fc;

		// Cylinder cooling:
		m_w = 6.18*Sp_bar; // Local average gas velocity in cylinder. C1=6.18, C2=0 in gas exchange process.;
		m_hcg = 3.26*pow(ppt->B,-0.2)*pow(m_pc/1000,0.8)*pow(m_Tg,-0.55)*pow(m_w,0.8);	// pc must be in kPa!
		m_hcgm = 3.26*pow(ppt->B,-0.2)*pow(m_pcm/1000,0.8)*pow(m_Tgm,-0.55)*pow(m_w,0.8);	// m_wm = m_w!
		//m_hcg = 0;
		//m_hcgm = 0;

		// Radiation terms???
		m_Twg = (m_hcg*m_Tg + m_coeffs*ppt->m_Tc)/(m_hcg + m_coeffs);
		m_Twgm = (m_hcgm*m_Tgm + m_coeffs*ppt->m_Tc)/(m_hcgm + m_coeffs);
	
		Qcoolrate=m_hcg*m_Fcc*(m_Tg - m_Twg);
		
		// Calculate new pressure rate change:
		dpcdt = (pow(ai,2.0)*(-dmidt) - pow(ae,2.0)*dmedt
					- (ppt->mGamma-1)*m_hcg*m_Fcc*(m_Tg - m_Twg)
					- (ppt->mGamma)*m_pc*dVdt)
					/V;

		// Motored version
		dpcmdt = (pow(aim,2.0)*(-dmimdt) - pow(aem,2.0)*dmemdt
					- (ppt->mGamma-1)*m_hcgm*m_Fcc*(m_Tgm - m_Twgm)
					- (ppt->mGamma)*m_pcm*dVdt)
					/V;

		dmdt=(-dmidt-dmedt);
		dmmdt=(-dmimdt-dmemdt);
		m_air+=dmdt*dt;
		m_airm+=dmmdt*dt;
		m_total=m_air;
		m_totalm=m_airm;
		m_pc+=dpcdt*dt;
		m_pcm+=dpcmdt*dt;
		m_Tg=m_pc*V/(m_air*ppt->mR);
		m_Tgm=m_pcm*V/(m_air*ppt->mR);

		fuel_status=0;
		openlastrun=true;
	}

	// Stuff needed to be done each iteration:
	
	if(m_Tg>m_Tg_max)m_Tg_max=m_Tg;		// Keep a record of maximum cylinder temperature achieved over whole simulation.
	if(m_Twg>m_Twg_max)m_Twg_max=m_Twg;		// Keep a record of maximum cylinder wall temperature achieved over whole simulation.
	
	if(m_pc>m_pc_max)m_pc_max=m_pc;		// Keep a record of maximum cylinder pressure achieved over whole cycle (resets at IVC).
	
	
	m_Torque=(m_pc*Fc)*sin(thetar);

//	FA_actual = m_fuel_unburnt/m_air;
//	phiFA = FA_actual/FA_stoich;
	m_Pinstant=m_pc*Fc*dssdt;
//	m_Pinstant2=pc*dVdt;
	m_CycleWork+=m_pc*dVdt*dt;
	m_dWork=m_pc*dVdt*dt;				// Work this iteration.
	m_CycleTime+=dt;
	//cout << "InF1 = " << InF1 << endl;
}
*/




/*
void CCylinder::HCyl(CProperties* pPpt, double DELZe, double DELZi)
{	
/*	
	if(ZREV>180.0*this->EnginePtr->CYCLE*rev_number)
	{
		if(!wait)
		{
			// Calculate mean entropy levels once per cycle

			// Model 1
		//	aA_1_mean = sum_mc_aA_1/sum_mc;
		//	cout << "Cylinder [" << this->ID << "]: aA_1_mean = " << aA_1_mean << endl;
		//	char pause;
		//	cin >> pause;
			// Reset entropy calculations every cycle
		//	sum_mc_aA_1 = 0;

			// Model 2
		//	aA_2_mean = sum_mc_aA_2/sum_mc;
		//	cout << "Cylinder [" << this->ID << "]: aA_2_mean = " << aA_2_mean << endl;
		//	cin >> pause;
			// Reset entropy calculations every cycle
		//	sum_mc_aA_2 = 0;

			// Model 3
			// Choose aA_mean from either model 1 or model 2
		//	aA_3_mean = aA_1_mean;
		//	aA_3_mean = aA_2_mean;

		//	K_entropy = rAREFe/aA_3_mean; 
			// This rAREFe should be the reference speed of sound from the
			// previous cycle

			// For all models
		//	sum_mc = 0;

			// Now select a model for AREFe or not at all
		//	rAREFe = aA_1_mean;
		//	rAREFe = aA_2_mean;
		//	rAREFe = aA_3_mean; // Also use K_entropy in the valve
		}

		//reset_release = true; 
		rev_number += 1; 
		wait=false;
	}
*/
/*
	// Average entropy levels
	sum_mc += DMEDT;

	// Model 1
	aA_1 = pow(PC/pPpt->PREF, (pPpt->gammaAir()-1)/(2*pPpt->gammaAir()))*pPpt->aref_e;
	sum_mc_aA_1 += DMEDT*aA_1;	// This will be 0 when exhaust valve is closed

	// Model 2
	P1 = TC*(pBN[exhaust_side]->p_dash*pPpt->PREF)/pBN[exhaust_side]->T 
			+//- 
			((pPpt->gammaAir()-1)/2)*pow(fabs(pBN[exhaust_side]->U)*pPpt->aref_e,2)*(pBN[exhaust_side]->p_dash*pPpt->PREF)/(R_air*pBN[exhaust_side]->T);
	P1 = pBN[exhaust_side]->p_dash*pPpt->PREF;
	aA_2 = pow(P1/pPpt->PREF, (pPpt->gammaAir()-1)/(2*pPpt->gammaAir()))*pPpt->aref_e;
	sum_mc_aA_2 += DMEDT*aA_2;	// This will be 0 when exhaust valve is closed

	// Model 3
	// Chooses aA_mean from either model 1 or model 2
*/
/*
	return;
}
*/

void CCylinder::LoadCylinderData(char* InputFile)
//--------------------------------------------------//
// Load cylinder data file							//
// -----------------------							//
// Loads cylinder pressure, temperature from file, 	//
// as a function of theta (degrees).				//
//--------------------------------------------------//
{
	//cout << InputFile << endl;
	int c = 0;
	//float temp;
	double temp;
	int num_columns = 3;
	int col, row;
	datapoints_cyl = 0;
	
	FILE* stream;
	stream = fopen(InputFile, "r");

	if(stream == NULL)
	{
		printf("Error opening valve lift arrays file\n");
		exit(1);
	}
	else
	{
		fseek(stream, 0L, SEEK_SET);	// Set pointer to beginning of file
		do
		{
			fscanf(stream, "%d", &c);	// Scans for the line number (ints), puts it in c
//cout << "c = " << c << endl;
			if(c>datapoints_cyl) datapoints_cyl = c;
			fscanf(stream, "%lf", &temp);		// Runs over the theta value
			fscanf(stream, "%lf", &temp);		// Runs over the pressure value
			fscanf(stream, "%lf", &temp);		// Runs over the temp. value
			//cin >> pause;
		}while(fscanf(stream, "%l")!=EOF);
//cout << "Cylinder datapoints = " << datapoints_cyl << endl;
		//cin >> pause;
		// Can now dimension arrays
		cyl_data = new double* [num_columns]; 
		for(col=0; col<num_columns; ++col) cyl_data[col] = new double [datapoints_cyl];

		fseek(stream, 0L, SEEK_SET);	// Reset pointer to beginning of file

		for (row=0; row<datapoints_cyl; row++)
		{
			fscanf(stream, "%d", &c);	// Scans past the line number
			for (col=0; col<num_columns; col++)
			{
				// Finds the angle and area values
				fscanf(stream, "%lf", &cyl_data[col][row]);	
//cout << "cyl_data[col][row] = " << cyl_data[col][row] << endl;
			}
		}
	}
	fclose(stream);
}

void CCylinder::Interpolate(CProperties* pPpt, double ca)
//--------------------------------------------------//
// Read cylinder pressure, temperature from file	//
// ---------------------------------------------	//
// Interpolates for pressure and temperature		//
// from a data file.								//
//--------------------------------------------------//
{
	int THETA = 0; // Labels
	int P = 1;
	int T = 2;
	int row;

	// Check whether the lookup asked for is valid
	row = 0;
	if(ca<cyl_data[THETA][0] || ca>cyl_data[THETA][datapoints_cyl-1])
	{
		//cout << "Cylinder data interpolation: lookup out of range." << endl;
		PC = 1.0;
		TC = 300.0;
	}
	else
	{
		// Find the two datapoints either side of lookup
		while(!(ca>=cyl_data[THETA][row] && ca<cyl_data[THETA][row+1])) ++row;
		// While lookup is not between consecutive datapoints, keep searching
		
		// Interpolate for PC, TC. cyl_data[THETA][row] is the lookup, cyl_data[P][row] is the pressure etc.
		PC = cyl_data[P][row] + (ca - cyl_data[THETA][row])/(cyl_data[THETA][row+1] - cyl_data[THETA][row])
								*(cyl_data[P][row+1] - cyl_data[P][row]);

		TC = cyl_data[T][row] + (ca - cyl_data[THETA][row])/(cyl_data[THETA][row+1] - cyl_data[THETA][row])
								*(cyl_data[T][row+1] - cyl_data[T][row]);
	}
}


void CCylinder::AirReceiver(CProperties* pPpt, double& rDMIDT, double PAIR, double TAIR, double PCYL, double TCYL, double AAIR, double ACYL, double F2)
{
	double B, C, D, RCR, FRM, T;
	double RM = PCYL/PAIR;
	if(RM>1.0) RM = 1.0/RM;
	if(RM==1.0){rDMIDT = 0.0; return;}	// Pressures balance so no flow
	if(PCYL>PAIR) T = TCYL; else T = TAIR;

	// Test for choked flow
	D = pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1.0);
	RCR = pow(2.0/(pPpt->gammaAir(T)+1.0), D);
	if(RM>RCR)	// Subsonic flow
	{
		B = 2.0/pPpt->gammaAir(T);
		C = (pPpt->gammaAir(T)-1.0)/pPpt->gammaAir(T);
		D = (2.0*pow(pPpt->gammaAir(T), 2))/(pPpt->gammaAir(T)-1.0);
		FRM = sqrt(D*pow(RM, B)*(1.0 - pow(RM, C)));
	}
	else		// Choked flow
	{
		B = (pPpt->gammaAir(T)+1.0)/(2.0*(pPpt->gammaAir(T)-1.0));
		FRM = pPpt->gammaAir(T)*pow(2.0/(pPpt->gammaAir(T)+1.0), B);
	}

	// Mass flow rate
	// Test flow direction
	if(PCYL>PAIR)	// Outflow
	{
		rDMIDT = -PCYL*F2*FRM*1.0E5/ACYL;
	}
	else			// Inflow
	{
		rDMIDT = PAIR*F2*FRM*1.0E5/AAIR;
	}
	return;
}

void CCylinder::PrintToScreen(CProperties* pPpt)
{
	int i;
	cout << Underline(Identify(), "=", "\t");	
	cout << "\tCrank angle relative to TDCF, THETA\t=\t" << THETA << Deg() << endl;
	cout << endl;

	for(i=0; i<EnginePtr->NEXVALVES; ++i)
	{
		ExhaustValve[i].PrintToScreen(pPpt);
		cout << "\t\tMass flow, DMEDT[" << i << "]\t\t=\t" << DMEDT[i] << " kg/s" << endl;
		cout << endl;
	}

//	if(EnginePtr->IAIR) // If there is an intake valve
//	{
		for(i=0; i<EnginePtr->NINVALVES; ++i)
		{
			cout << endl;
			IntakeValve[i].PrintToScreen(pPpt);
			cout << "\t\tMass flow (DMIDT[" << i << "])\t\t=\t" << DMIDT[i] << " kg/s" << endl;
			cout << endl;
		}
//	}
	
	if(EnginePtr->MODEL == WATSON)
	{
		cout << Underline("Combustion", "-", "\t");
		cout << "\tFuel mass inj. this cycle (m_finj)\t=\t" << m_finj << " kg\n";
		cout << "\tUnburnt fuel mass in cylinder (m_fub)\t=\t" << m_fub << " kg\n";
		cout << "\tMass of fuel burned this cycle (m_fb)\t=\t" << m_fb << " kg\n";
		cout << "\tCurrent heat release rate (QLHVrate)\t=\t" << QLHVrate << " J/s\n";
		cout << endl;
	}

	if(EnginePtr->HT_MODEL == WOSCHNI)
	{
		cout << Underline("Heat transfer", "-", "\t");

		cout << endl;
	}

	cout << Underline("Working variables", "-", "\t");
	cout << "\tVol. rate of change (DVCDT)\t\t=\t" << DVCDT*1000 << " l/s\n";
	cout << "\tPressure rate of change (DPCDT)\t\t=\t" << DPCDT << " bar/s\n";
	cout << "\tCyl. pressure (PC)\t\t\t=\t" << PC << " bar\n";
	cout << "\tCyl. temp. (TC)\t\t\t\t=\t" << TC << " K\n";
	cout << "\tCylinder volume, VC\t\t\t=\t" << VC*1000 << " litres" << endl;;
	cout << "\tCyl. mass (MC)\t\t\t\t=\t" << MC << " kg\n";
	cout << endl;

	if(EnginePtr->HT_MODEL==WOSCHNI)
	{
		cout << Underline("Motored values", "-", "\t");
		cout << "\tMot. pressure rate change (DPC_motDT)\t=\t" << DPC_motDT << " bar/s\n";
		cout << "\tMot. cyl. pressure (PC_mot)\t\t=\t" << PC_mot << " bar\n";
		cout << "\tMot. cyl. temp. (TC_mot)\t\t=\t" << TC_mot << " K\n";
		cout << "\tMot. cyl. mass (MC_mot)\t\t\t=\t" << MC_mot << " kg\n";
		cout << endl;
	}

	cout << Underline("Cycle values", "-", "\t");
	cout << "\tNo. of cycles completed (NCYCLES)\t=\t" << NCYCLES << endl;
	cout << endl;
	cout << "\tWk done previous cycle\n";
	cout << "\t- gross (W_cig_prev)\t\t\t=\t" << W_cig_prev << " J\n";
	cout << "\t- net (W_cin_prev)\t\t\t=\t" << W_cin_prev << " J\n";
	cout << endl;
	cout << "\tWk done so far this cycle\n";
	cout << "\t- gross (W_cig)\t\t\t\t=\t" << W_cig << " J\n";
	cout << "\t- net (W_cin)\t\t\t\t=\t" << W_cin << " J\n";
	cout << endl;
	cout << "\tMotored wk done previous cycle\n";
	cout << "\t- gross (W_cig_mot_prev)\t\t=\t" << W_cig_mot_prev << " J\n";
	cout << "\t- net (W_cin_mot_prev)\t\t\t=\t" << W_cin_mot_prev << " J\n";
	cout << endl;
	cout << "\tMotored wk done so far this cycle\n";
	cout << "\t- gross (W_cig_mot)\t\t\t=\t" << W_cig_mot << " J\n";
	cout << "\t- net (W_cin_mot)\t\t\t=\t" << W_cin_mot << " J\n";
	cout << endl;
	cout << endl;
}

void CCylinder::ReadInput(CProperties* pPpt, char *InputFile)
{
	int last_entry;
	CommonReadInput(InputFile, pPpt->NUM_PARAMS, this->labels, this->values, this->strings, last_entry);
	for(int r=0; r<last_entry+1; ++r)
	{
		// ====================================================================================================
		// Parameter file for Cylinder [0]
		// ====================================================================================================

		// Optional object description (max. 500 characters - use underscores for spaces)
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "strDesc") == 0) strDesc = strings[r];

		// Operation
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "offset") == 0) offset = values[r];

		// Measurements
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "PRINT_TO_FILE_MOV") == 0) PRINT_TO_FILE_MOV = DoubleToBool(values[r]);
		if(strcmp(labels[r], "mov_freq") == 0) mov_freq = int(values[r]);
		if(strcmp(labels[r], "print_from_cycle") == 0) print_from_cycle = int(values[r]);
		
		// Parameters to record
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "num_props_measured") == 0) num_props_measured = int(values[r]);
		if(strcmp(labels[r], "CYL_PRESSURE") == 0) CYL_PRESSURE = DoubleToBool(values[r]);
		if(strcmp(labels[r], "CYL_TEMPERATURE") == 0) CYL_TEMPERATURE = DoubleToBool(values[r]);
	}
}

void CCylinder::ListProperties(CProperties* pPpt)
{
	if (pPpt->SHOW_calls) { pPpt->Out(Identify()); pPpt->Out(".ListProperties\n"); }

	// ====================================================================================================
	// Parameter file for Exhaust End Environment []
	// ====================================================================================================

	pPpt->Out(Underline(Identify(), "=", "\t", strDesc)); pPpt->Out("\n");

	// Operation
	// ----------------------------------------------------------------------------------------------------
	pPpt->Out("\tCrank angle offset, offset\t\t\t=\t"); pPpt->Out(offset); pPpt->Out(Deg()); pPpt->Out("\n");
	
	// Measurements
	// ----------------------------------------------------------------------------------------------------
	if (PRINT_TO_FILE_MOV) {
		pPpt->Out("\tRecording data for animation, PRINT_TO_FILE_MOV\t=\t"); pPpt->Out(TrueOrFalse(PRINT_TO_FILE_MOV)); pPpt->Out("\n");
		pPpt->Out("\tMovie sampling frequency once per, mov_freq\t=\t"); pPpt->Out(mov_freq); pPpt->Out(" timesteps"); pPpt->Out("\n");
		pPpt->Out("\tPrint to file from of cycle, print_from_cycle\t=\t"); pPpt->Out(print_from_cycle); pPpt->Out("\n");
	}
	else pPpt->Out("\tNot recording movie data, PRINT_TO_FILE_MOV\t=\t"); pPpt->Out(TrueOrFalse(PRINT_TO_FILE_MOV)); pPpt->Out("\n");

	// Parameters to record
	// ----------------------------------------------------------------------------------------------------
	if (num_props_measured > 0) {
		pPpt->Out("\tRecording:\n");
		if (CYL_PRESSURE) pPpt->Out("\t- cylinder pressure"); pPpt->Out("\n");
		if (CYL_TEMPERATURE) pPpt->Out("\t- cylinder temperature"); pPpt->Out("\n");
	}
	pPpt->Out("\n");
	pPpt->Out("\n");
}

void CCylinder::PrintToFile(CProperties* pPpt, int timestep, double time, double ca)
{
	if(timestep==1)
	{
		fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR CYLINDER", "-", ID));
		fprintf(OUTPUT_FILE,"\n");
		fprintf(OUTPUT_FILE,"%s", "Time(s)");
		fprintf(OUTPUT_FILE,"\t%s%c%s\t%s%c%s", "Elapsed (", Deg(), "CA)", "Crank angle (", Deg(), "CA)");
		fprintf(OUTPUT_FILE,"\t\t%s\t%s\t%s\t%s\t%s", "PC (bar)", "TC (K)", "VC (m^3)", "MC (kg)", "QLHVrate (J/s)");
		fprintf(OUTPUT_FILE,"\n");
	}

	if(timestep%pPpt->freq==0) // Print data at the main default specified sampling frequency
	{
		fprintf(OUTPUT_FILE,"%f\t%f\t%f\t\t%f\t%f\t%f\t%f\t%f",
				time,
				EnginePtr->ca_elapsed,
				THETA,
				PC,
				TC,
				VC,
				MC,
				QLHVrate
				);
		fprintf(OUTPUT_FILE,"\n");
	}
}

void CCylinder::PrintToFileMovie(CProperties* pPpt, int timestep, double time, double ca)
{
	if(PRINT_TO_FILE_MOV)
	{
		fprintf(FILE_OVERALL_MOV,"%.6f\t", time);
		fprintf(FILE_OVERALL_MOV,"%.6f\t", ca);
		fprintf(FILE_OVERALL_MOV,"%f\t", PC);
		fprintf(FILE_OVERALL_MOV,"%f\t", TC);
		for(int i=0; i<EnginePtr->NEXVALVES; ++i)
			fprintf(FILE_OVERALL_MOV,"%f\t", ExhaustValve[i].Get_final_lift());
		fprintf(FILE_OVERALL_MOV,"\n");
	}
}

void CCylinder::SetupFiles(CProperties* pPpt, string parent_assy_res_dir)
{
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

	std::string bcname_str = "res_cyl";
	OUTPUT_FILE = fopen(ConstructString(pPpt, RES_DIR, bcname_str, ID), "w");

	bcname_str = "res_cyc_cyl";
	FILE_CYCLE = fopen(ConstructString(pPpt, RES_DIR, bcname_str, ID), "w");

	bcname_str = "res_cyl_mov";
	FILE_OVERALL_MOV = fopen(ConstructString(pPpt, RES_DIR, bcname_str, ID), "w");
}

void CCylinder::CloseFiles()
{
	fclose(OUTPUT_FILE);
	fclose(FILE_CYCLE);
	if(PRINT_TO_FILE_MOV) fclose(FILE_OVERALL_MOV);
}

char* CCylinder::Identify()
// ============================================================ //
// Returns identification of the current object					//
// ============================================================ //
{
	std::string sz;
	sz = "Assembly [";
	sz += IntToString(this->EnginePtr->AssyID);
	sz += "], ";
	sz += "Engine [";
	sz += IntToString(this->EnginePtr->ID);
	sz += "], ";
	sz += "Cylinder [";
	sz += IntToString(ID);
	sz += "]";

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str());
	return szz;
}