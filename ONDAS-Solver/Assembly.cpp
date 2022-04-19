// Assembly.cpp: implementation of the CAssembly class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Assembly.h"
#include "Properties.h"
#include "Tools.h"
#include <windows.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CAssembly::CAssembly()
{

}

CAssembly::~CAssembly()
{

}

void CAssembly::Initialise(CProperties* pPpt, int id)
{
	ID = id;
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Initialise\n");}
	string bcname_str = "ASSEMBLY";

	ReadAssyDir(pPpt, ConstructString(pPpt, pPpt->case_dir, bcname_str, ID)); // This reads and sets param_dir
	ReadInput(pPpt, ConstructString(pPpt, param_dir, bcname_str, true));

	pPpt->Out("\n");

	// Dimension pipe and boundary objects
	// ====================================================================================================
	Exhaust		= new CPipe[NEXPIPES];
	Intake		= new CPipe[NINPIPES];
	EX_ANEC		= new CAnechoic[NEXANEC];
	IN_ANEC		= new CAnechoic[NINANEC];
	EX_END		= new CEndEnvironment[NEXEND];
	IN_END		= new CEndEnvironment[NINEND];
	CE			= new CEndCap[NEXENDCAP];
	CI			= new CEndCap[NINENDCAP];
	TransmE		= new CTransmissive[NEXTRANSM];
	TransmI		= new CTransmissive[NINTRANSM];
	SE			= new CSudden[NEXSUD];
	SI			= new CSudden[NINSUD];
	JE			= new CJunction[NEXJUNCS];
	JI			= new CJunction[NINJUNCS];
	Eng			= new CEngine[NENGINES];
	Cyl			= new CCylinder[NCYLS];
	APLDevE		= new CAPLDev[NEXAPLDev];
	APLDevI		= new CAPLDev[NINAPLDev];
	Turbine		= new CTurbine[NTURBINE];
	ExVolume	= new CVolume[NEXVOLUMES];
	InVolume	= new CVolume[NEXVOLUMES];

	// Create a res_assembly directory for this assembly within numbered case results directory
	string assy_res_str;
	assy_res_str = "assembly";
	if(ID>=10) assy_res_str += int(ID/10) + 48;
	assy_res_str += (ID - int(ID/10)*10) + 48;
	assy_res_str += ".res";
	LPSECURITY_ATTRIBUTES attr; attr = NULL;
	CreateDirectory(StringToChar(pPpt->case_res_dir + assy_res_str), attr);
	this->assy_res_dir = pPpt->case_res_dir + assy_res_str + "\\"; // e.g. ../dat/case/case-name.001.res/assembly0.res/
}

void CAssembly::InitialiseBCs(CTime* pMyTime, CProperties* pPpt, CAssembly* &rAssy)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".InitialiseBCs\n");}

	// Initialise boundary conditions
	// ====================================================================================================
	int i, p;
	string calling_object_str;

	// End cap
	// ----------------------------------------------------------------------------------------------------
	for (i = 0; i < NEXENDCAP; ++i) {
		calling_object_str = "CE["; calling_object_str += to_string(i); calling_object_str += "]";
		CE[i].Initialise(pPpt, PipePtrs(EXEndCapPIPES_NUMS[i], rAssy, EXEndCapASSYS, EXEndCapPIPES, i, true),
			rAssy[EXEndCapASSYS[i][ONE_SIDE]].Exhaust, EXEndCapPIPES, EXEndCapPIPES_ENDS, EXEndCapENDCORR, i, true, 1, ID, this->assy_res_dir, calling_object_str);
	}

	for (i = 0; i < NINENDCAP; ++i) {
		calling_object_str = "CI["; calling_object_str += to_string(i); calling_object_str += "]";
		CI[i].Initialise(pPpt, PipePtrs(INEndCapPIPES_NUMS[i], rAssy, INEndCapASSYS, INEndCapPIPES, i, false),
			rAssy[INEndCapASSYS[i][ONE_SIDE]].Intake, INEndCapPIPES, INEndCapPIPES_ENDS, INEndCapENDCORR, i, false, 1, ID, this->assy_res_dir, calling_object_str);
	}

	// Anechoic ends
	// ----------------------------------------------------------------------------------------------------
	for (i = 0; i < NEXANEC; ++i) {
		calling_object_str = "EX_ANEC["; calling_object_str += to_string(i); calling_object_str += "]";
		EX_ANEC[i].Initialise(pPpt, PipePtrs(EXANECPIPES_NUMS[i], rAssy, EXANECASSYS, EXANECPIPES, i, true),
			rAssy[EXANECASSYS[i][ONE_SIDE]].Exhaust, EXANECPIPES, EXANECPIPES_ENDS, EXANECENDCORR, i, true, 1/*1*/, param_dir, ID, this->assy_res_dir, calling_object_str);
	}

	for (i = 0; i < NINANEC; ++i) {
		calling_object_str = "IN_ANEC["; calling_object_str += to_string(i); calling_object_str += "]";
		IN_ANEC[i].Initialise(pPpt, PipePtrs(INANECPIPES_NUMS[i], rAssy, INANECASSYS, INANECPIPES, i, false),
			rAssy[INANECASSYS[i][ONE_SIDE]].Intake, INANECPIPES, INANECPIPES_ENDS, INANECENDCORR, i, false, 1/*1*/, param_dir, ID, this->assy_res_dir, calling_object_str);
	}

	// End environments
	// ----------------------------------------------------------------------------------------------------
	for (i = 0; i < NEXEND; ++i) {
		calling_object_str = "EX_END["; calling_object_str += to_string(i); calling_object_str += "]";
		EX_END[i].Initialise(pPpt, PipePtrs(EXEPIPES_NUMS[i], rAssy, EXEASSYS, EXEPIPES, i, true),
			rAssy[EXEASSYS[i][ONE_SIDE]].Exhaust, EXEPIPES, EXEPIPES_ENDS, EXEENDCORR, i, true, 1, param_dir, ID, this->assy_res_dir, calling_object_str);
	}

	for (i = 0; i < NINEND; ++i) {
		calling_object_str = "IN_END["; calling_object_str += to_string(i); calling_object_str += "]";
		IN_END[i].Initialise(pPpt, PipePtrs(INEPIPES_NUMS[i], rAssy, INEASSYS, INEPIPES, i, false),
			rAssy[INEASSYS[i][ONE_SIDE]].Intake, INEPIPES, INEPIPES_ENDS, INEENDCORR, i, false, 1, param_dir, ID, this->assy_res_dir, calling_object_str);
	}

	// Transmissive boundary conditions
	// ----------------------------------------------------------------------------------------------------	
	for (i = 0; i < NEXTRANSM; ++i) {
		calling_object_str = "TransmE["; calling_object_str += to_string(i); calling_object_str += "]";
		TransmE[i].Initialise(pPpt, PipePtrs(EXTPIPES_NUMS[i], rAssy, EXTASSYS, EXTPIPES, i, true),
			rAssy[EXTASSYS[i][ONE_SIDE]].Exhaust, EXTPIPES, EXTPIPES_ENDS, EXTENDCORR, i, true, 1, param_dir, ID, this->assy_res_dir, NEXTRANSM, NEXPIPES, calling_object_str);
	}

	for (i = 0; i < NINTRANSM; ++i) {
		calling_object_str = "TransmI["; calling_object_str += to_string(i); calling_object_str += "]";
		TransmI[i].Initialise(pPpt, PipePtrs(INTPIPES_NUMS[i], rAssy, INTASSYS, INTPIPES, i, false),
			rAssy[INTASSYS[i][ONE_SIDE]].Intake, INTPIPES, INTPIPES_ENDS, INTENDCORR, i, false, 1, param_dir, ID, this->assy_res_dir, NINTRANSM, NINPIPES, calling_object_str);
	}

	// Sudden area changes
	// ----------------------------------------------------------------------------------------------------
	for (i = 0; i < NEXSUD; ++i) {
		calling_object_str = "SE["; calling_object_str += to_string(i); calling_object_str += "]";
		SE[i].Initialise(pPpt, PipePtrs(EXSPIPES_NUMS[i], rAssy, EXSASSYS, EXSPIPES, i, true),
			Exhaust, EXSPIPES, EXSPIPES_ENDS, EXSENDCORR, i, true, 2, ID, this->assy_res_dir, calling_object_str);
	}
	
	for (i = 0; i < NINSUD; ++i) {
		calling_object_str = "SI["; calling_object_str += to_string(i); calling_object_str += "]";
		SI[i].Initialise(pPpt, PipePtrs(INSPIPES_NUMS[i], rAssy, INSASSYS, INSPIPES, i, false),
			Intake, INSPIPES, INSPIPES_ENDS, INSENDCORR, i, false, 2, ID, this->assy_res_dir, calling_object_str);
	}

	// Junctions
	// ----------------------------------------------------------------------------------------------------
	for (i = 0; i < NEXJUNCS; ++i) {
		calling_object_str = "JE["; calling_object_str += to_string(i); calling_object_str += "]";
		JE[i].Initialise(pPpt, PipePtrs(EXJPIPES_NUMS[i], rAssy, EXJASSYS, EXJPIPES, i, true),
			Exhaust, EXJPIPES, EXJPIPES_ENDS, EXJENDCORR, i, true, EXJPIPES_NUMS[i], param_dir, res_dir, ID, this->assy_res_dir, calling_object_str);
	}
	
	for (i = 0; i < NINJUNCS; ++i) {
		calling_object_str = "JI["; calling_object_str += to_string(i); calling_object_str += "]";
		JI[i].Initialise(pPpt, PipePtrs(INJPIPES_NUMS[i], rAssy, INJASSYS, INJPIPES, i, false),
			Intake, INJPIPES, INJPIPES_ENDS, INJENDCORR, i, false, INJPIPES_NUMS[i], param_dir, res_dir, ID, this->assy_res_dir, calling_object_str);
	}

	// Engines
	// ----------------------------------------------------------------------------------------------------
	for (i = 0; i < NENGINES; ++i) {
		calling_object_str = "Eng["; calling_object_str += to_string(i); calling_object_str += "]";
		Eng[i].Initialise(pPpt, i, param_dir, ENG_FILE, this->NCYLS, this->NINPIPES, ID, this->assy_res_dir, calling_object_str);
	}

	// Cylinders
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NCYLS; ++i) {
		// Work out which connections are to exhaust pipes and which are to intake pipes
		bool* EX;
		EX = new bool [CYLPIPES_NUMS[i]];
		switch(CYLPIPES_NUMS[i])
		{
		case 0:
			break;
		case 1:			// One exhaust valve only
			EX[0] = true;
			break;
		case 2:			// Two exhaust valves, or one of each with exhaust first
			EX[0] = true;
			if(Eng[i].NEXVALVES==2) EX[1] = true;
			else EX[1] = false;
			break;
		case 3:			// Two exhaust, then one intake
			EX[0] = true;
			EX[1] = true;
			EX[2] = false;
			break;
		case 4:			// Two exhaust, then two intake
			EX[0] = true;
			EX[1] = true;
			EX[2] = false;
			EX[3] = false;
			break;
		default:
			cout << "Assembly->InitialiseBCs->Cylinders - problem with number of cylinders" << endl;
			exit(1);
		}
		//for(int j=0; j<CYLPIPES_NUMS[i]; ++j) cout << "EX[" << j << "] = " << TrueOrFalse(EX[j]) << endl;
		
		/*for (i = 0; i<NEXEND; ++i)
			EX_END[i].Initialise(pPpt, PipePtrs(EXEPIPES_NUMS[i], rAssy, EXEASSYS, EXEPIPES, i, true),
				rAssy[EXEASSYS[i][ONE_SIDE]].Exhaust, EXEPIPES, EXEPIPES_ENDS, EXEENDCORR, i, true, 1, param_dir, ID, this->assy_res_dir);
		
		for(i=0; i<NINEND; ++i)
			IN_END[i].Initialise(pPpt, PipePtrs(INEPIPES_NUMS[i], rAssy, INEASSYS, INEPIPES, i, false),
				rAssy[INEASSYS[i][ONE_SIDE]].Intake, INEPIPES, INEPIPES_ENDS, INEENDCORR, i, false, 1, param_dir, ID, this->assy_res_dir);

		*/

		/*Cyl[i].Initialise(pPpt,
						  PipePtrsCyl(CYLPIPES_NUMS[i], rAssy, CYLASSYS, CYLPIPES, i, EX /true?/),// false?),	// Exhaust side
						  // Intake side?
						  Exhaust, Intake, CYLPIPES, CYLPIPES_ENDS, CYLENDCORR, i, &Eng[0], param_dir, ID, this->assy_res_dir);
		*/
		
		//EXHAUST
		calling_object_str = "Cyl["; calling_object_str += to_string(i); calling_object_str += "]";
		Cyl[i].Initialise(pPpt,
			PipePtrsCyl(CYLPIPES_NUMS[i], rAssy, CYLASSYS, CYLPIPES, i, EX),
			Exhaust, CYLPIPES, CYLPIPES_ENDS, CYLENDCORR, i, EX[i], CYLPIPES_NUMS[i], &Eng[0], param_dir, ID, this->assy_res_dir, calling_object_str, NEXPIPES, NINPIPES);

		//INTAKE
		//Cyl[i].Initialise(pPpt,
		//	PipePtrsCyl(CYLPIPES_NUMS[i], rAssy, CYLASSYS, CYLPIPES, i, EX),
		//	rAssy[EXEASSYS[i][ONE_SIDE]].Exhaust, rAssy[INEASSYS[i][ONE_SIDE]].Intake, CYLPIPES, CYLPIPES_ENDS, CYLENDCORR, i, EX[i], CYLPIPES_NUMS[i], &Eng[0], param_dir, ID, this->assy_res_dir);
	}

	// Adiabatic pressure loss devices
	// ----------------------------------------------------------------------------------------------------
	for (i = 0; i < NEXAPLDev; ++i) {
		calling_object_str = "APLDevE["; calling_object_str += to_string(i); calling_object_str += "]";
		APLDevE[i].Initialise(pMyTime, pPpt, PipePtrs(EXAPLDevPIPES_NUMS[i], rAssy, EXAPLDevASSYS, EXAPLDevPIPES, i, true),
			Exhaust, EXAPLDevPIPES, EXAPLDevPIPES_ENDS, EXAPLDevENDCORR, i, true, 2, &TransmE[0], param_dir, ID,
			NEXTRANSM > 0 ? true : false, this->assy_res_dir, EX_END, TransmE, JE, calling_object_str);
	}
	
	for (i = 0; i < NINAPLDev; ++i) {
		calling_object_str = "APLDevI["; calling_object_str += to_string(i); calling_object_str += "]";
		APLDevI[i].Initialise(pMyTime, pPpt, PipePtrs(INAPLDevPIPES_NUMS[i], rAssy, INAPLDevASSYS, INAPLDevPIPES, i, false),
			Intake, INAPLDevPIPES, INAPLDevPIPES_ENDS, INAPLDevENDCORR, i, false, 2, &TransmI[0], param_dir, ID,
			NINTRANSM > 0 ? true : false, this->assy_res_dir, IN_END, TransmI, JI, calling_object_str);
	}

	// Turbines
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NTURBINE; ++i) {
		calling_object_str = "Turbine["; calling_object_str += to_string(i); calling_object_str += "]";
		Turbine[i].Initialise(pPpt, PipePtrs(TURBPIPES_NUMS[i], rAssy, TURBASSYS, TURBPIPES, i, true),
			Exhaust, TURBPIPES, TURBPIPES_ENDS, TURBPIPES_NUMS, TURBENDCORR, i, true, param_dir, ID, this->assy_res_dir, calling_object_str);
	}

	// Volumes
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXVOLUMES; ++i) {
		calling_object_str = "ExVolume["; calling_object_str += to_string(i); calling_object_str += "]";
		ExVolume[i].Initialise(pPpt, PipePtrs(EXVOLPIPES_NUMS[i], rAssy, EXVOLASSYS, EXVOLPIPES, i, true),
			Exhaust, EXVOLPIPES, EXVOLPIPES_ENDS, EXVOLENDCORR, i, true, EXVOLPIPES_NUMS[i], param_dir, res_dir, ID, this->assy_res_dir, calling_object_str);
	}
	for (i = 0; i < NINVOLUMES; ++i) {
		calling_object_str = "InVolume["; calling_object_str += to_string(i); calling_object_str += "]";
		InVolume[i].Initialise(pPpt, PipePtrs(INVOLPIPES_NUMS[i], rAssy, INVOLASSYS, INVOLPIPES, i, false),
			Intake, INVOLPIPES, INVOLPIPES_ENDS, INVOLENDCORR, i, false, INVOLPIPES_NUMS[i], param_dir, res_dir, ID, this->assy_res_dir, calling_object_str);
	}

	// Initialise pipes (must initialise boundary conditions first to get pipe end corrections)
	// ----------------------------------------------------------------------------------------------------
	exVolume = 0;
	for(p=0; p<NEXPIPES; ++p) {
		calling_object_str = "Exhaust["; calling_object_str += to_string(p); calling_object_str += "]";
		Exhaust[p].Initialise(pPpt, ID, p, true, &Eng[0], param_dir, this->assy_res_dir, calling_object_str);
		exVolume += Exhaust[p].vol;
	}
	inVolume = 0;
	for(p=0; p<NINPIPES; ++p) {
		calling_object_str = "Intake["; calling_object_str += to_string(p); calling_object_str += "]";
		Intake[p].Initialise(pPpt, ID, p, false, &Eng[0], param_dir, this->assy_res_dir, calling_object_str);
		inVolume += Intake[p].vol;
	}

	// POST pipe initialization
	// ----------------------------------------------------------------------------------------------------
///*
	// Adiabatic pressure loss devices
	// ----------------------------------------------------------------------------------------------------
	for (i = 0; i < NEXAPLDev; ++i) {
		calling_object_str = "APLDevE["; calling_object_str += to_string(i); calling_object_str += "]";
		APLDevE[i].InitialisePost(pMyTime, pPpt, PipePtrs(EXAPLDevPIPES_NUMS[i], rAssy, EXAPLDevASSYS, EXAPLDevPIPES, i, true),
			Exhaust, EXAPLDevPIPES, EXAPLDevPIPES_ENDS, EXAPLDevENDCORR, i, true, 2, &TransmE[0], param_dir, ID,
			NEXTRANSM > 0 ? true : false, this->assy_res_dir, EX_END, TransmE, JE, calling_object_str);
	}
	
	for (i = 0; i < NINAPLDev; ++i) {
		calling_object_str = "APLDevI["; calling_object_str += to_string(i); calling_object_str += "]";
		APLDevI[i].InitialisePost(pMyTime, pPpt, PipePtrs(INAPLDevPIPES_NUMS[i], rAssy, INAPLDevASSYS, INAPLDevPIPES, i, false),
			Intake, INAPLDevPIPES, INAPLDevPIPES_ENDS, INAPLDevENDCORR, i, false, 2, &TransmI[0], param_dir, ID,
			NINTRANSM > 0 ? true : false, this->assy_res_dir, IN_END, TransmI, JI, calling_object_str);
	}
//*/
}

CPipe** CAssembly::PipePtrs(int NPIPES, CAssembly* &rAssy, int** ASSYS, int** PIPES, int i, bool EX)
{
	CPipe** pPipes;
	pPipes = new CPipe* [NPIPES];	// Dimension as many pointers as there are connections
	for(int p = 0; p<NPIPES; ++p) pPipes[p] = (EX ? &(rAssy[ASSYS[i][p]].Exhaust[PIPES[i][p]]) : &(rAssy[ASSYS[i][p]].Intake[PIPES[i][p]]));
	return pPipes;					// Return actual array (so it can be deleted)
}

CPipe** CAssembly::PipePtrsCyl(int NPIPES, CAssembly* &rAssy, int** ASSYS, int** PIPES, int i, bool* EX)
{
	CPipe** pPipes;
	pPipes = new CPipe*[NPIPES];	// Dimension as many pointers as there are connections
	for (int p = 0; p < NPIPES; ++p) pPipes[p] = (EX[p] ? &(rAssy[ASSYS[i][p]].Exhaust[PIPES[i][p]]) : &(rAssy[ASSYS[i][p]].Intake[PIPES[i][p]]));
	return pPipes;
}

void CAssembly::ConfigureBCs(CProperties* pPpt)
// ============================================================ //
// Configures boundary nodes (attaches boundaries to pipe		//
// boundary nodes)												//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ConfigureBCs\n");}
	
	int i;

	// Endcaps
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXENDCAP; ++i) CE[i].Configure(pPpt);
	for(i=0; i<NINENDCAP; ++i) CI[i].Configure(pPpt);
	
	// Anechoic ends
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXANEC; ++i){EX_ANEC[i].Configure(pPpt);EX_ANEC[i].InitialiseBufferDamper(pPpt);}
	for(i=0; i<NINANEC; ++i){IN_ANEC[i].Configure(pPpt);IN_ANEC[i].InitialiseBufferDamper(pPpt);}

	// End environments
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXEND; ++i) EX_END[i].Configure(pPpt);
	for(i=0; i<NINEND; ++i) IN_END[i].Configure(pPpt);

	// End environments
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXEND; ++i)
	{
		EX_END[i].InitialiseAnechoic(pPpt);
	}
	for(i=0; i<NINEND; ++i)
	{
		IN_END[i].InitialiseAnechoic(pPpt);
	}

	// Transmissive boundaries
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXTRANSM; ++i) TransmE[i].Configure(pPpt);
	for(i=0; i<NINTRANSM; ++i) TransmI[i].Configure(pPpt);

	// Sudden area changes
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXSUD; ++i) SE[i].Configure(pPpt);
	for(i=0; i<NINSUD; ++i) SI[i].Configure(pPpt);

	// Junctions
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXJUNCS; ++i)
	{
		JE[i].Configure(pPpt);
		JE[i].ConfigureBranchAreas(pPpt);
	}
	for(i=0; i<NINJUNCS; ++i)
	{
		JI[i].Configure(pPpt);
		JI[i].ConfigureBranchAreas(pPpt);
	}

	// Cylinders
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NCYLS; ++i) Cyl[i].Configure(pPpt, &Exhaust, &Intake, NEXPIPES, NINPIPES);
	
	// Adiabatic pressure loss devices
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXAPLDev; ++i)
  {
    APLDevE[i].Configure(pPpt);
    APLDevE[i].InitialisePostConfig(pPpt);
  }

  for(i=0; i<NINAPLDev; ++i)
  {
    APLDevI[i].Configure(pPpt);
    APLDevI[i].InitialisePostConfig(pPpt);
  }

	// Turbines
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NTURBINE; ++i)
	{
		Turbine[i].Configure(pPpt);
		Turbine[i].InitialiseMap(pPpt); // Do this after Configure() as this requires pBN->f
	}

	// Volumes
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXVOLUMES; ++i)
	{
//		cout << "i = " << i << endl;
		ExVolume[i].Configure(pPpt, &Exhaust, &Intake, NEXPIPES, NINPIPES);
	}
	for(i=0; i<NINVOLUMES; ++i) InVolume[i].Configure(pPpt, &Exhaust, &Intake, NEXPIPES, NINPIPES);
}

void CAssembly::ConfigureAll(CProperties* pPpt, vector<string> &rConfStrs)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ConfigureAll\n");}
	if(pPpt->SHOW_config) {
		std::string temp;
		temp = "Assembly ["; temp += ID + 48; temp += "]: reading configuration files for boundary conditions";
		pPpt->Out(Underline(StringToChar(temp), "=", true, 80, 1));
	}
	
	// Cylinders
	// ----------------------------------------------------------------------------------------------------
	if(NCYLS>0) pPpt->Configure(ConstructString(pPpt, param_dir, "CYLINDER.txt"), CYLCONFP, CYLENDCORR, CYLASSYS, CYLPIPES, CYLPIPES_ENDS, CYLPIPES_NUMS, NCYLS, pPpt->STRLEN, "CYLINDER", 4, rConfStrs, ID);  // Maximum 4; 2 exhaust valves, 2 intake valves
	
	// Endcaps
	// ----------------------------------------------------------------------------------------------------
	if(NEXENDCAP>0) pPpt->Configure(ConstructString(pPpt, param_dir, "EX_CLOSED.txt"), EXEndCapCONFP, EXEndCapENDCORR, EXEndCapASSYS, EXEndCapPIPES, EXEndCapPIPES_ENDS, EXEndCapPIPES_NUMS, NEXENDCAP, pPpt->STRLEN, "EXHAUST CLOSED END", 1, rConfStrs, ID);
	if(NINENDCAP>0) pPpt->Configure(ConstructString(pPpt, param_dir, "IN_CLOSED.txt"), INEndCapCONFP, INEndCapENDCORR, INEndCapASSYS, INEndCapPIPES, INEndCapPIPES_ENDS, INEndCapPIPES_NUMS, NINENDCAP, pPpt->STRLEN, "INTAKE CLOSED END", 1, rConfStrs, ID);
	
	// Anechoic ends
	// ----------------------------------------------------------------------------------------------------
	if(NEXANEC>0) pPpt->Configure(ConstructString(pPpt, param_dir, "EX_ANEC.txt"), EXANECCONFP, EXANECENDCORR, EXANECASSYS, EXANECPIPES, EXANECPIPES_ENDS, EXANECPIPES_NUMS, NEXANEC, pPpt->STRLEN, "EXHAUST ANECHOIC END", 1, rConfStrs, ID);
	if(NINANEC>0) pPpt->Configure(ConstructString(pPpt, param_dir, "IN_ANEC.txt"), INANECCONFP, INANECENDCORR, INANECASSYS, INANECPIPES, INANECPIPES_ENDS, INANECPIPES_NUMS, NINANEC, pPpt->STRLEN, "INTAKE ANECHOIC END", 1, rConfStrs, ID);

	// End environments
	// ----------------------------------------------------------------------------------------------------
	if(NEXEND>0) pPpt->Configure(ConstructString(pPpt, param_dir, "EX_END.txt"), EXECONFP, EXEENDCORR, EXEASSYS, EXEPIPES, EXEPIPES_ENDS, EXEPIPES_NUMS, NEXEND, pPpt->STRLEN, "EXHAUST END ENVIRONMENT", 1, rConfStrs, ID);
	if(NINEND>0) pPpt->Configure(ConstructString(pPpt, param_dir, "IN_END.txt"), INECONFP, INEENDCORR, INEASSYS, INEPIPES, INEPIPES_ENDS, INEPIPES_NUMS, NINEND, pPpt->STRLEN, "INTAKE END ENVIRONMENT", 1, rConfStrs, ID);

	// Junctions
	// ----------------------------------------------------------------------------------------------------
	if(NEXJUNCS>0) pPpt->Configure(ConstructString(pPpt, param_dir, "EX_JUNCTION.txt"), EXJCONFP, EXJENDCORR, EXJASSYS, EXJPIPES, EXJPIPES_ENDS, EXJPIPES_NUMS, NEXJUNCS, pPpt->STRLEN, "EXHAUST JUNCTION", pPpt->max_branches, rConfStrs, ID);
	if(NINJUNCS>0) pPpt->Configure(ConstructString(pPpt, param_dir, "IN_JUNCTION.txt"), INJCONFP, INJENDCORR, INJASSYS, INJPIPES, INJPIPES_ENDS, INJPIPES_NUMS, NINJUNCS, pPpt->STRLEN, "INTAKE JUNCTION", pPpt->max_branches, rConfStrs, ID);

	// Adiabatic pressure loss devices
	// ----------------------------------------------------------------------------------------------------
	if(NEXAPLDev>0) pPpt->Configure(ConstructString(pPpt, param_dir, "EX_APLDev.txt"), EXAPLDevCONFP, EXAPLDevENDCORR, EXAPLDevASSYS, EXAPLDevPIPES, EXAPLDevPIPES_ENDS, EXAPLDevPIPES_NUMS, NEXAPLDev, pPpt->STRLEN, "EXHAUST ADIABATIC PRESSURE LOSS DEVICE", 2, rConfStrs, ID);
	if(NINAPLDev>0) pPpt->Configure(ConstructString(pPpt, param_dir, "IN_APLDev.txt"), INAPLDevCONFP, INAPLDevENDCORR, INAPLDevASSYS, INAPLDevPIPES, INAPLDevPIPES_ENDS, INAPLDevPIPES_NUMS, NINAPLDev, pPpt->STRLEN, "INTAKE ADIABATIC PRESSURE LOSS DEVICE", 2, rConfStrs, ID);

	// Sudden area changes
	// ----------------------------------------------------------------------------------------------------
	if(NEXSUD>0) pPpt->Configure(ConstructString(pPpt, param_dir, "EX_SUDDEN.txt"), EXSCONFP, EXSENDCORR, EXSASSYS, EXSPIPES, EXSPIPES_ENDS, EXSPIPES_NUMS, NEXSUD, pPpt->STRLEN, "EXHAUST SUDDEN AREA CHANGE", 2, rConfStrs, ID);
	if(NINSUD>0) pPpt->Configure(ConstructString(pPpt, param_dir, "IN_SUDDEN.txt"), INSCONFP, INSENDCORR, INSASSYS, INSPIPES, INSPIPES_ENDS, INSPIPES_NUMS, NINSUD, pPpt->STRLEN, "INTAKE SUDDEN AREA CHANGE", 2, rConfStrs, ID);

	// Transmissive boundaries
	// ----------------------------------------------------------------------------------------------------
	if(NEXTRANSM>0) pPpt->Configure(ConstructString(pPpt, param_dir, "EX_TRANSM.txt"), EXTCONFP, EXTENDCORR, EXTASSYS, EXTPIPES, EXTPIPES_ENDS, EXTPIPES_NUMS, NEXTRANSM, pPpt->STRLEN, "EXHAUST TRANSMISSIVE BOUNDARY", 1, rConfStrs, ID); 
	if(NINTRANSM>0) pPpt->Configure(ConstructString(pPpt, param_dir, "IN_TRANSM.txt"), INTCONFP, INTENDCORR, INTASSYS, INTPIPES, INTPIPES_ENDS, INTPIPES_NUMS, NINTRANSM, pPpt->STRLEN, "INTAKE TRANSMISSIVE BOUNDARY", 1, rConfStrs, ID); 

	// Turbines
	// ----------------------------------------------------------------------------------------------------
	if(NTURBINE>0) pPpt->Configure(ConstructString(pPpt, param_dir, "TURBINE.txt"), TURBCONFP, TURBENDCORR, TURBASSYS, TURBPIPES, TURBPIPES_ENDS, TURBPIPES_NUMS, NTURBINE, pPpt->STRLEN, "TURBINE", 2, rConfStrs, ID);

	// Volumes
	// ----------------------------------------------------------------------------------------------------
	if(NEXVOLUMES>0) pPpt->Configure(ConstructString(pPpt, param_dir, "EX_VOLUME.txt"), EXVOLCONFP, EXVOLENDCORR, EXVOLASSYS, EXVOLPIPES, EXVOLPIPES_ENDS, EXVOLPIPES_NUMS, NEXVOLUMES, pPpt->STRLEN, "EXHAUST VOLUME", 4, rConfStrs, ID);  // Maximum 4 valves/pipes join the volume
	if(NINVOLUMES>0) pPpt->Configure(ConstructString(pPpt, param_dir, "IN_VOLUME.txt"), INVOLCONFP, INVOLENDCORR, INVOLASSYS, INVOLPIPES, INVOLPIPES_ENDS, INVOLPIPES_NUMS, NINVOLUMES, pPpt->STRLEN, "INTAKE VOLUME", 4, rConfStrs, ID);  // Maximum 4 valves/pipes join the volume
	
/*
	// Recursive call
	for(int i=0; i<NASSEMBLY; ++i) Assembly[i].ConfigureAll(pPpt, rConfStrs);
*/
}

void CAssembly::ListPropertiesAll(CProperties* pPpt)
// ============================================================ //
// Prints individual object properties to screen				//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify());pPpt->Out(".ListPropertiesAll\n");}
	if(pPpt->SHOW_params) 
	{
		pPpt->Out("\n");
		ListProperties(pPpt);
		pPpt->Out("\n");
		//Line(80, 22);
		//pPpt->Out("\n");
		std::string temp;
		temp = "Assembly ["; temp += ID + 48; temp += "] object parameters";
//pPpt->Out(temp);
//char pause;
//cin >>pause;
		pPpt->Out(Underline(StringToChar(temp), "="));
		pPpt->Out("\n");
		
		int i, j;

		// Pipe properties
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<NEXPIPES; ++i) Exhaust[i].ListProperties(pPpt);
		for(i=0; i<NINPIPES; ++i) Intake[i].ListProperties(pPpt);

		// Volume properties
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<NEXVOLUMES; ++i) ExVolume[i].ListProperties(pPpt);
		for(i=0; i<NINVOLUMES; ++i) InVolume[i].ListProperties(pPpt);
		
		// Engine properties
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<NENGINES; ++i) Eng[i].ListProperties(pPpt);
		
		// Cylinder properties
		// ----------------------------------------------------------------------------------------------------
		for (i = 0; i < NCYLS; ++i) {
			Cyl[i].ListProperties(pPpt);
			for (j = 0; j < Cyl[i].EnginePtr->NEXVALVES; ++j) {
				Cyl[i].ExhaustValve[j].ListProperties(pPpt, "\t");
			}
			for (j = 0; j < Cyl[i].EnginePtr->NINVALVES; ++j) {
				Cyl[i].IntakeValve[j].ListProperties(pPpt, "\t");
			}
		}
		
		// Anechoic end properties
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<NEXANEC; ++i) EX_ANEC[i].ListProperties(pPpt);
		for(i=0; i<NINANEC; ++i) IN_ANEC[i].ListProperties(pPpt);
		
		// End environent properties
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<NEXEND; ++i) EX_END[i].ListProperties(pPpt);
		for(i=0; i<NINEND; ++i) IN_END[i].ListProperties(pPpt);
		
		// Junction properties
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<NEXJUNCS; ++i) JE[i].ListProperties(pPpt);
		for(i=0; i<NINJUNCS; ++i) JI[i].ListProperties(pPpt);
		
		// Adiabatic pressure loss device properties
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<NEXAPLDev; ++i) APLDevE[i].ListProperties(pPpt);
		for(i=0; i<NINAPLDev; ++i) APLDevI[i].ListProperties(pPpt);
		
		// Turbine properties
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<NTURBINE; ++i) Turbine[i].ListProperties(pPpt);
		
		// Transmissive boundary properties
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<NEXTRANSM; ++i) TransmE[i].ListProperties(pPpt);
		for(i=0; i<NINTRANSM; ++i) TransmI[i].ListProperties(pPpt);
		//pPpt->Out("\n");
	}
}

void CAssembly::ReadAssyDir(CProperties* pPpt, char *InputFile)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ReadAssyDir\n");}
	int last_entry;
	CommonReadInput(InputFile, pPpt->NUM_PARAMS, labels, values, strings, last_entry);
	for(int r=0; r<last_entry+1; ++r)
	{
		if(strcmp(labels[r], "assy_dir") == 0)
		{
			assy_dir = strings[r];
		}
	}

	std::string s;
	s = assy_dir;
	
	std::string sz;
	for(int i=0; i<int(s.length()); ++i) sz += s[i]; 
	sz += "\\"; // Add directory slash

	assy_dir_slash = new char[sz.length() + 1];
	strcpy(assy_dir_slash, sz.c_str());

	param_dir = ConstructString(pPpt, pPpt->assemblies_dir, assy_dir_slash);
	//res_dir = ConstructString(pPpt, pPpt->case_dir, assy_dir_slash);
/*			
	//cout << assy_dir << endl;
	//cout << assy_dir_slash << endl;
	//cout << param_dir << endl;
	//cout << res_dir << endl;
	//exit(1);
//*/
}

void CAssembly::ReadInput(CProperties* pPpt, char *InputFile)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ReadInput\n");}

	int last_entry;
	CommonReadInput(InputFile, pPpt->NUM_PARAMS, labels, values, strings, last_entry);
	for(int r=0; r<last_entry+1; ++r)
	{
		// Pipe parameters
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NEXPIPES") == 0) NEXPIPES = int(values[r]);
		if(strcmp(labels[r], "NINPIPES") == 0) NINPIPES = int(values[r]);

		// Engines, cylinders & valves
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "ENG_FILE") == 0) ENG_FILE = strings[r];
		if(strcmp(labels[r], "NENGINES") == 0) NENGINES = int(values[r]);
		if(strcmp(labels[r], "NCYLS") == 0) NCYLS = int(values[r]);
//		if(strcmp(labels[r], "CYL_DIR") == 0) CYL_DIR = strings[r];
		//if(strcmp(labels[r], "POPPET_H_TOL1") == 0) POPPET_H_TOL1 = values[r];
		//if(strcmp(labels[r], "POPPET_H_TOL2") == 0) POPPET_H_TOL2 = values[r];
//		if(strcmp(labels[r], "VT_DIR") == 0) VT_DIR = strings[r];

		// End cap parameters
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NEXENDCAP") == 0) NEXENDCAP = int(values[r]);
		if(strcmp(labels[r], "NINENDCAP") == 0) NINENDCAP = int(values[r]);

		// Anechoic end parameters
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NEXANEC") == 0) NEXANEC = int(values[r]);
		if(strcmp(labels[r], "NINANEC") == 0) NINANEC = int(values[r]);

		// End environment parameters
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NEXEND") == 0) NEXEND = int(values[r]);
		if(strcmp(labels[r], "NINEND") == 0) NINEND = int(values[r]);

		// Junction parameters
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NEXJUNCS") == 0) NEXJUNCS = int(values[r]);
		if(strcmp(labels[r], "NINJUNCS") == 0) NINJUNCS = int(values[r]);

		// Adiabatic pressure loss device parameters
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NEXAPLDev") == 0) NEXAPLDev = int(values[r]);
		if(strcmp(labels[r], "NINAPLDev") == 0) NINAPLDev = int(values[r]);

		// Sudden area change parameters
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NEXSUD") == 0) NEXSUD = int(values[r]);
		if(strcmp(labels[r], "NINSUD") == 0) NINSUD = int(values[r]);

		// Transmissive boundary parameters
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NEXTRANSM") == 0) NEXTRANSM = int(values[r]);
		if(strcmp(labels[r], "NINTRANSM") == 0) NINTRANSM = int(values[r]);
		
		// Turbine parameters
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NTURBINE") == 0) NTURBINE = int(values[r]);

		// Volume parameters
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NEXVOLUMES") == 0) NEXVOLUMES = int(values[r]);
		if(strcmp(labels[r], "NINVOLUMES") == 0) NINVOLUMES = int(values[r]);
	}

	
	// If simulation type is set to periodic, but there are no cylinders, reset it to continuous
	if(NCYLS==0)
	{
		pPpt->CONTINUOUS = true;
	}
}


void CAssembly::ListProperties(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ListProperties\n");}
	string temp;
	temp = "Assembly [";
	temp += IntToString(ID);
	temp += "] parameters";
	pPpt->Out(Underline(StringToChar(temp), "="));
	pPpt->Out("\n");
	pPpt->Out("\tAssembly data file directory, assy_dir\t\t=\t"); pPpt->Out(assy_dir); pPpt->Out("\n");
	pPpt->Out("\n");

	// Pipe parameters
	// ----------------------------------------------------------------------------------------------------
	if(NEXPIPES>0 || NINPIPES>0)
	{
		pPpt->Out(Underline("Pipe parameters", "-", "\t"));
		pPpt->Out("\tNo. of pipes, exhaust system, NEXPIPES\t\t=\t"); pPpt->Out(NEXPIPES); pPpt->Out("\n");
		if(exVolume>0){pPpt->Out("\tTotal volume of exhaust system, exVolume\t=\t"); pPpt->Out(exVolume*1000); pPpt->Out(" litres\n");}
		pPpt->Out("\tNo. of pipes, intake system, NINPIPES\t\t=\t"); pPpt->Out(NINPIPES); pPpt->Out("\n");
		if(inVolume>0){pPpt->Out("\tTotal volume of intake system, inVolume\t=\t"); pPpt->Out(inVolume*1000); pPpt->Out(" litres\n");}
		pPpt->Out("\n");
	}

	// Engines, cylinders & valves
	// ----------------------------------------------------------------------------------------------------
	if(NENGINES>0)
	{
		pPpt->Out(Underline("Engine parameters", "-", "\t"));
		pPpt->Out("\tNo. of engines, NENGINES\t\t\t=\t"); pPpt->Out(NENGINES); pPpt->Out("\n");
		pPpt->Out("\tEngine parameter file, ENG_FILE\t\t\t=\t"); pPpt->Out(ENG_FILE); pPpt->Out("\n");
		pPpt->Out("\tNo. of cylinders, NCYLS\t\t\t\t=\t"); pPpt->Out(NCYLS); pPpt->Out("\n");
		/*	
		if(!pPpt->HOMENTROPIC)
		{
		cout << endl;
		cout << Underline("Valve parameters", "-", "\t\t");
		//		cout << "\t\tValve timing files directory (VT_DIR)\t=\t" << VT_DIR << "\n";
		}
		else
		{
		cout << endl;
		cout << Underline("Valve parameters", "-", "\tzt");
		//cout << "\tHom. tol. 1 (POPPET_H_TOL1)\t\t\t=\t" << POPPET_H_TOL1 << " \n";
		//cout << "\tHom. tol. 2 (POPPET_H_TOL2)\t\t\t=\t" << POPPET_H_TOL2 << " \n";
		//		cout << "\t\tValve timing files directory (VT_DIR)\t\t\t=\t" << VT_DIR << "\n";
		}
		*/
		cout << endl;
		pPpt->Out("\n");
	}

	// End cap parameters
	// ----------------------------------------------------------------------------------------------------
	if(NEXENDCAP>0 || NINENDCAP>0)
	{
		pPpt->Out(Underline("End cap parameters", "-", "\t"));
		if(NEXENDCAP>0)
		{
			pPpt->Out("\tNo. closed ends, exhaust sys., NEXENDCAP\t=\t"); pPpt->Out(NEXENDCAP); pPpt->Out("\n");
		}
		if(NINENDCAP>0)
		{
			pPpt->Out("\tNo. closed ends, intake sys., NINENDCAP\t=\t"); pPpt->Out(NINENDCAP); pPpt->Out("\n");
		}
		pPpt->Out("\n");
	}

	// Anechoic end parameters
	// ----------------------------------------------------------------------------------------------------
	if(NEXANEC>0 || NINANEC>0)
	{
		pPpt->Out(Underline("Anechoic end parameters", "-", "\t"));
		if(NEXANEC>0)
		{
			pPpt->Out("\tNo. anechoic ends, exhaust sys., NEXANEC\t=\t"); pPpt->Out(NEXANEC); pPpt->Out("\n");
		}
		if(NINANEC>0) 
		{
			pPpt->Out("\tNo. anechoic ends, intake sys., NINANEC\t=\t"); pPpt->Out(NINANEC); pPpt->Out("\n");
		}
		pPpt->Out("\n");
	}

	// End environment parameters
	// ----------------------------------------------------------------------------------------------------
	if(NEXEND>0 || NINEND>0)
	{
		pPpt->Out(Underline("End environment parameters", "-", "\t"));
		if(NEXEND>0)
		{
			pPpt->Out("\tNo. end environments, exhaust sys., NEXEND\t=\t"); pPpt->Out(NEXEND); pPpt->Out("\n");
		}
		if(NINEND>0) 
		{
			pPpt->Out("\tNo. end environments, intake sys., NINEND\t=\t"); pPpt->Out(NINEND); pPpt->Out("\n");
		}
		pPpt->Out("\n");
	}

	// Junction parameters
	// ----------------------------------------------------------------------------------------------------
	if(NEXJUNCS>0 || NINJUNCS>0)
	{
		pPpt->Out(Underline("Junction parameters", "-", "\t"));
		if(NEXJUNCS>0)
		{
			pPpt->Out("\tNo. junctions, exhaust sys., NEXJUNCS\t\t=\t"); pPpt->Out(NEXJUNCS); pPpt->Out("\n");
		}
		if(NINJUNCS>0)
		{
			pPpt->Out("\tNo. junctions, intake sys., NINJUNCS\t\t=\t"); pPpt->Out(NINJUNCS); pPpt->Out("\n");
		}
		pPpt->Out("\n");
	}

	// Adiabatic pressure loss device parameters
	// ----------------------------------------------------------------------------------------------------
	if(NEXAPLDev>0 || NINAPLDev>0)
	{
		pPpt->Out(Underline("Adiabatic pressure loss device parameters", "-", "\t"));
		if(NEXAPLDev>0)
		{
			pPpt->Out("\tNo. a.p.l. devs, exhaust sys., NEXAPLDev\t=\t"); pPpt->Out(NEXAPLDev); pPpt->Out("\n");
		}
		if(NINAPLDev>0)
		{
			pPpt->Out("\tNo. a.p.l. devs, intake sys., NINAPLDev\t=\t"); pPpt->Out(NINAPLDev); pPpt->Out("\n");
		}
		pPpt->Out("\n");
	}

	// Sudden area change parameters
	// ----------------------------------------------------------------------------------------------------
	if(NEXSUD>0 || NINSUD>0)
	{
		pPpt->Out(Underline("Sudden area change parameters", "-", "\t"));
		if(NEXSUD>0)
		{
			pPpt->Out("\tNo. sudden area changes, exhaust sys., NEXSUD\t=\t"); pPpt->Out(NEXSUD); pPpt->Out("\n");
		}
		if(NINSUD>0)
		{
			pPpt->Out("\tNo. sudden area changes, intake sys., NINSUD\t=\t"); pPpt->Out(NINSUD); pPpt->Out("\n");
		}
		pPpt->Out("\n");
	}

	// Transmissive boundary parameters
	// ----------------------------------------------------------------------------------------------------
	if(NEXTRANSM>0 || NINTRANSM>0)
	{
		pPpt->Out(Underline("Transmissive boundary parameters", "-", "\t"));
		if(NEXTRANSM>0)
		{
			pPpt->Out("\tNo. transmissive bcs, exhaust sys., NEXTRANSM\t=\t"); pPpt->Out(NEXTRANSM); pPpt->Out("\n");
		}
		if(NINTRANSM>0)
		{
			pPpt->Out("\tNo. transmissive bcs, intake sys., NINTRANSM\t=\t"); pPpt->Out(NINTRANSM); pPpt->Out("\n");
		}
		pPpt->Out("\n");
	}

	// Turbine parameters
	// ----------------------------------------------------------------------------------------------------
	if(NTURBINE>0)
	{
		pPpt->Out(Underline("Turbine parameters", "-", "\t"));
		pPpt->Out("\tNo. of turbines, NTURBINE\t\t\t=\t"); pPpt->Out(NTURBINE); pPpt->Out("\n");
		pPpt->Out("\n");
	}

	// Volume boundary parameters
	// ----------------------------------------------------------------------------------------------------
	if(NEXVOLUMES>0 || NINVOLUMES>0)
	{
		pPpt->Out(Underline("Volume boundary parameters", "-", "\t"));
		if(NEXVOLUMES>0) 
		{
			pPpt->Out("\tNo. volumes, exhaust sys., NEXVOLUMES\t\t=\t"); pPpt->Out(NEXVOLUMES); pPpt->Out("\n");
		}
		if(NINVOLUMES>0)
		{
			pPpt->Out("\tNo. volumes, intake sys., NINVOLUMES\t\t=\t"); pPpt->Out(NINVOLUMES); pPpt->Out("\n");
		}
		pPpt->Out("\n");
	}
/*
	Line(80, 22);
	cout << endl;
	cout << Underline("Simulation object parameters", "=");
	cout << endl;

	// Pipe properties
	for(int i=0; i<NEXPIPES; ++i) Exhaust[i].ListProperties();
	for(i=0; i<NINPIPES; ++i) Intake[i].ListProperties();

	// Engine properties
	for(i=0; i<NENGINES; ++i) Eng[i].ListProperties(pPpt);

	// Cylinder properties
	for(i=0; i<NCYLS; ++i) Cyl[i].ListProperties();

	// End environent properties
	for(i=0; i<NEXEND; ++i) EE[i].ListProperties();
	for(i=0; i<NINEND; ++i) EI[i].ListProperties();

	// Junction properties
	for(i=0; i<NEXJUNCS; ++i) JE[i].ListProperties(pPpt);
	for(i=0; i<NINJUNCS; ++i) JI[i].ListProperties(pPpt);
	
	// Adiabatic pressure loss device properties
	for(i=0; i<NEXAPLDev; ++i) APLDevE[i].ListProperties();
	for(i=0; i<NINAPLDev; ++i) APLDevI[i].ListProperties();

	// Turbine properties
	for(i=0; i<NTURBINE; ++i) Turbine[i].ListProperties();
	
	// Transmissive boundary properties
	for(i=0; i<NEXTRANSM; ++i) TransmE[i].ListProperties();
	for(i=0; i<NINTRANSM; ++i) TransmI[i].ListProperties();
*/
}

void CAssembly::RunBoundary(CProperties* pPpt, CTime MyTime, double DELZe, double DELZi, double* TIMEe, double* TIMEi, int timestep, bool RESTORE)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".RunBoundary\n");}

	int i, p;
	int R=0;
		
	// Engines, Cylinders
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NENGINES; ++i)
	{
		Eng[i].RunBoundary(pPpt, MyTime, DELZe, DELZi, TIMEe[R+1]);
		if(Eng[i].NEW_CYCLE)
		{
			// For all the turbines linked to this engine
			for(int t=0; t<NTURBINE; ++t) Turbine[t].CycleWorkAndMassFlow(pPpt, Eng[i].reveng/2);
		}
	}
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NCYLS; ++i) Cyl[i].RunBoundary(pPpt, MyTime, DELZe, DELZi, Eng[0].ca, fabs(Eng[0].ca - Eng[0].ca_old), timestep, TIMEe[R+1]);
	
	// Endcaps
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXENDCAP; ++i) CE[i].RunBoundary();
	for(i=0; i<NINENDCAP; ++i) CI[i].RunBoundary();
	
	// Anechoic end
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXANEC; ++i) EX_ANEC[i].RunBoundary(pPpt, timestep, TIMEe[R+1]);
	for(i=0; i<NINANEC; ++i) IN_ANEC[i].RunBoundary(pPpt, timestep, TIMEe[R+1]);

	// End environments
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXEND; ++i) EX_END[i].RunBoundary(pPpt, timestep, TIMEe[R+1]);
	for(i=0; i<NINEND; ++i) IN_END[i].RunBoundary(pPpt, timestep, TIMEe[R+1]);

	// Junctions
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXJUNCS; ++i) JE[i].RunBoundary(pPpt, TIMEe[R+1], RESTORE, timestep);
	for(i=0; i<NINJUNCS; ++i) JI[i].RunBoundary(pPpt, TIMEe[R+1], RESTORE, timestep);
	
	// Sudden area changes
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXSUD; ++i) SE[i].RunBoundary(pPpt, TIMEe[R+1], timestep);
	for(i=0; i<NINSUD; ++i) SI[i].RunBoundary(pPpt, TIMEe[R+1], timestep);

	// Turbines
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NTURBINE; ++i) Turbine[i].RunBoundary(pPpt, timestep, TIMEe[R+1], DELZe);

	// Transmissive boundaries
	// ----------------------------------------------------------------------------------------------------
	bool EX_STEADY, IN_STEADY;
	EX_STEADY = true;
	for(p=0; p<NEXPIPES; ++p)
	{
		// Exclude pipes which are joiner pipes (method=4) from steady test
		//if(Exhaust[p].METHOD!=4 && !Exhaust[p].Steady(pPpt)) EX_STEADY = false;

		if(!Exhaust[p].Steady(pPpt)) EX_STEADY = false;
	}
	IN_STEADY = true;
	for(p=0; p<NINPIPES; ++p)
	{
		if(!Intake[p].Steady(pPpt)) IN_STEADY = false;
	}
	for(i=0; i<NEXTRANSM; ++i) TransmE[i].RunBoundary(pPpt, MyTime, DELZe, DELZi, Eng[0].ca, fabs(Eng[0].ca - Eng[0].ca_old), TIMEe[R+1], TIMEi[R+1], timestep, EX_STEADY);
	for(i=0; i<NINTRANSM; ++i) TransmI[i].RunBoundary(pPpt, MyTime, DELZe, DELZi, Eng[0].ca, fabs(Eng[0].ca - Eng[0].ca_old), TIMEe[R+1], TIMEi[R+1], timestep, IN_STEADY);

	double vel_grad = TransmE[0].vel_grad;
	double mfr_under_est = TransmE[0].mfr_under_est;

//if(timestep>=3369) cout << "bottom Timestep = " << timestep << endl;	
	// Adiabatic pressure loss devices
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXAPLDev; ++i) APLDevE[i].RunBoundary(pPpt, TIMEe[R+1], vel_grad, mfr_under_est, MyTime, DELZe, DELZi, Eng[0].ca, fabs(Eng[0].ca - Eng[0].ca_old), TIMEe[R+1], TIMEi[R+1], timestep, EX_STEADY);
	for(i=0; i<NINAPLDev; ++i) APLDevI[i].RunBoundary(pPpt, TIMEe[R+1], vel_grad, mfr_under_est, MyTime, DELZe, DELZi, Eng[0].ca, fabs(Eng[0].ca - Eng[0].ca_old), TIMEe[R+1], TIMEi[R+1], timestep, IN_STEADY);

	// Volumes
	// ----------------------------------------------------------------------------------------------------
	for(i=0; i<NEXVOLUMES; ++i) ExVolume[i].RunBoundary(pPpt, DELZe, TIMEe[R+1], timestep);
	for(i=0; i<NINVOLUMES; ++i) InVolume[i].RunBoundary(pPpt, DELZi, TIMEe[R+1], timestep);
}

void CAssembly::PrintBoundaryConnections(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".PrintBoundaryConnections\n");}
	if(pPpt->SHOW_conn && ((NEXPIPES>0 || NINPIPES>0)))
	{
		int p;
		std::string temp;
		temp = "Assembly ["; temp += ID + 48; temp += "] boundary conditions";
		pPpt->Out(Underline(StringToChar(temp), "="));
		if(NEXPIPES>0 || NINPIPES>0)
		{
			for(p=0; p<NEXPIPES; ++p){cout << endl;	Exhaust[p].PrintBoundaries(pPpt);}
			for(p=0; p<NINPIPES; ++p){cout << endl;	Intake[p].PrintBoundaries(pPpt);}
			pPpt->Out("\n");
		}
	}
}

char* CAssembly::Identify()
// ============================================================ //
// Returns identification of the current object					//
// ============================================================ //
{
	std::string sz;
	sz = "Assembly [";
	sz += IntToString(ID);
	sz += "]";

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str());
	return szz;
}