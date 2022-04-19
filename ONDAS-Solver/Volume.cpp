// Volume.cpp: implementation of the CVolume class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Volume.h"
#include "Tools.h"
#include "Valve.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CVolume::CVolume()
{

}

CVolume::~CVolume()
{

}

//void CCylinder::Initialise(CProperties* pPpt, CPipe** pPipes,
//						   CPipe* &rExhaustPipe, CPipe* &rIntakePipe, int** &rCYLPIPES, int** &rCYLPIPES_ENDS, double* &rENDCORR,
//							int id, CEngine* EngPtr, std::string param_dir/*, char* cyl_dir*//*, char* vt_dir*/, int assyid)

void CVolume::Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rVOLPIPES, int** &rVOLPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, std::string param_dir, std::string res_dir, int assyid, string parent_assy_res_dir, string calling_object_str)
{
	if (pPpt->SHOW_calls) { pPpt->Out("CVolume.Initialise ("); pPpt->Out(calling_object_str); pPpt->Out(")\n"); }

	ID = id;
	AssyID = assyid;				// ID of assembly on which boundary belongs
	EX = ex;

	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Initialise\n");}

	NVOLUMEVALVES = npipes;

//	cout << "NVOLUMEVALVES = " << NVOLUMEVALVES << endl;
//	pPpt->Out("\n");

	std::string bcname_str = "VOLUME";
	ReadInput(pPpt, ConstructString(pPpt, param_dir, bcname_str, EX, ID));

	// Setup volume valve(s)
	// ---------------------
	int i;
	
	CPipe*** pPipesVolume;
	int** VVPIPES; int** VVPIPES_ENDS; double* VVENDCORR;
	
	pPipesVolume = new CPipe** [NVOLUMEVALVES];
	for(i=0; i<NVOLUMEVALVES; ++i)
		pPipesVolume[i] = new CPipe* [1]; // Only one pipe connection per valve


	VVPIPES = new int* [NVOLUMEVALVES]; VVPIPES_ENDS = new int* [NVOLUMEVALVES]; VVENDCORR = new double [NVOLUMEVALVES];	
	
	
	for(i=0; i<NVOLUMEVALVES; ++i)
	{VVPIPES[i] = new int [1]; VVPIPES_ENDS[i] = new int [1];}
	for(i=0; i<NVOLUMEVALVES; ++i)
	{
		pPipesVolume[i][ONE_SIDE] = pPipes[i];
		VVPIPES[i][ONE_SIDE] = rVOLPIPES[ID][i];
		VVPIPES_ENDS[i][ONE_SIDE] = rVOLPIPES_ENDS[ID][i];
		VVENDCORR[i] = rENDCORR[ID];
/*
		cout << "VolumeValve[" << i << "]:" << endl;
		cout << "pPipesVolume[" << i << "][ONE_SIDE] = " << pPipesVolume[i][ONE_SIDE] << endl;
		cout << "VVPIPES[" << i << "][ONE_SIDE] = " << VVPIPES[i][ONE_SIDE] << endl;
		cout << "VVPIPES_ENDS[" << i << "][ONE_SIDE] = " << VVPIPES_ENDS[i][ONE_SIDE] << endl;
		cout << "VVENDCORR[" << i << "] = " << VVENDCORR[i] << endl;
		pPpt->Out("\n");
*/
	}

	VolumeValve = new CValve [NVOLUMEVALVES];
	for(i=0; i<NVOLUMEVALVES; ++i)
		VolumeValve[i].InitialiseVolumeValve(pPpt, pPipesVolume[i], rPipe, VVPIPES, VVPIPES_ENDS, VVENDCORR, i, EX, 1, param_dir, AssyID, ID, parent_assy_res_dir);

	dmdt = new double [NVOLUMEVALVES];
	for(i=0; i<NVOLUMEVALVES; ++i) dmdt[i] = 0.0;
	
	av = new double [this->NVOLUMEVALVES];
	for(i=0; i<NVOLUMEVALVES; ++i) av[i] = (EX ? pPpt->AREFe : pPpt->AREFi);

	if(USE_V_INIT) V = V_INIT;		// V_INIT is the volume specified by the user
	else V = vol;					// vol is calculated from pipe geometry

	p0 = p0_INIT;
	T0 = T0_INIT;
	m = 1.0E5*p0*V/(287.0*T0);
	a0 = sqrt(pPpt->gammaAir(T0)*287.0*T0);
	dVdt = 0;
	dp0dt = -pPpt->gammaAir(T0)*dVdt*p0/V;

	pX_temp = new double [NVOLUMEVALVES];
	TX_temp = new double [NVOLUMEVALVES];
	regime = new double [NVOLUMEVALVES];
	FORWARD = new bool [NVOLUMEVALVES];
	for(i=0; i<NVOLUMEVALVES; ++i)
	{
		pX_temp[i] = p0;
		TX_temp[i] = T0;
		regime[i] = 0;

		FORWARD[i] = true;
	}

	SetupFiles(pPpt, parent_assy_res_dir);
}

void CVolume::Configure(CProperties* pPpt, CPipe** pExhaustSystem, CPipe** pIntakeSystem, int nExPipes, int nInPipes)
{
	int i;
	// Configure valves
//	cout << "NVOLUMEVALVES = " << NVOLUMEVALVES << endl;
	for(i=0; i<NVOLUMEVALVES; ++i)
		VolumeValve[i].Configure(pPpt);

	for(i=0; i<NVOLUMEVALVES; ++i)
		VolumeValve[i].ConfigureExtra(pPpt, pExhaustSystem, pIntakeSystem, nExPipes, nInPipes); // Sets FP, the area of the adjoining pipe

	for(i=0; i<NVOLUMEVALVES; ++i)
	{
		VolumeValve[i].Set_USE_PIPE_DIA(this->USE_PIPE_DIA);
		VolumeValve[i].Set_ref_dia(this->ref_dia); // Sets valve ref_dia
		VolumeValve[i].EffectiveArea(pPpt); // Sets an up-to-date area
		VolumeValve[i].Set_psi_val(VolumeValve[i].Get_eff_area()/VolumeValve[i].Get_FP());
		// Sets an up-to-date psi_val
	}
}

void CVolume::RunBoundary(CProperties* pPpt, double DELZ, double time, int timestep)
{
	int i;
	double dt;
	dt = DELZ*pPpt->xref/(EX ? pPpt->AREFe : pPpt->AREFi); // Set up time step

	// Calculate volume pressure
	// =========================
	double dQdt = 0;
/*
	switch(this->MODEL)
	{
	case CONST_P:
		p0 = this->p0_INIT; // T0 must be allowed to vary (especilly if variable volume cylinder)
		break;
//	case READ_P_T:
//		Interpolate(pPpt, THETA);	// Sets PC, TC
//		break;
	case RESET_P_T:
		p0 += dp0dt*dt; // For both RESET_P and WATSON cylinder models
		break;

	default:
		this->MODEL = RESET_P_T;
		cout << "CVolume::Update [" << ID << "]: Unknown volume model; selecting RESET_P_T\n";
		break;
	}
*/
	p0 += dp0dt*dt;

	// Volume heat transfer
	// ====================
//	switch(this->HT_MODEL)
//	{
//	case NONE:
//		dQdt += 0;
//		break;
//	default:
//		this->HT_MODEL = NONE;
//		cout << "CVolume::Update [" << ID << "]: Unknown heat transfer model; selecting NONE\n";
//		break;
//	}

	// Valve mass fluxes
	// =================
	bool UPDATE_PIPE = true; // Update pipe values with valve calls (as opposed to motored valve calls)
	
	// Get valve area
	for(i=0; i<this->NVOLUMEVALVES; ++i)
	{
		VolumeValve[i].EffectiveArea(pPpt); // Sets an up-to-date area

		if(VolumeValve[i].Get_eff_area()<=0.0) {dmdt[i]=0.0; VolumeValve[i].Set_open(false);}
		else VolumeValve[i].Set_open(true);
	}

	// Enter exhaust valve boundary
	double* aval; aval = new double [this->NVOLUMEVALVES];
	for(i=0; i<this->NVOLUMEVALVES; ++i)
	{
		VolumeValve[i].Set_psi_val(VolumeValve[i].Get_eff_area()/VolumeValve[i].Get_FP());

		pX_temp[i] = p0;
		TX_temp[i] = T0;
/*	
		pX_temp[i] = VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF;	
		TX_temp[i] = VolumeValve[i].Get_pBN_T();
		TX_temp[i] = T0;
*/
		FORWARD[i] = true;
	
		double p0_temp;
		//p0_temp = VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF + VolumeValve[i].Get_dynamic_head(pPpt)/1e5;
		//p0_temp = VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF;
		p0_temp = TotalPressureBar(pPpt, VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF/*ps, bar*/, VolumeValve[i].Get_pBN_T(), VolumeValve[i].Get_pBN_U()*(EX ? pPpt->AREFe : pPpt->AREFi));

///*		
		if(dmdt[i] > 0)
		{
/*
			cout << "Time = " << time << endl;
			cout << "p0 = " << p0 << endl;
			cout << "node ps = " << VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF << endl;
			pPpt->Out("\n");
//*/
			double Re = VolumeValve[i].Get_pBN_Re();
			
			double del_p_ideal = VolumeValve[i].Get_dynamic_head(pPpt)/1e5;
			//double del_p_ideal = TotalPressureBar(pPpt, VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF/*ps, bar*/, VolumeValve[i].Get_pBN_T(), VolumeValve[i].Get_pBN_U*(EX ? pPpt->AREFe : pPpt->AREFi)); 

			

			//double beta = Re/1e6;
			//double beta = Re/1.2e6;
			double beta = Re/1.3e6;

			beta = 1;

			if(beta<1) beta = 1;

		
			pX_temp[i] = p0 - beta*del_p_ideal;


			//pX_temp[i] = (pX_temp[i] + VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF)/2;

			//pX_temp[i] = (1-beta)*p0 + beta*VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF;
			//pX_temp[i] = beta*p0 + (1-beta)*VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF;
			//pX_temp[i] = VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF;
			//pX_temp[i] = (VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF + p0)/2;
			//pX_temp[i] = p0 + VolumeValve[i].Get_dynamic_head(pPpt)/1e5;
			//pX_temp[i] = p0_temp;
			
			// Inflow - temperature has no effect
			regime[i] = 1.8;
			FORWARD[i] = true;
		}
		else
		{

			//double del_p_ideal = VolumeValve[i].Get_dynamic_head(pPpt)/1e5;
			//pX_temp[i] = p0 + del_p_ideal;

			//pX_temp[i] = p0;
			//pX_temp[i] = VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF + VolumeValve[i].Get_dynamic_head(pPpt)/1e5;
			//pX_temp[i] = VolumeValve[i].Get_pBN_p_dash()*pPpt->PREF;
			//TX_temp[i] = T0;
			//TX_temp[i] = VolumeValve[i].Get_pBN_T();
			//TX_temp[i] = T0 + pow(VolumeValve[i].Get_pBN_U()*(EX ? pPpt->AREFe : pPpt->AREFi),2)/(2*pPpt->cpAir());
			//TX_temp[i] = T0 - pow(VolumeValve[i].Get_pBN_U()*(EX ? pPpt->AREFe : pPpt->AREFi),2)/(2*pPpt->cpAir());
			regime[i] = 0.2;

			FORWARD[i] = false;
		}
//*/

		if(!pPpt->HOMENTROPIC) 
			dmdt[i] = VolumeValve[i].Poppet_NH
			(pPpt, VolumeValve[i].Get_psi_val(), /*p0*/pX_temp[i], /*T0*/TX_temp[i], (EX ? pPpt->AREFe : pPpt->AREFi), ID, UPDATE_PIPE, timestep, time);
		else 
			VolumeValve[i].Poppet_H(pPpt, timestep, time, dmdt[i], VolumeValve[i].Get_psi_val(), pPpt->gammaAir(/*T0*/TX_temp[i]), (EX ? pPpt->AREFe : pPpt->AREFi), pPpt->PREF, p0, VolumeValve[i].Get_FP());

		aval[i] = ((VolumeValve[i].Get_CLOUT() + VolumeValve[i].Get_CLIN())/2)*(EX ? pPpt->AREFe : pPpt->AREFi);
	}

	// Sum mass fluxes for all valves
	double dmdt_total = 0;	for(i=0; i<this->NVOLUMEVALVES; ++i) dmdt_total += dmdt[i];
	
	m += dmdt_total*dt;	// New mass of gas in volume - here dmdt_total is positive when entering volume
	

	// Calculate cylinder volume and dVdt
	dVdt = 0.0;			// Constant volume volume


	T0 = 1.0E5*p0*V/(287.0*m); // Do not calculate T0 if read from file
	a0 = sqrt(pPpt->gammaAir(T0)*287.0*T0);


	// Check direction of flow for exhaust and select appropriate speed of sound
	for(i=0; i<this->NVOLUMEVALVES; ++i)
	{

//		if(FORWARD[i]) av[i] = aval[i];
//		else av[i] = a0;

///*
		if(dmdt[i]<0.0) av[i] = a0;
		if(dmdt[i]>0.0) av[i] = aval[i];
		if(dmdt[i]==0.0) av[i] = 0.0;
//*/
//av[i] = aval[i];
//av[i] = 300;
//av[i] = 0.1*av[i];
//av[i] = 0;
	}

	// Calculate pressure rates of change due to individual mass fluxes
	double dp0dt_dmdt = 0; for(i=0; i<this->NVOLUMEVALVES; ++i) dp0dt_dmdt += dmdt[i]*pow(av[i], 2)/(V*1.0E5);

	// Calculate overall pressure rate of change
	dp0dt = -pPpt->gammaAir(T0)*dVdt*p0/V					// Volume change component
				+ dp0dt_dmdt								// Valve mass flux component
				+ ((pPpt->gammaAir(T0)-1)/(V*1.0e5))*dQdt		// Heat transfer component (includes heat release due to combustion) 
				;
	return;
}


void CVolume::ReadInput(CProperties* pPpt, char *InputFile)
{
	int last_entry;
	CommonReadInput(InputFile, pPpt->NUM_PARAMS, this->labels, this->values, this->strings, last_entry);
	for(int r=0; r<last_entry+1; ++r)
	{
		// Volume geometry
		// --------------------------------------------------------------------------------
		if(strcmp(labels[r], "V_INIT") == 0) V_INIT = values[r];
		if(strcmp(labels[r], "USE_V_INIT") == 0) USE_V_INIT = DoubleToBool(values[r]);
		if(strcmp(labels[r], "length") == 0) length = values[r]/1000; // Conversion from mm to m
		if(strcmp(labels[r], "LINEAR") == 0) LINEAR = DoubleToBool(values[r]);
		if(strcmp(labels[r], "D_odd") == 0) D_odd = values[r]/1000;
		if(strcmp(labels[r], "D_int") == 0) D_int = values[r]/1000;
		if(strcmp(labels[r], "x_int") == 0) x_int = values[r]/1000;
		if(strcmp(labels[r], "D_even") == 0) D_even = values[r]/1000;

		// Initial conditions
		// --------------------------------------------------------------------------------
		if(strcmp(labels[r], "p0_INIT") == 0) p0_INIT = values[r];
		if(strcmp(labels[r], "T0_INIT") == 0) T0_INIT = values[r];

		// Valve parameters
		// --------------------------------------------------------------------------------
		if(strcmp(labels[r], "USE_PIPE_DIA") == 0) USE_PIPE_DIA = DoubleToBool(values[r]);
		if(strcmp(labels[r], "ref_dia") == 0) ref_dia = values[r];

		// Measurements
		// --------------------------------------------------------------------------------
		if(strcmp(labels[r], "USE_DEF_FREQ") == 0) USE_DEF_FREQ = DoubleToBool(values[r]);
		if(strcmp(labels[r], "freq") == 0) freq = int(values[r]);
	}

	// Set some derived parameters
	if(USE_DEF_FREQ) freq = pPpt->freq;		// Use the default sampling rate

	// Pipe gradient calculation - don't need to consider end corrections as we want the physical gradient of the pipe
	C = (D_even - D_odd)/(length/pPpt->xref); // For linear variation
	if(D_even==D_odd) vol = (PI/4)*pow(D_odd,2)*length;
	else vol = (PI/12)*(length/(D_even - D_odd))*(pow(D_even,3) - pow(D_odd,3));
	if(!LINEAR)
	{
		// Pipe diameter variation
		if(x_int>length || x_int<0)
		{
			cout << "Interior point must lie between odd and even ends of pipe\n" << endl;
			exit(1);
		}
		else
		{			
			// Caclulate quadratic variation of pipe diameter given 3 points: odd, even and interior
			// dia = aX^2 + bX + c
			double X1, y1, X2, y2, X3, y3;
			X1 = 0/pPpt->xref; y1 = D_odd;
			X2 = x_int/pPpt->xref; y2 = D_int;
			X3 = length/pPpt->xref; y3 = D_even;
			a = ((y3-y2)/((X3-X2)*(X3-X1))) - ((y1-y2)/((X1-X2)*(X3-X1)));
			//b = ((y3-y2) + a*(pow(X2,2)-pow(X3,2)))/(X2-X3);
			b = ((y1-y2) + a*(pow(X2,2)-pow(X1,2)))/(X1-X2);
			//c = y1 - a*pow(X1,2) - b*X1;
			c = y2 - a*pow(X2,2) - b*X2;
			//c = y3 - a*pow(X3,2) - b*X3;

			double xtemp = length;
			vol = (PI/4)*(
				  (pow(a,2)*pow(xtemp,5))/(5*pow(pPpt->xref,4))
				+ (a*b*pow(xtemp,4))/(2*pow(pPpt->xref,3))
				+ ((2*a*c + pow(b,2))*pow(xtemp,3))/(3*pow(pPpt->xref,2))
				+ (b*c*pow(xtemp,2))/pPpt->xref
				+ pow(c,2)*xtemp );
		}
	}
}

void CVolume::ListProperties(CProperties* pPpt)
{
	pPpt->Out(Underline(Identify(), "="));
	if(USE_V_INIT)
	{
		pPpt->Out("\tUsing V_INIT directly as volume:\n");
		pPpt->Out("\n");
		pPpt->Out("\t(Initial) volume, V_INIT\t\t\t=\t"); pPpt->Out(V_INIT); pPpt->Out(" m^3\n");
	}
	else
	{
		pPpt->Out("\tVolume derived from the following pipe geometry:\n");
		pPpt->Out("\n");
		pPpt->Out("\tPhysical length\t\t\t\t\t=\t"); pPpt->Out(length*1000); pPpt->Out(" mm\n");
		if(LINEAR)
		{
			pPpt->Out("\tLinear pipe diameter variation:\n");
			pPpt->Out("\t- dia. odd end, D_odd\t\t\t\t=\t"); pPpt->Out(D_odd*1000); pPpt->Out(" mm\n");
			pPpt->Out("\t- dia. even end, D_even\t\t\t=\t"); pPpt->Out(D_even*1000); pPpt->Out(" mm\n");
			pPpt->Out("\t- D(X) = "); pPpt->Out(D_odd); pPpt->Out(" + "); pPpt->Out(C); pPpt->Out("X\n");
		}
		else
		{
			pPpt->Out("\tQuadratic pipe diameter variation:\n");
			pPpt->Out("\t- dia. odd end, D_odd\t\t\t\t=\t"); pPpt->Out(D_odd*1000); pPpt->Out(" mm\n");
			pPpt->Out("\t- dia. interior point, D_int\t\t\t=\t"); pPpt->Out(D_int*1000); pPpt->Out(" mm\n");
			pPpt->Out("\t- distance of int. point from odd end, x_int\t=\t"); pPpt->Out(x_int*1000); pPpt->Out(" mm\n");
			pPpt->Out("\t- dia. even end, D_even\t\t\t=\t"); pPpt->Out(D_even*1000); pPpt->Out(" mm\n");
			pPpt->Out("\t- D(X) = "); pPpt->Out(a); pPpt->Out("X^2 + "); pPpt->Out(b); pPpt->Out("X + "); pPpt->Out(c); pPpt->Out("\n");
		}
		pPpt->Out("\t- volume, vol\t\t\t\t\t=\t"); pPpt->Out(vol*1000); pPpt->Out(" l^3\n");
	}
	pPpt->Out("\n");
	pPpt->Out("\tInitial pressure, p0_INIT\t\t\t=\t"); pPpt->Out(p0_INIT); pPpt->Out(" bar\n");
	pPpt->Out("\tInitial temperature, T0_INIT\t\t\t=\t"); pPpt->Out(T0_INIT); pPpt->Out(" K\n");
	pPpt->Out("\n");
	pPpt->Out("\tNo. pipes joining this volume, NVOLUMEVALVES\t=\t"); pPpt->Out(NVOLUMEVALVES); pPpt->Out(":\n");
	pPpt->Out("\n");
	for(int i=0; i<NVOLUMEVALVES; ++i)
	{
		VolumeValve[i].ListProperties(pPpt, "\t");
	}
	pPpt->Out("\n");
	if(USE_DEF_FREQ)
	{
		if(freq==1) pPpt->Out("\tUsing default sampling rate, pPpt->freq\t=\tonce per timestep\n");
		else
		{
			pPpt->Out("\tUsing default sampling rate, pPpt->freq\t=\tonce per "); pPpt->Out(freq); pPpt->Out(" timesteps\n");
		}
	}
	else
	{
		if(freq==1) pPpt->Out("\tUsing local sampling rate (freq)\t\t=\tonce per timestep\n");
		else
		{
			pPpt->Out("\tUsing local sampling rate (freq)\t\t=\tonce per "); pPpt->Out(freq); pPpt->Out(" timesteps\n");
		}
	}
	pPpt->Out("\n");
	pPpt->Out("\n");
}

void CVolume::PrintToScreen(CProperties* pPpt)
{
	if(EX) cout << Underline("Exhaust Volume", "-", ID);
	else cout << Underline("Intake Volume", "-", ID);
	// -------------------------------------------------
	cout << "Current volume (V)\t\t\t\t\t=\t" << V << " m^3\n";
	cout << "Stagnation pressure (p0)\t\t\t\t=\t" << p0 << " bar\n";
	cout << "Stagnation temperature (T0)\t\t\t\t=\t" << T0 << " K\n";
	pPpt->Out("\n");
	for(int i=0; i<NVOLUMEVALVES; ++i)
	{
		VolumeValve[i].PrintToScreenVolume(pPpt);
	}

	pPpt->Out("\n");
}

void CVolume::PrintToFile(CProperties* pPpt, int timestep, double time, double ca)
{
	int i;
	if(timestep==0)
	{
		fprintf(OUTPUT_FILE,"%s", Underline("RESULTS FILE FOR VOLUME", "-", ID));
		fprintf(OUTPUT_FILE,"\n");
		fprintf(OUTPUT_FILE,"%s", "Time(s)");
		fprintf(OUTPUT_FILE,"\t\t%s\t%s\t%s\t%s", "p0 (bar)", "T0 (K)", "V (m^3)", "m (kg)");
		
		for(i=0; i<NVOLUMEVALVES; ++i)
		{
			fprintf(OUTPUT_FILE,"\t\t%s%i%s\t%s%i%s\t%s%i%s\t%s%i%s", "pX_temp[" , i ,"]", "TX_temp[" , i ,"]", "regime[", i, "]", "dmdt[" , i ,"]");

		}		
		fprintf(OUTPUT_FILE,"\n");
	}

	if(timestep%freq==0) // Print data at the specified sampling frequency
	{
		fprintf(OUTPUT_FILE,"%f\t\t%f\t%f\t%f\t%f", time, p0, T0, V, m);
		
		for(i=0; i<NVOLUMEVALVES; ++i)
		{
			fprintf(OUTPUT_FILE,"\t\t%f\t%f\t%f\t%f",
				pX_temp[i],
				TX_temp[i],
				regime[i],
				dmdt[i]
				);
		}
		fprintf(OUTPUT_FILE,"\n");
	}
}

void CVolume::SetupFiles(CProperties* pPpt, string parent_assy_res_dir)
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
	cout << RES_DIR << endl;
*/
	RES_DIR = StringToChar(parent_assy_res_dir);

	std::string bcname_str;
	if(EX) bcname_str = "res_ex_vol"; else  bcname_str = "res_in_vol";
	OUTPUT_FILE = fopen(ConstructString(pPpt, RES_DIR, bcname_str, ID), "w");
}

void CVolume::CloseFiles()
{
	fclose(OUTPUT_FILE);
}

char* CVolume::Identify()
// ============================================================ //
// Returns identification of the current object					//
// ============================================================ //
{
	std::string sz;
	sz = "Assembly [";
	sz += IntToString(AssyID);
	sz += "], ";
	sz += "Volume [";
	sz += IntToString(ID);
	sz += "]";

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str());
	return szz;
}