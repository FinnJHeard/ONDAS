// Anechoic.cpp: implementation of the CAnechoic class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Anechoic.h"
#include "Tools.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CAnechoic::CAnechoic()
{

}

CAnechoic::~CAnechoic()
{

}

void CAnechoic::Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rEPIPES, int** &rEPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, std::string param_dir, int assyid, string parent_assy_res_dir, string calling_object_str)
{
	if (pPpt->SHOW_calls) { pPpt->Out("CAnechoic.Initialise ("); pPpt->Out(calling_object_str); pPpt->Out(")\n"); }

	// Standard initialization
	// ====================================================================================================
	InitialiseGen(pPpt, pPipes, rPipe, rEPIPES, rEPIPES_ENDS, rENDCORR, id, ex, 1/*npipes*/, assyid, calling_object_str, parent_assy_res_dir);
	// Note only want to configure it with the real pipe; the buffer pipe is internal, hence npipes==1 above
	
	std::string bcname_str = "ANEC";

	Buffer.BUFFER = true;
	Buffer.DAMPER = false;

	Damper.DAMPER = true;
	Damper.BUFFER = false;

	ReadInput(pPpt, ConstructString(pPpt, param_dir, bcname_str, EX, ID));

	// Boundary node labels
	REAL_PIPE = 0;
	BUFFER_PIPE = 1;
	
	NAME = new int [NPIPES];
	NAME[REAL_PIPE] = ANECHOIC;
	//NAME[BUFFER_PIPE] = ANECHOIC;

	SONIC = false;
}

void CAnechoic::InitialiseBufferDamper(CProperties* pPpt)
{
	//cout << "HERE" << endl;
	//cout << Buffer.METHOD << endl;
	//exit(1);
	// Buffer initialization
	// ====================================================================================================
	int K, S;
	double k;
	// Set NULL values:
	Buffer.bend_angle = 0;
	Buffer.C_p = 0;
	Damper.bend_angle = 0;
	Damper.C_p = 0;
	
	Buffer.odd_end_flow = NOFLOW; Buffer.even_end_flow = NOFLOW;
	Buffer.MMOC=1; Buffer.W_ALPHA_BETA=2; Buffer.FandE = 3; Buffer.JOINER=4;	// Numerate METHOD labels
	Damper.odd_end_flow = NOFLOW; Damper.even_end_flow = NOFLOW;
	Damper.MMOC=1; Damper.W_ALPHA_BETA=2; Damper.FandE = 3; Damper.JOINER=4;	// Numerate METHOD labels

	// Need to put in anechoic file:
	//Buffer.length = this->pPipe[this->REAL_PIPE]->length;
	//Buffer.min_meshes = this->pPipe[this->REAL_PIPE]->min_meshes;

	//Buffer.discret = this->pPipe[this->REAL_PIPE]->discret;//10; //mm
	Buffer.end_corr_odd_p = 0;
	Buffer.end_corr_even_p = 0;
	Damper.end_corr_odd_p = 0;
	Damper.end_corr_even_p = 0;
	//Buffer.nsections = 1;

	if(pBN[REAL_PIPE]->side == ODD){
		Buffer.d_odd = pPipe[REAL_PIPE]->d_odd;
		Buffer.d_even = Buffer.d_odd;
	}
	else{
		Buffer.d_even = pPipe[REAL_PIPE]->d_even;
		Buffer.d_odd = Buffer.d_even;
	}
	Damper.d_odd = Buffer.d_even;
	Damper.d_even = Damper.d_odd;

	// Set end corrections
	Buffer.end_corr_odd		= Buffer.end_corr_odd_p*(Buffer.d_odd);
	Buffer.end_corr_even	= Buffer.end_corr_even_p*(Buffer.d_even);
	Buffer.eff_length		= Buffer.length + (Buffer.end_corr_odd + Buffer.end_corr_even);
	Damper.end_corr_odd		= Damper.end_corr_odd_p*(Damper.d_odd);
	Damper.end_corr_even	= Damper.end_corr_even_p*(Damper.d_even);
	Damper.eff_length		= Damper.length + (Damper.end_corr_odd + Damper.end_corr_even);

	// Determine correct number of meshes based on the desired discretisation length
	if((Buffer.eff_length/Buffer.discret) - int(Buffer.eff_length/Buffer.discret)>1e-6) 
			Buffer.meshes = int(Buffer.eff_length/Buffer.discret) + 1;
	else	Buffer.meshes = int(Buffer.eff_length/Buffer.discret);
	if((Damper.eff_length/Damper.discret) - int(Damper.eff_length/Damper.discret)>1e-6) 
			Damper.meshes = int(Damper.eff_length/Damper.discret) + 1;
	else	Damper.meshes = int(Damper.eff_length/Damper.discret);


	// Enforce a minimum number of meshes
	if(Buffer.meshes<Buffer.min_meshes) Buffer.meshes = Buffer.min_meshes; 
	if(Damper.meshes<Damper.min_meshes) Damper.meshes = Damper.min_meshes; 

	Buffer.XPIPE				= Buffer.length/pPpt->xref;		// N.D. pipe physical length ()
	Buffer.N					= Buffer.meshes + 1;			// For a pipe with seperate ends
	Buffer.Node					= new CNode[Buffer.N];			// Dimension Node vector
	Buffer.Node_Backup			= new CNode[Buffer.N];			// Dimension backup Node vector
	Damper.XPIPE				= Damper.length/pPpt->xref;		// N.D. pipe physical length ()
	Damper.N					= Damper.meshes + 1;			// For a pipe with seperate ends
	Damper.Node					= new CNode[Damper.N];			// Dimension Node vector
	Damper.Node_Backup			= new CNode[Damper.N];			// Dimension backup Node vector

	if(Buffer.meshes>0) Buffer.xmesh	= Buffer.eff_length/Buffer.meshes;	// Mesh effective length (m)
	else				Buffer.xmesh	= Buffer.eff_length;					// Special case for a single node
	Buffer.XMESH						= Buffer.xmesh/pPpt->xref;			// N.D. mesh effective length ()
	Buffer.XPIPE_EFF					= Buffer.eff_length/pPpt->xref;		// N.D. pipe effective length ()

	if(Damper.meshes>0) Damper.xmesh	= Damper.eff_length/Damper.meshes;	// Mesh effective length (m)
	else				Damper.xmesh	= Damper.eff_length;					// Special case for a single node
	Damper.XMESH						= Damper.xmesh/pPpt->xref;			// N.D. mesh effective length ()
	Damper.XPIPE_EFF					= Damper.eff_length/pPpt->xref;		// N.D. pipe effective length ()


	for(S=0; S<Buffer.N; ++S)
	{
		Buffer.Node[S].bc = INTERIOR;
		Buffer.Node[S].side = INSIDE;
	}
	Buffer.Node[0].side = ODD;
	Buffer.Node[Buffer.N-1].side = EVEN;
	for(S=0; S<Damper.N; ++S)
	{
		Damper.Node[S].bc = INTERIOR;
		Damper.Node[S].side = INSIDE;
	}
	Damper.Node[0].side = ODD;
	Damper.Node[Damper.N-1].side = EVEN;

	if(pPpt->NUM_PATH_MULT*(2*Buffer.N - 1) > 2*Buffer.N - 1) Buffer.num_pathlines = pPpt->NUM_PATH_MULT*(2*Buffer.N - 1);
	else Buffer.num_pathlines = 2*Buffer.N - 1;
	if(pPpt->NUM_PATH_MULT*(2*Damper.N - 1) > 2*Damper.N - 1) Damper.num_pathlines = pPpt->NUM_PATH_MULT*(2*Damper.N - 1);
	else Damper.num_pathlines = 2*Damper.N - 1;

	if(!pPpt->HOMENTROPIC) // Pathlines exist only for non-homentropic cases
	{
		Buffer.PathLine = new CPathLine[Buffer.num_pathlines];
		Buffer.PathLine_Backup = new CPathLine[Buffer.num_pathlines];
	}
	if(!pPpt->HOMENTROPIC) // Pathlines exist only for non-homentropic cases
	{
		Damper.PathLine = new CPathLine[Damper.num_pathlines];
		Damper.PathLine_Backup = new CPathLine[Damper.num_pathlines];
	}

	// Need to define x (m), X () for all nodes
	for(S=0; S<Buffer.N; ++S)
	{
		Buffer.Node[S].x = -Buffer.end_corr_odd + S*Buffer.xmesh;
		Buffer.Node[S].X = Buffer.Node[S].x/pPpt->xref;
	}
	for(S=0; S<Damper.N; ++S)
	{
		Damper.Node[S].x = -Damper.end_corr_odd + S*Damper.xmesh;
		Damper.Node[S].X = Damper.Node[S].x/pPpt->xref;
	}

	// Calculate quadratic variation of pipe diameter given 3 points: odd, even and interior
	Buffer.x_int = Buffer.length/2;
	Buffer.d_int = (Buffer.d_odd + Buffer.d_even)/2;
	Damper.x_int = Damper.length/2;
	Damper.d_int = (Damper.d_odd + Damper.d_even)/2;

	double y1, y2, y3;
	y1 = Buffer.d_odd;
	y2 = Buffer.d_int;
	y3 = Buffer.d_even;
	// dia = ax^2 + bx + c
	double x1, x2, x3;
	x1 = 0; 
	x2 = Buffer.x_int; 
	x3 = Buffer.length;
/*
	Buffer.a = ((y3-y2)/((x3-x2)*(x3-x1))) - ((y1-y2)/((x1-x2)*(x3-x1)));
	Buffer.b = ((y1-y2) + Buffer.a*(pow(x2,2)-pow(x1,2)))/(x1-x2);
	Buffer.c = y2 - Buffer.a*pow(x2,2) - Buffer.b*x2;
*/	
/*
	//double y1, y2, y3;
	y1 = Damper.d_odd;
	y2 = Damper.d_int;
	y3 = Damper.d_even;
	//double x1, x2, x3;
	x1 = 0; 
	x2 = Damper.x_int; 
	x3 = Damper.length;
	Damper.a = ((y3-y2)/((x3-x2)*(x3-x1))) - ((y1-y2)/((x1-x2)*(x3-x1)));
	Damper.b = ((y1-y2) + Damper.a*(pow(x2,2)-pow(x1,2)))/(x1-x2);
	Damper.c = y2 - Damper.a*pow(x2,2) - Damper.b*x2;
	double x = Buffer.length;
	Buffer.vol = (PI/4)*( (pow(Buffer.a,2)/5)*pow(x,5) + (Buffer.a*Buffer.b/2)*pow(x,4) 
							+ ((2*Buffer.a*Buffer.c + pow(Buffer.b,2))/3)*pow(x,3) + Buffer.b*Buffer.c*pow(x,2) + pow(Buffer.c,2)*x );

	x = Damper.length;
	Damper.vol = (PI/4)*( (pow(Damper.a,2)/5)*pow(x,5) + (Damper.a*Damper.b/2)*pow(x,4) 
							+ ((2*Damper.a*Damper.c + pow(Damper.b,2))/3)*pow(x,3) + Damper.b*Damper.c*pow(x,2) + pow(Damper.c,2)*x );
*/

/*
	// Still required for MOC:
	// dia = AX^2 + BX + C
	{
		double X1, X2, X3;
		X1 = 0/pPpt->xref; 
		X2 = Buffer.x_int/pPpt->xref; 
		X3 = Buffer.length/pPpt->xref;
		Buffer.A = ((y3-y2)/((X3-X2)*(X3-X1))) - ((y1-y2)/((X1-X2)*(X3-X1)));
		Buffer.B = ((y1-y2) + Buffer.A*(pow(X2,2)-pow(X1,2)))/(X1-X2);
		Buffer.C = y2 - Buffer.A*pow(X2,2) - Buffer.B*X2;

		//double X1, X2, X3;
		X1 = 0/pPpt->xref; 
		X2 = Damper.x_int/pPpt->xref; 
		X3 = Damper.length/pPpt->xref;
		Damper.A = ((y3-y2)/((X3-X2)*(X3-X1))) - ((y1-y2)/((X1-X2)*(X3-X1)));
		Damper.B = ((y1-y2) + Damper.A*(pow(X2,2)-pow(X1,2)))/(X1-X2);
		Damper.C = y2 - Damper.A*pow(X2,2) - Damper.B*X2;
	}
//cout << "dia = " << a << "x^2 + " << b << "x + " << c << endl << endl;	
//*/

	for(S=0; S<Buffer.N; ++S)
	{
		double x = Buffer.Node[S].x;
		Buffer.Node[S].d = Buffer.d(x);
		Buffer.Node[S].f = Buffer.f(x);
		Buffer.Node[S].f_dash = Buffer.Node[S].f/pPpt->fref;
		Buffer.Node[S].dfdx = Buffer.dfdx(x);
	}
	for(S=0; S<Damper.N; ++S)
	{
		double x = Damper.Node[S].x;
		Damper.Node[S].d = Damper.d(x);
		Damper.Node[S].f = Damper.f(x);
		Damper.Node[S].f_dash = Damper.Node[S].f/pPpt->fref;
		Damper.Node[S].dfdx = Damper.dfdx(x);
	}



	// Check that x value of last node equals (physical length of pipe plus even end correction)
	if(fabs(Buffer.Node[Buffer.N-1].x - (Buffer.length + Buffer.end_corr_even)) > 1e-9)
	{
		pPpt->Out(Identify());  pPpt->Out(": geometry mismatch:\n");
		pPpt->Out("Buffer.Node[Buffer.N-1].x = "); pPpt->Out(Buffer.Node[Buffer.N-1].x); pPpt->Out(" m\n");
		pPpt->Out("(Buffer.length + Buffer.end_corr_even) = "); pPpt->Out(Buffer.length + Buffer.end_corr_even); pPpt->Out(" m\n");
		exit(1);
	}
	// Check that x value of last node equals (physical length of pipe plus even end correction)
	if(fabs(Damper.Node[Damper.N-1].x - (Damper.length + Damper.end_corr_even)) > 1e-9)
	{
		pPpt->Out(Identify()); pPpt->Out(": geometry mismatch:\n");
		pPpt->Out("Damper.Node[Damper.N-1].x = "); pPpt->Out(Damper.Node[Damper.N-1].x); pPpt->Out(" m\n");
		pPpt->Out("(Damper.length + Damper.end_corr_even) = "); pPpt->Out(Damper.length + Damper.end_corr_even); pPpt->Out(" m\n");
		exit(1);
	}

	// Dimension intialisation parameters using nsections
	Buffer.nsections = 1;
	Buffer.xri = new double [Buffer.nsections];
	Buffer.pri = new double [Buffer.nsections];
	Buffer.Tri = new double [Buffer.nsections];
	Buffer.vri = new double [Buffer.nsections];
			
	Buffer.xri[Buffer.nsections-1] = 1.0; // Always add an imaginary division at the end of the pipe
	
	Buffer.pri[0] = this->pBN[this->REAL_PIPE]->p_dash;
	Buffer.Tri[0] = this->pBN[this->REAL_PIPE]->T;
	Buffer.vri[0] = this->pBN[this->REAL_PIPE]->U/this->pPipe[this->REAL_PIPE]->AREF;

	Damper.nsections = 1;
	Damper.xri = new double [Damper.nsections];
	Damper.pri = new double [Damper.nsections];
	Damper.Tri = new double [Damper.nsections];
	Damper.vri = new double [Damper.nsections];
			
	Damper.xri[Damper.nsections-1] = 1.0; // Always add an imaginary division at the end of the pipe
	
	Damper.pri[0] = this->pBN[this->REAL_PIPE]->p_dash;
	Damper.Tri[0] = this->pBN[this->REAL_PIPE]->T;
	Damper.vri[0] = this->pBN[this->REAL_PIPE]->U/this->pPipe[this->REAL_PIPE]->AREF;




	// Initial Riemann values
	double CLI;
	int section;
	section = 0;

	for(S=0; S<Buffer.N; ++S)
	{
		// Node spacing
		if(S<Buffer.N-1) Buffer.Node[S].DELX_R = Buffer.Node[S+1].X - Buffer.Node[S].X;	else Buffer.Node[S].DELX_R = 0;		
		if(S>0) Buffer.Node[S].DELX_L = Buffer.Node[S].X - Buffer.Node[S-1].X; else Buffer.Node[S].DELX_L = 0;
	
		if(Buffer.Node[S].DELX_R < 0) 
		{
			pPpt->Out(Identify()); pPpt->Out(": Buffer.DELX_R negative!\n");
			pPpt->Out("Buffer.end_corr_odd_p = "); pPpt->Out(Buffer.end_corr_odd_p); pPpt->Out("\n");
			pPpt->Out("Buffer.end_corr_even_p = "); pPpt->Out(Buffer.end_corr_even_p); pPpt->Out("\n");
			pPpt->Out("\n");
			exit(1);
		}
		if(Buffer.Node[S].DELX_L < 0) 
		{
			pPpt->Out(Identify()); pPpt->Out(": DELX_L negative!\n"); 
			pPpt->Out("Buffer.end_corr_odd_p = "); pPpt->Out(Buffer.end_corr_odd_p); pPpt->Out("\n");
			pPpt->Out("Buffer.end_corr_even_p = "); pPpt->Out(Buffer.end_corr_even_p); pPpt->Out("\n");
			pPpt->Out("\n");
			exit(1);
		}
			
		// Determine which section this node is in
//		section = 0;
//		while(Buffer.eff_length*Buffer.xri[section] < Buffer.Node[S].x && section < Buffer.nsections-1)
//		{
//			//cout << "eff_length*xri[section=" << section << "] = " << eff_length*xri[section] << endl;
//			//cout << "Node[S=" << S << "].x = " << Node[S].x << endl;
//			++section;
//		}

		if(pPpt->HOMENTROPIC)
		{
			//CLI = pow(Buffer.pri[section]/pPpt->PREF, pPpt->Q);
			CLI = pow(this->pBN[this->REAL_PIPE]->p_dash, pPpt->Q);

			// Sets pressure == PRI; see p271 Simple Homentropic Program

			Buffer.Node[S].CL1[0] = CLI;
			Buffer.Node[S].CL2[0] = CLI;
			Buffer.Node[S].CL1[1] = Buffer.Node[S].CL1[0];
			Buffer.Node[S].CL2[1] = Buffer.Node[S].CL2[0];

			// For homentropic, need to set TREF in order to show correct TRI; see p271 Simple Homentropic Program
			// Cannot have multiple section temps under homentropic flow though due to there only being one TREF
			double TREF_temp;
			if(pPpt->CONTINUOUS){
				k = pPpt->gammaAir(Buffer.Tri[section]);
				TREF_temp = pow(pPpt->PREF/Buffer.pri[section], (k-1)/k)*Buffer.Tri[section];
			}
			else{
				TREF_temp = 300;
			}
			
			// Need to have one AREF for entire exhaust system, and another separate one for entire intake system
			if(Buffer.EX)
			{
				k = pPpt->gammaAir(pPpt->TREFe);
				pPpt->TREFe = TREF_temp;
				pPpt->AREFe = sqrt(k*pPpt->R_air*pPpt->TREFe);
				Buffer.AREF = pPpt->AREFe;
			}
			else
			{
				k = pPpt->gammaAir(pPpt->TREFi);
				pPpt->TREFi = TREF_temp;
				pPpt->AREFi = sqrt(k*pPpt->R_air*pPpt->TREFi);
				Buffer.AREF = pPpt->AREFi;
			}
			/*
			cout << "TREF_temp = " << TREF_temp << endl;
			cout << "pPpt->TREFe = " << pPpt->TREFe << endl;
			cout << "pPpt->TREFi = " << pPpt->TREFi << endl;
			cout << "pPpt->AREFe = " << pPpt->AREFe << endl;
			cout << "pPpt->AREFi = " << pPpt->AREFi << endl;
			*/
		}
		else
		{			
			// Choose TREF as the initial temperature in the first pipe (first section) of that system and calculate AREF		
			if(Buffer.EX)
			{
				if(Buffer.ID==0 && Buffer.AssyID==0) // If this is the first exhaust pipe in the first assembly
				{
					k = pPpt->gammaAir(pPpt->TREFe);
					pPpt->TREFe = Buffer.Tri[/*section=*/0];
					//pPpt->TREFe = this->pBN[this->REAL_PIPE]->T;
					pPpt->AREFe = sqrt(k*pPpt->R_air*pPpt->TREFe);
				}
				Buffer.AREF = pPpt->AREFe;
			}
			else
			{
				if(Buffer.ID==0 && Buffer.AssyID==0) // If this is the first intake pipe in the first assembly
				{
					k = pPpt->gammaAir(pPpt->TREFi);
					pPpt->TREFi = Buffer.Tri[/*section=*/0];
					//pPpt->TREFe = this->pBN[this->REAL_PIPE]->T;
					pPpt->AREFi = sqrt(k*pPpt->R_air*pPpt->TREFi);
				}
				Buffer.AREF = pPpt->AREFi;
			}

			// For variable initial velocity
			k = pPpt->gammaAir(pBN[REAL_PIPE]->T);
			Buffer.Node[S].CL1[0] = sqrt(pBN[REAL_PIPE]->T/(Buffer.EX ? pPpt->TREFe : pPpt->TREFi)) + ((k-1)/2)*((pBN[REAL_PIPE]->U*Buffer.AREF)/(Buffer.EX ? pPpt->AREFe : pPpt->AREFi));
			Buffer.Node[S].CL2[0] = sqrt(pBN[REAL_PIPE]->T/(Buffer.EX ? pPpt->TREFe : pPpt->TREFi)) - ((k-1)/2)*((pBN[REAL_PIPE]->U*Buffer.AREF)/(Buffer.EX ? pPpt->AREFe : pPpt->AREFi));	
			Buffer.Node[S].CL1[1] = Buffer.Node[S].CL1[0];
			Buffer.Node[S].CL2[1] = Buffer.Node[S].CL2[0];
		}

		// Generate initial A and U
		Buffer.Node[S].A = (Buffer.Node[S].CL1[0] + Buffer.Node[S].CL2[0])/2;
		Buffer.Node[S].T = pow(Buffer.Node[S].A,2)*(Buffer.EX ? pPpt->TREFe : pPpt->TREFi);
		k = pPpt->gammaAir(Buffer.Node[S].T);
		Buffer.Node[S].U = (Buffer.Node[S].CL1[0] - Buffer.Node[S].CL2[0])/(k-1);
		//Buffer.Node[S].p_dash = Buffer.pri[section]/pPpt->PREF;
		Buffer.Node[S].p_dash = this->pBN[this->REAL_PIPE]->p_dash;
//if(ID==1) cout << "pri[section] = " << pri[section] << endl;
//cout << "Node[S=" << S << "].p_dash = " << Node[S].p_dash << endl;
		Buffer.Node[S].rho = (Buffer.Node[S].p_dash*pPpt->PREF*1e5)/(pPpt->R_air*Buffer.Node[S].T);
		Buffer.Node[S].rho_prev = Buffer.Node[S].rho;
		Buffer.Node[S].mdot = Buffer.Node[S].rho*(Buffer.Node[S].U*Buffer.AREF)*Buffer.Node[S].f;
		Buffer.Node[S].CHOKED = false;
	}
	for(S=0; S<Damper.N; ++S)
	{
		// Node spacing
		if(S<Damper.N-1) Damper.Node[S].DELX_R = Damper.Node[S+1].X - Damper.Node[S].X;	else Damper.Node[S].DELX_R = 0;		
		if(S>0) Damper.Node[S].DELX_L = Damper.Node[S].X - Damper.Node[S-1].X; else Damper.Node[S].DELX_L = 0;
	
		if(Damper.Node[S].DELX_R < 0) 
		{
			cout << Identify() << ": Damper.DELX_R negative!" << endl;
			cout << "Damper.end_corr_odd_p = " << Damper.end_corr_odd_p << endl;
			cout << "Damper.end_corr_even_p = " << Damper.end_corr_even_p << endl;
			cout << endl;
			exit(1);
		}
		if(Damper.Node[S].DELX_L < 0) 
		{
			cout << Identify() << ": DELX_L negative!" << endl; 
			cout << "Damper.end_corr_odd_p = " << Damper.end_corr_odd_p << endl;
			cout << "Damper.end_corr_even_p = " << Damper.end_corr_even_p << endl;
			cout << endl;
			exit(1);
		}
			
		// Determine which section this node is in
//		section = 0;
//		while(Damper.eff_length*Damper.xri[section] < Damper.Node[S].x && section < Damper.nsections-1)
//		{
//			//cout << "eff_length*xri[section=" << section << "] = " << eff_length*xri[section] << endl;
//			//cout << "Node[S=" << S << "].x = " << Node[S].x << endl;
//			++section;
//		}

		if(pPpt->HOMENTROPIC){
			//CLI = pow(Damper.pri[section]/pPpt->PREF, pPpt->Q);
			CLI = pow(this->pBN[this->REAL_PIPE]->p_dash, pPpt->Q);
			// Sets pressure == PRI; see p271 Simple Homentropic Program

			Damper.Node[S].CL1[0] = CLI;
			Damper.Node[S].CL2[0] = CLI;
			Damper.Node[S].CL1[1] = Damper.Node[S].CL1[0];
			Damper.Node[S].CL2[1] = Damper.Node[S].CL2[0];

			// For homentropic, need to set TREF in order to show correct TRI; see p271 Simple Homentropic Program
			// Cannot have multiple section temps under homentropic flow though due to there only being one TREF
			double TREF_temp;
			if(pPpt->CONTINUOUS){
				k = pPpt->gammaAir(Damper.Tri[section]);
				TREF_temp = pow(pPpt->PREF/Damper.pri[section], (k-1)/k)*Damper.Tri[section];
			}
			else{
				TREF_temp = 300;
			}
			
			// Need to have one AREF for entire exhaust system, and another separate one for entire intake system
			if(Damper.EX){
				k = pPpt->gammaAir(pPpt->TREFe);
				pPpt->TREFe = TREF_temp;
				pPpt->AREFe = sqrt(k*pPpt->R_air*pPpt->TREFe);
				Damper.AREF = pPpt->AREFe;
			}
			else{
				k = pPpt->gammaAir(pPpt->TREFi);
				pPpt->TREFi = TREF_temp;
				pPpt->AREFi = sqrt(k*pPpt->R_air*pPpt->TREFi);
				Damper.AREF = pPpt->AREFi;
			}
		}
		else{			
			// Choose TREF as the initial temperature in the first pipe (first section) of that system and calculate AREF		
			if(Damper.EX){
				if(Damper.ID==0 && Damper.AssyID==0){ // If this is the first exhaust pipe in the first assembly
					k = pPpt->gammaAir(pPpt->TREFe);
					pPpt->TREFe = Damper.Tri[/*section=*/0];
					//pPpt->TREFe = this->pBN[this->REAL_PIPE]->T;
					pPpt->AREFe = sqrt(k*pPpt->R_air*pPpt->TREFe);
				}
				Damper.AREF = pPpt->AREFe;
			}
			else
			{
				if(Damper.ID==0 && Damper.AssyID==0) // If this is the first intake pipe in the first assembly
				{
					k = pPpt->gammaAir(pPpt->TREFi);
					pPpt->TREFi = Damper.Tri[/*section=*/0];
					//pPpt->TREFe = this->pBN[this->REAL_PIPE]->T;
					pPpt->AREFi = sqrt(k*pPpt->R_air*pPpt->TREFi);
				}
				Damper.AREF = pPpt->AREFi;
			}

			// For variable initial velocity
			k = pPpt->gammaAir(pBN[REAL_PIPE]->T);
			Damper.Node[S].CL1[0] = sqrt(pBN[REAL_PIPE]->T/(Damper.EX ? pPpt->TREFe : pPpt->TREFi)) + ((k-1)/2)*((pBN[REAL_PIPE]->U*Damper.AREF)/(Damper.EX ? pPpt->AREFe : pPpt->AREFi));
			Damper.Node[S].CL2[0] = sqrt(pBN[REAL_PIPE]->T/(Damper.EX ? pPpt->TREFe : pPpt->TREFi)) - ((k-1)/2)*((pBN[REAL_PIPE]->U*Damper.AREF)/(Damper.EX ? pPpt->AREFe : pPpt->AREFi));	
			Damper.Node[S].CL1[1] = Damper.Node[S].CL1[0];
			Damper.Node[S].CL2[1] = Damper.Node[S].CL2[0];
		}

		// Generate initial A and U
		Damper.Node[S].A = (Damper.Node[S].CL1[0] + Damper.Node[S].CL2[0])/2;
		Damper.Node[S].T = pow(Damper.Node[S].A,2)*(Damper.EX ? pPpt->TREFe : pPpt->TREFi);
		k = pPpt->gammaAir(Damper.Node[S].T);
		Damper.Node[S].U = (Damper.Node[S].CL1[0] - Damper.Node[S].CL2[0])/(k-1);
		//Damper.Node[S].p_dash = Damper.pri[section]/pPpt->PREF;
		Damper.Node[S].p_dash = this->pBN[this->REAL_PIPE]->p_dash;
//if(ID==1) cout << "pri[section] = " << pri[section] << endl;
//cout << "Node[S=" << S << "].p_dash = " << Node[S].p_dash << endl;
		Damper.Node[S].rho = (Damper.Node[S].p_dash*pPpt->PREF*1e5)/(pPpt->R_air*Damper.Node[S].T);
		Damper.Node[S].rho_prev = Damper.Node[S].rho;
		Damper.Node[S].mdot = Damper.Node[S].rho*(Damper.Node[S].U*Damper.AREF)*Damper.Node[S].f;
		Damper.Node[S].CHOKED = false;
	}

	if(!pPpt->HOMENTROPIC){ // Requirements for non-homentropic flow: initial entropy levels
		for(S=0; S<Buffer.N; ++S){
			k = pPpt->gammaAir(Buffer.Node[S].T);
			Buffer.Node[S].AA[0] = Buffer.Node[S].A/pow(Buffer.Node[S].p_dash, (k-1)/(2*k));
			Buffer.Node[S].AA[1] = Buffer.Node[S].AA[0];
		}
		for(S=0; S<Damper.N; ++S){
			k = pPpt->gammaAir(Damper.Node[S].T);
			Damper.Node[S].AA[0] = Damper.Node[S].A/pow(Damper.Node[S].p_dash, (k-1)/(2*k));
			Damper.Node[S].AA[1] = Damper.Node[S].AA[0];
		}

		// Pathlines
		for(K=0; K<Buffer.num_pathlines; ++K){
			// Just need to position the pathlines initially and provide an entropy level value
			// PathLines() will calculate the remaining member variables
			Buffer.PathLine[K].ID = K;
			Buffer.PathLine[K].XK = K*(Buffer.XPIPE_EFF/double(Buffer.num_pathlines)) 
									+ 0.5*(Buffer.XPIPE_EFF/double(Buffer.num_pathlines)); // if XK has dimension (m)
			Buffer.PathLine[K].XK_old = Buffer.PathLine[K].XK;
			// Must select pathline AA from the surrounding nodes
			S=0;
			while(Buffer.Node[S].X<Buffer.PathLine[K].XK && S+1!=Buffer.N) ++S;

			if(Buffer.N>1){
				if(S==0) S=1; // Else S-1 undefined
				if(S==Buffer.N) S=Buffer.N-1;	// Else S undefined
				
				Buffer.PathLine[K].AA = ((Buffer.PathLine[K].XK - Buffer.Node[S-1].X)/(Buffer.Node[S].X - Buffer.Node[S-1].X))
					*(Buffer.Node[S].AA[0] - Buffer.Node[S-1].AA[0])
					+ Buffer.Node[S-1].AA[0];
			}
			else{ // For single node, joine pipes
				S=0;
				Buffer.PathLine[K].AA = Buffer.Node[S].AA[0];
			}
		}
		for(K=0; K<Damper.num_pathlines; ++K){
			// Just need to position the pathlines initially and provide an entropy level value
			// PathLines() will calculate the remaining member variables
			Damper.PathLine[K].ID = K;
			Damper.PathLine[K].XK = K*(Damper.XPIPE_EFF/double(Damper.num_pathlines)) 
									+ 0.5*(Damper.XPIPE_EFF/double(Damper.num_pathlines)); // if XK has dimension (m)
			Damper.PathLine[K].XK_old = Damper.PathLine[K].XK;
			// Must select pathline AA from the surrounding nodes
			S=0;
			while(Damper.Node[S].X<Damper.PathLine[K].XK && S+1!=Damper.N) ++S;

			if(Damper.N>1){
				if(S==0) S=1; // Else S-1 undefined
				if(S==Damper.N) S=Damper.N-1;	// Else S undefined
				
				Damper.PathLine[K].AA = ((Damper.PathLine[K].XK - Damper.Node[S-1].X)/(Damper.Node[S].X - Damper.Node[S-1].X))
					*(Damper.Node[S].AA[0] - Damper.Node[S-1].AA[0])
					+ Damper.Node[S-1].AA[0];
			}
			else{ // For single node, joine pipes
				S=0;
				Damper.PathLine[K].AA = Damper.Node[S].AA[0];
			}
		}
	}

	// W_alpha_beta schemes
	// ====================
	for(S=0; S<Buffer.N; ++S){
		// Allocate space for solution, flux, source vectors
		Buffer.Node[S].W = new double [3];
		Buffer.Node[S].F = new double [3];
		Buffer.Node[S].C = new double [3];

		Buffer.Node[S].W_pred = new double [3];
		Buffer.Node[S].F_pred = new double [3];
		Buffer.Node[S].C_pred = new double [3];

		Buffer.Node[S].S = new double [3];
		Buffer.Node[S].S_pred = new double [3];

		Buffer.Node[S].W_prev = new double [3];

		// Generate remaining primitives
//		Buffer.Node[S].rho = ((Buffer.Node[S].p_dash*pPpt->PREF)*1e5)/(pPpt->R_air*Buffer.Node[S].T);
//		Buffer.Node[S].rho_prev = Buffer.Node[S].rho;

		// Generate initial solution vector
		// --------------------------------
		k = pPpt->gammaAir(Buffer.Node[S].T);
		Buffer.Node[S].W[0] = Buffer.Node[S].rho*Buffer.Node[S].f;
		Buffer.Node[S].W[1] = Buffer.Node[S].rho*Buffer.Node[S].f*(Buffer.Node[S].U*Buffer.AREF);
		Buffer.Node[S].W[2] = Buffer.Node[S].rho * (((Buffer.Node[S].p_dash*pPpt->PREF)*1e5)
													/(Buffer.Node[S].rho*(k-1)) + 0.5*pow((Buffer.Node[S].U*Buffer.AREF), 2)) * Buffer.Node[S].f;
	}
	for(S=0; S<Damper.N; ++S)
	{
		// Allocate space for solution, flux, source vectors
		Damper.Node[S].W = new double [3];
		Damper.Node[S].F = new double [3];
		Damper.Node[S].C = new double [3];

		Damper.Node[S].W_pred = new double [3];
		Damper.Node[S].F_pred = new double [3];
		Damper.Node[S].C_pred = new double [3];

		Damper.Node[S].S = new double [3];
		Damper.Node[S].S_pred = new double [3];

		Damper.Node[S].W_prev = new double [3];

		// Generate remaining primitives
//		Damper.Node[S].rho = ((Damper.Node[S].p_dash*pPpt->PREF)*1e5)/(pPpt->R_air*Damper.Node[S].T);
//		Damper.Node[S].rho_prev = Damper.Node[S].rho;

		// Generate initial solution vector
		// --------------------------------
		k = pPpt->gammaAir(Damper.Node[S].T);
		Damper.Node[S].W[0] = Damper.Node[S].rho*Damper.Node[S].f;
		Damper.Node[S].W[1] = Damper.Node[S].rho*Damper.Node[S].f*(Damper.Node[S].U*Damper.AREF);
		Damper.Node[S].W[2] = Damper.Node[S].rho * (((Damper.Node[S].p_dash*pPpt->PREF)*1e5)
													/(Damper.Node[S].rho*(k-1)) + 0.5*pow((Damper.Node[S].U*Damper.AREF), 2)) * Damper.Node[S].f;
	}
}

/*
void CEndEnvironment::InitialiseAnechoic(CProperties* pPpt)
{
//	if(ANECHOIC)
	{
		p_anechoic = pBN[ONE_SIDE]->p_dash*pPpt->PREF;
		T_anechoic = pBN[ONE_SIDE]->T;

		// Zero velocity, so lambda_in = lambda_out = A
		lambda_in_an = sqrt(pPpt->gammaAir()*pPpt->R_air*T_anechoic)/pPipe[ONE_SIDE]->AREF;
		lambda_out_an = sqrt(pPpt->gammaAir()*pPpt->R_air*T_anechoic)/pPipe[ONE_SIDE]->AREF;
	}
	
} 
*/	

void CAnechoic::ReadInput(CProperties* pPpt, char *InputFile)
{
	int last_entry;
	CommonReadInput(InputFile, pPpt->NUM_PARAMS, this->labels, this->values, this->strings, last_entry);
	for(int r=0; r<last_entry+1; ++r)
	{
		// ====================================================================================================
		// Parameter file for Anechoic End [0]
		// ====================================================================================================
		
		// Pipe propagation method
		// --------------------------------------------------------------------------------
		if(strcmp(labels[r], "METHOD") == 0)
		{
			Buffer.METHOD = int(values[r]);
			Damper.METHOD = Buffer.METHOD;
		}

		// Pipe geometry
		// --------------------------------------------------------------------------------
		if(strcmp(labels[r], "length") == 0)
		{
			Buffer.length = values[r]/1000; // Conversion from mm to m
			Damper.length = 10*Buffer.length;
		}
		//if(strcmp(labels[r], "d_odd") == 0) d_odd = values[r]/1000;
		//if(strcmp(labels[r], "d_even") == 0) d_even = values[r]/1000;
		//if(strcmp(labels[r], "LINEAR") == 0) Buffer.LINEAR = bool(values[r]);
		//if(strcmp(labels[r], "LINEAR_F") == 0) LINEAR_F = bool(values[r]);
		//if(strcmp(labels[r], "x_int") == 0) x_int = values[r]/1000;
		//if(strcmp(labels[r], "d_int") == 0) d_int = values[r]/1000;
		if(strcmp(labels[r], "bend_angle") == 0)
		{
			Buffer.bend_angle = values[r];
			Damper.bend_angle = Buffer.bend_angle;
		}
		if(strcmp(labels[r], "discret") == 0)
		{
			Buffer.discret = values[r]/1000;
			Damper.discret = Buffer.discret;
		}
		if(strcmp(labels[r], "min_meshes") == 0)
		{
			Buffer.min_meshes = int(values[r]);
			Damper.min_meshes = Buffer.min_meshes;
		}
		if(strcmp(labels[r], "epsilon") == 0)
		{
			Buffer.epsilon = values[r]/1000;
			Damper.epsilon = Buffer.epsilon;
		}
		
		// Pipe initialisation
		// --------------------------------------------------------------------------------
		if(strcmp(labels[r], "Tw") == 0)
		{
			Buffer.Tw = values[r];
			Damper.Tw = Buffer.Tw;
		}
		if(strcmp(labels[r], "GLOBAL_EF") == 0)
		{
			Buffer.GLOBAL_EF = DoubleToBool(values[r]);
			Damper.GLOBAL_EF = Buffer.GLOBAL_EF;
		}
		if(strcmp(labels[r], "CFTRANS") == 0)
		{
			Buffer.CFTRANS = values[r];
			Damper.CFTRANS = Buffer.CFTRANS;
		}
		if(strcmp(labels[r], "HGTRANS") == 0)
		{
			Buffer.HGTRANS = values[r];
			Damper.HGTRANS = Buffer.HGTRANS;
		}
		if(strcmp(labels[r], "CPTRANS") == 0)
		{
			Buffer.CPTRANS = values[r];
			Damper.CPTRANS = Buffer.CPTRANS;
		}
	
		// Pipe measurements
		// --------------------------------------------------------------------------------
		if(strcmp(labels[r], "USE_DEF_FREQ") == 0)
		{
			Buffer.USE_DEF_FREQ = DoubleToBool(values[r]);
			Damper.USE_DEF_FREQ = Buffer.USE_DEF_FREQ;
		}
		if(strcmp(labels[r], "freq") == 0)
		{
			Buffer.freq = int(values[r]);	
			Damper.freq = Buffer.freq;
		}
		if(strcmp(labels[r], "print_from_time") == 0)
		{
			Buffer.print_from_time = values[r];
			Damper.print_from_time = Buffer.print_from_time;
		}		
		if(strcmp(labels[r], "ntappings") == 0)
		{
			Buffer.ntappings = int(values[r]);
			Damper.ntappings = Buffer.ntappings;
			Buffer.loc_measure = new double [Buffer.ntappings];
			Damper.loc_measure = new double [Damper.ntappings];
		}
		
		if(strcmp(labels[r], "loc_measure0") == 0 && Buffer.ntappings>0) Buffer.loc_measure[0] = values[r];
		if(strcmp(labels[r], "loc_measure1") == 0 && Buffer.ntappings>1) Buffer.loc_measure[1] = values[r];
		if(strcmp(labels[r], "loc_measure2") == 0 && Buffer.ntappings>2) Buffer.loc_measure[2] = values[r];
		if(strcmp(labels[r], "loc_measure3") == 0 && Buffer.ntappings>3) Buffer.loc_measure[3] = values[r];
		if(strcmp(labels[r], "loc_measure4") == 0 && Buffer.ntappings>4) Buffer.loc_measure[4] = values[r];
		if(strcmp(labels[r], "loc_measure5") == 0 && Buffer.ntappings>5) Buffer.loc_measure[5] = values[r];
		if(strcmp(labels[r], "loc_measure6") == 0 && Buffer.ntappings>6) Buffer.loc_measure[6] = values[r];
		if(strcmp(labels[r], "loc_measure7") == 0 && Buffer.ntappings>7) Buffer.loc_measure[7] = values[r];
		if(strcmp(labels[r], "loc_measure8") == 0 && Buffer.ntappings>8) Buffer.loc_measure[8] = values[r];
		if(strcmp(labels[r], "loc_measure9") == 0 && Buffer.ntappings>9) Buffer.loc_measure[9] = values[r];
		if(strcmp(labels[r], "loc_measure10") == 0 && Buffer.ntappings>10) Buffer.loc_measure[10] = values[r];

		if(strcmp(labels[r], "loc_measure0") == 0 && Damper.ntappings>0) Damper.loc_measure[0] = values[r];
		if(strcmp(labels[r], "loc_measure1") == 0 && Damper.ntappings>1) Damper.loc_measure[1] = values[r];
		if(strcmp(labels[r], "loc_measure2") == 0 && Damper.ntappings>2) Damper.loc_measure[2] = values[r];
		if(strcmp(labels[r], "loc_measure3") == 0 && Damper.ntappings>3) Damper.loc_measure[3] = values[r];
		if(strcmp(labels[r], "loc_measure4") == 0 && Damper.ntappings>4) Damper.loc_measure[4] = values[r];
		if(strcmp(labels[r], "loc_measure5") == 0 && Damper.ntappings>5) Damper.loc_measure[5] = values[r];
		if(strcmp(labels[r], "loc_measure6") == 0 && Damper.ntappings>6) Damper.loc_measure[6] = values[r];
		if(strcmp(labels[r], "loc_measure7") == 0 && Damper.ntappings>7) Damper.loc_measure[7] = values[r];
		if(strcmp(labels[r], "loc_measure8") == 0 && Damper.ntappings>8) Damper.loc_measure[8] = values[r];
		if(strcmp(labels[r], "loc_measure9") == 0 && Damper.ntappings>9) Damper.loc_measure[9] = values[r];
		if(strcmp(labels[r], "loc_measure10") == 0 && Damper.ntappings>10) Damper.loc_measure[10] = values[r];

		if(strcmp(labels[r], "max_pts") == 0)
		{
			Buffer.max_pts = int(values[r]);
			Damper.max_pts = Buffer.max_pts;
		}

		if (strcmp(labels[r], "DIAMETER") == 0) {
			Buffer.DIAMETER = DoubleToBool(values[r]);
			Damper.DIAMETER = Buffer.DIAMETER;
		}
		if (strcmp(labels[r], "AREA") == 0) {
			Buffer.AREA = DoubleToBool(values[r]);
			Damper.AREA = Buffer.AREA;
		}
		if(strcmp(labels[r], "STATIC_PRESSURE") == 0) {
			Buffer.STATIC_PRESSURE = DoubleToBool(values[r]);
			Damper.STATIC_PRESSURE = Buffer.STATIC_PRESSURE;
		}
		if(strcmp(labels[r], "TEMPERATURE") == 0) {
			Buffer.TEMPERATURE = DoubleToBool(values[r]);
			Damper.TEMPERATURE = Buffer.TEMPERATURE;
		}
		if(strcmp(labels[r], "VELOCITY") == 0) {
			Buffer.VELOCITY = DoubleToBool(values[r]);
			Damper.VELOCITY = Buffer.VELOCITY;
		}
		if(strcmp(labels[r], "MACH_NUMBER") == 0) {
			Buffer.MACH_NUMBER = DoubleToBool(values[r]);
			Damper.MACH_NUMBER = Buffer.MACH_NUMBER;
		}
		if(strcmp(labels[r], "MASS_FLOW_RATE") == 0) {
			Buffer.MASS_FLOW_RATE = DoubleToBool(values[r]);
			Damper.MASS_FLOW_RATE = Buffer.MASS_FLOW_RATE;
		}
		if(strcmp(labels[r], "REYNOLDS_NO") == 0) {
			Buffer.REYNOLDS_NO = DoubleToBool(values[r]);
			Damper.REYNOLDS_NO = Buffer.REYNOLDS_NO;
		}

		// Screen output
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "SHOW_DATA") == 0)
		{
			Buffer.SHOW_DATA = DoubleToBool(values[r]);
			Damper.SHOW_DATA = Buffer.SHOW_DATA;
		}

		// ====================================================================================================
		// End of file
		// ====================================================================================================
	}

	// Set some derived parameters
	// ---------------------------
	if(Buffer.METHOD==0) Buffer.METHOD=pPpt->DEF_METHOD;	// Use the default method for this pipe
	if(Damper.METHOD==0) Damper.METHOD=pPpt->DEF_METHOD;	// Use the default method for this pipe
	if(Buffer.GLOBAL_EF) // Apply global enhancement factors
	{
		Buffer.CFTRANS = pPpt->CFTRANS;
		Buffer.HGTRANS = pPpt->HGTRANS;
		Buffer.CPTRANS = pPpt->CPTRANS;
	}
	if(Damper.GLOBAL_EF) // Apply global enhancement factors
	{
		Damper.CFTRANS = pPpt->CFTRANS;
		Damper.HGTRANS = pPpt->HGTRANS;
		Damper.CPTRANS = pPpt->CPTRANS;
	}
	if(Buffer.USE_DEF_FREQ) Buffer.freq = pPpt->freq;		// Use the default sampling rate for this pipe
	if(Damper.USE_DEF_FREQ) Damper.freq = pPpt->freq;		// Use the default sampling rate for this pipe

	Buffer.num_props_measured = int(Buffer.DIAMETER) + int(Buffer.AREA) + int(Buffer.STATIC_PRESSURE) + int(Buffer.TEMPERATURE) + int(Buffer.VELOCITY) + int(Buffer.MACH_NUMBER) + int(Buffer.MASS_FLOW_RATE) + int(Buffer.REYNOLDS_NO);
	Damper.num_props_measured = int(Damper.DIAMETER) + int(Damper.AREA) + int(Damper.STATIC_PRESSURE) + int(Damper.TEMPERATURE) + int(Damper.VELOCITY) + int(Damper.MACH_NUMBER) + int(Damper.MASS_FLOW_RATE) + int(Damper.REYNOLDS_NO);
}

void CAnechoic::ListProperties(CProperties* pPpt)
{
	pPpt->Out(Underline(Identify(), "=", "\t"));
	pPpt->Out("\n");
	Buffer.ListProperties(pPpt);
	pPpt->Out("\n");
	Damper.ListProperties(pPpt);

//	if(ANECHOIC)
//	{
//		pPpt->Out("\tAnechoic termination, ANECHOIC\t\t\t=\t"); pPpt->Out(TrueOrFalse(ANECHOIC)); pPpt->Out("\n");
//		pPpt->Out("\tAnechoic ref. pressure, p_anechoic\t\t=\t"); pPpt->Out(p_anechoic); pPpt->Out(" bar\n");
//		pPpt->Out("\tAnechoic ref. temperature, T_anechoic\t\t=\t"); pPpt->Out(T_anechoic); pPpt->Out(" K\n");
//	}
//	else
/*
	{
		if(VAR_P0)
		{
			pPpt->Out("\tVariable stagnation conditions");
			if(P0_COS)
			{
				pPpt->Out(" - sinusoidal variation:\n");
				pPpt->Out("\tSinusoidal frequency, FREQ\t\t\t=\t"); pPpt->Out(FREQ); pPpt->Out(" Hz\n");
				pPpt->Out("\tPulse duty cycle, PHI\t\t\t\t=\t"); pPpt->Out(PHI); pPpt->Out("\n");
				pPpt->Out("\tLowest stagnation pressure, P0_LOW\t\t=\t"); pPpt->Out(P0_LOW); pPpt->Out(" bar\n");
				pPpt->Out("\tHighest stagnation pressure, P0_HIGH\t\t=\t"); pPpt->Out(P0_HIGH); pPpt->Out(" bar\n");
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
		if(VAR_PHI)
		{
			pPpt->Out("\n");
			pPpt->Out("\tNozzle variable operation");
			if(!SQUARE)
			{
				pPpt->Out(" - linear ramp:\n");
				pPpt->Out("\tStart area ratio, PHI_START\t\t\t=\t"); pPpt->Out(PHI_START); pPpt->Out("\n");
				pPpt->Out("\tEnd area ratio, PHI_END\t\t\t=\t"); pPpt->Out(PHI_END); pPpt->Out("\n");
				pPpt->Out("\tVariation starts at, START_TIME\t\t=\t"); pPpt->Out(START_TIME); pPpt->Out(" s\n");
				pPpt->Out("\tTime taken, VAR_TIME\t\t\t\t=\t"); pPpt->Out(VAR_TIME); pPpt->Out(" s\n");
			}
			else
			{
				pPpt->Out(" - step change:\n");
				pPpt->Out("\tStart area ratio, PHI_START\t\t\t=\t"); pPpt->Out(PHI_START); pPpt->Out("\n");
				pPpt->Out("\tStep area ratio, PHI_END\t\t\t=\t"); pPpt->Out(PHI_END); pPpt->Out("\n");
				pPpt->Out("\tStep change starts at, START_TIME\t\t=\t"); pPpt->Out(START_TIME); pPpt->Out(" s\n");
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
	pPpt->Out("\n");

	//if(NOTIFY) pPpt->Out("\tNozzle flow direction changes will be notified\n");
	//else pPpt->Out("\tNozzle flow direction changes will not be notified\n");
	//if(PRINT_MOVIE_FILE) pPpt->Out("\tPrinting movie file\n");
	pPpt->Out("\n");
	pPpt->Out("\n");
*/
}

void CAnechoic::RunBoundary(CProperties* pPpt, int timestep, double time)
{
	Anechoic(pPpt, timestep, time);
}

void CAnechoic::Anechoic(CProperties* pPpt, int timestep, double time)
//--------------------------------------------------//
// Non-homentropic anechoic termination				//
// ------------------------------------				//
//													//
//--------------------------------------------------//
{
	double k;
	double lambda_in_n_real, lambda_out_n_real, AA_n_real;
	double lambda_in_n_buffer, lambda_out_n_buffer, AA_n_buffer;
	double *lambda_in_c, *lambda_out_c, *AA_c;
	int *pipe_flow;
	lambda_in_c = new double [NPIPES]; lambda_out_c = new double [NPIPES]; AA_c = new double [NPIPES];	pipe_flow = new int [NPIPES];
	bool* CHOKED; CHOKED = new bool [NPIPES]; 
	CHOKED[ONE_SIDE] = false; // Reset here for all cases

	// Real pipe values
	lambda_in_n_real = (*(pCLIN[ONE_SIDE]))[R+1];
	lambda_out_n_real = (*(pCLOUT[ONE_SIDE]))[R+1];
	AA_n_real = pBN[ONE_SIDE]->AA[R+1];
	
	// ODD end of buffer always joins the anechoic boundary with the real pipe
	lambda_in_n_buffer = Buffer.Node[0].CL2[R+1];
	lambda_out_n_buffer = Buffer.Node[0].CL1[R+1];
	AA_n_buffer = Buffer.Node[0].AA[R+1];

	// Test flow direction
	k = pPpt->gammaAir(pBN[ONE_SIDE]->T);
	double PIp_real = pow(lambda_in_n_real/AA_n_real, (2*k)/(k-1));
	k = pPpt->gammaAir(Buffer.Node[0].T);
	double PIp_buffer = pow(lambda_in_n_buffer/AA_n_buffer, (2*k)/(k-1));

/*
	if(fabs(PIp_real - PIp_buffer) < 1e-12) 
	{
		// NOFLOW
		pipe_flow[ONE_SIDE] = NOFLOW;
		Buffer.odd_end_flow = NOFLOW;

		lambda_in_c[ONE_SIDE] = lambda_in_n_real;		// No change
		lambda_out_c[ONE_SIDE] = lambda_in_n_buffer;	// Passes through boundary unchanged
		AA_c[ONE_SIDE] = AA_n_real;						// No change

		Buffer.Node[0].CL2[R+1] = lambda_in_n_buffer;	// No change	
		Buffer.Node[0].CL1[R+1] = lambda_in_n_real;		// Passes through boundary unchanged
		Buffer.Node[0].AA[R+1] = AA_n_buffer;			// No change		
	}
	else
	{
		if(PIp_real < PIp_buffer)	
		{
			// INFLOW into real pipe
			pipe_flow[ONE_SIDE] = INFLOW;
			Buffer.odd_end_flow = OUTFLOW;

			lambda_out_c[ONE_SIDE] = lambda_in_n_buffer;		// Passes through boundary unchanged
			AA_c[ONE_SIDE] = AA_n_buffer;						// Since inflow from buffer
			lambda_in_c[ONE_SIDE] = 2*lambda_in_n_real*(AA_c[ONE_SIDE]/(AA_c[ONE_SIDE] + AA_n_real)) 
									+ lambda_out_c[ONE_SIDE]*((AA_c[ONE_SIDE] - AA_n_real)/(AA_c[ONE_SIDE] + AA_n_real));
																// Corrected due to inflow
			
			Buffer.Node[0].CL2[R+1] = lambda_in_n_buffer;		// No change	
			Buffer.Node[0].CL1[R+1] = lambda_in_c[ONE_SIDE];	// The corrected version
			Buffer.Node[0].AA[R+1] = AA_n_buffer;				// No change


			//	lambda_in_c = 2*lambda_in_n*(AA_c/(AA_c + AA_n)) + lambda_out_c*((AA_c - AA_n)/(AA_c + AA_n));
		}
		else 
		{
			// OUTFLOW from real pipe
			pipe_flow[ONE_SIDE] = OUTFLOW;
			Buffer.odd_end_flow = INFLOW;

			Buffer.Node[0].CL1[R+1] = lambda_in_n_real;			// Passes through boundary unchanged
			Buffer.Node[0].AA[R+1] = AA_n_real;					// Since inflow into buffer from real
			Buffer.Node[0].CL2[R+1] = 2*lambda_in_n_buffer*(Buffer.Node[0].AA[R+1]/(Buffer.Node[0].AA[R+1] + AA_n_buffer)) 
									+ Buffer.Node[0].CL1[R+1]*((Buffer.Node[0].AA[R+1] - AA_n_buffer)/(Buffer.Node[0].AA[R+1] + AA_n_buffer));
																// Corrected due to inflow into buffer

			lambda_in_c[ONE_SIDE] = lambda_in_n_real;			// No change
			lambda_out_c[ONE_SIDE] = Buffer.Node[0].CL2[R+1];	// The corrected version
			AA_c[ONE_SIDE] = AA_n_real;							// No change
		}
	}

	// Buffer-Damper exchange
	double PIp_buffer_even = pow(Buffer.Node[Buffer.N-1].CL1[R+1]/Buffer.Node[Buffer.N-1].AA[R+1], (2*pPpt->gammaAir())/(pPpt->gammaAir()-1));
	double PIp_damper = pow(Damper.Node[0].CL2[R+1]/Buffer.Node[0].AA[R+1], (2*pPpt->gammaAir())/(pPpt->gammaAir()-1));

	if(fabs(PIp_buffer_even - PIp_damper) < 1e-12) 
	{
		// NOFLOW
		Buffer.even_end_flow = NOFLOW;
		Damper.odd_end_flow = NOFLOW;

		Buffer.Node[Buffer.N-1].CL1[R+1] = Buffer.Node[Buffer.N-1].CL1[R+1];	// No change	
		Buffer.Node[Buffer.N-1].CL2[R+1] = Damper.Node[0].CL2[R+1];				// Passes through boundary unchanged
		Buffer.Node[Buffer.N-1].AA[R+1] = Buffer.Node[Buffer.N-1].AA[R+1];		// No change		

		Damper.Node[0].CL2[R+1] = Damper.Node[0].CL2[R+1];						// No change	
		Damper.Node[0].CL1[R+1] = Buffer.Node[Buffer.N-1].CL1[R+1];				// Passes through boundary unchanged
		Damper.Node[0].AA[R+1] = Damper.Node[0].AA[R+1];						// No change		
	}
	else
	{
		if(PIp_buffer_even < PIp_damper)	
		{
			// INFLOW into buffer pipe from damper
			Buffer.even_end_flow = INFLOW;
			Damper.odd_end_flow = OUTFLOW;

			Buffer.Node[Buffer.N-1].CL2[R+1] = Damper.Node[0].CL2[R+1];;		// Passes through boundary unchanged
			Buffer.Node[Buffer.N-1].CL1[R+1] = 2*Buffer.Node[Buffer.N-1].CL1[R+1]*(Damper.Node[0].AA[R+1]/(Damper.Node[0].AA[R+1] + Buffer.Node[Buffer.N-1].AA[R+1])) 
									+ Buffer.Node[Buffer.N-1].CL2[R+1]*((Damper.Node[0].AA[R+1] - Buffer.Node[Buffer.N-1].AA[R+1])/(Damper.Node[0].AA[R+1] + Buffer.Node[Buffer.N-1].AA[R+1]));
																// Corrected due to inflow
			Buffer.Node[Buffer.N-1].AA[R+1] = Damper.Node[0].AA[R+1];			// Since inflow from damper


			Damper.Node[0].CL2[R+1] = Damper.Node[0].CL2[R+1];					// No change	
			Damper.Node[0].CL1[R+1] = Buffer.Node[Buffer.N-1].CL1[R+1];			// The corrected version
			Damper.Node[0].AA[R+1] = Damper.Node[0].AA[R+1];					// No change

			//	lambda_in_c = 2*lambda_in_n*(AA_c/(AA_c + AA_n)) + lambda_out_c*((AA_c - AA_n)/(AA_c + AA_n));
		}
		else 
		{
			// OUTFLOW from buffer pipe into damper
			Buffer.even_end_flow = OUTFLOW;
			Damper.odd_end_flow = INFLOW;

			Damper.Node[0].CL1[R+1] = Buffer.Node[Buffer.N-1].CL1[R+1];			// Passes through boundary unchanged
			
			Damper.Node[0].CL2[R+1] = 2*Damper.Node[0].CL2[R+1]*(Buffer.Node[Buffer.N-1].AA[R+1]/(Buffer.Node[Buffer.N-1].AA[R+1] + Damper.Node[0].AA[R+1])) 
									+ Damper.Node[0].CL1[R+1]*((Buffer.Node[Buffer.N-1].AA[R+1] - Damper.Node[0].AA[R+1])/(Buffer.Node[Buffer.N-1].AA[R+1] + Damper.Node[0].AA[R+1]));
																// Corrected due to inflow into buffer
			Damper.Node[0].AA[R+1] = Buffer.Node[Buffer.N-1].AA[R+1];			// Since inflow into buffer from real
			

			Buffer.Node[Buffer.N-1].CL1[R+1] = Buffer.Node[Buffer.N-1].CL1[R+1];	// No change	
			Buffer.Node[Buffer.N-1].CL2[R+1] = Damper.Node[0].CL2[R+1];				// Passes through boundary unchanged
			Buffer.Node[Buffer.N-1].AA[R+1] = Buffer.Node[Buffer.N-1].AA[R+1];		// No change	
		}
	}
//*/

///*
	// Buffer pipe only (no damper):
	if(fabs(PIp_real - PIp_buffer) < 1e-12) 
	{
		// NOFLOW
		pipe_flow[ONE_SIDE] = NOFLOW;
		Buffer.odd_end_flow = NOFLOW;

		lambda_in_c[ONE_SIDE] = lambda_in_n_real;		// No change
		lambda_out_c[ONE_SIDE] = lambda_in_n_buffer;	// Passes through boundary unchanged
		AA_c[ONE_SIDE] = AA_n_real;						// No change

		Buffer.Node[0].CL2[R+1] = lambda_in_n_buffer;	// No change	
		Buffer.Node[0].CL1[R+1] = lambda_in_n_real;		// Passes through boundary unchanged
		Buffer.Node[0].AA[R+1] = AA_n_buffer;			// No change
	}
	else
	{
		if(PIp_real < PIp_buffer)	
		{
			// INFLOW into real pipe
			pipe_flow[ONE_SIDE] = INFLOW;
			Buffer.odd_end_flow = OUTFLOW;

			lambda_out_c[ONE_SIDE] = lambda_in_n_buffer;		// Passes through boundary unchanged
			AA_c[ONE_SIDE] = AA_n_buffer;						// Since inflow from buffer
			lambda_in_c[ONE_SIDE] = 2*lambda_in_n_real*(AA_c[ONE_SIDE]/(AA_c[ONE_SIDE] + AA_n_real)) 
									+ lambda_out_c[ONE_SIDE]*((AA_c[ONE_SIDE] - AA_n_real)/(AA_c[ONE_SIDE] + AA_n_real));
																// Corrected due to inflow
			
			Buffer.Node[0].CL2[R+1] = lambda_in_n_buffer;		// No change	
			Buffer.Node[0].CL1[R+1] = lambda_in_c[ONE_SIDE];	// The corrected version
			Buffer.Node[0].AA[R+1] = AA_n_buffer;				// No change


			//	lambda_in_c = 2*lambda_in_n*(AA_c/(AA_c + AA_n)) + lambda_out_c*((AA_c - AA_n)/(AA_c + AA_n));
		}
		else 
		{
			// OUTFLOW from real pipe
			pipe_flow[ONE_SIDE] = OUTFLOW;
			Buffer.odd_end_flow = INFLOW;

			Buffer.Node[0].CL1[R+1] = lambda_in_n_real;			// Passes through boundary unchanged
			Buffer.Node[0].AA[R+1] = AA_n_real;					// Since inflow into buffer from real
			Buffer.Node[0].CL2[R+1] = 2*lambda_in_n_buffer*(Buffer.Node[0].AA[R+1]/(Buffer.Node[0].AA[R+1] + AA_n_buffer)) 
									+ Buffer.Node[0].CL1[R+1]*((Buffer.Node[0].AA[R+1] - AA_n_buffer)/(Buffer.Node[0].AA[R+1] + AA_n_buffer));
																// Corrected due to inflow into buffer

			lambda_in_c[ONE_SIDE] = lambda_in_n_real;			// No change
			lambda_out_c[ONE_SIDE] = Buffer.Node[0].CL2[R+1];	// The corrected version
			AA_c[ONE_SIDE] = AA_n_real;							// No change
		}
	}
//*/
/*
	// Simple no update without buffer pipe:
	// OUTFLOW from real pipe
	pipe_flow[ONE_SIDE] = OUTFLOW;

	lambda_in_c[ONE_SIDE] = lambda_in_n_real;			// No change
	lambda_out_c[ONE_SIDE] = lambda_out_n_real;			// No update
	AA_c[ONE_SIDE] = AA_n_real;							// No change
//*/

	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out_c, AA_c, pipe_flow, CHOKED);


	// Other end of buffer pipe (considered to be at stagnant conditions, i.e., U=0, but shouldn't reflect, so not closed end)
	// No need to change lambda_out from its initial (stagnant conditions)
}

void CAnechoic::AnechoicOld(CProperties* pPpt, int timestep, double time)
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

void CAnechoic::PrintToScreen(CProperties* pPpt)
{
	pPpt->Out(Underline(Identify(), "="));

	Buffer.PrintToScreen(pPpt);

	Damper.PrintToScreen(pPpt);

	//pPpt->Out("\n");
}