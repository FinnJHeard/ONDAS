// Junction.cpp: implementation of the CJunction class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Junction.h"
#include "Node.h"
#include "Pipe.h"
#include "Properties.h"
#include "Tools.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CJunction::CJunction(void)
{


}

CJunction::~CJunction()
{
	// Though pCLIN, pCLOUT are ***, memory was only allocated to the top 2 levels, so no
	// need to delete at all levels
	for(int p=0; p<NPIPES; ++p)
	{
		delete [] pPipe[p];
		delete [] pBN[p];
/*		
		for(int i=0; i<2; ++i)
		{
			delete [] pCLIN[p][i];
			delete [] pCLOUT[p][i];
		}
*/
		delete [] pCLIN[p];
		delete [] pCLOUT[p];
		delete [] pend_flow[p];
		delete [] pPathLine[p];
	}
	delete [] pPipe;
	delete [] pBN;
	delete [] pCLIN;
	delete [] pCLOUT;
	delete [] pend_flow;
	delete [] pPathLine;

	delete [] end;
	delete [] Fb;
	delete [] lambda_in_star;
	delete [] lambda_in_star_n;
	delete [] AA;		// Working entropy level
	delete [] AA_n;	// Uncorrected entropy levels
	delete [] A_star;		// Starred N.D. speed of sound
	delete [] U_star;
	delete [] pipe_flow;
	delete [] lambda_out;
	delete [] pipe_flow_old;
	delete [] del_A_star;
	delete [] del_A_star_old;


	// BWPLJ
	// -----
	delete [] del_star;


	delete [] m_dot;
	delete [] W;	// Mass flow ratio, other pipe 'a' and 'b'
	delete [] L;	// Loss coefficients
	delete [] G1;
	delete [] G2;
	delete [] coeff;	// Which loss coefficient to lookup [La] [Lb]
	delete [] x_star;
	delete [] x_star_old;
	delete [] Kq;
}

// Copy constructor
CJunction::CJunction(const CJunction& inJunc)
{
	int p;
///*
	flow_type_winterbone = inJunc.flow_type_winterbone;
	flow_type_winterbone_old = inJunc.flow_type_winterbone_old;
	flow_type_winterbone_orig = inJunc.flow_type_winterbone_orig;
	flow_type_benson = inJunc.flow_type_benson;

	switch_counter = inJunc.switch_counter;
//*/
//	rollback = inJunc.rollback;
	rolledback = inJunc.rolledback;
	stoprollback = inJunc.stoprollback;

/*
	// Do not restore these!
	ALTERNATING = inJunc.ALTERNATING;
	FIX_FLOW_TYPE = inJunc.FIX_FLOW_TYPE;
	fix_counter = inJunc.fix_counter;
	OPPOSITE = inJunc.OPPOSITE; FIRST_OPPOSITE = inJunc.FIRST_OPPOSITE;
	OTHER = inJunc.OTHER;
	other_counter = inJunc.other_counter;
	ZERO_FLOW = inJunc.ZERO_FLOW;

	flow_type_actual = inJunc.flow_type_actual;

  	TESTED_ACTUAL = inJunc.TESTED_ACTUAL;
	TESTED_ORIGINAL = inJunc.TESTED_ORIGINAL;
*/ 

	normal_counter_total = inJunc.normal_counter_total;
	fix_counter_total = inJunc.fix_counter_total;
	opposite_counter_total = inJunc.opposite_counter_total;
	other_counter_total = inJunc.other_counter_total;
	zero_counter_total = inJunc.zero_counter_total;

	ID = inJunc.ID;
	EX = inJunc.EX;
	NPIPES = inJunc.NPIPES;
//	AREF = inJunc.AREF;
//	angle_deg = inJunc.angle_deg;		
//	angle_rad = inJunc.angle_rad;
	psiT = inJunc.psiT;

	pipe1 = inJunc.pipe1;
	pipe2 = inJunc.pipe2;
	pipe3 = inJunc.pipe3;

	bpipe1 = inJunc.bpipe1;
	bpipe2 = inJunc.bpipe2;
	bpipe1 = inJunc.bpipe3;

	a = inJunc.a;
	b = inJunc.b;

//	tol_cont = inJunc.tol_cont;
//	loop_limit_cont = inJunc.loop_limit_cont;
//	tol_main = inJunc.tol_main;

	loop_limit_main = inJunc.loop_limit_main;
//	loop_limit_same = inJunc.loop_limit_same;

	com = inJunc.com;
	i2 = inJunc.i2;
	i3 = inJunc.i3;

	C1 = inJunc.C1; C2 = inJunc.C2; C3 = inJunc.C3; C4 = inJunc.C4; C5 = inJunc.C5; C6 = inJunc.C6;
	
	Mn = inJunc.Mn;
//	type = inJunc.type;

	first_loop = inJunc.first_loop;
	first_loop_ever = inJunc.first_loop_ever;
	hcpj_first_loop = inJunc.hcpj_first_loop;
	nhplj_first_loop = inJunc.nhplj_first_loop;
	FT = inJunc.FT;
	
	pPipe = new CPipe* [NPIPES];
	pBN = new CNode* [NPIPES];
	pCLIN = new double** [NPIPES];
	pCLOUT = new double** [NPIPES];
	pend_flow = new int* [NPIPES];
	pPathLine = new CPathLine* [NPIPES];
	end = new int [NPIPES];

	Fb = new double [NPIPES];		
	G1 = new double [NPIPES];
	G2 = new double [NPIPES];

	lambda_in_star = new double [NPIPES];
	lambda_in_star_n = new double [NPIPES];
	AA = new double [NPIPES];
	AA_n = new double [NPIPES];

	pipe_flow = new int [NPIPES];
	pipe_flow_old = new int [NPIPES];

	lambda_out = new double [NPIPES];
	A_star = new double [NPIPES];
	U_star = new double [NPIPES];

	del_A_star = new double [NPIPES];
	del_A_star_old = new double [NPIPES];

	del_star = new double [NPIPES];

	m_dot = new double [NPIPES];
	x_star = new double [NPIPES];
	x_star_old = new double [NPIPES];
	Kq = new double [NPIPES];

	for(p=0; p<NPIPES; ++p)
	{
		pPipe[p] = inJunc.pPipe[p]; // Yes, copy pointers, this is what is desired
		pBN[p] = inJunc.pBN[p]; // Yes, copy pointers, this is what is desired
		pCLIN[p] = inJunc.pCLIN[p]; // Yes, copy pointers, this is what is desired
		pCLOUT[p] = inJunc.pCLOUT[p]; // Yes, copy pointers, this is what is desired
		pend_flow[p] = inJunc.pend_flow[p]; // Yes, copy pointers, this is what is desired
		pPathLine[p] = inJunc.pPathLine[p]; // Yes, copy pointers, this is what is desired
		end[p] = inJunc.end[p];

		Fb[p] = inJunc.Fb[p];
		G1[p] = inJunc.G1[p];
		G2[p] = inJunc.G2[p];

		lambda_in_star[p] = inJunc.lambda_in_star[p];
		lambda_in_star_n[p] = inJunc.lambda_in_star_n[p];
		AA[p] = inJunc.AA[p];
		AA_n[p] = inJunc.AA_n[p];

		pipe_flow[p] = inJunc.pipe_flow[p];
		pipe_flow_old[p] = inJunc.pipe_flow_old[p];

		lambda_out[p] = inJunc.lambda_out[p];
		A_star[p] = inJunc.A_star[p];
		U_star[p] = inJunc.U_star[p];

		del_A_star[p] = inJunc.del_A_star[p];
		del_A_star_old[p] = inJunc.del_A_star_old[p];

		del_star[p] = inJunc.del_A_star[p];

		m_dot[p] = inJunc.m_dot[p];
		x_star[p] = inJunc.x_star[p];
		x_star_old[p] = inJunc.x_star_old[p];
		Kq[p] = inJunc.Kq[p];
	}

	W = new double [2];
	L = new double [2];
	coeff = new int [2];	
	for(p=0; p<2; ++p)
	{
		W[p] = inJunc.W[p];
		L[p] = inJunc.L[p];
		coeff[p] = inJunc.coeff[p];
	}
}

CJunction& CJunction::operator=(const CJunction& inJunc)
{
	int p;
	if(this != &inJunc)
	{
///*
		flow_type_winterbone = inJunc.flow_type_winterbone;
		flow_type_winterbone_old = inJunc.flow_type_winterbone_old;
		flow_type_winterbone_orig = inJunc.flow_type_winterbone_orig;
		flow_type_benson = inJunc.flow_type_benson;

		switch_counter = inJunc.switch_counter;
//*/

//		rollback = inJunc.rollback;
		rolledback = inJunc.rolledback;
		stoprollback = inJunc.stoprollback;
/*
		// Do not restore these!
		ALTERNATING = inJunc.ALTERNATING;
		FIX_FLOW_TYPE = inJunc.FIX_FLOW_TYPE;
		fix_counter = inJunc.fix_counter;
		OPPOSITE = inJunc.OPPOSITE; FIRST_OPPOSITE = inJunc.FIRST_OPPOSITE;
		OTHER = inJunc.OTHER;
		other_counter = inJunc.other_counter;
		ZERO_FLOW = inJunc.ZERO_FLOW;

		flow_type_actual = inJunc.flow_type_actual;

		TESTED_ACTUAL = inJunc.TESTED_ACTUAL;
		TESTED_ORIGINAL = inJunc.TESTED_ORIGINAL;
*/
		normal_counter_total = inJunc.normal_counter_total;
		fix_counter_total = inJunc.fix_counter_total;
		opposite_counter_total = inJunc.opposite_counter_total;
		other_counter_total = inJunc.other_counter_total;
		zero_counter_total = inJunc.zero_counter_total;

		ID = inJunc.ID;
		EX = inJunc.EX;
		NPIPES = inJunc.NPIPES;
//		AREF = inJunc.AREF;
//		angle_deg = inJunc.angle_deg;		
//		angle_rad = inJunc.angle_rad;
		psiT = inJunc.psiT;
		
		pipe1 = inJunc.pipe1;
		pipe2 = inJunc.pipe2;
		pipe3 = inJunc.pipe3;
	
		bpipe1 = inJunc.bpipe1;
		bpipe2 = inJunc.bpipe2;
		bpipe3 = inJunc.bpipe3;
	
		a = inJunc.a;
		b = inJunc.b;
	
//		tol_cont = inJunc.tol_cont;
//		loop_limit_cont = inJunc.loop_limit_cont;
//		tol_main = inJunc.tol_main;
	
		loop_limit_main = inJunc.loop_limit_main;
//		loop_limit_same = inJunc.loop_limit_same;
	
		com = inJunc.com;
		i2 = inJunc.i2;
		i3 = inJunc.i3;
	
		C1 = inJunc.C1; C2 = inJunc.C2; C3 = inJunc.C3; C4 = inJunc.C4; C5 = inJunc.C5; C6 = inJunc.C6;
		
		Mn = inJunc.Mn;
//		type = inJunc.type;
	
		first_loop = inJunc.first_loop;
		first_loop_ever = inJunc.first_loop_ever;
		hcpj_first_loop = inJunc.hcpj_first_loop;
		nhplj_first_loop = inJunc.nhplj_first_loop;
		FT = inJunc.FT;
///*		
		pPipe = new CPipe* [NPIPES];
		pBN = new CNode* [NPIPES];
		pCLIN = new double** [NPIPES];
		pCLOUT = new double** [NPIPES];
		pend_flow = new int* [NPIPES];
		pPathLine = new CPathLine* [NPIPES];
		end = new int [NPIPES];
	
		Fb = new double [NPIPES];		
		G1 = new double [NPIPES];
		G2 = new double [NPIPES];

		lambda_in_star = new double [NPIPES];
		lambda_in_star_n = new double [NPIPES];
		AA = new double [NPIPES];
		AA_n = new double [NPIPES];

		pipe_flow = new int [NPIPES];
		pipe_flow_old = new int [NPIPES];

		lambda_out = new double [NPIPES];
		A_star = new double [NPIPES];
		U_star = new double [NPIPES];

		del_A_star = new double [NPIPES];
		del_A_star_old = new double [NPIPES];

		C = new double [NPIPES];
		del_star = new double [NPIPES];

		m_dot = new double [NPIPES];
		x_star = new double [NPIPES];
		x_star_old = new double [NPIPES];
		Kq = new double [NPIPES];
//*/	
		for(p=0; p<NPIPES; ++p)
		{
			pPipe[p] = inJunc.pPipe[p];
			end[p] = inJunc.end[p];
			pBN[p] = inJunc.pBN[p];
			pCLIN[p] = inJunc.pCLIN[p];
			pCLOUT[p] = inJunc.pCLOUT[p];
			pend_flow[p] = inJunc.pend_flow[p];
			pPathLine[p] = inJunc.pPathLine[p];
	
			Fb[p] = inJunc.Fb[p];
			G1[p] = inJunc.G1[p];
			G2[p] = inJunc.G2[p];

			lambda_in_star[p] = inJunc.lambda_in_star[p];
			lambda_in_star_n[p] = inJunc.lambda_in_star_n[p];
			AA[p] = inJunc.AA[p];
			AA_n[p] = inJunc.AA_n[p];

			pipe_flow[p] = inJunc.pipe_flow[p];
			pipe_flow_old[p] = inJunc.pipe_flow_old[p];

			lambda_out[p] = inJunc.lambda_out[p];
			A_star[p] = inJunc.A_star[p];
			U_star[p] = inJunc.U_star[p];

			del_A_star[p] = inJunc.del_A_star[p];
			del_A_star_old[p] = inJunc.del_A_star_old[p];

			C[p] = inJunc.C[p];
			del_star[p] = inJunc.del_star[p];

			m_dot[p] = inJunc.m_dot[p];
			x_star[p] = inJunc.x_star[p];
			x_star_old[p] = inJunc.x_star_old[p];
			Kq[p] = inJunc.Kq[p];

		}
///*
		W = new double [2];
		L = new double [2];
		coeff = new int [2];
//*/	
		for(p=0; p<2; ++p)
		{
			W[p] = inJunc.W[p];
			L[p] = inJunc.L[p];
			coeff[p] = inJunc.coeff[p];
		}
	}
	return *this;
}

void CJunction::Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, /*int** &rJASSYS,*/ int** &rJPIPES, int** &rJPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, std::string param_dir, std::string res_dir, int assyid, string parent_assy_res_dir, string calling_object_str)
{
	if (pPpt->SHOW_calls) { pPpt->Out("CJunction.Initialise ("); pPpt->Out(calling_object_str); pPpt->Out(")\n"); }

	int p;
	NPIPES = npipes;
	
	// Enumerate array identifers (all notations)
	pipe1 = 0; // Pipe '1'=[0]
	pipe2 = 1; // Pipe '2'=[1]
	pipe3 = 2; // Pipe '3'=[2]

	bpipe1 = 0; // Benson '1' = my [0]
	bpipe2 = 2; // Benson '2' = my [2]
	bpipe3 = 1; // Benson '3' = my [1]

	opposite = new int* [7]; // For flowtypes 0 through 7
	for(int type=0; type<7; ++type) opposite[type] = new int [3]; // One for each branch

	opposite[0][0] = 0;
	opposite[0][1] = 0;
	opposite[0][2] = 0;

	opposite[1][0] = 13; // N/A in this case, must consider next_lowest_branch instead
	opposite[1][1] = 6;
	opposite[1][2] = 5;

	opposite[2][0] = 6; 
	opposite[2][1] = 13; // N/A in this case, must consider next_lowest_branch instead
	opposite[2][2] = 4;

	opposite[3][0] = 5; 
	opposite[3][1] = 4;
	opposite[3][2] = 13; // N/A in this case, must consider next_lowest_branch instead

	opposite[4][0] = 13; // N/A in this case, must consider next_lowest_branch instead
	opposite[4][1] = 3;
	opposite[4][2] = 2;

	opposite[5][0] = 3; 
	opposite[5][1] = 13; // N/A in this case, must consider next_lowest_branch instead
	opposite[5][2] = 1;

	opposite[6][0] = 2; 
	opposite[6][1] = 1;
	opposite[6][2] = 13; // N/A in this case, must consider next_lowest_branch instead
/*
	for(type=0; type<7; ++type)
	{
		for(int branch=0; branch<3; ++branch)
		{
			cout << "opposite[" << type << "][" << branch << "] = " << opposite[type][branch];
			cout << endl;
		}
		cout << endl;
	}
*/
	low_flow_branch = 0; // Arbitrary - just to fill the int
	next_lowest_branch = 1; // Arbitrary - just to fill the int

	flow_type_winterbone = 0;
	flow_type_winterbone_old = 0;
	flow_type_winterbone_next = 0;
	flow_type_actual = 0;
	
	// Enumerate array identifiers
	a = 0; // e.g. L[a] = L[0]
	b = 1;

	label = 0;
	value = 1;

	// Set arbitrary flowtype
	com = 0;
	i2 = 1;
	i3 = 2;

	// Loop limits
//	loop_limit_same = 2000;	// Number of main loops before continuing if stuck in one flowtype
	no_flow = 1e-12;

	switch_counter = 0;

	// NHPLJ:

	ALTERNATING = false;
	FIX_FLOW_TYPE = false;
	OPPOSITE = false; FIRST_OPPOSITE = true;
	OTHER = false;
	ZERO_FLOW = false;
	
	TESTED_ACTUAL = false;
	TESTED_ORIGINAL = false;

	fix_counter = 0;
	other_counter = 0;
	zero_counter = 0;

	normal_counter_total = 0;
	fix_counter_total = 0;
	opposite_counter_total = 0;
	other_counter_total = 0;
	zero_counter_total = 0;
	
	// Join the pipe pointer(s) to the relevant pipe(s)
	InitialiseGen(pPpt, pPipes, rPipe, rJPIPES, rJPIPES_ENDS, rENDCORR, id, ex, NPIPES, assyid, "CJunction", parent_assy_res_dir);

//cout << NPIPES;
//exit(1);
	// Boundary name
	NAME = new int [NPIPES];
	for(p=0; p<NPIPES; ++p)
	{
		NAME[p] = JUNCTION;
	}

	std::string bcname_str = "JUNCTION";
	ReadInput(pPpt, ConstructString(pPpt, param_dir, bcname_str, EX, ID));

	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Initialise\n");}

//	std::string res_str;
//	if(EX) res_str = "res_ex_junc"; else res_str = "res_in_junc";
//	std::cout << ConstructString(pPpt, res_dir, res_str, ID) << std::endl;
//	OUTPUT_FILE = fopen(ConstructString(pPpt, res_dir, res_str, ID), "w");

	// Set up files
	PRINT_MOVIE_FILE = false; // Default

	std::string res_str;
	if(EX) res_str = "res_ex_junc"; else res_str = "res_in_junc";
	SetupFiles(pPpt, res_str);

/*
	if(this->type==NHPLJID)
	{
		fprintf(OUTPUT_FILE,"%s\t%c%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
			"Time (s)", Deg(), "CA", "angle_rad", "psiT", "flow_type_winterbone", 
			"Mn", "W[a]", "W[b]", "K_loss[a][value]", "K_loss[b][value]", "L[a]", "L[b]",
			"x_star[pipe2]", "x_star[pipe3]");
	}
*/

	// NHPLTJ - three-way branch with pressure loss, 'equal-area' tee-junction (right-angle)
	G1 = new double [3];
	G2 = new double [3];
	// Empirical loss coefficients
	C1 = 0.3, C2 = 0.6; C3 = 0.75; C4 = 0.9; C5 = 0.9; C6 = 0.85;

	// NHPLJ stuff:
	lambda_in_star = new double [NPIPES];
	lambda_in_star_n = new double [NPIPES];
	AA = new double [NPIPES];
	AA_n = new double [NPIPES];
	pipe_flow = new int [NPIPES];
	pipe_flow_old = new int [NPIPES];
	lambda_out = new double [NPIPES];
	
	A_star = new double [NPIPES];
	U_star = new double [NPIPES];
	del_A_star = new double [NPIPES];
	del_A_star_old = new double [NPIPES];

	C = new double [NPIPES];
	del_star = new double [NPIPES];
	DATUM = 0;

	m_dot = new double [NPIPES];

	W = new double [2];
	L = new double [2];
	K_loss = new double* [2];
	for(int k=0; k<2; ++k) K_loss[k] = new double [2];

	coeff = new int [2];
	x_star = new double [3]; // Different from x
	x_star_old = new double [3]; // Different from x

	first_loop = true;
	first_loop_ever = true;

	// HCPJ stuff:
	Kq = new double [NPIPES];
	hcpj_first_loop = true;

	// NHPLJ stuff:
	nhplj_first_loop = true;

	// Pressure loss models
	m_dot_total = 0;
	for(p=0; p<NPIPES; ++p)
	{
		A_star[p] = 1.0;
		del_A_star[p] = 0;
		del_A_star_old[p] = del_A_star[p];
		del_star[p] = 0;
		C[p] = 0;
		m_dot[p] = 0;
	}

	error1 = 1000;
	error2 = 1000;

	MAX_FLOW_BRANCH = 0;

//	rollback = false;
}

void CJunction::ConfigureBranchAreas(CProperties* pPpt) // Call after boundary configure function
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ConfigureBranchAreas\n");}
	// Enter branch areas, Fb[NPIPES] 
	Fb = new double [NPIPES];			// Pipe actual areas
	for(int p=0; p<NPIPES; ++p) Fb[p] = pBN[p]->f_dash*pPpt->fref;
	
	this->psiT = Fb[pipe3]/Fb[pipe2]; // com area/branch area, pipe2 is always the branch pipe
}

void CJunction::RunBoundary(CProperties* pPpt, double time, bool &rRESTORE, int timestep)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".RunBoundary\n");}
//if(timestep>=291) cout << "timestep = " << timestep << " (1)\n"; 
	if(pPpt->HOMENTROPIC) HCPJ(pPpt);
	else
	{
		if(CONSTP || time < constp_until_time)
		{
			if(time < constp_until_time){pPpt->Out("Forcing constant pressure junction model since time ("); pPpt->Out(time); pPpt->Out(") < constp_until_time ("); pPpt->Out(constp_until_time); pPpt->Out(")\n");}
			NHCPJ(pPpt);
		}
		else
		{
			if(!BWPLJ(pPpt, timestep))
			{
        pPpt->Out("BWPLJ model solution didn't converge; running constant pressure model\n"); 
				NHCPJ(pPpt); // Run constant pressure model if solution didn't converge
			} 
		}
	}
//if(timestep>=291) cout << "timestep = " << timestep << " (2)\n"; 
	// Measure junction flow conditions
	// ----------------------------------------------------------------------------------------------------

/*
	// Identify MAX_FLOW_BRANCH branch, i.e., that with greatest mass flow TOWARDS junction
	// ----------------------------------------------------------------------------------------------------
	double temp_mfr, greatest_positive_mfr; greatest_positive_mfr = 0;
	int j;
	for(j=0; j<NPIPES; ++j)
	{
		temp_mfr = (pow(A_star[j], 2/(pPpt->gammaAir(pBN[j]->T) - 1))/AA[j])*U_star[j]*Fb[j];
		if(temp_mfr > greatest_positive_mfr)
		{
			greatest_positive_mfr = temp_mfr;
			MAX_FLOW_BRANCH = j;
		}
	}
*/
	// Resultant mass flow rate (positive towards junction)
	// ----------------------------------------------------------------------------------------------------
//	m_dot_total = 0;
//	int j;
//	Ws_total = 0;
//	for(j=0; j<NPIPES; ++j)
//	{
//		m_dot[j] = (pow(A_star[j], 2/(pPpt->gammaAir(pBN[j]->T) - 1))/AA[j])*U_star[j]*Fb[j]; // Expression for mass flow used in BWPLJ
		
//cout << "m_dot[j=" << j << "] = " << m_dot[j] << " kg/s\n\n"; 

//		m_dot[j] = 
//		pPpt->gammaAir(pBN[j]->T) * ((pPpt->PREF*1e5)/pow(pPipe[0]->AREF,2)) * (pow(A_star[j],(2*pPpt->gammaAir(pBN[j]->T))/(pPpt->gammaAir(pBN[j]->T)-1))/pow(A_star[j]*AA[j],2))									// rho (kg/m^3)
//		*(U_star[j]*AA[j]*pPipe[0]->AREF)	// u (m/s)
//		*Fb[j];								// F (m^2)


//cout << "(U_star[j]*AA[j]*pPipe[0]->AREF) = " << (U_star[j]*AA[j]*pPipe[0]->AREF) << " m/s\n";
//cout << "Fb[j] = " << Fb[j] << " m^2\n";
//cout << "m_dot[j=" << j << "] = " << m_dot[j] << " kg/s\n\n";

//cout << "rho pBN[j=" << j << "] = " << pBN[j]->rho << " kg/m^3\n";
//cout << "u pBN[j=" << j << "] = " << pBN[j]->U*pPipe[0]->AREF << " m/s\n";
//cout << "F pBN[j=" << j << "] = " << pBN[j]->f_dash*pPpt->fref << " m^2\n";
//cout << "m_dot pBN[j=" << j << "] = " << pBN[j]->mdot << " kg/s\n";

//		m_dot_total += m_dot[j];

		// Sum isentropic powers in all branches with flow towards the junction

//		if(m_dot[j]>0)
//		{
//			Ws_total += m_dot[j]*pow(sqrt(2*pPpt->cpAir(pBN[j]->T)*T0[limbID][ROTOR]*(1-pow(p2_shaft/p0[limbID][ROTOR],(pPpt->gammaAir(pBN[MAX_FLOW_BRANCH]->T)-1)/pPpt->gammaAir(pBN[MAX_FLOW_BRANCH]->T)))),2)/2;
//	}
//cout << "\n";
	//Ws_max_flow_branch = m_dot[MAX_FLOW_BRANCH]*pow(sqrt(2*pPpt->cpAir(pBN[MAX_FLOW_BRANCH]->T)*T0[limbID][ROTOR]*(1-pow(p2_shaft/p0[limbID][ROTOR],(pPpt->gammaAir(pBN[MAX_FLOW_BRANCH]->T)-1)/pPpt->gammaAir(pBN[MAX_FLOW_BRANCH]->T)))),2)/2;
}

void CJunction::HCPJ(CProperties* pPpt)
//--------------------------------------------------//
// Homentropic, constant pressure junction			//
// -------------------------------------------		//
// NPIPES-way branch with no pressure loss.			//
// Based on algorithm by Benson.					//
//--------------------------------------------------//
{
	double *lambda_in, *lambda_out;
	lambda_in = new double [NPIPES];
	lambda_out = new double [NPIPES];
	int * pipe_flow;
	pipe_flow = new int [NPIPES];

	int p;
	double KK;

	for(p=0; p<NPIPES; ++p) lambda_in[p] = (*pCLIN[p])[R+1];

	if(this->hcpj_first_loop) // Only need to set FT and K[] once per simulation
	{
		this->FT = 0;
		for(p=0; p<NPIPES; ++p)
			this->FT += this->Fb[p]; // Total area (m^2)
		
		for(p=0; p<NPIPES; ++p)
			this->Kq[p] = 2*this->Fb[p]/this->FT;

		this->hcpj_first_loop = false;
	}

	KK = 0;
	for(p=0; p<NPIPES; ++p)
			KK += this->Kq[p]*lambda_in[p];

	for(p=0; p<NPIPES; ++p)
	{
		lambda_out[p] = KK - lambda_in[p];
		
		double U_temp;
		U_temp = (lambda_in[p] - lambda_out[p])/(pPpt->gammaAir(pBN[p]->T)-1);
		// This U will be +ve toward the junction, which is OUTFLOW
		if(U_temp>0) pipe_flow[p] = OUTFLOW;
		else
		{ 
			if(U_temp<0) pipe_flow[p] = INFLOW;
			else pipe_flow[p] = NOFLOW;
		}
	}

	for(p=0; p<NPIPES; ++p) (*pCLOUT[p])[R+1] = lambda_out[p];
	// Note there is no update to lambda_in required in HCPJ
	for(p=0; p<NPIPES; ++p) *(pend_flow[p]) = pipe_flow[p];

	delete [] lambda_in;
	delete [] lambda_out;
	delete [] pipe_flow;
}

void CJunction::NHCPJ(CProperties* pPpt)
//--------------------------------------------------//
// Non-homentropic, constant pressure junction      //
// -------------------------------------------		//
// NPIPES-way branch with no pressure loss.				//
// Based on algorithm by Benson.					//
//--------------------------------------------------//
{
	int p;

	double E, E_old;
	double A_star_N;
	double sum_num, sum_den;
	double sum_num_ent, sum_den_ent;
	double AA_NS;

	E = 0;
	E_old = 0;
	A_star_N = 0;

	// Store uncorrected values and intialise working variables
	// --------------------------------------------------------
	CommonSetup();
	
	int counter=0;
	int loop_limit=1000;
	bool stop=false;
	do
	{
		++counter;
		if(counter>loop_limit)
		{
			cout << this->ID << endl;
			if(this->EX) cout << "Exhaust" << endl;
			else cout << "Intake" << endl;
			cout << counter << endl;
			cout << "E_old = " << E_old << endl;
			cout << "E = " << E << endl;
			cout << "A_star_N = " << A_star_N << endl;
			cout << "(A_star_N - E) = " << (A_star_N - E) << endl;
			cout << "fabs(A_star_N - E) = " << fabs(A_star_N - E) << endl;
			cout << "Setting A_star_N as the average of the oscillatory values..." << endl;
		//	A_star_N = 0.5*(A_star_N + E);
			stop = true;
		}

		E_old = E; // Record old value for oscillation test
		E = A_star_N; // Record old value for convergence test

		sum_num = 0;
		sum_den = 0;
		for(p=0; p<NPIPES; ++p) sum_num += (lambda_in_star[p]/AA[p])*Fb[p];
		for(p=0; p<NPIPES; ++p) sum_den += Fb[p]/AA[p];
		A_star_N = sum_num/sum_den;
		// Eq. 8.132 of Benson p439
	
		for(p=0; p<NPIPES; ++p) U_star[p] = (2/(pPpt->gammaAir(pBN[p]->T)-1))*(lambda_in_star[p] - A_star_N);
		// Eq. 8.133 of Benson p439

		// Test for joining flow only in each pipe
		sum_num_ent = 0;
		sum_den_ent = 0;
		AA_NS =  0;
		for(p=0; p<NPIPES; ++p)
		{
			if(U_star[p]>0) // Joining/towards the junction
			{
				// Calculate the weighted mean of the entropy levels
				// of the JOINING flows
				sum_num_ent += U_star[p]*Fb[p]*AA[p];
				sum_den_ent += U_star[p]*Fb[p];
			}
		}

		// Calculate the weighted mean of the entropy levels
		// of the JOINING flows
		AA_NS = sum_num_ent/sum_den_ent;
		// Eq. 8.135 of Benson p439

		// Test for joining/separating flow in each pipe, 
		// and give it the calculated weighted entropy in 
		// the case of separating flow
		double limit = 1e-6;
		for(p=0; p<NPIPES; ++p)
		{
			if(U_star[p]>limit/*0*/) // Joining/towards the junction
			{
				pipe_flow[p] = OUTFLOW;
			//	AA[p] = AA[p];	// i.e., remains unchanged
			}
			else
			{
				if(U_star[p]<(-1*limit)/*0*/) // Away from the junction
				{
					pipe_flow[p] = INFLOW;
					AA[p] = AA_NS;
				}
				else
				{
					pipe_flow[p] = NOFLOW;
				//	AA[p] = AA[p];	// i.e., remains unchanged ????? or
				//	AA[p] = AA_NS; // ????
				}
			}
		}

		// Apply the correction to lambda_in_star
		for(p=0; p<NPIPES; ++p)
		{
			//change described 15-Oct-03 
			//lambda_in_star[p] = A_star_N + (AA_n[p]/AA[p])*(lambda_in_star[p] - A_star_N);
			lambda_in_star[p] = A_star_N + (AA_n[p]/AA[p])*(lambda_in_star_n[p] - A_star_N);
			// Eq. 8.137 of Benson p440
		
			lambda_out[p] = AA[p]*(2*A_star_N - lambda_in_star[p]);
		}

	}while(fabs(A_star_N - E)>tol_A_star/*1e-6*/ && !stop);

	// Converged, evaluate flow conditions
	// ----------------------------------------------------------------------------------------------------
	m_dot_total = 0;
	double greatest_positive_mfr; greatest_positive_mfr = 0;
	for(p=0; p<NPIPES; ++p)
	{
		m_dot[p] = 
		pPpt->gammaAir(pBN[p]->T) * ((pPpt->PREF*1e5)/pow(pPipe[0]->AREF,2)) * (pow(A_star_N,(2*pPpt->gammaAir(pBN[p]->T))/(pPpt->gammaAir(pBN[p]->T)-1))/pow(A_star_N*AA[p],2))						// rho (kg/m^3)
		*(U_star[p]*AA[p]*pPipe[0]->AREF)	// u (m/s)
		*Fb[p];								// F (m^2)

//cout << "NHCPJ m_dot[p=" << p << "] = " << m_dot[p] << " kg/s\n";

		if(m_dot[p] > greatest_positive_mfr)
		{
			greatest_positive_mfr = m_dot[p];
			MAX_FLOW_BRANCH = p;
		}
		m_dot_total += m_dot[p];
	}

	// Converged, pass back the new values by pointers
	// ----------------------------------------------------------------------------------------------------
	CommonUpdate(pPpt);
/*
	//////////////////////////////////////////
	double *lambda_in_c; lambda_in_c = new double [NPIPES];
	for(p=0; p<NPIPES; ++p)	lambda_in_c[p] = lambda_in_star[p]*AA[p];

	bool* CHOKED; CHOKED = new bool [NPIPES];
	for(p=0; p<NPIPES; ++p) CHOKED[p] = false;

	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out, AA, pipe_flow, CHOKED);
	//////////////////////////////////////////
*/
	// Delete any local arrays
	// -----------------------
}

bool CJunction::BWPLJ(CProperties* pPpt, int timestep)
// ====================================================================================================
// Generalized Bassett-Winterbone pressure loss junction
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".BWPLJ\n");}

	int j;
	double k;
  double max_del_star_error;

	// Store uncorrected values and intialise working variables
	// ----------------------------------------------------------------------------------------------------
	CommonSetup();

	// Dimension temporary variables
	// ----------------------------------------------------------------------------------------------------
	double* A_star_prev;			A_star_prev = new double [NPIPES];
	double* U_star_prev;			U_star_prev = new double [NPIPES];
	int* pipe_flow_prev;			pipe_flow_prev = new int [NPIPES];
	double* AA_prev;				AA_prev = new double [NPIPES];
	double* del_star_prev;			del_star_prev = new double [NPIPES];
	double* lambda_in_star_prev;	lambda_in_star_prev = new double [NPIPES];

	// Estimate initial values of A* and U* using lambda_in_star and old time level lambda_out_star
	// ----------------------------------------------------------------------------------------------------
	for(j=0; j<NPIPES; ++j)
	{
		k = pPpt->gammaAir(pBN[j]->T);
		A_star[j] = (lambda_in_star[j] + ((*pCLOUT[j])[R+1]/AA[j]))/2;
		U_star[j] = (lambda_in_star[j] - ((*pCLOUT[j])[R+1]/AA[j]))/(k - 1);
	}
	
	// Initialize DATUM and pressure loss terms - set zero values on first time step, else use previous solution values
	// ----------------------------------------------------------------------------------------------------
	if(timestep==1)
	{
		for(j=0; j<NPIPES; ++j){DATUM = 0; C[j] = 0; del_star[j] = 0; del_star_prev[j] = 1000;}
		h_0_sep = 0;
	}

	// Initially identify datum branch and determine pressure loss terms
	// ----------------------------------------------------------------------------------------------------
	int counter, lambda_counter;
	bool SOLUTION_CONVERGED, ESCAPE, LAMBDA_CONVERGED, USTAR_CONVERGED;
	//double tol_entropy = 1e-15;//0.1; // Starting value
	//double tol_U_star = 1e-15;//0.1; // Starting value
	double enthalpy_joining, enthalpy_separating;
	double h_0_sep_prev, h_0_sep_num, h_0_sep_den, h_0_sep_old;

	counter = 0;
	ESCAPE = false;
	BWPLJLossTerms(pPpt, timestep, counter, ESCAPE);

//if(timestep>=17980 && counter>=153) cout << "1: before outer loop\n\n";

	do{ // Outer loop

		++counter;
//if(timestep>=17980 && counter>=153) cout << "\t2: top of outer loop, counter = " << counter << endl << endl;
		lambda_counter = 0;
		do{ // Lambda loop
//if(timestep>=17980 && counter>=153) cout << "\t\t3: top of lambda loop, lambda_counter = " << lambda_counter << endl << endl;
			++lambda_counter;

//			entropy_counter = 0;
//			do{ // Entropy loop
//				++entropy_counter;

				//flow_type_counter = 0;
				//do{ // Flow type loop
					//++ flow_type_counter;

//					U_star_counter = 0;
//					do{ // Continuity, converging on A_star, U_star loop	
//						++U_star_counter;

						// Solve continuity equation for A_star[DATUM] by Newton-Raphson
						// ----------------------------------------------------------------------------------------------------
//if(timestep>=291) cout << "timestep = " << timestep << " (1.0)\n";
						BWPLJContinuityNR(pPpt, timestep, counter);
						//BWPLJContinuityFP(pPpt, timestep);
//if(timestep>=17980 && counter>=153) cout << "\t\t5: after BWPLJContinuityNR\n\n";
						// Solving continuity for A_star[DATUM] provides updated values for A_star[j] and U_star[j]
//if(timestep>=291) cout << "timestep = " << timestep << " (1.01)\n";
						// ----------------------------------------------------------------------------------------------------
						USTAR_CONVERGED = true;
						for(j=0; j<NPIPES; ++j)
						{
							k = pPpt->gammaAir(pBN[j]->T);
							A_star_prev[j] = A_star[j];
							U_star_prev[j] = U_star[j];
							A_star[j] = A_star[DATUM] - del_star[j]; // From definition of del_star[j] = A_star[DATUM] - A_star[j];
							U_star[j] = (2/(k - 1))*(lambda_in_star[j] - A_star[j]); // From definiton of lambda_in_star																   

							// Check for convergence of A_star, U_star
							// ----------------------------------------------------------------------------------------------------
						//	if(fabs(A_star[j] - A_star_prev[j]) > tol_U_star) USTAR_CONVERGED = false;
						//	if(fabs(U_star[j] - U_star_prev[j]) > tol_U_star) USTAR_CONVERGED = false;
						}
/*
if(timestep>=171 && counter>=1)
{
cout << "main counter = " << counter << endl;
cout << "lambda_counter = " << lambda_counter << endl;
//cout << "entropy_counter = " << entropy_counter << endl;
cout << "U_star_counter = " << U_star_counter << endl;
//cout << "entropy_counter_inner = " << entropy_counter_inner << endl;
cout << "USTAR_CONVERGED = " << TrueOrFalse(USTAR_CONVERGED) << endl;
//cout << "tol_entropy = " << tol_entropy << endl;
for(j=0; j<NPIPES; ++j)
{
//cout << "lambda_in_star_prev[j=" << j << "] = " << lambda_in_star_prev[j] << endl;
cout << "A_star_prev[j=" << j << "] = " << A_star_prev[j] << endl;
cout << "A_star[j=" << j << "] = " << A_star[j] << endl;
cout << "fabs(A_star - A_star_prev) = " << fabs(A_star[j] - A_star_prev[j]) << endl;
cout << "U_star_prev[j=" << j << "] = " << U_star_prev[j] << endl;
cout << "U_star[j=" << j << "] = " << U_star[j] << endl;
cout << "fabs(U_star - U_star_prev) = " << fabs(U_star[j] - U_star_prev[j]) << endl;
//cout << "lambda_in_star[j=" << j << "] = " << lambda_in_star[j] << endl;
//cout << "fabs(lambda_in_star - lambda_in_star_prev) = " << fabs(lambda_in_star[j] - lambda_in_star_prev[j]) << endl;
}
cout << endl;
cin >> pause;
}
//*/

/*
if(U_star_counter>1000)
{
	cout << "timestep = " << timestep << endl;
	cout << "main counter = " << counter << endl;
	cout << "lambda_counter = " << lambda_counter << endl;
	//cout << "entropy_counter_inner = " << entropy_counter_inner << endl;
	cout << "U_star_counter = " << U_star_counter << endl;
	exit(1);
}
//*/
					//}while(!USTAR_CONVERGED);

					// Record pipe flow directions
					// ----------------------------------------------------------------------------------------------------
					//FLOW_TYPE_CONVERGED = true;
					for(j=0; j<NPIPES; ++j)
					{
						pipe_flow_prev[j] = pipe_flow[j];
						if(U_star[j] > 0) pipe_flow[j] = OUTFLOW;	// Flow is towards junction, i.e., out of pipe
						else
						{
							if(U_star[j] < 0) pipe_flow[j] = INFLOW;// Flow is away from junction, i.e., into pipe
							else pipe_flow[j] = NOFLOW;
						}
					}

					//entropy_counter_inner = 0;
					//do{
					//	++entropy_counter_inner;

						// Calculate mass-averaged stagnation enthalpy of joining flows
						// ----------------------------------------------------------------------------------------------------
						h_0_sep_prev = h_0_sep;
						//h_0_sep = 0; 
						h_0_sep_num = 0; h_0_sep_den = 0;
						h_0_sep_old = 0;
						for(j=0; j<NPIPES; ++j){
							k = pPpt->gammaAir(pBN[j]->T);
							if(U_star[j] > 0) // Sum joining flows only
							{
								h_0_sep_num += (pow(pPipe[j]->AREF,2)/(k - 1))*
												( pow(A_star[j],2/(k-1)) 
													* U_star[j] * Fb[j] * AA[j] 
													* (pow(A_star[j],2) + ((k - 1)/2)*pow(U_star[j],2)) 
												);

								h_0_sep_den += (pow(A_star[j],2/(k-1))/AA[j]) * U_star[j] * Fb[j];

								h_0_sep_old += ((pow(pPipe[j]->AREF,2)/(k - 1))*
												( pow(A_star[j],2/(k - 1)) 
													* U_star[j] * Fb[j] * AA[j] 
													* (pow(A_star[j],2) + ((k - 1)/2)*pow(U_star[j],2)) 
												))/
												((pow(A_star[j],2/(k - 1))/AA[j]) * U_star[j] * Fb[j]);
							}
						}
						if(h_0_sep_den!=0) h_0_sep = h_0_sep_num/h_0_sep_den; else h_0_sep = 0;
				
						// Calculate entropy of separating flows
						// ----------------------------------------------------------------------------------------------------
					//	ENTROPY_CONVERGED_INNER = true;
						for(j=0; j<NPIPES; ++j){
							k = pPpt->gammaAir(pBN[j]->T);
							AA_prev[j] = AA[j];
							if(U_star[j] < 0) // Only need to calculate for separating flows; joining entropies unchanged
							{
								if(h_0_sep > 0) // Only calculate for non-zero enthalpies, else AA[j] = 0 and causing div by 0 in lambda_in_star
									AA[j] = sqrt( h_0_sep*((k - 1)/pow(pPipe[j]->AREF,2))
												 /
												  (pow(A_star[j],2) + ((k - 1)/2)*pow(U_star[j],2))
												);
							}
			
							// Check for convergence of entropy at each pipe end
							// ----------------------------------------------------------------------------------------------------
					//		if(fabs(AA[j] - AA_prev[j]) > tol_entropy) ENTROPY_CONVERGED_INNER = false;
						}


						// Calculate energy balance
						// ----------------------------------------------------------------------------------------------------
						enthalpy_joining = 0;
						enthalpy_separating = 0;
						for(j=0; j<NPIPES; ++j){
							k = pPpt->gammaAir(pBN[j]->T);
							if(U_star[j] > 0){ // Sum joining flows
								//enthalpy_joining += ((k*pPpt->PREF*pPipe[j]->AREF)/(k-1))*
								//				    ( pow(A_star[j],2/(k-1)) * U_star[j] * Fb[j] * AA[j] * (pow(A_star[j],2) + ((k-1)/2)*pow(U_star[j],2)) );
								enthalpy_joining += ((k*(pPpt->PREF*1e5)*pPipe[j]->AREF)/(k-1))*
												    ( pow(A_star[j],2/(k-1)) * U_star[j] * Fb[j] * AA[j] * (pow(A_star[j],2) + ((k-1)/2)*pow(U_star[j],2)) );
							}
							if(U_star[j] < 0){ // Sum separating flows
								//enthalpy_separating += ((k*pPpt->PREF*pPipe[j]->AREF)/(k-1))*
								//					   ( pow(A_star[j],2/(k-1)) * U_star[j] * Fb[j] * AA[j] * (pow(A_star[j],2) + ((k-1)/2)*pow(U_star[j],2)) );
								enthalpy_separating += ((k*(pPpt->PREF*1e5)*pPipe[j]->AREF)/(k-1))*
													   ( pow(A_star[j],2/(k-1)) * U_star[j] * Fb[j] * AA[j] * (pow(A_star[j],2) + ((k-1)/2)*pow(U_star[j],2)) );
							}
						}
/*
if(entropy_counter_inner>1000)
{
	cout << "timestep = " << timestep << endl;
	cout << "main counter = " << counter << endl;
	cout << "lambda_counter = " << lambda_counter << endl;
	cout << "entropy_counter_inner = " << entropy_counter_inner << endl;
	exit(1);
}
//*/
				//	}while(!ENTROPY_CONVERGED_INNER);				

/*
if(h_0_sep_old != h_0_sep)
{
cout << "main counter = " << counter << endl;
cout << "h_0_sep_old = " << h_0_sep_old << endl;
cout << "h_0_sep = " << h_0_sep << endl;
//cout << "entropy_counter = " << entropy_counter << endl;
//cout << "flow_type_counter = " << flow_type_counter << endl;
//cout << "U_star_counter = " << U_star_counter << endl;
//cout << "entropy_counter_inner = " << entropy_counter_inner << endl;
//cout << "FLOW_TYPE_CONVERGED = " << TrueOrFalse(FLOW_TYPE_CONVERGED) << endl;
//cout << "tol_entropy = " << tol_entropy << endl;
for(j=0; j<NPIPES; ++j)
{

//cout << "AA_prev[j=" << j << "] = " << AA_prev[j] << endl;
//cout << "AA[j=" << j << "] = " << AA[j] << endl;
//cout << "fabs(AA - AA_prev) = " << fabs(AA[j] - AA_prev[j]) << endl;
}
cout << endl;
cin >> pause;
}
//*/

		//	}while(!ENTROPY_CONVERGED);

/*
if(timestep>=82 && counter>=1)
{
cout << "main counter = " << counter << endl;
cout << "lambda_counter = " << lambda_counter << endl;
//cout << "entropy_counter = " << entropy_counter << endl;
cout << "U_star_counter = " << U_star_counter << endl;
cout << "entropy_counter_inner = " << entropy_counter_inner << endl;
//cout << "ENTROPY_CONVERGED = " << TrueOrFalse(ENTROPY_CONVERGED) << endl;
//cout << "tol_entropy = " << tol_entropy << endl;
for(j=0; j<NPIPES; ++j)
{
cout << "lambda_in_star_prev[j=" << j << "] = " << lambda_in_star_prev[j] << endl;
cout << "A_star[j=" << j << "] = " << A_star[j] << endl;
cout << "U_star[j=" << j << "] = " << U_star[j] << endl;
cout << "lambda_in_star[j=" << j << "] = " << lambda_in_star[j] << endl;
cout << "fabs(lambda_in_star - lambda_in_star_prev) = " << fabs(lambda_in_star[j] - lambda_in_star_prev[j]) << endl;
}
cout << endl;
cin >> pause;
}
//*/

/*
if(lambda_counter>1000)
{
	cout << "timestep = " << timestep << endl;
	cout << "main counter = " << counter << endl;
	cout << "lambda_counter = " << lambda_counter << endl;
	exit(1);
}
//*/

			// Apply converged entropy correction to lambda_in_star, to use in loop next time round
			// ----------------------------------------------------------------------------------------------------
			LAMBDA_CONVERGED = true;
			for(j=0; j<NPIPES; ++j){
				lambda_in_star_prev[j] = lambda_in_star[j];
				lambda_in_star[j] = lambda_in_star_n[j]*(AA_n[j]/AA[j]) + A_star[j]*(1 - (AA_n[j]/AA[j]));
				
				// Check for convergence of entropy at each pipe end
				// ----------------------------------------------------------------------------------------------------
				//if(fabs(lambda_in_star[j] - lambda_in_star_prev[j]) > tol_del_star) LAMBDA_CONVERGED = false;

				if(fabs(lambda_in_star[j] - lambda_in_star_prev[j]) > tol_del_star){
					LAMBDA_CONVERGED = false;
///*
					if(lambda_counter>loop_limit_main)
          //if(lambda_counter>100)
					{
            int temp_precision = cout.precision(); cout << setprecision(12);
						cout << "lambda_counter = " << lambda_counter << endl;
            cout << "lambda_in_star[j=" << j << "] = " << lambda_in_star[j] << endl;
            cout << "lambda_in_star_prev[j=" << j << "] = " << lambda_in_star_prev[j] << endl;
            cout << "error [j=" << j << "] = " << lambda_in_star[j] - lambda_in_star_prev[j] << endl;
            cout << endl;
            cout << setprecision(temp_precision);
						// Flow type is switching; wait until flow type is same as previous time step, then allow convergence
						//if(pipe_flow[j] == pipe_flow_old[j]) LAMBDA_CONVERGED = true;
					}
//*/
				}
			}

      //if(timestep>=291) cout << "timestep = " << timestep << " (1.1)\n";

		}while(!LAMBDA_CONVERGED && lambda_counter<loop_limit_main);
		//}while(!LAMBDA_CONVERGED);

/*
if(timestep>=82)
{
cout << "main counter = " << counter << endl;
cout << "lambda_counter = " << lambda_counter << endl;
cout << "entropy_counter = " << entropy_counter << endl;
cout << "entropy_counter_inner = " << entropy_counter_inner << endl;
cout << "U_star_counter = " << U_star_counter << endl;
//cout << "entropy_counter_inner = " << entropy_counter_inner << endl;
cout << "LAMBDA_CONVERGED = " << TrueOrFalse(LAMBDA_CONVERGED) << endl;
//cout << "tol_entropy = " << tol_entropy << endl;
cout << "enthalpy_joining = " << enthalpy_joining << endl;
cout << "enthalpy_separating = " << enthalpy_separating << endl;
for(j=0; j<NPIPES; ++j)
{
cout << "A_star[j=" << j << "] = " << A_star[j] << endl;
cout << "U_star[j=" << j << "] = " << U_star[j] << endl;
cout << "C[j=" << j << "] = " << C[j] << endl;
cout << "del_star_prev[j=" << j << "] = " << del_star_prev[j] << endl;
cout << "del_star[j=" << j << "] = " << del_star[j] << endl;
cout << "fabs(del_star - del_star_prev) = " << fabs(del_star[j] - del_star_prev[j]) << endl;
//cout << "h_0_sep = " << h_0_sep << endl;
//cout << "AA_prev[j=" << j << "] = " << AA_prev[j] << endl;
//cout << "AA[j=" << j << "] = " << AA[j] << endl;
//cout << "fabs(AA - AA_prev) = " << fabs(AA[j] - AA_prev[j]) << endl;
cout << endl;
}
cout << endl;
cin >> pause;
}
//*/
//if(timestep>=291) cout << "timestep = " << timestep << " (1.2)\n";
		// Determine new pressure loss terms, once everything else has converged
		// ----------------------------------------------------------------------------------------------------
		for(j=0; j<NPIPES; ++j) del_star_prev[j] = del_star[j];	// Save previous values	
		BWPLJLossTerms(pPpt, timestep, counter, ESCAPE);
//if(timestep>=291) cout << "timestep = " << timestep << " (1.3)\n";
		// Check for convergence of pressure loss terms
		// ----------------------------------------------------------------------------------------------------
		SOLUTION_CONVERGED = true;
    max_del_star_error = 0;
		for(j=0; j<NPIPES; ++j)
    {
      if(fabs(del_star[j] - del_star_prev[j]) > max_del_star_error) max_del_star_error = fabs(del_star[j] - del_star_prev[j]);

      if(fabs(del_star[j] - del_star_prev[j]) > tol_del_star)
      {
        SOLUTION_CONVERGED = false;
/*
        if(counter>9000)
        {
          cout << "counter = " << counter << endl;
          cout << "del_star[j=" << j << "] = " << del_star[j] << endl;
          cout << "del_star_prev[j=" << j << "] = " << del_star_prev[j] << endl;
          cout << "error = " << del_star[j] - del_star_prev[j] << endl;
          cout << endl;
        }
*/
      }
    }

		if(counter <= 1) SOLUTION_CONVERGED = false; // Must repeat process at least once to calculate new pressure loss terms
/*
if(counter>1000)
{
	cout << "timestep = " << timestep << endl;;
	exit(1);
}
//*/
//if(timestep>=291) cout << "timestep = " << timestep << " (1.4)\n";
	}while(!SOLUTION_CONVERGED && counter<loop_limit_main && !ESCAPE);
	//}while(!SOLUTION_CONVERGED);

if(ESCAPE) cout << "ESCAPE = " << TrueOrFalse(ESCAPE) << endl;

//if(timestep>=291) cout << "timestep = " << timestep << " (1.5)\n";
/*
cout << "main counter = " << counter << endl;
cout << "entropy_counter = " << entropy_counter << endl;
cout << "U_star_counter = " << U_star_counter << endl;
cout << "SOLUTION_CONVERGED = " << TrueOrFalse(SOLUTION_CONVERGED) << endl;
cout << "tol_del_star = " << tol_del_star << endl;
for(j=0; j<NPIPES; ++j)
{
;
}
cout << endl;
cin >> pause;
//*/


	// Calculate lambda_out at new time level
	// ----------------------------------------------------------------------------------------------------
	for(j=0; j<NPIPES; ++j) lambda_out[j] = (2*A_star[j] - lambda_in_star[j])*AA[j];

	// Converged, pass back the new values by pointers
	// ----------------------------------------------------------------------------------------------------
	if(SOLUTION_CONVERGED) 
    CommonUpdate(pPpt); // Only update pipes if solution converged, else run constant pressure model
     
  //if(!SOLUTION_CONVERGED) cout << Identify() << ": BWPLJ not converged, max_del_star_error = " << max_del_star_error << endl;

	// Delete any local arrays
	// -----------------------
	delete [] A_star_prev;
	delete [] U_star_prev;
	delete [] pipe_flow_prev;
	delete [] AA_prev;
	delete [] del_star_prev;
	delete [] lambda_in_star_prev;


	//if(counter>100) cout << "counter = " << counter << endl;
	//if(lambda_counter>10) cout << "lambda_counter = " << lambda_counter << endl << endl;
	//cout << "counter = " << counter << endl;
	//cout << endl;


	return SOLUTION_CONVERGED;
}

void CJunction::BWPLJContinuityNR(CProperties* pPpt, int timestep, int counter)
{
	int j;
	// Solve continuity equation for A_star[DATUM] by Newton-Raphson
	// ----------------------------------------------------------------------------------------------------
	double k;
	double mass_sum, mass_sum_prev, mass_sum_prev_prev, grad_mass_sum, grad_mass_sum_prev, A_star_prev, factor;
	double mass_sum_left, mass_sum_left_prev, mass_sum_left_prev_prev;
	double mass_sum_right, mass_sum_right_prev, mass_sum_right_prev_prev;
	double A_star_left, A_star_right;

	bool CONTINUITY_CONVERGED = false;
	int continuity_counter = 0;
/*
	double remember_A_star = A_star[DATUM];
	int a_max = 100;
	double min_A_star = -0.5;
	double max_A_star = 2.5;
	for(int a=0; a<a_max; ++a)
	{
		mass_sum = 0;
		A_star[DATUM] = min_A_star + a*((max_A_star - min_A_star)/a_max);

		for(j=0; j<NPIPES; ++j)
		{
			// Function value at x = A_star[DATUM]
			mass_sum += ( pow(A_star[DATUM] - del_star[j], 2/(pPpt->gammaAir() - 1))/AA[j] )
						* (lambda_in_star[j] - A_star[DATUM] + del_star[j]) 
						* Fb[j];
		}
		pPpt->Out(A_star[DATUM]); pPpt->Out("\t"); pPpt->Out(mass_sum); pPpt->Out("\n");
	}
	A_star[DATUM] = remember_A_star;
exit(1);
*/

	A_star[DATUM] = 1.5;	// False position method
	//double sign = 1.0;		// Direction

	mass_sum_prev = 1000;
	mass_sum = 1000;

	mass_sum_left_prev = 1000;
	mass_sum_left = 1000;

	mass_sum_right_prev = 1000;
	mass_sum_right = 1000;

	grad_mass_sum = 0;
	//A_star_prev = A_star[DATUM];

	factor = 1;
	do{
		mass_sum_prev_prev = mass_sum_prev; // Save previous value
		mass_sum_prev = mass_sum;			// Save previous value
		mass_sum = 0;						// Reset

		mass_sum_left_prev_prev = mass_sum_left_prev;	// Save previous value
		mass_sum_left_prev = mass_sum_left;				// Save previous value
		mass_sum_left = 0;								// Reset

		mass_sum_right_prev_prev = mass_sum_right_prev;	// Save previous value
		mass_sum_right_prev = mass_sum_right;				// Save previous value
		mass_sum_right = 0;								// Reset

		grad_mass_sum_prev = grad_mass_sum; // Save previous value
		grad_mass_sum = 0;					// Reset

		A_star_left = A_star[DATUM]*0.9999;
		A_star_right = A_star[DATUM]*1.0001;

//if(timestep>=17980 && counter>=153) cout << "\t\t\t5: A_star[DATUM] = " << A_star[DATUM] << endl;

		for(j=0; j<NPIPES; ++j)
		{
			k = pPpt->gammaAir(pBN[j]->T);
			// Function value at x = A_star[DATUM]
			mass_sum += ( pow(A_star[DATUM] - del_star[j], 2/(k - 1))/AA[j] )
						* (lambda_in_star[j] - A_star[DATUM] + del_star[j]) 
						* Fb[j];
/*
if(timestep>=21631) // && counter>=153)
{
	cout << endl;
	cout << "\t\t\t5: A_star[DATUM] = " << A_star[DATUM] << endl;
	cout << "\t\t\t5: del_star[j=" << j << "] = " << del_star[j] << endl;
	cout << "\t\t\t5: lambda_in_star[j=" << j << "] = " << lambda_in_star[j] << endl;
	cout << "\t\t\t5: AA[j=" << j << "] = " << AA[j] << endl;
	
	
	cout << "\t\t\t5: (A_star[DATUM] - del_star[j]) = " << (A_star[DATUM] - del_star[j]) << endl;
	

	cout << "\t\t\t5: 2/(pPpt->gammaAir(pBN[j]->T) - 1) = " << 2/(pPpt->gammaAir(pBN[j]->T) - 1) << endl;
	
	cout << "\t\t\t5: ( pow(A_star[DATUM] - del_star[j], 2/(pPpt->gammaAir(pBN[j]->T) - 1))/AA[j] ) = " << ( pow(A_star[DATUM] - del_star[j], 2/(pPpt->gammaAir(pBN[j]->T) - 1))/AA[j] ) << endl;
	
	
	
	cout << "\t\t\t5: (lambda_in_star[j] - A_star[DATUM] + del_star[j])  = " << (lambda_in_star[j] - A_star[DATUM] + del_star[j])  << endl;
	cout << "\t\t\t5: mass_sum = " << mass_sum << endl;
	cout << endl;
	cout << endl;
}
//*/
			mass_sum_left += ( pow(A_star_left - del_star[j], 2/(k - 1))/AA[j] )
						* (lambda_in_star[j] - A_star_left + del_star[j]) 
						* Fb[j];

			mass_sum_right += ( pow(A_star_right - del_star[j], 2/(k - 1))/AA[j] )
						* (lambda_in_star[j] - A_star_right + del_star[j]) 
						* Fb[j];


	
			// Gradient of function at x = A_star[DATUM] required for N-R iteration

			
/*
			grad_mass_sum += (Fb[j]/AA[j]) 
							*(
								(lambda_in_star[j] - A_star[DATUM] + del_star[j])*
								( 
									(2/(pPpt->gammaAir()-1)) * pow(A_star[DATUM] - del_star[j], (3 - pPpt->gammaAir())/(pPpt->gammaAir() - 1))
									-
									pow(A_star[DATUM] - del_star[j], 2/(pPpt->gammaAir() - 1))
								)
							 );
*/


			

//cout << "mass_sum up to [j=" << j << "] = " << mass_sum << endl;
//cout << "del_star[j=" << j << "] = " << del_star[j] << endl;
//cout << "AA[j=" << j << "] = " << AA[j] << endl;
//cout << "del_star[j=" << j << "] = " << del_star[j] << endl;
//cout << "lambda_in_star[j=" << j << "] = " << lambda_in_star[j] << endl;
		}
//cout << endl;

		grad_mass_sum = (mass_sum_right - mass_sum_left)/(A_star_right - A_star_left);



/*	
//if(this->ID==0)
{
cout << "A_star[DATUM=" << DATUM << "] = " << A_star[DATUM] << endl;
cout << "grad_mass_sum[iteration=" << continuity_counter << "] = " << grad_mass_sum << endl;
cout << "mass_sum[iteration=" << continuity_counter << "] = " << mass_sum << endl;
cout << endl;
}
//exit(1);
//*/

		// Check for sign change in mass sum; if there is, go back and try again with smaller step
		// ----------------------------------------------------------------------------------------------------	
		if(mass_sum/mass_sum_prev < 0 && continuity_counter > 0) // There has been a sign change
		{
/*
			// Sign changes indicates solution lies between A_star[DATUM] and A_star_previous
			if(fabs(mass_sum) < fabs(mass_sum_prev))
			{
//cout << "SIGN CHANGE in mass sum; ALLOW!" << endl; // Allow sign change since result is closest to solution so far
				sign *= -1;			// Change sign since we are on other side of solution and must step in other direction
				factor *= 0.5;				// Also halve factor else will simply return to previous A_star[DATUM]
			}
			else	// Disallows sign change as result is further away from solution
			{
*/

/*
if(this->ID==0)
{
cout << "SIGN CHANGE in mass sum; restore previous values and reduce step size" << endl;
//cout << "Current values:" << endl;
//cout << "A_star[DATUM=" << DATUM << "] = " << A_star[DATUM] << endl;
//cout << "grad_mass_sum[iteration=" << continuity_counter << "] = " << grad_mass_sum << endl;
//cout << "mass_sum[iteration=" << continuity_counter << "] = " << mass_sum << endl;
}
//*/
				//A_star[DATUM] = A_star_prev;		// Restore previous value
				//mass_sum = mass_sum_prev;			// Restore previous value	
				//grad_mass_sum = grad_mass_sum_prev;	// Restore previous value
				factor *= 0.5;						// Apply only a fraction of the calculated step size by N-R method

/*
if(this->ID==0)
{
//cout << "Restored values:" << endl;
//cout << "A_star[DATUM=" << DATUM << "] = " << A_star[DATUM] << endl;
//cout << "grad_mass_sum[iteration=" << continuity_counter << "] = " << grad_mass_sum << endl;
//cout << "mass_sum[iteration=" << continuity_counter << "] = " << mass_sum << endl;
cout << "factor = " << factor << endl;
cout << endl;
}
//*/

//			}

		}
		//else factor = 1;  // Reset once there is progress without a sign change

		// Check for continuity convergence or apply change to A_star[DATUM] calculated by N-R method
		// ----------------------------------------------------------------------------------------------------	
		if(fabs(mass_sum) < tol_mass_sum) CONTINUITY_CONVERGED = true;
		else
		{
			A_star_prev = A_star[DATUM];						// Save previous value
			A_star[DATUM] -= /*sign**/factor*(mass_sum/grad_mass_sum);	// Determine A_star[DATUM] closer to solution by N-R method
		}

		++continuity_counter;

if(continuity_counter>1000)
{
cout << "continuity_counter = " << continuity_counter << endl;
cout << "A_star[DATUM] = " << A_star[DATUM] << endl;
cout << "A_star_prev = " << A_star_prev << endl;
cout << "mass_sum = " << mass_sum << endl;
cout << endl;
}

/*
if(timestep>=18950)
{
	cout << "\t\t\t5: bottom of BWPLJContinuityNR\n";
	cout << "\t\t\t5: continuity_counter = " << continuity_counter << endl;
	cout << "\t\t\t5: mass_sum = " << mass_sum << endl;
	cout << "\t\t\t5: tol_mass_sum = " << tol_mass_sum << endl;
	cout << endl;
	cout << endl;
}
*/
	}while(!CONTINUITY_CONVERGED);

//cout << endl;
//if(continuity_counter>10) cout << "continuity_counter = " << continuity_counter << endl;
//cout << "continuity converged!" << endl;

/*	
		for(j=0; j<NPIPES; ++j)
		{
			cout << "mass_sum[j=" << j << "] = " << ( pow(A_star[DATUM] - del_star[j], 2/(pPpt->gammaAir() - 1))/AA[j] )
							* (lambda_in_star[j] - A_star[DATUM] + del_star[j]) 
							* Fb[j] << endl;
		}
*/
}

void CJunction::BWPLJContinuityFP(CProperties* pPpt, int timestep)
{
	int j;
	// Solve continuity equation for A_star[DATUM] by false position
	// ----------------------------------------------------------------------------------------------------
	double mass_sum, mass_sum_prev, mass_sum_prev_prev;//, grad_mass_sum, grad_mass_sum_prev;
	double A_star_prev;//, factor;
	bool CONTINUITY_CONVERGED = false;
	int continuity_counter = 0;

	A_star[DATUM] = 1.5;	// Starting value
	double step_size = 0.1;	// Initial step size
	double sign = 1.0;		// Step direction

/*
	if(timestep == 100)
	{
	double remember_A_star = A_star[DATUM];
	int a_max = 1000;
	double min_A_star = 0;
	double max_A_star = 2.5;
	for(int a=0; a<a_max; ++a)
	{
		mass_sum = 0;
		A_star[DATUM] = min_A_star + a*((max_A_star - min_A_star)/a_max);

		for(j=0; j<NPIPES; ++j)
		{
			// Function value at x = A_star[DATUM]
			mass_sum += ( pow(A_star[DATUM] - del_star[j], 2/(pPpt->gammaAir() - 1))/AA[j] )
						* (lambda_in_star[j] - A_star[DATUM] + del_star[j]) 
						* Fb[j];
		}
		pPpt->Out(A_star[DATUM]); pPpt->Out("\t"); pPpt->Out(mass_sum); pPpt->Out("\n");
	}
	A_star[DATUM] = remember_A_star;
exit(1);
	}
//*/


	mass_sum_prev = 1000;
	mass_sum = 1000;
	//grad_mass_sum = 0;
	//A_star_prev = A_star[DATUM];
	//factor = 1;
	do{
		mass_sum_prev_prev = mass_sum_prev; // Save previous value
		mass_sum_prev = mass_sum;			// Save previous value
		mass_sum = 0;						// Reset
		//grad_mass_sum_prev = grad_mass_sum; // Save previous value
		//grad_mass_sum = 0;					// Reset
		for(j=0; j<NPIPES; ++j)
		{
			// Function value at x = A_star[DATUM]
			mass_sum += ( pow(A_star[DATUM] - del_star[j], 2/(pPpt->gammaAir(pBN[j]->T) - 1))/AA[j] )
						* (lambda_in_star[j] - A_star[DATUM] + del_star[j]) 
						* Fb[j];		

//cout << "mass_sum up to [j=" << j << "] = " << mass_sum << endl;
//cout << "del_star[j=" << j << "] = " << del_star[j] << endl;
//cout << "AA[j=" << j << "] = " << AA[j] << endl;
//cout << "del_star[j=" << j << "] = " << del_star[j] << endl;
//cout << "lambda_in_star[j=" << j << "] = " << lambda_in_star[j] << endl;
		}
//cout << endl;

/*	
if(timestep >= 333)
//if(this->ID==0)
{
cout << endl;
cout << "continuity_counter = " << continuity_counter << endl;
cout << "A_star[DATUM=" << DATUM << "] = " << A_star[DATUM] << endl;
//cout << "grad_mass_sum[iteration=" << continuity_counter << "] = " << grad_mass_sum << endl;
cout << "AA[j=" << j << "] = " << AA[j] << endl;
cout << "mass_sum[iteration=" << continuity_counter << "] = " << mass_sum << endl;
cout << "sign*step_size = " << sign*step_size << endl;
}
//exit(1);
//*/
/*
		// Check for non-improvement in mass sum; if there none, go back and try again with smaller step
		// ----------------------------------------------------------------------------------------------------	
		if(mass_sum==mass_sum_prev)		
		{
cout << "NO CHANGE in mass sum; restore previous values and reduce step size" << endl;
cout << "step_size = " << step_size << endl;
cout << endl;
			A_star[DATUM] = A_star_prev;		// Restore previous value
			mass_sum = mass_sum_prev;			// Restore previous value	
			step_size *= 3.;					// Halve the step size
		}
*/

		// Check for sign change in mass sum; if there is, go back and try again with smaller step
		// ----------------------------------------------------------------------------------------------------	
		if(mass_sum/mass_sum_prev < 0 && continuity_counter>0)			// There has been a sign change
		{
			// Sign changes indicates solution lies between A_star[DATUM] and A_star_previous
			if(fabs(mass_sum) < fabs(mass_sum_prev))
			{
//cout << "SIGN CHANGE in mass sum; ALLOW!" << endl; // Allow sign change since result is closest to solution so far
				sign *= -1;			// Change sign since we are on other side of solution and must step in other direction
				step_size *= 0.5;				// Also halve factor else will simply return to previous A_star[DATUM]
			}
			else	// Disallows sign change as result is further away from solution
			{
/*
if(timestep >= 333)
//if(this->ID==0)
{
cout << "SIGN CHANGE in mass sum; restore previous values and reduce step size" << endl;
cout << "Current values:" << endl;
cout << "A_star[DATUM=" << DATUM << "] = " << A_star[DATUM] << endl;
//cout << "grad_mass_sum[iteration=" << continuity_counter << "] = " << grad_mass_sum << endl;
cout << "mass_sum[iteration=" << continuity_counter << "] = " << mass_sum << endl;
}
//*/
				A_star[DATUM] = A_star_prev;		// Restore previous value
				mass_sum = mass_sum_prev;			// Restore previous value	
				step_size *= 0.5;					// Halve the step size

/*
if(timestep >= 333)
//if(this->ID==0)
{
cout << "Restored values:" << endl;
cout << "A_star[DATUM=" << DATUM << "] = " << A_star[DATUM] << endl;
//cout << "grad_mass_sum[iteration=" << continuity_counter << "] = " << grad_mass_sum << endl;
cout << "mass_sum[iteration=" << continuity_counter << "] = " << mass_sum << endl;
//cout << "factor = " << factor << endl;
//cout << "step_size = " << step_size << endl;
cout << "sign*step_size = " << sign*step_size << endl;
cout << endl;
cout << endl;
}
//*/
			}
		}
		//else factor = 1;  // Reset once there is progress without a sign change

		// Check for continuity convergence or apply change to A_star[DATUM] calculated by N-R method
		// ----------------------------------------------------------------------------------------------------	
		if(fabs(mass_sum) < tol_mass_sum) CONTINUITY_CONVERGED = true;
		else
		{
			A_star_prev = A_star[DATUM];						// Save previous value
			A_star[DATUM] -= sign*step_size;					// Determine A_star[DATUM] closer to solution by halving step size
		}

		++continuity_counter;
/*
if(continuity_counter>10)
{
//	cout << "timestep = " << timestep << endl;
	cout << "continuity_counter = " << continuity_counter << endl;
//	exit(1);
}
//*/

//cout << "mass_sum = " << mass_sum << endl;
//cout << "continuity_counter = " << continuity_counter << endl;
	}while(!CONTINUITY_CONVERGED);
//cout << endl;

//cout << "continuity_counter = " << continuity_counter << endl;


/*	
		for(j=0; j<NPIPES; ++j)
		{
			cout << "mass_sum[j=" << j << "] = " << ( pow(A_star[DATUM] - del_star[j], 2/(pPpt->gammaAir() - 1))/AA[j] )
							* (lambda_in_star[j] - A_star[DATUM] + del_star[j]) 
							* Fb[j] << endl;
		}
*/
}	

void CJunction::BWPLJLossTerms(CProperties* pPpt, int timestep, int counter, bool &rESCAPE)
{
	int j;
	double k;
	double theta, qpsi;//, C;
	double C_prev;

	// Identify datum branch, i.e., that with greatest mass flow TOWARDS junction
	// ----------------------------------------------------------------------------------------------------
	//int DATUM_PREV = DATUM; // Save number of previous datum branch; if first time step, DATUM will have been set to [0]
	double temp_mfr, greatest_positive_mfr; greatest_positive_mfr = 0;
	for(j=0; j<NPIPES; ++j)
	{
		k = pPpt->gammaAir(pBN[j]->T);
		temp_mfr = (pow(A_star[j], 2/(k - 1))/AA[j])*U_star[j]*Fb[j];
		if(temp_mfr > greatest_positive_mfr)
		{
			greatest_positive_mfr = temp_mfr;
			DATUM = j;
		}
	}

	// Determine pressure loss terms
	// ----------------------------------------------------------------------------------------------------
	for(j=0; j<NPIPES; ++j)
	{
		k = pPpt->gammaAir(pBN[j]->T);
		C_prev = C[j];

		// Obtain appropriate angle (in radians) between datum and current branch j
		theta = ref_datum_branch[DATUM][j];

//cout << "j = " << j << endl;
//cout << "id_branch[j=" << j << "] = " << id_branch[j] << endl;
//cout << "theta = ref_datum_branch[DATUM=" << DATUM << "=[" << id_branch[DATUM] << "]][j=" << j << "=[" << id_branch[j] << "]] = " << theta*180/PI << endl;
//cout << endl;

		if(U_star[j] > 0)
		{
			del_star[j] = 0;	// No pressure loss if joining (INFLOW) branch
		}
		else // Separating (OUTFLOW) branch
		{
			if(U_star[DATUM]==0) // If this is zero, then all flows zero
			{
				qpsi = 1;
				C[j] = 1 - (1/(qpsi))*cos(0.75*(PI - theta));
			}
			else
			{
				//qpsi = (U_star[j]*AA[j])/(U_star[DATUM]*AA[DATUM]);	
				qpsi = fabs((U_star[j]*AA[j])/(U_star[DATUM]*AA[DATUM]));	
				// q_j = m_dot_j/m_dot_dat = (rho*u_j*F_j)/(rho*u_dat*F_dat)
				// psi_j = F_dat/F_j
				// Hence for constant rho, q_j*psi_j = u_j/u_dat = (U_star_j*AA_j)/(U_star_dat*AA_dat)
				
				if(qpsi == 0) C[j] = 1;
				else C[j] = 1 - (1/(qpsi))*cos(0.75*(PI - theta));
			}

			// Relax C
      //cout << "counter = " << counter << endl;
			C[j] = ((loop_limit_main - counter)/loop_limit_main)*C[j] + (counter/loop_limit_main)*C_prev;

			// Before calculating del_star make sure pPpt->gammaAir(pBN[j]->T)*C[j]*pow(U_star[j]/A_star[j], 2) >= -1 else will get sqrt(-ve)
			if(k*C[j]*pow(U_star[j]/A_star[j], 2) < -1)
			{
				pPpt->Out(Identify()); pPpt->Out("BWPLJLossTerms: qpsi = "); pPpt->Out(qpsi); pPpt->Out(" at timestep = "); pPpt->Out(timestep); pPpt->Out("\n\n");
				
				pPpt->Out(Identify()); pPpt->Out("BWPLJLossTerms: original C[j="); pPpt->Out(j); pPpt->Out("] = "); pPpt->Out(C[j]); pPpt->Out(" at timestep = "); pPpt->Out(timestep); pPpt->Out("\n\n");
				
				C[j] = -1/(k*pow(U_star[j]/A_star[j], 2)); // Recalculate C[j] to ensure calculation proceeds 
				//if(timestep>=18950)	
				pPpt->Out(Identify()); pPpt->Out("BWPLJLossTerms: fixing C[j="); pPpt->Out(j); pPpt->Out("] = "); pPpt->Out(C[j]); pPpt->Out(" at timestep = "); pPpt->Out(timestep); pPpt->Out("\n\n");
				//rESCAPE = true;
			}

			del_star[j] = A_star[j]*
				( 
					pow(1 + k*C[j]*pow(U_star[j]/A_star[j], 2), (k - 1)/(2*k)) 
			- 1 );
/*
if(timestep>=13043 || counter>=100)
{
	char pause;
	//cin >> pause;
	cout << endl;
	cout << "\t\t\t5: C[j=" << j << "] = " << C[j] << endl;
	cout << "\t\t\t5: U_star[j=" << j << "] = " << U_star[j] << endl;
	cout << "\t\t\t5: A_star[j=" << j << "] = " << A_star[j] << endl;
	//cout << "\t\t\t5: pow(U_star[j]/A_star[j], 2) = " << pow(U_star[j]/A_star[j], 2) << endl;
	cout << "\t\t\t5: (1 + ... ) = " << (1 + pPpt->gammaAir(pBN[j]->T)*C[j]*pow(U_star[j]/A_star[j], 2)) << endl;
	cout << "\t\t\t5: exp = " << (pPpt->gammaAir(pBN[j]->T) - 1)/(2*pPpt->gammaAir(pBN[j]->T)) << endl;
	cout << "\t\t\t5: pow(1 + ... ) = " << pow(1 + pPpt->gammaAir(pBN[j]->T)*C[j]*pow(U_star[j]/A_star[j], 2), (pPpt->gammaAir(pBN[j]->T) - 1)/(2*pPpt->gammaAir(pBN[j]->T))) << endl;
	cout << "\t\t\t5: GIVING del_star[j=" << j << "] = " << del_star[j] << endl;
	cout << endl;
	cout << endl;
}
//*/
		}

//del_star[j] = 0.0;

/*
if(timestep >= 600)
{
//cout << "j = " << j << endl;
cout << "id_branch[j=" << j << "] = " << id_branch[j] << endl;
cout << "theta = ref_datum_branch[DATUM=" << DATUM << "=[" << id_branch[DATUM] << "]][j=" << j << "=[" << id_branch[j] << "]] = " << theta*180/PI << endl;
//cout << endl;
//cout << "theta = ref_datum_branch[DATUM=" << DATUM << "][j=" << j << "] = " << theta << endl;
cout << "cos((3/4)*(PI - theta)) = " << cos(0.75*(PI - theta)) << endl;
cout << "qpsi[j=" << j << "] = " << qpsi << endl;
cout << "C[j=" << j << "] = " << C[j] << endl;				
cout << "del_star[j=" << j << "] = " << del_star[j] << endl;
cout << endl;
}
//*/

	}
//cout << endl;
}
/*
void CJunction::NHPLTJ(CProperties* pPpt, double time, bool &rRESTORE)
//--------------------------------------------------//
// Non-homentropic, pressure loss 'T' junction      //
// -------------------------------------------		//		
// Three-way branch with pressure loss,				//
// 'equal-area', 'T' junction (right-angle).		//
// Based on algorithm by Benson.					//
//--------------------------------------------------//
{
	if(NPIPES!=3){cout << "Attempting to apply NHPLTJ to a non 3-way junction!\n"; exit(1);}
	int p;

	int counter;
	bool first_loop_local = true;
	bool stop;

	// Newton-Raphson variables
	// ------------------------
	int x_counter;
	bool converged;
	double* sum; sum = new double [3];
	double* grad; grad = new double [3];
	double sum_old;
	double factor;
	bool changed_sign;

	// Declare local working arrays
	// ----------------------------

	// Dimension local working arrays
	// ------------------------------

	// Store uncorrected values and intialise working variables
	// --------------------------------------------------------
	CommonSetup();
	flow_type_winterbone_orig = flow_type_winterbone; // Record flowtype upon entering

	// Main loop
	// ==================================================
	stop = false;
	counter = 0;
	do
	{
		++ counter;
		// Continuity loop
		// --------------------------------------------------
		// Obtain A_star[pipe1] via continuity - eq. 8.124 of Benson p430
		CommonContinuityFP(pPpt, rRESTORE);
		// --------------------------------------------------

		// Continuity has converged, so have a value for A_star[pipe1]
		for(p=0; p<3; ++p)
		{
			// Calculate A_star[N] - eq. 8.122 of Benson p430
			this->A_star[p] = this->A_star[0] - this->del_A_star[p];
			// Calculate U_star[N] - eq. 8.110 of Benson p427
			U_star[p] = (2/(pPpt->gammaAir(pBN[p]->T)-1))*(lambda_in_star[p] - this->A_star[p]);
		}

		// Test for overall convergence
		// ----------------------------
		if(	fabs(this->del_A_star_old[pipe2] - this->del_A_star[pipe2])<tol_main &&
			fabs(this->del_A_star_old[pipe3] - this->del_A_star[pipe3])<tol_main
			&& !first_loop_local)
		{
			stop=true;
			CommonLambdaOut(pPpt); // Now calculate lambda_out
		}
		else
		{
			first_loop_local = false;

			// Determine flow type and common branch
			CommonFlowType();

			// Given flow type set up loss parameters
			CommonSetBranches(pPpt, flow_type_winterbone);

			// Evaluate mass flow ratios and common branch Mach no.
			CommonMassFlowMachNo(pPpt);

			// Calculate A_star[N]/A_star1 - eq. 8.119 of Benson p429
			// Uses Newton-Raphson method
			// --------------------------------------------------
			// Newton-Raphson convergence upon A_star[pipe2, pipe3]/A_star[pipe1] variables
			changed_sign = false;
			for(p=1; p<=2; ++p)
			{
				converged = false;
				x_counter = 0;
				x_star[p] = 1.25;
				factor = 1;
				do
				{
					++ x_counter;
					x_star_old[p] = x_star[p];

					sum[p] = pow(x_star[p], 2*pPpt->gammaAir(pBN[p]->T)/(pPpt->gammaAir(pBN[p]->T)-1)) + this->G1[p]*pow(x_star[p], 2/(pPpt->gammaAir(pBN[p]->T)-1)) - this->G2[p];
					sum_old = sum[p];
					grad[p] = (2*pPpt->gammaAir(pBN[p]->T))/(pPpt->gammaAir(pBN[p]->T)-1)*pow(x_star[p], (pPpt->gammaAir(pBN[p]->T)+1)/(pPpt->gammaAir(pBN[p]->T)-1)) + (2/(pPpt->gammaAir(pBN[p]->T)-1))*this->G1[p]*pow(x_star[p], (3-pPpt->gammaAir(pBN[p]->T))/(pPpt->gammaAir(pBN[p]->T)-1));
					x_star[p] = x_star[p] - factor*(sum[p]/grad[p]);
					sum[p] = pow(x_star[p], 2*pPpt->gammaAir(pBN[p]->T)/(pPpt->gammaAir(pBN[p]->T)-1)) + this->G1[p]*pow(x_star[p], 2/(pPpt->gammaAir(pBN[p]->T)-1)) - this->G2[p]; // Get new sum for new x
					// If sign of sum has changed, reduce the distance moved by half
					if( (sum[p]>0 && sum_old<0) || (sum[p]<0 && sum_old>0) ) changed_sign = true;
					if(changed_sign)
					{ 
						factor *= 0.5;
						changed_sign = false;
					}
					if( fabs(sum[p]) < 1e-12 && fabs(x_star[p] - x_star_old[p]) < 1e-12 ) converged = true;
					else converged = false;
				}while(!converged);
			}
			// A_star_2/A_star_1 and A_star_3/A_star_1 are then stored in x_star[pipe2] and [pipe3]

			// Now calculate del_A_star[N] (8.123) N = 2,3
			for(p=1; p<=2; ++p)
			{
				this->del_A_star_old[p] = this->del_A_star[p];
				this->del_A_star[p] = A_star[pipe1]*(1 - x_star[p]);
			}

			// --------------------------------------------------------------------
			// Calculate new entropy levels depending on joining or separating flow
			CommonEntropyLevels(pPpt, flow_type_benson); 
			// Though flow_type_benson & flow_type_winterbone differ, they fall into the same
			// statements for joining or separating flow. The com pipe does matter though.
			// --------------------------------------------------------------------

			// Apply the correction to lambda_in_star
			CommonLambdaInStarCorrection();
		}
	}while(!stop);

	// Converged, pass back the new values by pointers
	// -----------------------------------------------
//	CommonUpdate(pPpt);

	//////////////////////////////////////////
	double *lambda_in_c; lambda_in_c = new double [NPIPES];
	for(p=0; p<NPIPES; ++p)	lambda_in_c[p] = lambda_in_star[p]*AA[p];

	bool* CHOKED; CHOKED = new bool [NPIPES];
	for(p=0; p<NPIPES; ++p) CHOKED[p] = false;

	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out, AA, pipe_flow, CHOKED);
	//////////////////////////////////////////

	// Delete any local arrays
	// -----------------------
	delete [] sum;
	delete [] grad;
}
*/
/*
void CJunction::NHPLJ(CProperties* pPpt, double time, bool &rRESTORE, int timestep)
//--------------------------------------------------//
// Non-homentropic, pressure loss junction			//
// ---------------------------------------			//		
// Three-way branch with pressure loss,				//
// and user-specified branch angle.					//
//													//
// Algorithm structure derived from Benson's		//
// 'T' junction (see above) coupled with			//
// loss coefficient equations of Winterbone.		//
//--------------------------------------------------//
{
	if(NPIPES!=3){cout << "Attempting to apply NHPLJ to a non 3-way junction!\n"; exit(1);}
	int p;

	int counter;
	bool first_loop_local = true;
	bool stop;

	// Declare local working arrays
	// ----------------------------

	// Dimension local working arrays
	//-------------------------------

	// Store uncorrected values and intialise working variables
	// --------------------------------------------------------
	CommonSetup();
	flow_type_winterbone_orig = flow_type_winterbone; // Record flowtype upon entering

	error1 = 1e9; error2 = 1e9;
	error1_old = error1; error2_old = error2;

	bool SUCCESS = false;

	// Main loop
	// ==================================================
	stop = false;
	counter = 0;
	bool CONVERGING = false;
	do
	{
		++ counter;

		// Test for exceeding maximum number of loops to be carried out
		if((counter>loop_limit_main && !CONVERGING) || counter>1000)
		{
			// Always stop and restore whenever loop limit exceeded
			stop=true;

			if(CONVERGING && error1<(tol_main*10) && error2<(tol_main*10))
				SUCCESS=true;
			else
				SUCCESS=false;
		}
		else
		{
			// Continuity loop
			// --------------------------------------------------
			// Obtain A_star[pipe1] via continuity - eq. 8.124 of Benson p430
			CommonContinuityFP(pPpt, rRESTORE);
			// --------------------------------------------------
				
			// Continuity has converged, so have a value for A_star[pipe1]
			for(p=0; p<3; ++p)
			{
				// Calculate A_star[N] - eq. 8.122 of Benson p430
				A_star[p] = A_star[0] - del_A_star[p];
				// Calculate U_star[N] - eq. 8.110 of Benson p427
				U_star[p] = (2/(pPpt->gammaAir(pBN[p]->T)-1))*(lambda_in_star[p] - A_star[p]);
			}

			if(_isnan(U_star[0])||_isnan(U_star[1])||_isnan(U_star[2])) // Tests for -1#IND errors
			{
				stop = true;
				SUCCESS = false;
			}
			else
			{
				// Test for overall convergence
				// ----------------------------
				error1_old = error1; error2_old = error2;
				error1 = fabs(this->del_A_star_old[pipe2] - this->del_A_star[pipe2]);
				error2 = fabs(this->del_A_star_old[pipe3] - this->del_A_star[pipe3]);

				if(error1<=error1_old && error2<=error2_old) CONVERGING = true;
				else CONVERGING = false;

				if(error1<tol_main && error2<tol_main && !first_loop_local)
				{
					// Function has converged successfully
					stop = true;									
					SUCCESS = true;
				}
				else
				{
					first_loop_local = false;

					// Determine flow type
					// -------------------

					if(FIX_FLOW_TYPE && OPPOSITE && !OTHER && !ZERO_FLOW) 
					// Having tried switching try opposite before trying other flowtypes
					{
						// Find the low flow pipe, and set the flowtype that would occur if flow in the
						// low flow branch was reversed (the "opposite" flowtype)
						
						if(false//counter<10
								)
						{
						cout << endl;
						cout << "timestep = " << timestep << endl;
						cout << "counter = " << counter << endl;
						cout << "TRY_OPPOSITE:\n";
						if(CONVERGING) cout << "Converging:\n"; else cout << "Not converging:\n";
						cout << "error1_old\t= " << error1_old << "   \terror2_old\t= " << error2_old << endl;
						cout << "error1\t\t= " << error1 << "\terror2\t\t= " << error2 << endl << endl;

						cout << "U_star[0] = " << U_star[0] << "; U_star[1] = " << U_star[1] << "; U_star[2] = " << U_star[2] << endl;
						cout << "flow_type_winterbone = " << flow_type_winterbone << endl;
						cout << "flow_type_actual = " << flow_type_actual << endl;
						}

						if(FIRST_OPPOSITE) // Must do this the first time only, else alternates
						{
							FIRST_OPPOSITE = false;

							// Determine the low flow branch based on U_star values
							CommonLowFlowBranch();

							// Set the opposite according to the initial fresh flow_type_actual and the low flow branch
							if(opposite[flow_type_actual][low_flow_branch] >=0 && opposite[flow_type_actual][low_flow_branch] <=6)
							{
								// Set the flowtype to the "opposite" type for this actual flowtype low_flow_branch
								flow_type_winterbone = opposite[flow_type_actual][low_flow_branch];
							}
							else
							{
								// Use next_lowest_branch since low_flow_branch "opposite" would give rise to impossible flowtype
								flow_type_winterbone = opposite[flow_type_actual][next_lowest_branch];
							}
						}
						
						if(false //counter<10
							)
						{
						cout << "low_flow_branch = " << low_flow_branch << endl;
						cout << "next_lowest_branch = " << next_lowest_branch << endl;
						cout << "Setting flowtype to " << flow_type_winterbone << endl;
						cin >> pause;
						}
					}
					else
					{
						if(FIX_FLOW_TYPE && OPPOSITE && OTHER && !ZERO_FLOW) 
						// Have tried switching & opposite already, try other flowtypes
						{
							// Set other flowtypes
							if(false
								//counter<10
								)
							{
							cout << endl;
							cout << "timestep = " << timestep << endl;
							cout << "counter = " << counter << endl;
							cout << "TRY_OTHER:\n";
							if(CONVERGING) cout << "Converging:\n"; else cout << "Not converging:\n";
							cout << "error1_old\t= " << error1_old << "   \terror2_old\t= " << error2_old << endl;
							cout << "error1\t\t= " << error1 << "\terror2\t\t= " << error2 << endl << endl;

							cout << "flow_type_winterbone_orig = " << flow_type_winterbone_orig << endl;
							cout << "flow_type_winterbone = " << flow_type_winterbone << endl;
							cout << "flow_type_actual = " << flow_type_actual << endl;
							
							//cout << "U_star[0] = " << U_star[0] << "; U_star[1] = " << U_star[1] << "; U_star[2] = " << U_star[2] << endl;
							//cout << "low_flow_branch = " << low_flow_branch << endl;
							//cout << "next_lowest_branch = " << next_lowest_branch << endl;
							//cout << "Setting flowtype to " << flow_type_winterbone << endl;
							
							cin >> pause;
							}
						}
						else
						{
							if(FIX_FLOW_TYPE && OPPOSITE && OTHER && ZERO_FLOW) 
							// Have tried switching, opposite and other already, try setting zero flow
							{
								// Run with zero flow with the actual flowtype

								// Determine the low flow branch based on U_star values
								CommonLowFlowBranch();
				
								lambda_in_star[low_flow_branch] = A_star[low_flow_branch]; // Will cause U_star[p] = 0
								U_star[low_flow_branch] = (2/(pPpt->gammaAir(pBN[low_flow_branch]->T)-1))*(lambda_in_star[low_flow_branch] - A_star[low_flow_branch]);
		
								// Then with zero flow and the opposite flowtype
								// After one go with SET_ZERO_FLOW, start trying other flowtypes
								if(zero_counter>1 && zero_counter<=7)
								{
									flow_type_winterbone = zero_counter - 1;
									cout << endl;
									cout << "timestep = " << timestep << endl;
									cout << "ZERO_FLOW: setting artificial flowtype " << flow_type_winterbone << endl;
								}
								else
									CommonFlowType();

								cout << endl;
								cout << "timestep = " << timestep << endl;
								cout << "ZERO_FLOW:\n";
								
								//cout << "U_star[0] = " << U_star[0] << "; U_star[1] = " << U_star[1] << "; U_star[2] = " << U_star[2] << endl;
								//cout << "low_flow_branch = " << low_flow_branch << endl;
								//cout << "next_lowest_branch = " << next_lowest_branch << endl;
								//cout << "Setting flowtype to " << flow_type_winterbone << endl;
								
						//		cin >> pause;
							}
							else
							{
								if(FIX_FLOW_TYPE && !OPPOSITE && !OTHER && !ZERO_FLOW)
								{
									// Has been restored, try switching prevention from the start
									
									// Fix the flowtype by setting it to that used the previous iteration
									flow_type_winterbone = flow_type_winterbone_orig;
									if(false
										//counter<10
										)
									{
									cout << endl;
									cout << "timestep = " << timestep << endl;
									cout << "counter = " << counter << endl;
									cout << "SWITCHING = true\n";
									cout << "TRY_SWITCHING\n";
									if(CONVERGING) cout << "Converging:\n"; else cout << "Not converging:\n";
									cout << "error1_old\t= " << error1_old << "   \terror2_old\t= " << error2_old << endl;
									cout << "error1\t\t= " << error1 << "\terror2\t\t= " << error2 << endl << endl;

									cout << "flow_type_winterbone_orig = " << flow_type_winterbone_orig << endl;
									cout << "flow_type_winterbone = " << flow_type_winterbone << endl;
									cout << "flow_type_actual = " << flow_type_actual << endl;
									cout << "Use flow_type_winterbone_orig; setting flowtype to " << flow_type_winterbone << endl;
									cin >> pause;
									}
								}
								else
								{
									if(!FIX_FLOW_TYPE && !OPPOSITE && !OTHER && !ZERO_FLOW)
									{
										if(!ALTERNATING)
										{
											// Will come here first at the start of every fresh iteration (non-restored)

											// The normal operation
											CommonFlowType();

											// Record actual flowtype in case of restores
											flow_type_actual = flow_type_winterbone;

											// Test for switching
											if(flow_type_winterbone!=flow_type_winterbone_old)
											{
												++switch_counter;
												if(switch_counter>loop_limit_switch)
												{
													ALTERNATING=true;
												
													// Determine which flowtype is being switched to
													if(flow_type_winterbone!=flow_type_winterbone_orig)
													{
														flow_type_winterbone_next = flow_type_winterbone;
													}
													else flow_type_winterbone_next = flow_type_winterbone_old;
												}
											}
										}
										else
										{
											// Prevent switching within the first ever loop

											// If switching fix the flowtype by setting it to that used the previous iteration
											flow_type_winterbone = flow_type_winterbone_orig;

											if(false
												//counter<10
												)
											{
											cout << endl;
											cout << "timestep = " << timestep << endl;
											cout << "counter = " << counter << endl;
											cout << "SWITCHING = true\n";
											cout << "not TRY_SWITCHING yet\n";
											if(CONVERGING) cout << "Converging:\n"; else cout << "Not converging:\n";
											cout << "error1_old\t= " << error1_old << "   \terror2_old\t= " << error2_old << endl;
											cout << "error1\t\t= " << error1 << "\terror2\t\t= " << error2 << endl << endl;
											cout << "flow_type_winterbone_orig = " << flow_type_winterbone_orig << endl;
											cout << "flow_type_winterbone = " << flow_type_winterbone << endl;
											cout << "flow_type_actual = " << flow_type_actual << endl;
											cout << "Use flow_type_winterbone_orig; setting flowtype to " << flow_type_winterbone << endl;
											cin >> pause;
											}
										}
									}
									else
									{
										cout << "UNCLEAR WHAT ACTION TO TAKE\n";
										cin >> pause;
									}
								}
							}
						}
					}

					// Given flow type set up loss parameters
					CommonSetBranches(pPpt, flow_type_winterbone);
	
					// Evaluate mass flow ratios and common branch Mach no.
					CommonMassFlowMachNo(pPpt);

					// Calculate x_star[pipe2] and x_star[pipe3] depending upon flow type
					LossCoeffs_Equations(pPpt, flow_type_winterbone); // Use equations rather than empirical data

					// Now calculate del_A_star[N] (8.123) N = 2,3
					for(p=0; p<=2; ++p)
					//for(p=1; p<=2; ++p)
					{
						this->del_A_star_old[p] = this->del_A_star[p];
						this->del_A_star[p] = A_star[pipe1]*(1 - x_star[p]);
					}

					// --------------------------------------------------------------------
					// Calculate new entropy levels depending on joining or separating flow
					CommonEntropyLevels(pPpt, flow_type_winterbone); 
					// Though flow_type_benson & flow_type_winterbone differ, they fall into the same
					// statements for joining or separating flow. The com pipe does matter though.
					// --------------------------------------------------------------------

					// Apply the correction to lambda_in_star
					CommonLambdaInStarCorrection();
				}
			}
		}	
	}while(!stop);

	if(!SUCCESS)
	{
		// Simulation will be restored - now decide appropriate order of action on the replay step
		if(!FIX_FLOW_TYPE && ALTERNATING) // Only try switching solution from start if it wasswitching to begin with
		{
			// Loop limit exceeded during first ever iteration while attempting to converge using flow_type_winterbone_orig only
			// Replay step - will use flow_type_winterbone_orig from the start since SWITCHING flag
			// will remain true (not reset by the Restore())

			FIX_FLOW_TYPE = true;
			OPPOSITE = false;
			OTHER = false;
			ZERO_FLOW = false;
		}
		else
		{
			if(!OPPOSITE)
			{
				// Loop limit exceeded but no switching
				// Replay step with "opposite" flowtype
				FIX_FLOW_TYPE = true;
				OPPOSITE = true; FIRST_OPPOSITE = true;
				OTHER = false;
				ZERO_FLOW = false;
			}
			else
			{
				if(!OTHER)
				{
					// Loop limit exceeded while trying opposite flowtype
					// Start trying other flowtypes
					FIX_FLOW_TYPE = true;
					OPPOSITE = true;
					OTHER = true;
					ZERO_FLOW = false;
				}
				else
				{
					FIX_FLOW_TYPE = true;
					OPPOSITE = true;
					OTHER = true;
					ZERO_FLOW = true;
				}
			}
		}

		rRESTORE=true;


//		this->SET_ZERO_FLOW = true;
//		++this->set_zero_flow_count;

	
//		if(set_zero_flow_count<=7) rRESTORE=true; // Restore only a certain number of times
//		else
//		{
//			cout << "NEVER SUCCESS\n";
//	//		cin >> pause;
//		}

	}
	else
	{
		// Reset flags upon a successful completion
		ALTERNATING = false;
		FIX_FLOW_TYPE = false;
		OPPOSITE = false;
		OTHER = false;
		ZERO_FLOW = false; 
		zero_counter = 0;

		CommonLambdaOut(pPpt); // Now calculate lambda_out

		// Converged, pass back the new values by pointers
		// -----------------------------------------------
//		CommonUpdate(pPpt);

		//////////////////////////////////////////
		double *lambda_in_c; lambda_in_c = new double [NPIPES];
		for(p=0; p<NPIPES; ++p)	lambda_in_c[p] = lambda_in_star[p]*AA[p];

		bool* CHOKED; CHOKED = new bool [NPIPES];
		for(p=0; p<NPIPES; ++p) CHOKED[p] = false;

		common_UPDATE_NH(pPpt, lambda_in_c, lambda_out, AA, pipe_flow, CHOKED);
		//////////////////////////////////////////
	}

	// Delete any local arrays
	// -----------------------
}
*/
/*
void CJunction::NHPLJ2(CProperties* pPpt, double time, bool &rRESTORE, int timestep)
//--------------------------------------------------//
// Non-homentropic, pressure loss junction			//
// ---------------------------------------			//		
// Three-way branch with pressure loss,				//
// and user-specified branch angle.					//
//													//
// Algorithm structure derived from Benson's		//
// 'T' junction (see above) coupled with			//
// loss coefficient equations of Winterbone.		//
//--------------------------------------------------//
{
	if(NPIPES!=3){cout << "Attempting to apply NHPLJ to a non 3-way junction!\n"; exit(1);}
	int p;

	int counter;
	bool first_loop_local = true;
	bool stop;

	// Declare local working arrays
	// ----------------------------
	//none

	// Dimension local working arrays
	//-------------------------------
	//none

	// Store uncorrected values and intialise working variables
	// --------------------------------------------------------
	CommonSetup();
	flow_type_winterbone_orig = flow_type_winterbone; // Record flowtype upon entering

	// Reset errors
	error1 = 1e9; 
	error2 = 1e9;
	error1_old = error1; 
	error2_old = error2;

	bool SUCCESS = false;

	// Main loop
	// ==================================================
	stop = false;
	counter = 0;
	bool CONVERGING = false;
	do
	{
		++ counter;

		// Test for exceeding maximum number of loops to be carried out
		if((counter>loop_limit_main && !CONVERGING) || counter>loop_limit_main*10)//1000)
		{
			// Always stop and restore whenever loop limit exceeded
			stop=true;

			if(CONVERGING && error1<(tol_main*10) && error2<(tol_main*10))
				SUCCESS=true;
			else
				SUCCESS=false;
		}
		else
		{
			// Continuity loop
			// --------------------------------------------------
			// Obtain A_star[pipe1] via continuity - eq. 8.124 of Benson p430
			CommonContinuityFP(pPpt, rRESTORE);
			if(rRESTORE) {stop=true; SUCCESS=false;} // If continuity fails
			// --------------------------------------------------
				
			// Continuity has converged, so have a value for A_star[pipe1]
			for(p=0; p<3; ++p)
			{
				// Calculate A_star[N] - eq. 8.122 of Benson p430
				A_star[p] = A_star[0] - del_A_star[p];
				// Calculate U_star[N] - eq. 8.110 of Benson p427
				U_star[p] = (2/(pPpt->gammaAir(pBN[p]->T)-1))*(lambda_in_star[p] - A_star[p]);
			}

			if(_isnan(U_star[0])||_isnan(U_star[1])||_isnan(U_star[2])) // Tests for -1#IND errors
			{
				stop = true;
				SUCCESS = false;
			}
			else
			{
				// Test for overall convergence
				// ----------------------------
				error1_old = error1; 
				error2_old = error2;
				error1 = fabs(this->del_A_star_old[pipe2] - this->del_A_star[pipe2]);
				error2 = fabs(this->del_A_star_old[pipe3] - this->del_A_star[pipe3]);

				if(error1<=error1_old && error2<=error2_old) CONVERGING = true;
				else CONVERGING = false;

				if(error1<tol_main && error2<tol_main && !first_loop_local)
				{
					// Function has converged successfully
					stop = true;									
					SUCCESS = true;
				}
				else
				{
					first_loop_local = false;

					// Determine flow type
					// -------------------
					if(ZERO_FLOW) // Simulation restored; set a zero flow in the low_flow_branch
					{
						// Run with zero flow with the actual flowtype

						// Determine the low flow branch based on U_star values
						CommonLowFlowBranch();
				
						lambda_in_star[low_flow_branch] = A_star[low_flow_branch]; // Will cause U_star[p] = 0
						U_star[low_flow_branch] = (2/(pPpt->gammaAir(pBN[low_flow_branch]->T)-1))*(lambda_in_star[low_flow_branch] - A_star[low_flow_branch]);
		
						// Then with zero flow and the opposite flowtype
						// After one go with SET_ZERO_FLOW, start trying other flowtypes
						if(zero_counter>1 && zero_counter<=7)
						{
							flow_type_winterbone = zero_counter - 1;
							cout << endl;
							cout << "timestep = " << timestep << endl;
							cout << "ZERO_FLOW: setting artificial flowtype " << flow_type_winterbone << endl;
						}
						else
							CommonFlowType();

						cout << endl;
						cout << "timestep = " << timestep << endl;
						cout << "ZERO_FLOW:\n";
					}
					else
					{
						if(OTHER) // Simulation restored; try all other flow types that haven't been tested yet
						{
							while( (other_counter==flow_type_actual && TESTED_ACTUAL) 
								|| (other_counter==flow_type_winterbone_orig && TESTED_ORIGINAL) && other_counter<6)
							{
								++other_counter;
							}
							
							flow_type_winterbone = other_counter;
						}
						else
						{
							if(OPPOSITE) // Simulation restored; set the type if flow in the low flow branch was reversed (the 'opposite')
							{
								if(FIRST_OPPOSITE) // Do this the first time only else flow type will alternate
								{
									FIRST_OPPOSITE = false;

									// Determine the low flow branch based on U_star values
									CommonLowFlowBranch();

									// Set the 'opposite' flow type according to the non-restored flow_type_actual and the low flow branch
									if(opposite[flow_type_actual][low_flow_branch] >=0 && opposite[flow_type_actual][low_flow_branch] <=6)
									{
										// Set the 'opposite' flow type for this actual flowtype low_flow_branch
										flow_type_winterbone = opposite[flow_type_actual][low_flow_branch];
									}
									else
									{
										// Use next_lowest_branch since low_flow_branch 'opposite' would give rise to impossible situation
										flow_type_winterbone = opposite[flow_type_actual][next_lowest_branch];
									}
								}
							}
							else
							{
								if(FIX_FLOW_TYPE) // Simulation restored; prevents alternating flow types from the very first loop
								{
									if(fix_counter==0)
									{
										// Fix the flow type to that found by the first loop in this iteration
										flow_type_winterbone = flow_type_actual;
										TESTED_ACTUAL = true;
									}
									else
									{
										// Fix the flow type to that found by the solution to the previous iteration
										flow_type_winterbone = flow_type_winterbone_orig;
										TESTED_ORIGINAL = true;
									}
								}
								else
								{
									if(ALTERNATING) // Prevents alternating flow types in the non-restored/first iteration
									{
										// Fix the flow type to that found by the solution to the previous iteration
										flow_type_winterbone = flow_type_winterbone_orig;
									}
									else // The normal operation
									{
										CommonFlowType();

										// Record the actual new flow type in case of restores
										flow_type_actual = flow_type_winterbone;

										// Reset flags
										TESTED_ACTUAL = false;
										TESTED_ORIGINAL = false;

										// Test for alternating flow types
										if(flow_type_winterbone!=flow_type_winterbone_old)
										{
											++switch_counter;
											if(switch_counter>loop_limit_switch)
											{
												ALTERNATING = true;
												
												// Determine the alternating flow types
												if(flow_type_winterbone!=flow_type_winterbone_orig)
													flow_type_winterbone_next = flow_type_winterbone;
												else flow_type_winterbone_next = flow_type_winterbone_old;
											}
										}
									}
								}
							}
						}	
					}			
				
					// Given flow type set up loss parameters
					CommonSetBranches(pPpt, flow_type_winterbone);
	
					// Evaluate mass flow ratios and common branch Mach no.
					CommonMassFlowMachNo(pPpt);

					// Calculate x_star[pipe2] and x_star[pipe3] depending upon flow type
					LossCoeffs_Equations(pPpt, flow_type_winterbone); // Use equations rather than empirical data

					// Now calculate del_A_star[N] (8.123) N = 2,3
					for(p=0; p<=2; ++p)
					//for(p=1; p<=2; ++p)
					{
						this->del_A_star_old[p] = this->del_A_star[p];
						this->del_A_star[p] = A_star[pipe1]*(1 - x_star[p]);
					}

					// --------------------------------------------------------------------
					// Calculate new entropy levels depending on joining or separating flow
					CommonEntropyLevels(pPpt, flow_type_winterbone); 
					// Though flow_type_benson & flow_type_winterbone differ, they fall into the same
					// statements for joining or separating flow. The com pipe does matter though.
					// --------------------------------------------------------------------

					// Apply the correction to lambda_in_star
					CommonLambdaInStarCorrection();
				}
			}
		}	
	}while(!stop);

	if(!SUCCESS) // Simulation will be restored to the solution of the previous iteration
	{
		rRESTORE=true;
		
		if(!FIX_FLOW_TYPE && ALTERNATING) // Fix flow type from first loop only if it was alternating to begin with
		{
			// Loop limit exceeded during non-restored/first iteration under alternating flow types
			FIX_FLOW_TYPE = true;
		}
		else
		{
			if(FIX_FLOW_TYPE && fix_counter==0) ++fix_counter;
			else
			{
				if(!OPPOSITE)
				{
					// Loop limit exceeded but no switching - replay iteration with 'opposite' flow type
					OPPOSITE = true; FIRST_OPPOSITE = true;
				}
				else
				{
					if(!OTHER)
					{
						// Loop limit exceeded while under 'opposite' flow type - try others
						OTHER = true;
						++other_counter; // Will now start at flow type 1
					}
					else
					{
						if(OTHER && other_counter<6)
						{
							++other_counter;
						}
						else
						{
							if(!ZERO_FLOW)
							{
								ZERO_FLOW = true;
								// Start with zero_counter = 0
							}
							else
							{
								++zero_counter;
							}
						}
					}
				}
			}
		}
	}
	else
	{
		if(ZERO_FLOW)
		{
			++zero_counter_total;
		}
		else
		{
			if(OTHER)
			{
				++other_counter_total;
			}
			else
			{
				if(OPPOSITE)
				{
					++opposite_counter_total;
				}
				else
				{
					if(FIX_FLOW_TYPE)
					{
						++fix_counter_total;
					}
					else
					{
						++normal_counter_total;
					}
				}
			}
		}

		// Reset flags upon a successful completion
		ALTERNATING = false;
		FIX_FLOW_TYPE = false;
		OPPOSITE = false;
		OTHER = false;
		ZERO_FLOW = false; 

		fix_counter = 0;
		other_counter = 0;
		zero_counter = 0;

		CommonLambdaOut(pPpt); // Now calculate lambda_out

		// Converged, pass back the new values by pointers
		// -----------------------------------------------
//		CommonUpdate(pPpt);

		//////////////////////////////////////////
		double *lambda_in_c; lambda_in_c = new double [NPIPES];
		for(int p=0; p<NPIPES; ++p)	lambda_in_c[p] = lambda_in_star[p]*AA[p];

		bool* CHOKED; CHOKED = new bool [NPIPES];
		for(p=0; p<NPIPES; ++p) CHOKED[p] = false;

		common_UPDATE_NH(pPpt, lambda_in_c, lambda_out, AA, pipe_flow, CHOKED);
		//////////////////////////////////////////
	}

	// Delete any local arrays
	// -----------------------
	//none
}
*/
/*
void CJunction::NHPLPC(CProperties* pPpt, double time, bool &rRESTORE)
//--------------------------------------------------//
// Non-homentropic, pressure loss pulse converter	//
// ----------------------------------------------	//		
// Three-way branch with pressure loss.				//
// Based on a particular pulse converter geometry.	//
//													//
// Algorithm based on Benson.						//
//                                                  //
// [0] and [1] must be identical, [2] is the		//
// entry to the turbine.							//
//													//	
//  Notation:										//
//													//
//													//
//	[pipe1==0]-----									//
//					>-----[pipe3==2]				//
//	[pipe2==1]-----									//
//													//
//--------------------------------------------------//
{
	if(NPIPES!=3){cout << "Attempting to apply NHPLPC to a non 3-way junction!\n"; exit(1);}
	int p;
	
	// Loop counters
	int counter;
	bool first_loop_local = true;

	// Convergence flags
	bool stop, diverging;
	diverging = false;
	
	int same_counter = 0;

	// Declare local working variables
	double* del_diff;
	double* del_diff_old;

	// Dimension local arrays
	del_diff = new double [3];
	del_diff_old = new double [3];

	// Set large initial values
	del_diff[0] = 1000;
	del_diff[1] = 1000;
	del_diff[2] = 1000;

	// Store uncorrected values and intialise working variables
	CommonSetup();

	// Main loop
	// --------------------------------------------------
	// --------------------------------------------------
	stop = false;
	counter = 0;
	same_counter=0;
	do
	{
		++ counter;
		// --------------------------------------------------
		CommonContinuityFP(pPpt, rRESTORE);// Obtain A_star[pipe1] via continuity - eq. 8.124 of Benson p430 (false position)
		// --------------------------------------------------

		// Continuity has converged, so have a value for A_star[pipe1]
		for(p=0; p<3; ++p)
		{
			// Calculate A_star[N] - eq. 8.122 of Benson p430
			this->A_star[p] = this->A_star[0] - this->del_A_star[p];
			// Calculate U_star[N] - eq. 8.110 of Benson p427
			U_star[p] = (2/(pPpt->gammaAir(pBN[p]->T)-1))*(lambda_in_star[p] - this->A_star[p]);
		}

		// Test for overall convergence
		// --------------------------------------------------
		if(	(fabs(this->del_A_star_old[1] - this->del_A_star[1])<tol_main &&
			 fabs(this->del_A_star_old[2] - this->del_A_star[2])<tol_main &&
			 !first_loop_local) || counter>loop_limit_main 
			// || same_counter>loop_limit_same 
			 || diverging )
		{
			stop=true;
			CommonLambdaOut(pPpt);
		}
		else
		{
			// Determine flow type and common branch
			// Test flow type (Fig. 8.45)
			CommonFlowType();
			CommonSetBranches(pPpt, flow_type_winterbone);
	
			if(flow_type_winterbone != flow_type_winterbone_old)
					same_counter=0; // Reset this counter
			else ++same_counter;

			// Evaluate mass flow ratio and common branch Mach no.
			CommonMassFlowMachNo(pPpt);
			
			// Determine La, Lb by linear interpolation
			// According to the flow type calculate x2, x3 from Eq (8.198) to (8.210)
			LossCoeffs_Interpolation(pPpt, flow_type_winterbone, time);

			// Calculate del_A_star[N] (8.214) N = 2,3
			for(p=pipe2; p<=pipe3; ++p)
			{
				this->del_A_star_old[p] = this->del_A_star[p];
				this->del_A_star[p] = A_star[pipe1]*(1 - x_star[p]);

				del_diff_old[p] = del_diff[p];
				del_diff[p] = fabs(this->del_A_star_old[p] - this->del_A_star[p]);
				
				if( del_diff[p]>del_diff_old[p] && counter>5 )
				{
					cout << "del_diff[" << p << "] = " << del_diff[p] << endl;
					cout << "del_diff_old[" << p << "] = " << del_diff_old[p] << endl;
					cout << "diverging!!!" << endl;
					diverging = true;
				}
				// Stops the loop if one or more del_A_star diffs is diverging after 5 loops
			}

			// --------------------------------------------------------------------
			// Calculate new entropy levels depending on joining or separating flow
			CommonEntropyLevels(pPpt, flow_type_winterbone);
			// --------------------------------------------------------------------

			// Apply the correction to lambda_in_star
			CommonLambdaInStarCorrection();
							
			first_loop_local = false; // Convergence is possible now that at least one main loop is complete
		}
	}while(!stop);

	// Converged, pass back the new values by pointers
	// -----------------------------------------------
//	CommonUpdate(pPpt);

	//////////////////////////////////////////
	double *lambda_in_c; lambda_in_c = new double [NPIPES];
	for(p=0; p<NPIPES; ++p)	lambda_in_c[p] = lambda_in_star[p]*AA[p];

	bool* CHOKED; CHOKED = new bool [NPIPES];
	for(p=0; p<NPIPES; ++p) CHOKED[p] = false;

	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out, AA, pipe_flow, CHOKED);
	//////////////////////////////////////////
	
	// Delete any local arrays
	// -----------------------
	delete [] del_diff;
	delete [] del_diff_old;
}
*/
/*
// --------------------------------------------------
void CJunction::NHPLJ_Winterbone(CProperties* pPpt, double time, bool &rRESTORE)
// --------------------------------------------------
// --------------------------------------------------
// "Non-homentropic, pressure loss junction"
// The most general pressure loss junction
// based on Winterbone's three-way branch with pressure 
// loss for tee and Y junctions at angle theta
//
// Notation:
//
//	[0]-------[2]		[pipe1]-------[pipe3]		
//  theta |				    theta |			
//		  |			OR			  |			
//		  |						  |			
//	     [1]				   [pipe2]
// --------------------------------------------------
// --------------------------------------------------
{
	if(NPIPES!=3){cout << "Attempting to apply NHPLJ to a non 3-way junction!\n"; exit(1);}
	int p;
	
	// Loop counters
	int pressure_counter, entropy_counter, flow_type_loop_counter;
	int loop_limit_entropy = 100;//1000;
	int loop_limit_flow_type = 100;//1000;
	bool first_loop_local = true;

	// Convergence flags
	bool pressure_converged, entropy_converged, flow_type_converged, diverging, entropy_diverging;
	diverging = false;
	entropy_diverging = false;
	
	int same_counter = 0;

	double limit_entropy = 1e-6; // Entropy loop specific to NHPLJ

	// Declare local working variables
	double* AA_old;
	double* del_diff;
	double* del_diff_old;
	double* entropy_diff;
	double* entropy_diff_old;

	// Dimension local arrays
	AA_old = new double [3];
	del_diff = new double [3];
	del_diff_old = new double [3];
	entropy_diff = new double [3];
	entropy_diff_old = new double [3];

	// Set large initial values
	del_diff[0] = 1000;
	del_diff[1] = 1000;
	del_diff[2] = 1000;	
	entropy_diff[0] = 1000;
	entropy_diff[1] = 1000;
	entropy_diff[2] = 1000;

	// Store uncorrected values and intialise working variables
	CommonSetup();

	// Main loop
	// --------------------------------------------------
	// --------------------------------------------------
	pressure_converged = false;
	pressure_counter = 0;
	do // Pressure
	{
		++ pressure_counter;

		if(!first_loop_local) // Only calculate following after having done one loop
		{
			// Evaluate mass flow ratio and common branch Mach no.
			CommonMassFlowMachNo(pPpt);

//			// Interpolate for appropriate values of the loss coefficients
//			// &
//			// Calculate x_star[pipe2] and x_star[pipe3] depending upon flow type
//			LossCoeffs_Equations(pPpt, flow_type_winterbone); // Use equations rather than empirical data

			// Calculate A_star[N]/A_star1 - eq. 8.119 of Benson p429
			// Uses Newton-Raphson method
			// --------------------------------------------------
			// Newton-Raphson convergence upon A_star[pipe2, pipe3]/A_star[pipe1] variables
			int x_counter;
			bool converged;
			double sum [3];
			double grad [3];
			double sum_old;
			double factor;
			bool changed_sign = false;
			for(p=1; p<=2; ++p)
			{
				converged = false;
				x_counter = 0;
				x_star[p] = 1.25;
				factor = 1;
				do
				{
					++ x_counter;
					x_star_old[p] = x_star[p];

					sum[p] = pow(x_star[p], 2*pPpt->gammaAir(pBN[p]->T)/(pPpt->gammaAir(pBN[p]->T)-1)) + this->G1[p]*pow(x_star[p], 2/(pPpt->gammaAir(pBN[p]->T)-1)) - this->G2[p];
					sum_old = sum[p];
					grad[p] = (2*pPpt->gammaAir(pBN[p]->T))/(pPpt->gammaAir(pBN[p]->T)-1)*pow(x_star[p], (pPpt->gammaAir(pBN[p]->T)+1)/(pPpt->gammaAir(pBN[p]->T)-1)) + (2/(pPpt->gammaAir(pBN[p]->T)-1))*this->G1[p]*pow(x_star[p], (3-pPpt->gammaAir(pBN[p]->T))/(pPpt->gammaAir(pBN[p]->T)-1));
					x_star[p] = x_star[p] - factor*(sum[p]/grad[p]);
					sum[p] = pow(x_star[p], 2*pPpt->gammaAir(pBN[p]->T)/(pPpt->gammaAir(pBN[p]->T)-1)) + this->G1[p]*pow(x_star[p], 2/(pPpt->gammaAir(pBN[p]->T)-1)) - this->G2[p]; // Get new sum for new x
					// If sign of sum has changed, reduce the distamce moved by half
					if( (sum[p]>0 && sum_old<0) || (sum[p]<0 && sum_old>0) ) changed_sign = true;
					if(changed_sign)
					{ 
						factor *= 0.5;
						changed_sign = false;
					}
					if( fabs(sum[p]) < 1e-12 && fabs(x_star[p] - x_star_old[p]) < 1e-12 ) converged = true;
					else converged = false;
				}while(!converged);
			}
			// A_star_2/A_star_1 and A_star_3/A_star_1 are then stored in x_star[pipe2] and [pipe3]


			for(p=pipe2; p<=pipe3; ++p)
			{
				// Set previous iteration values of pressure terms equal to new values
				this->del_A_star_old[p] = this->del_A_star[p];
				// Calculate new pressure difference terms del_A_star[pipe2] and del_A_star[pipe3]
				this->del_A_star[p] = A_star[pipe1]*(1 - x_star[p]); // 6.108
				
				del_diff_old[p] = del_diff[p];
				del_diff[p] = fabs(this->del_A_star_old[p] - this->del_A_star[p]);
				
				if( del_diff[p]>del_diff_old[p] && pressure_counter>5 )
				{
					//cout << "del_diff[" << p << "] = " << del_diff[p] << endl;
					//cout << "del_diff_old[" << p << "] = " << del_diff_old[p] << endl;
					//cout << "diverging!!!" << endl;
					diverging = true;
				}
				// Stops the loop if one or more del_A_star diffs is diverging after 5 loops
			}
		}


		entropy_diverging = false;

		entropy_converged = false;
		entropy_counter = 0;
		do // Entropy
		{
			++ entropy_counter;
			cout << "> entropy_counter = " << entropy_counter << endl;
			flow_type_converged = false;
			flow_type_loop_counter = 0;
			do // Flow type
			{
				++ flow_type_loop_counter;
				//cout << ">> flow_type_loop_counter = " << flow_type_loop_counter << endl;

				// --------------------------------------------------
				// Solve continuity using method of false position
				CommonContinuityFP(pPpt, rRESTORE);// Obtain A_star[pipe1] via continuity
				// --------------------------------------------------
		
				// Continuity has converged, so have a value for A_star[pipe1]
				for(p=0; p<3; ++p)
				{
					// Calculate A_star[N] - eq. 8.122 of Benson p430
					this->A_star[p] = A_star[0] - this->del_A_star[p];
					// Calculate U_star[N] - eq. 8.110 of Benson p427
					U_star[p] = (2/(pPpt->gammaAir(pBN[p]->T)-1))*(lambda_in_star[p] - this->A_star[p]);
				}
		
				// Determine flow type and common branch
				CommonFlowType();
				
				// Given flow type set up loss parameters
				CommonSetBranches(pPpt, flow_type_winterbone);

				// --------------------------------------------------------------------
				// Calculate new entropy levels depending on joining or separating flow
				for(p=0; p<3; ++p) AA_old[p] = AA[p]; // Record old entropy levels
				CommonEntropyLevels(pPpt, flow_type_winterbone);
				// --------------------------------------------------------------------
		
				//cout << endl;
				//cout << "Entropy loop: " << entropy_counter << endl;
				for(p=0; p<3; ++p)
				{
					//cout << "fabs(AA_old[p] - AA[p]) = " << fabs(AA_old[p] - AA[p]) << endl;
			
					entropy_diff_old[p] = entropy_diff[p];
					entropy_diff[p] = fabs(AA_old[p] - AA[p]);

					if( entropy_diff[p]>entropy_diff_old[p] && entropy_counter>995 )
					{
						cout << "entropy_diff[" << p << "] = " << entropy_diff[p] << endl;
						cout << "entropy_diff_old[" << p << "] = " << entropy_diff_old[p] << endl;
						cout << "entropy diverging!!!" << endl;
						entropy_diverging = true;
					}
				}

				// Apply the correction to lambda_in_star
				CommonLambdaInStarCorrection();
				
				//if(flow_type_winterbone == flow_type_winterbone_old || 
				//   flow_type_loop_counter > loop_limit_flow_type)
				//{
				//	++same_counter;
				//	flow_type_converged = true;
				//}
				//else
				//{
				//	same_counter=0;
				//	flow_type_converged = false;
				}

//				flow_type_converged = true;
			}
			while(!flow_type_converged);
			//cout << "flow_type loop converged after " << flow_type_loop_counter << " loops" << endl;
	
			// Entropy levels converged at all junction pipe ends?
			if((fabs(AA_old[0] - AA[0]) < limit_entropy &&
				fabs(AA_old[1] - AA[1]) < limit_entropy &&
				fabs(AA_old[2] - AA[2]) < limit_entropy) 
				|| entropy_counter > loop_limit_entropy
				|| entropy_diverging
				) 

				entropy_converged = true;
			
			else entropy_converged = false;
//////////
entropy_converged = true;			
		}
		while(!entropy_converged);// && entropy_counter<loop_entropy);
		//cout << "entropy loop converged after " << entropy_counter << " loops" << endl;

		// Pressure levels (del_A_star) converged at all junction pipes ends?
		if(	(fabs(this->del_A_star_old[pipe2] - this->del_A_star[pipe2])<tol_main &&
			 fabs(this->del_A_star_old[pipe3] - this->del_A_star[pipe3])<tol_main
			 && !first_loop_local)
			 
		//	 || pressure_counter>loop_limit_main 
		//	 || pressure_counter>loop_limit_same 
		//	 || diverging
			 )
		{
			pressure_converged = true;
			CommonLambdaOut(pPpt); // Now calculate lambda_out	
		}
		else pressure_converged = false;

		first_loop_local = false; // Convergence is possible now that at least one main loop is complete	
	}
	while(!pressure_converged);
	//cout << "pressure loop converged after " << pressure_counter << " loops" << endl;

	// Converged, pass back the new values by pointers
	// ------------------------------------------------------------------------
//	CommonUpdate(pPpt);

	//////////////////////////////////////////
	double *lambda_in_c; lambda_in_c = new double [NPIPES];
	for(p=0; p<NPIPES; ++p)	lambda_in_c[p] = lambda_in_star[p]*AA[p];

	bool* CHOKED; CHOKED = new bool [NPIPES];
	for(p=0; p<NPIPES; ++p) CHOKED[p] = false;

	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out, AA, pipe_flow, CHOKED);
	//////////////////////////////////////////

	// Delete any local arrays
	// ------------------------------------------------------------------------
	delete [] AA_old;
	delete [] del_diff;
	delete [] del_diff_old;
	delete [] entropy_diff;
	delete [] entropy_diff_old;
}
*/
void CJunction::CommonSetup()
// -----------------------------------
// Code common to pressure loss models
// -----------------------------------
{
	for(int p=0; p<NPIPES; ++p)
	{
		// Store uncorrected values
		// ------------------------
		// Calculate initial lambda_in_star[N]
		lambda_in_star_n[p] = (*pCLIN[p])[R+1]/pBN[p]->AA[R+1];
		// Enter initial entropy levels, AA_n[N]
		AA_n[p] = pBN[p]->AA[R+1];

		// Initialise working values
		// -------------------------
		lambda_in_star[p] = lambda_in_star_n[p];
		AA[p] = AA_n[p];
	}
}
/*
void CJunction::CommonContinuityFP(CProperties* pPpt, bool &rRESTORE)
// --------------------------------------------------
// Solves continuity equation by method of false
// position
// --------------------------------------------------
{
	int counter;
	bool converged;

	double A_star_pipe1_0, A_star_pipe1_old;
	double f0, fn;
	double cont_sum;

	// Obtain A_star[pipe1] via continuity - eq. 8.124 of Benson p430 (false position)
	A_star_pipe1_0 = 1.25;
	
	f0 =	pow(A_star_pipe1_0, 2/(pPpt->gammaAir(pBN[pipe1]->T)-1))/AA[pipe1] * (lambda_in_star[pipe1] - A_star_pipe1_0)															
			  +	pow(A_star_pipe1_0 - del_A_star[pipe2], 2/(pPpt->gammaAir(pBN[pipe2]->T)-1))/AA[pipe2] * (lambda_in_star[pipe2] - A_star_pipe1_0 + del_A_star[pipe2])*(Fb[pipe2]/Fb[pipe1])
			  +	pow(A_star_pipe1_0 - del_A_star[pipe3], 2/(pPpt->gammaAir(pBN[pipe3]->T)-1))/AA[pipe3] * (lambda_in_star[pipe3] - A_star_pipe1_0 + del_A_star[pipe3])*(Fb[pipe3]/Fb[pipe1]);
		
	counter = 0;
	converged = false;
	do																																		
	{	
		++counter;

		fn = (pow(A_star[pipe1], 2/(pPpt->gammaAir(pBN[pipe1]->T)-1))/AA[pipe1]) * (lambda_in_star[pipe1] - A_star[pipe1])
			+(pow(A_star[pipe1] - del_A_star[pipe2], 2/(pPpt->gammaAir(pBN[pipe2]->T)-1))/AA[pipe2]) * (lambda_in_star[pipe2] - A_star[pipe1] + del_A_star[pipe2])*(Fb[pipe2]/Fb[pipe1])
			+(pow(A_star[pipe1] - del_A_star[pipe3], 2/(pPpt->gammaAir(pBN[pipe3]->T)-1))/AA[pipe3]) * (lambda_in_star[pipe3] - A_star[pipe1] + del_A_star[pipe3])*(Fb[pipe3]/Fb[pipe1]);

		A_star_pipe1_old = A_star[pipe1];
		A_star[pipe1] = A_star[pipe1] - (A_star[pipe1] - A_star_pipe1_0)*(fn/(fn-f0));

//		if(counter>loop_limit_cont-5) cout << "CJunction::CommonContinuityFP() counter = " << counter << endl;

		if(counter>loop_limit_cont)
		{
			converged = true;
			rRESTORE = true;
		}

		if( (fabs(A_star_pipe1_old - this->A_star[pipe1])<tol_cont && fabs(fn)<1e-12) || counter>loop_limit_cont) 
			converged = true;

	}while(!converged);	
		
	// Check continuity (fn is continuity/Fb[pipe1]); both should be very close to 0, depending on tol_cont
	cont_sum = (pow(A_star[pipe1], 2/(pPpt->gammaAir(pBN[pipe1]->T)-1))/AA[pipe1])*(lambda_in_star[pipe1] - A_star[pipe1])*(Fb[pipe1])
				+ (pow(A_star[pipe1] - del_A_star[pipe2], 2/(pPpt->gammaAir(pBN[pipe2]->T)-1))/AA[pipe2])*(lambda_in_star[pipe2] - A_star[pipe1] + del_A_star[pipe2])*(Fb[pipe2])
				+ (pow(A_star[pipe1] - del_A_star[pipe3], 2/(pPpt->gammaAir(pBN[pipe3]->T)-1))/AA[pipe3])*(lambda_in_star[pipe3] - A_star[pipe1] + del_A_star[pipe3])*(Fb[pipe3]);
	
	//cout << "CJunction::CommonContinuityFP() counter = " << counter << ", fn = " << fn << endl;
}
*/
/*
void CJunction::CommonFlowType()
{
	flow_type_winterbone_old = flow_type_winterbone; // Record previous flow type
				
	// Determine flow type
	if(U_star[bpipe1]>0) // I, V or VI
	{
		if(U_star[bpipe2]>0) // VI
		{
			flow_type_winterbone = 5; // JOINING
		}
		else	 // I or V	
		{
			if(U_star[bpipe3]>0)
			{
				flow_type_winterbone = 6; // JOINING
			}
			else
			{
				flow_type_winterbone = 1; // SEPARATING
			}
		}
	}
	else // II, III or IV
	{
		if(U_star[bpipe2]<0) // III
		{
			flow_type_winterbone = 2; // SEPARATING
		}
		else	 // II or IV	
		{
			if(U_star[bpipe3]<0)
			{
				flow_type_winterbone = 3; // SEPARATING
			}
			else
			{
				flow_type_winterbone = 4; // JOINING
			}
		}
	}

	// If flows are 'exactly' 0, the above will still select a finite flow type, hence:
	if(fabs(U_star[bpipe1])<no_flow && fabs(U_star[bpipe2])<no_flow && fabs(U_star[bpipe3])<no_flow)
	{
		flow_type_winterbone = 0;
	}
}
*/
/*
void CJunction::CommonLowFlowBranch()
{
	// Determine which branch has the lowest flow, and that with the second lowest flow
	if(fabs(U_star[bpipe1])<=fabs(U_star[bpipe2])) // 1 or 3
	{
		if(fabs(U_star[bpipe1])<=fabs(U_star[bpipe3]))
		{
			low_flow_branch = bpipe1;
			if(fabs(U_star[bpipe2])<=fabs(U_star[bpipe3]))
			{
				next_lowest_branch = bpipe2;
			}
			else next_lowest_branch = bpipe3;
		}
		else
		{
			low_flow_branch = bpipe3;
			next_lowest_branch = bpipe1;
		}
	}
	else // 2 or 3
	{
		if(fabs(U_star[bpipe2])<=fabs(U_star[bpipe3]))
		{
			low_flow_branch = bpipe2;
			if(fabs(U_star[bpipe1])<=fabs(U_star[bpipe3]))
			{
				next_lowest_branch = bpipe1;
			}
			else next_lowest_branch = bpipe3;
		}
		else
		{
			low_flow_branch = bpipe3;
			next_lowest_branch = bpipe2;
		}
	}
}
*/
/*
void CJunction::CommonFlowTypeOpposite(bool &rRESTORE)
{
	// Determine which branch has the lowest flow
	if(fabs(U_star[bpipe1])<=fabs(U_star[bpipe2])) // 1 or 3
	{
		if(fabs(U_star[bpipe1])<=fabs(U_star[bpipe3]))
		{
			bpipe1;
			switch(flow_type_orig)
			{
			case 1:

		}
		else
		{
			bpipe3;
		}
	}
	else // 2 or 3
	{
		if(fabs(U_star[bpipe2])<=fabs(U_star[bpipe3]))
		{
			bpipe2;
		}
		else
		{
			bpipe3;
		}
	}




	switch(flow_type)
	{
	case 1:
		return 6;
		break;
	case 2:
		return 4;
		break;
	case 3:
		return 4;
		break;
	flow_type_winterbone_old = flow_type_winterbone; // Record previous flow type
				
	// Determine flow type
	if(U_star[bpipe1]>0) // I, V or VI
	{
		if(U_star[bpipe2]>0) // VI
		{
			flow_type_winterbone = 5; // JOINING
		}
		else	 // I or V	
		{
			if(U_star[bpipe3]>0)
			{
				flow_type_winterbone = 6; // JOINING
			}
			else
			{
				flow_type_winterbone = 1; // SEPARATING
			}
		}
	}
	else // II, III or IV
	{
		if(U_star[bpipe2]<0) // III
		{
			flow_type_winterbone = 2; // SEPARATING
		}
		else	 // II or IV	
		{
			if(U_star[bpipe3]<0)
			{
				flow_type_winterbone = 3; // SEPARATING
			}
			else
			{
				flow_type_winterbone = 4; // JOINING
			}
		}
	}

	// If flows are 'exactly' 0, the above will still select a finite flow type, hence:
	if(fabs(U_star[bpipe1])<no_flow && fabs(U_star[bpipe2])<no_flow && fabs(U_star[bpipe3])<no_flow)
	{
		flow_type_winterbone = 0;

		// If we get -1#IND at junction, flowtype will be zero also, so:
//		rRESTORE = true;
	}

	// Test for switching
	if(flow_type_winterbone!=flow_type_winterbone_old)
	{
		++switch_counter;
	}
	else switch_counter=0;
}
*/
/*
void CJunction::CommonSetBranches(CProperties* pPpt, int local_flow_type)
{	
	switch (local_flow_type) // Below is for Winterbone flow type numbers
	{
	case 0:
		com = pipe1; i2 = pipe2; i3 = pipe3; // Arbitrary
		pipe_flow[pipe1] = NOFLOW; pipe_flow[pipe2] = NOFLOW; pipe_flow[pipe3] = NOFLOW;
		coeff[a] = 1; // Arbitrary
		coeff[b] = 2; // Arbitrary
		
		// G loss values only appropriate if NHPLTJ is being used
		if(type==NHPLTJID)
		{
			flow_type_benson = 0; // Winterbone 0 is Benson 0
			G1[bpipe2] = 0;
			G2[bpipe2] = 1;
			G1[bpipe3] = 0;
			G2[bpipe3] = 1;
		}
		break;																											//
	case 1:								
		com = pipe1; i2 = pipe2; i3 = pipe3;
		pipe_flow[pipe1] = OUTFLOW; pipe_flow[pipe2] = INFLOW; pipe_flow[pipe3] = INFLOW;
		coeff[a] = 1;	// Lookup L1 for La/i2
		coeff[b] = 2;	// Lookup L2 for Lb/i3

		// G loss values only appropriate if NHPLTJ is being used
		if(type==NHPLTJID)
		{
			flow_type_benson = 1; // Winterbone 1 is Benson 1
			G1[bpipe2] = C1*pPpt->gammaAir(pBN[bpipe2]->T)*pow(U_star[bpipe2]/A_star[bpipe1],2);
			G2[bpipe2] = C1*pPpt->gammaAir(pBN[bpipe1]->T)*pow(U_star[bpipe1]/A_star[bpipe1],2) + 1;
			G1[bpipe3] = C2*pPpt->gammaAir(pBN[bpipe3]->T)*pow(U_star[bpipe3]/A_star[bpipe1],2);
			G2[bpipe3] = 1;
		}
		break;																																
	case 2:
		com = pipe2; i2 = pipe1; i3 = pipe3;
		pipe_flow[pipe1] = INFLOW; pipe_flow[pipe2] = OUTFLOW; pipe_flow[pipe3] = INFLOW;
		coeff[a] = 1;	// Lookup L1 for La/i2
		coeff[b] = 2;	// Lookup L2 for Lb/i3

		// G loss values only appropriate if NHPLTJ is being used
		if(type==NHPLTJID)
		{
			flow_type_benson = 3; // Winterbone 2 is Benson 3
			G1[bpipe2] = C3*pPpt->gammaAir(pBN[bpipe2]->T)*pow(U_star[bpipe2]/A_star[bpipe1],2);
			G2[bpipe2] = C3*pPpt->gammaAir(pBN[bpipe1]->T)*pow(U_star[bpipe1]/A_star[bpipe1],2) + 1;
			G1[bpipe3] = 0;
			G2[bpipe3] = C3*pPpt->gammaAir(pBN[bpipe1]->T)*pow(U_star[bpipe1]/A_star[bpipe1],2) + 1;
		}
		break;																																	
	case 3:
		com = pipe3; i2 = pipe1; i3 = pipe2;
		pipe_flow[pipe1] = INFLOW; pipe_flow[pipe2] = INFLOW; pipe_flow[pipe3] = OUTFLOW;		
		coeff[a] = 3;	// Lookup L3 for La/i2
		coeff[b] = 3;	// Lookup L4=L3 for Lb/i3

		// G loss values only appropriate if NHPLTJ is being used
		if(type==NHPLTJID)
		{
			flow_type_benson = 2; // Winterbone 3 is Benson 2
			G1[bpipe2] = C1*pPpt->gammaAir(pBN[bpipe2]->T)*pow(U_star[bpipe2]/A_star[bpipe1],2);
			G2[bpipe2] = C1*pPpt->gammaAir(pBN[bpipe1]->T)*pow(U_star[bpipe1]/A_star[bpipe1],2) + 1;
			G1[bpipe3] = C2*pPpt->gammaAir(pBN[bpipe3]->T)*pow(U_star[bpipe3]/A_star[bpipe1],2);
			G2[bpipe3] = pow(A_star[bpipe2]/A_star[bpipe1], (2*pPpt->gammaAir(pBN[bpipe2]->T))/(pPpt->gammaAir(pBN[bpipe2]->T)-1));
		}
		break;																																
	case 4:	
		com = pipe1; i2 = pipe2; i3 = pipe3;
		pipe_flow[pipe1] = INFLOW; pipe_flow[pipe2] = OUTFLOW; pipe_flow[pipe3] = OUTFLOW;
		coeff[a] = 5;	// Lookup L5 for La/i2
		coeff[b] = 6;	// Lookup L6 for Lb/i3
		
		// G loss values only appropriate if NHPLTJ is being used
		if(type==NHPLTJID)
		{
			flow_type_benson = 4; // Winterbone 4 is Benson 4
			G1[bpipe2] = C4*pPpt->gammaAir(pBN[bpipe2]->T)*pow(U_star[bpipe2]/A_star[bpipe1],2);
			G2[bpipe2] = C4*pPpt->gammaAir(pBN[bpipe1]->T)*pow(U_star[bpipe1]/A_star[bpipe1],2) + 1;
			G1[bpipe3] = 0;
			G2[bpipe3] = -C5*pPpt->gammaAir(pBN[bpipe2]->T)*pow(U_star[bpipe2]/A_star[bpipe1],2)*pow(A_star[bpipe2]/A_star[bpipe1], 2/(pPpt->gammaAir(pBN[bpipe2]->T)-1))
							+ 1 + C5*pPpt->gammaAir(pBN[bpipe1]->T)*pow(U_star[bpipe1]/A_star[bpipe1],2);
		}
		break;																																	
	case 5:
		com = pipe2; i2 = pipe1; i3 = pipe3;
		pipe_flow[pipe1] = OUTFLOW; pipe_flow[pipe2] = INFLOW; pipe_flow[pipe3] = OUTFLOW;
		coeff[a] = 5;	// Lookup L5 for La/i2
		coeff[b] = 6;	// Lookup L6 for Lb/i3

		// G loss values only appropriate if NHPLTJ is being used
		if(type==NHPLTJID)
		{
			flow_type_benson = 6; // Winterbone 5 is Benson 6
			G1[bpipe2] = 0;
			G2[bpipe2] = 1;
			G1[bpipe3] = C6*pPpt->gammaAir(pBN[bpipe3]->T)*pow(U_star[bpipe3]/A_star[bpipe1],2);
			G2[bpipe3] = 1;
		}
		break;																																
	case 6:
		com = pipe3; i2 = pipe1; i3 = pipe2;
		pipe_flow[pipe1] = OUTFLOW; pipe_flow[pipe2] = OUTFLOW; pipe_flow[pipe3] = INFLOW;
		coeff[a] = 7;	// Lookup L7 for La/i2
		coeff[b] = 7;	// Lookup L8=L7 for Lb/i3
		
		// G loss values only appropriate if NHPLTJ is being used
		if(type==NHPLTJID)
		{
			flow_type_benson = 5; // Winterbone 6 is Benson 5
			G1[bpipe2] = C4*pPpt->gammaAir(pBN[bpipe2]->T)*pow(U_star[bpipe2]/A_star[bpipe1],2);
			G2[bpipe2] = C4*pPpt->gammaAir(pBN[bpipe1]->T)*pow(U_star[bpipe1]/A_star[bpipe1],2) + 1;
			G1[bpipe3] = 0;
			G2[bpipe3] = 1;
		}
		break;
	default:
		cout << "For some reason, we have come to the default section of a switch!\n";
		break;
	}
}				
*/
/*
void CJunction::CommonMassFlowMachNo(CProperties* pPpt)
{
	// Eq. 8.112
	m_dot[com] = (pPpt->PREF/(EX ? pPpt->AREFe : pPpt->AREFi))*pPpt->gammaAir(pBN[com]->T)*U_star[com]*(pow(A_star[com], 2/(pPpt->gammaAir(pBN[com]->T)-1))/AA[com])*Fb[com];
	m_dot[i2] = (pPpt->PREF/(EX ? pPpt->AREFe : pPpt->AREFi))*pPpt->gammaAir(pBN[i2]->T)*U_star[i2]*(pow(A_star[i2], 2/(pPpt->gammaAir(pBN[i2]->T)-1))/AA[i2])*Fb[i2];
	m_dot[i3] = (pPpt->PREF/(EX ? pPpt->AREFe : pPpt->AREFi))*pPpt->gammaAir(pBN[i3]->T)*U_star[i3]*(pow(A_star[i3], 2/(pPpt->gammaAir(pBN[i3]->T)-1))/AA[i3])*Fb[i3];

	// Evaluate mass flow ratio and common branch Mach no.
	Mn = fabs(U_star[com]/A_star[com]);	// Evaluate common branch Mach No.
	W[a] = fabs(m_dot[i2]/m_dot[com]);	// Mass flow ratio, other pipe 'i2'
	W[b] = fabs(m_dot[i3]/m_dot[com]);	// Mass flow ratio, other pipe 'i3'
}
*/
/*
void CJunction::LossCoeffs_Interpolation(CProperties* pPpt, int local_flow_type, double time)
{
	int graph, curve, curve_above, curve_below, point_curve_above, point_curve_below, 
		point_curve_above_above, point_curve_above_below, 
		point_curve_below_above, point_curve_below_below;
	
	double value_curve_above, value_curve_below;

	for(int branch=a; branch<=b; ++branch)
	{
		// Find the appropriate loss coefficient graph
		graph=0;
		while(pPpt->Loss[graph][0][0][0] != this->coeff[branch]) ++graph;
		//cout << "Using graph[" << graph << "]\n";

		// Find the right curve
		curve=0;
		while(pPpt->Loss[graph][curve][0][1] < this->W[branch]) ++curve;
		// Mass flow ratio between curve (curve_above) and curve-1 (curve_below)
		curve_above = curve;
		curve_below = curve-1;
		if(curve==0) curve_below = curve_above;
		//cout << "Interploating between curves [" << curve_above << "] and [" << curve_below << "]\n";

		// Find the two points on each of the curves that Mn lies between and interpolate each
		// 'curve above'
		point_curve_above=0;
		while(pPpt->Loss[graph][curve_above][point_curve_above][2] < this->Mn) ++point_curve_above;
		// Mn lies between point_curve_above and point_curve_above-1
		point_curve_above_above = point_curve_above;
		point_curve_above_below = point_curve_above-1;
		if(point_curve_above==0) point_curve_above_below = point_curve_above_above;
		// Interpolate these points
		if(fabs(pPpt->Loss[graph][curve_above][point_curve_above_above][2] - pPpt->Loss[graph][curve_above][point_curve_above_below][2]) < 1e-6)
		{
			value_curve_above 
			= pPpt->Loss[graph][curve_above][point_curve_above_below][3] + 0;
		}
		else
		{
			value_curve_above 
				= pPpt->Loss[graph][curve_above][point_curve_above_below][3]
				+ 
				( ((this->Mn - pPpt->Loss[graph][curve_above][point_curve_above_below][2])
				 /(pPpt->Loss[graph][curve_above][point_curve_above_above][2] - pPpt->Loss[graph][curve_above][point_curve_above_below][2]))
				*(pPpt->Loss[graph][curve_above][point_curve_above_above][3] - pPpt->Loss[graph][curve_above][point_curve_above_below][3]) );
		}

		// 'curve below'
		point_curve_below=0;
		while(pPpt->Loss[graph][curve_below][point_curve_below][2] < this->Mn) ++point_curve_below;
		// Mn lies between point_curve_below and point_curve_below-1
		point_curve_below_above = point_curve_below;
		point_curve_below_below = point_curve_below-1;
		if(point_curve_below==0) point_curve_below_below = point_curve_below_above;
		// Interpolate these points
		if(fabs(pPpt->Loss[graph][curve_below][point_curve_below_above][2] - pPpt->Loss[graph][curve_below][point_curve_below_below][2]) < 1e-6)
		{
			value_curve_below 
			= pPpt->Loss[graph][curve_below][point_curve_below_below][3] + 0;
		}
		else
		{
			value_curve_below 
				= pPpt->Loss[graph][curve_below][point_curve_below_below][3]
				+ 
				( ((this->Mn - pPpt->Loss[graph][curve_below][point_curve_below_below][2])
				/(pPpt->Loss[graph][curve_below][point_curve_below_above][2] - pPpt->Loss[graph][curve_below][point_curve_below_below][2]))
				*(pPpt->Loss[graph][curve_below][point_curve_below_above][3] - pPpt->Loss[graph][curve_below][point_curve_below_below][3]) );
		}
	
		// Now interpolate the values based on the curves
		if(fabs(pPpt->Loss[graph][curve_above][0][1] - pPpt->Loss[graph][curve_below][0][1]) < 1e-6)
		{
			this->L[branch]
			= value_curve_below + 0;
		}
		else
		{
			this->L[branch]
			= value_curve_below
			+
			( ((this->W[branch] - pPpt->Loss[graph][curve_below][0][1])
			  /(pPpt->Loss[graph][curve_above][0][1] - pPpt->Loss[graph][curve_below][0][1]))
			 *(value_curve_above - value_curve_below) );
		}	
	}

	// Now calculate x_star
	switch(local_flow_type)
	{
	case 0:
		x_star[pipe1] = 1;
		x_star[pipe2] = 1;
		x_star[pipe3] = 1;
		break;
	case 1: // case 'A':
		x_star[pipe1] = 1;
		x_star[pipe2] = pow(1 - L[a], (pPpt->gammaAir(pBN[pipe2]->T)-1)/(2*pPpt->gammaAir(pBN[pipe2]->T)));
		x_star[pipe3] = pow(1 - L[b], (pPpt->gammaAir(pBN[pipe3]->T)-1)/(2*pPpt->gammaAir(pBN[pipe3]->T)));
		break;
	case 2: // case 'B':
		x_star[pipe1] = 1;
		x_star[pipe2] = pow(1 - L[a], -(pPpt->gammaAir(pBN[pipe2]->T)-1)/(2*pPpt->gammaAir(pBN[pipe2]->T)));
		x_star[pipe3] = pow((1 - L[b])/(1 - L[a]), (pPpt->gammaAir(pBN[pipe3]->T)-1)/(2*pPpt->gammaAir(pBN[pipe3]->T)));
		break;
	case 3: // case 'C':
		x_star[pipe1] = 1;
		x_star[pipe2] = pow((1 - L[b])/(1 - L[a]), (pPpt->gammaAir(pBN[pipe2]->T)-1)/(2*pPpt->gammaAir(pBN[pipe2]->T)));
		x_star[pipe3] = pow((1 - L[a]), -(pPpt->gammaAir(pBN[pipe3]->T)-1)/(2*pPpt->gammaAir(pBN[pipe3]->T)));
		break;
	case 4: // case 'D':
		x_star[pipe1] = 1;
		x_star[pipe2] = pow(1 + L[a], -(pPpt->gammaAir(pBN[pipe2]->T)-1)/(2*pPpt->gammaAir(pBN[pipe2]->T)));
		x_star[pipe3] = pow(1 + L[b], (pPpt->gammaAir(pBN[pipe3]->T)-1)/(2*pPpt->gammaAir(pBN[pipe3]->T)));
		break;
	case 5: // case 'E':
		x_star[pipe1] = 1;
		x_star[pipe2] = pow(1 + L[a], -(pPpt->gammaAir(pBN[pipe2]->T)-1)/(2*pPpt->gammaAir(pBN[pipe2]->T)));
		x_star[pipe3] = pow((1 + L[b])/(1 + L[a]), (pPpt->gammaAir(pBN[pipe3]->T)-1)/(2*pPpt->gammaAir(pBN[pipe3]->T)));
		break;
	case 6: // case 'F':
		x_star[pipe1] = 1;
		x_star[pipe2] = pow((1 + L[b])/(1 + L[a]), (pPpt->gammaAir(pBN[pipe2]->T)-1)/(2*pPpt->gammaAir(pBN[pipe2]->T)));
		x_star[pipe3] = pow((1 + L[a]), -(pPpt->gammaAir(pBN[pipe3]->T)-1)/(2*pPpt->gammaAir(pBN[pipe3]->T)));
		break;
	default:
		cout << "Pulse converter - unknown flow type!\n";
		exit(1);
		break;
	}
}
*/
/*
void CJunction::LossCoeffs_Equations(CProperties* pPpt, int local_flow_type)
{
	double alpha, beta;
	switch (local_flow_type)
	{
	case 0:
		K_loss[a][label] = 0; K_loss[b][label] = 0;
		K_loss[a][value] = 0; K_loss[b][value] = 0;
		L[a] = 0; L[b] = 0;
		x_star[pipe1] = 1; x_star[pipe2] = 1; x_star[pipe3] = 1;
		break;
	case 1: //A
		K_loss[a][label] = 1; K_loss[b][label] = 2;
		K_loss[a][value] = pow(W[a]*psiT, 2) + 1 - 2*W[a]*psiT*cos(3/4*(PI - this->angle_rad)); // K1
		K_loss[b][value] = pow(W[b], 2) - 3/2*W[b] + 1/2; // K2	
		L[a] = 0.5*pPpt->gammaAir(pBN[i2]->T)*this->Mn*this->Mn*(K_loss[a][value] + pow( (U_star[i2]*AA[i2])/(U_star[com]*AA[com]), 2 ) - 1); // L1
		L[b] = 0.5*pPpt->gammaAir(pBN[i3]->T)*this->Mn*this->Mn*(K_loss[b][value] + pow( (U_star[i3]*AA[i3])/(U_star[com]*AA[com]), 2 ) - 1); // L2
		x_star[pipe1] = 1;
		x_star[pipe2] = pow(1 - L[a], (pPpt->gammaAir(pBN[pipe2]->T)-1)/(2*pPpt->gammaAir(pBN[pipe2]->T)));
		x_star[pipe3] = pow(1 - L[b], (pPpt->gammaAir(pBN[pipe3]->T)-1)/(2*pPpt->gammaAir(pBN[pipe3]->T)));
		break;
	case 2: //B
		K_loss[a][label] = 3; K_loss[b][label] = 4;
		K_loss[a][value] = 1 + pow(W[a]/psiT, 2) - 2*W[a]/psiT*cos(3/4*(PI - this->angle_rad)); // K3
		K_loss[b][value] = 1 + pow(W[b]/psiT, 2) - 2*W[b]/psiT*cos(3/4*this->angle_rad); // K4
		L[a] = 0.5*pPpt->gammaAir(pBN[i2]->T)*this->Mn*this->Mn*(K_loss[a][value] + pow( (U_star[i2]*AA[i2])/(U_star[com]*AA[com]), 2 ) - 1); // L3
		L[b] = 0.5*pPpt->gammaAir(pBN[i3]->T)*this->Mn*this->Mn*(K_loss[b][value] + pow( (U_star[i3]*AA[i3])/(U_star[com]*AA[com]), 2 ) - 1); // L4
		x_star[pipe1] = 1;
		x_star[pipe2] = pow(1 - L[a], -(pPpt->gammaAir(pBN[pipe2]->T)-1)/(2*pPpt->gammaAir(pBN[pipe2]->T)));
		x_star[pipe3] = pow((1 - L[b])/(1 - L[a]), (pPpt->gammaAir(pBN[pipe3]->T)-1)/(2*pPpt->gammaAir(pBN[pipe3]->T)));
		break;
	case 3: //C
		K_loss[a][label] = 5; K_loss[b][label] = 6;
		K_loss[a][value] = pow(W[a],2) - 3/2*W[a] + 1/2; // K5
		K_loss[b][value] = pow(W[b]*psiT,2) + 1 - 2*W[b]*psiT*cos(3/4*this->angle_rad); // K6
		L[a] = 0.5*pPpt->gammaAir(pBN[i2]->T)*this->Mn*this->Mn*(K_loss[a][value] + pow( (U_star[i2]*AA[i2])/(U_star[com]*AA[com]), 2 ) - 1); // L5
		L[b] = 0.5*pPpt->gammaAir(pBN[i3]->T)*this->Mn*this->Mn*(K_loss[b][value] + pow( (U_star[i3]*AA[i3])/(U_star[com]*AA[com]), 2 ) - 1); // L6
		x_star[pipe1] = 1;
		x_star[pipe2] = pow((1 - L[b])/(1 - L[a]), (pPpt->gammaAir(pBN[pipe2]->T)-1)/(2*pPpt->gammaAir(pBN[pipe2]->T)));
		x_star[pipe3] = pow((1 - L[a]), -(pPpt->gammaAir(pBN[pipe3]->T)-1)/(2*pPpt->gammaAir(pBN[pipe3]->T)));
		break;
	case 4: //D
		K_loss[a][label] = 7; K_loss[b][label] = 8;
		K_loss[a][value] = 4*W[a] - 1 + W[a]*W[a]*(psiT*psiT - 2 +2*psiT*cos(this->angle_rad)); // K7 
		K_loss[b][value] = 1 - pow(W[b], 2) + 2*pow(1 - W[b], 2)*psiT*cos(this->angle_rad); // K8
		L[a] = 0.5*pPpt->gammaAir(pBN[i2]->T)*this->Mn*this->Mn*(K_loss[a][value] + 1 - pow( (U_star[i2]*AA[i2])/(U_star[com]*AA[com]), 2 )); // L7			
		L[b] = 0.5*pPpt->gammaAir(pBN[i3]->T)*this->Mn*this->Mn*(K_loss[b][value] + 1 - pow( (U_star[i3]*AA[i3])/(U_star[com]*AA[com]), 2 )); // L8	
		x_star[pipe1] = 1;
		x_star[pipe2] = pow(1 + L[a], -(pPpt->gammaAir(pBN[pipe2]->T)-1)/(2*pPpt->gammaAir(pBN[pipe2]->T)));
		x_star[pipe3] = pow(1 + L[b], (pPpt->gammaAir(pBN[pipe3]->T)-1)/(2*pPpt->gammaAir(pBN[pipe3]->T)));
		break;
	case 5: //E
		alpha = (PI - angle_rad)/4;
		beta = angle_rad/4;

		K_loss[a][label] = 9; K_loss[b][label] = 10;
		K_loss[a][value] = 1 + (4*W[a] - 2)/psiT*cos(this->angle_rad) + pow(W[a]/psiT, 2); // K9
				
//		K_loss[a][value] = 2*pow(W[a], 2)*(1/psiT)*cos(angle_rad + alpha)
//							- 2*pow((1 - W[a]), 2)*(1/psiT)*cos(angle_rad - beta)
//							+ pow(W[a], 2)*(1/psiT) + 1; // Alternative K9 with angle correction


		K_loss[b][value] = 1 + (2 - 4*W[b])/psiT*cos(this->angle_rad) + pow(W[b]/psiT, 2); // K10

//		K_loss[b][value] = 2*pow(1 - W[b], 2)*(1/psiT)*cos(angle_rad + alpha)
//							- 2*pow(W[b], 2)*(1/psiT)*cos(angle_rad - beta)
//							+ pow(W[b], 2)*(1/psiT) + 1; // Alternative K10 with angle correction
		
		
		L[a] = 0.5*pPpt->gammaAir(pBN[i2]->T)*this->Mn*this->Mn*(K_loss[a][value] + 1 - pow( (U_star[i2]*AA[i2])/(U_star[com]*AA[com]), 2 )); // L9	
		L[b] = 0.5*pPpt->gammaAir(pBN[i3]->T)*this->Mn*this->Mn*(K_loss[b][value] + 1 - pow( (U_star[i3]*AA[i3])/(U_star[com]*AA[com]), 2 )); // L10
		x_star[pipe1] = 1;
		x_star[pipe2] = pow(1 + L[a], -(pPpt->gammaAir(pBN[pipe2]->T)-1)/(2*pPpt->gammaAir(pBN[pipe2]->T)));
		x_star[pipe3] = pow((1 + L[b])/(1 + L[a]), (pPpt->gammaAir(pBN[pipe3]->T)-1)/(2*pPpt->gammaAir(pBN[pipe3]->T)));
		break;
	case 6: //F
		K_loss[a][label] = 11; K_loss[b][label] = 12;
		K_loss[a][value] = 2*psiT/(psiT + cos(this->angle_rad)/2)
							*(1 - pow(W[a], 2) - pow(1 - W[a], 2)*psiT*cos(this->angle_rad)) 
							+ pow(W[a], 2) - 1; // K11
		K_loss[b][value] = 2*psiT/(psiT + cos(this->angle_rad)/2)
							*(1 - pow(1 - W[b], 2) - pow(W[b], 2)*psiT*cos(this->angle_rad)) 
							+ pow(W[b]*psiT, 2) - 1; // K12
		L[a] = 0.5*pPpt->gammaAir(pBN[i2]->T)*this->Mn*this->Mn*(K_loss[a][value] + 1 - pow( (U_star[i2]*AA[i2])/(U_star[com]*AA[com]), 2 )); // L11
		L[b] = 0.5*pPpt->gammaAir(pBN[i3]->T)*this->Mn*this->Mn*(K_loss[b][value] + 1 - pow( (U_star[i3]*AA[i3])/(U_star[com]*AA[com]), 2 )); // L12
		x_star[pipe1] = 1;
		x_star[pipe2] = pow((1 + L[b])/(1 + L[a]), (pPpt->gammaAir(pBN[pipe2]->T)-1)/(2*pPpt->gammaAir(pBN[pipe2]->T)));
		x_star[pipe3] = pow((1 + L[a]), -(pPpt->gammaAir(pBN[pipe3]->T)-1)/(2*pPpt->gammaAir(pBN[pipe3]->T)));
		break;
	default:
		cout << "Winterbone junction - unknown flow type!\n";
		exit(1);
		break;
	}
}
*/
/*
void CJunction::CommonEntropyLevels(CProperties* pPpt, int flow_type)
{
	// Based on Winterbone flow types

	double* A0_star_sq;
	A0_star_sq = new double [3];
	
	// Some intermediate values used in the separating flow calculations:
	for(int p=0; p<3; ++p) A0_star_sq[p] = pow(this->A_star[p],2) + ((pPpt->gammaAir(pBN[p]->T)-1)/2)*pow(this->U_star[p],2);
		
	if(flow_type >= 4 && flow_type <= 6) // Joining flow
	{
		// Calculate AA[N]c (8.114) for joining flow types
		// Entropy levels in the two branches in which flow is toward the junction
		// are equal to their uncorrected values
		this->AA[i2] = this->AA_n[i2]; AA[i3] = AA_n[i3];
		AA[com] = 
				-(AA[i2]*pow(A_star[i2],2/(pPpt->gammaAir(pBN[i2]->T)-1))*U_star[i2]*(pow(A_star[i2],2) + ((pPpt->gammaAir(pBN[i2]->T)-1)/2)*pow(U_star[i2],2))*Fb[i2]
				+ AA[i3]*pow(A_star[i3],2/(pPpt->gammaAir(pBN[i3]->T)-1))*U_star[i3]*(pow(A_star[i3],2) + ((pPpt->gammaAir(pBN[i3]->T)-1)/2)*pow(U_star[i3],2))*Fb[i3])
				/ (pow(A_star[com],2/(pPpt->gammaAir(pBN[com]->T)-1))*U_star[com]*(pow(A_star[com],2) + ((pPpt->gammaAir(pBN[com]->T)-1)/2)*pow(U_star[com],2))*Fb[com]);
	}
	else
	{
		if(flow_type >= 1 && flow_type <= 3) // Separating flow
		{
			// Calculate AA[N]c (8.115) for separating flow types
			// AA[com] is equal to the uncorrected value AAn[com]
			AA[com] = AA_n[com];
			AA[i2] = AA[com]*sqrt(A0_star_sq[com]/A0_star_sq[i2]);
			AA[i3] = AA[com]*sqrt(A0_star_sq[com]/A0_star_sq[i3]);
		}
		else // Corner flows and no flow
		{
			if(flow_type == 0)
			{
				AA[com] = AA_n[com];
				AA[i2] = AA_n[i2];
				AA[i3] = AA_n[i3];
			}
			else cout << "Unknown flow type!" << endl;
		}
	}
}
*/
/*
void CJunction::CommonLambdaInStarCorrection()
{
	// ------------------------------------------
	// Calculate lambda_in[N]c (8.125) N=1,2,3
	// this is the same as for constant pressure model
	// ------------------------------------------
	// Apply the correction to lambda_in_star - Eq. 8.137 of Benson p440
	for(int p=0; p<NPIPES; ++p)
		lambda_in_star[p] = A_star[p] 
							+ (AA_n[p]/AA[p])*(lambda_in_star_n[p] - A_star[p]);
}
*/
/*
void CJunction::CommonLambdaOut(CProperties* pPpt)
{
	// Calculate lambda_out (8.126) N =1,2,3
	for(int p=0; p<NPIPES; ++p)
		lambda_out[p] = (A_star[p] - ((pPpt->gammaAir(pBN[p]->T)-1)/2)*U_star[p])*AA[p];
}
*/

void CJunction::CommonUpdate(CProperties* pPpt)
{
/*
	//////////////////////////////////////////
	double *lambda_in_c;
	lambda_in_c = new double [NPIPES];
	
	for(int p=0; p<NPIPES; ++p)
		lambda_in_c[p] = lambda_in_star[p]*AA[p];

	
	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out, AA, pipe_flow);
	//////////////////////////////////////////
*/
	
	for(int p=0; p<NPIPES; ++p)
	{
		(*pCLOUT[p])[R+1] = lambda_out[p];	// Always update lambda_out
		pipe_flow_old[p] = *pend_flow[p];	// Record previous flow directions
		*pend_flow[p] = pipe_flow[p];

		if(*pend_flow[p]==INFLOW) // Only INFLOW creates pathlines
		{
			(*pCLIN[p])[R+1] = lambda_in_star[p]*AA[p]; // Only need to update on INFLOW, because OUTFLOWs will maintain the uncorrected value
			if((pPipe[p]->METHOD==pPipe[p]->MMOC && !pPpt->HOMENTROPIC) || 
			   (pPipe[p]->METHOD==pPipe[p]->W_ALPHA_BETA && pPpt->COMBINED_WAB_MOC))
				(*pPathLine[p]).AA = AA[p]; // Only INFLOW will have a new pathline at the end of the pipe ready to go

//???		
(*pBN[p]).AA[R+1] = AA[p];
// This is done in all other boundary methods, so do it here too
//???
			if(pipe_flow_old[p] != *pend_flow[p])
			// If old flow direction was not INFLOW then there was no new path line created
			// so adjust existing one by setting its XK value to the appropriate end
			{
				if((pPipe[p]->METHOD==pPipe[p]->MMOC && !pPpt->HOMENTROPIC) || 
				   (pPipe[p]->METHOD==pPipe[p]->W_ALPHA_BETA && pPpt->COMBINED_WAB_MOC))
				{
					if(this->end[p]==ODD){(*pPathLine[p]).XK = 0; (*pPathLine[p]).XK_old = (*pPathLine[p]).XK;}
					else{(*pPathLine[p]).XK = pPipe[p]->XPIPE; (*pPathLine[p]).XK_old = (*pPathLine[p]).XK;}
				}
			}
			else 
			// rpipe_flow_old == *pend_flow
			// Check that a new path line HAS been created correctly if INFLOW
			{
				if((pPipe[p]->METHOD==pPipe[p]->MMOC && !pPpt->HOMENTROPIC) || 
				   (pPipe[p]->METHOD==pPipe[p]->W_ALPHA_BETA && pPpt->COMBINED_WAB_MOC))
				{
          if(*pend_flow[p]==INFLOW && ((end[p]==ODD && (*pPathLine[p]).XK!=0) || (end[p]==EVEN && ((*pPathLine[p]).XK - pPipe[p]->XPIPE > pPpt->ZERO_TOL))))
					{
						if(EX) cout << "Exhaust "; else cout << "Intake ";
						cout << "Junction [" << ID << "], Connection (" << p << "): ";
						if(pPipe[p]->EX) cout << "Exhaust "; else cout << "Intake ";
						cout << "Pipe [" << pPipe[p]->ID << "]:" << endl;
						cout << "Should've created new path line but incorrect XK: ";
            int temp_precision = cout.precision(); cout << setprecision(12);
						if(this->end[p]==ODD) cout << "(*pPathLine[p]).XK = " << (*pPathLine[p]).XK << ", but it should be 0" << endl;
						else cout << "(*pPathLine[p]).XK = " << (*pPathLine[p]).XK << ", but it should be " << pPipe[p]->XPIPE << endl;
						cout << endl;
            cout << setprecision(temp_precision);
					}
				}
			}
		}
	}
}


void CJunction::PrintToScreen(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".PrintToScreen\n");}
	pPpt->Out(Underline(Identify(), "=", "", strDesc)); pPpt->Out("\n");
	if(CONSTP)
	{
		//pPpt->Out("Constant pressure junction - no important output parameters"); pPpt->Out("\n");
		pPpt->Out("Constant pressure junction"); pPpt->Out("\n");
	}
	else
	{
		pPpt->Out("Pressure loss junction"); pPpt->Out("\n");
	}
	pPpt->Out("\n");
	pPpt->Out(Underline("Mass flow rates toward junction", "-"));
	int j;
	for(j=0; j<NPIPES; ++j)
	{
		//pPpt->Out("Branch ["); pPpt->Out(j); pPpt->Out("] of [0--"); pPpt->Out(NPIPES-1); pPpt->Out("]:\n");
		//pPpt->Out("Mass flow rate towards junction, m_dot[j="); pPpt->Out(j); pPpt->Out("]\t\t=\t"); pPpt->Out(m_dot[j]); pPpt->Out(" kg/s\n");
		pPpt->Out("Branch ["); pPpt->Out(j); pPpt->Out("] of [0--"); pPpt->Out(NPIPES-1); pPpt->Out("], m_dot[j="); pPpt->Out(j); pPpt->Out("]\t\t\t=\t"); pPpt->Out(m_dot[j]); pPpt->Out(" kg/s\n");
		//pPpt->Out("\n");
	}
	pPpt->Out("Check resultant mass flow is about zero, m_dot_total\t=\t"); pPpt->Out(m_dot_total); pPpt->Out(" kg/s\n");
	pPpt->Out("\n");
	pPpt->Out("Max +ve mass flow is in branch ["); pPpt->Out(MAX_FLOW_BRANCH); pPpt->Out("], m_dot[j="); pPpt->Out(MAX_FLOW_BRANCH); pPpt->Out("]\t\t=\t"); pPpt->Out(m_dot[MAX_FLOW_BRANCH]); pPpt->Out(" kg/s\n");
	pPpt->Out("\n");
/*
	switch(this->type)
	{
	case HCPJID:
		pPpt->Out("HCPJ:\n");
		break;
	case NHCPJID:
		pPpt->Out("NHCPJ:\n");
		break;
	case NHPLTJID:
		pPpt->Out("NHPLTJ:\n");
		break;
	case NHPLPCID:
		pPpt->Out("NHPLPC:\n");
		break;
	case NHPLJID:
		pPpt->Out("NHPLJ:\n");
		break;
	default:
		pPpt->Out("Unknown junction type\n");
		break;
	}
*/
/*
	cout << "flow_type_winterbone = " << flow_type_winterbone << endl;
	cout << "flow_type_winterbone_old = " << flow_type_winterbone_old << endl;
	cout << "flow_type_winterbone_orig = " << flow_type_winterbone_orig << endl;
//	cout << "flow_type_benson = " << flow_type_benson << endl;
	cout << "switch_counter = " << switch_counter << endl;
	if(this->type==NHPLJID)
	{
		cout << "Iterations completed under:\n";
		cout << "---------------------------\n";
		cout << " Normal operation\t" << normal_counter_total << endl;
		cout << " FIX_FLOW_TYPE regime\t" << fix_counter_total << endl;
		cout << " OPPOSITE regime\t" << opposite_counter_total << endl;
		cout << " OTHER regime\t\t" << other_counter_total << endl;
		cout << " ZERO_FLOW regime\t" << zero_counter_total << endl;
		cout << "\t\t\t__________" << endl;
		cout << "\t\t\t" << normal_counter_total + fix_counter_total + opposite_counter_total + other_counter_total + zero_counter_total << endl;
	}
//*/
/*
	cout << "com = " << com << endl;
	cout << "i2 = " << i2 << endl;
	cout << "i3 = " << i3 << endl;
	cout << "angle = " << this->angle_deg << Deg() << endl;
	cout << "psiT = " << this->psiT << endl;
	cout << "coeff[a] = " << this->coeff[0] << endl;
	cout << "coeff[b] = " << this->coeff[1] << endl;
	cout << "m_dot[com] = " << m_dot[com] << endl;
	cout << "m_dot[i2] = " << m_dot[i2] << endl;
	cout << "m_dot[i3] = " << m_dot[i3] << endl;
//*/
/*
	cout << "Mn = " << this->Mn << endl;
	cout << "W[a] = " << this->W[0] << endl;
	cout << "W[b] = " << this->W[1] << endl;
	cout << "K" << K_loss[a][label] << " = " << K_loss[a][value] << endl;
	cout << "K" << K_loss[b][label] << " = " << K_loss[b][value] << endl;
	cout << "L[a] = " << this->L[0] << endl;
	cout << "L[b] = " << this->L[1] << endl;
	cout << "x_star[pipe1] = " << this->x_star[pipe1] << endl;
	cout << "x_star[pipe2] = " << this->x_star[pipe2] << endl;
	cout << "x_star[pipe3] = " << this->x_star[pipe3] << endl;
*/
/*
	cout << "G1[1] = " << G1[1] << endl;
	cout << "G2[1] = " << G2[1] << endl;
	cout << "G1[2] = " << G1[2] << endl;
	cout << "G2[2] = " << G2[2] << endl;
*/
}

void CJunction::PrintToFile(double time, double ca)
// ============================================================ //
// Prints instantaneous junction data to file.					//
// This function is called from the main function.				//
// ============================================================ //
{
/*
	if(type==NHPLJID)
	{
		fprintf(OUTPUT_FILE,"%f\t%f\t%f\t%f\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
			time, ca, this->angle_rad, psiT, flow_type_winterbone, 
			Mn, W[a], W[b], K_loss[a][value], K_loss[b][value], L[a], L[b],
			x_star[pipe2], x_star[pipe3]);
	}
*/
}

void CJunction::ReadInput(CProperties* pPpt, char *InputFile)
{
	int j, jj;
	int last_entry;
	CommonReadInput(InputFile, pPpt->NUM_PARAMS, labels, values, strings, last_entry);
	for(int r=0; r<last_entry+1; ++r)
	{
		// ====================================================================================================
		// Parameter file for Exhaust Junction [0]
		// ====================================================================================================

		// Optional object description (max. 500 characters - use underscores for spaces)
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "strDesc") == 0) strDesc = strings[r];

		if(strcmp(labels[r], "num_branches") == 0)
		{
			num_branches = int(values[r]); // Only time num_branches is used. Everywhere else use NPIPES.

			if(num_branches < 2 || NPIPES < 2)
			{
				pPpt->Out("Error: ");
				pPpt->Out(Identify());
				pPpt->Out(" requires at least two pipes to be joined there. Exiting program.\n\n");
				exit(1);
			}

			if(num_branches != NPIPES)
			{
				pPpt->Out("Warning: number of branches specified in input file for ");
				pPpt->Out(Identify());
				pPpt->Out(" (num_branches = "); pPpt->Out(num_branches);
				pPpt->Out(") does no match number of pipe connections (NPIPES = "); pPpt->Out(NPIPES);
				pPpt->Out(") at this boundary.\n");	
				if(NPIPES < num_branches)
				{
					pPpt->Out("\nCan continue, but will only use the first ");
					pPpt->Out(NPIPES); pPpt->Out(" branches specified in the input file.\n");
					pPpt->Out("\nContinue? (y/n):");
					char cont;
					cin >> cont;
					pPpt->Out("\n");
					if(cont=='n' || cont=='N')
					{
						pPpt->Out("User answers NO. Exiting program.\n\n");
						exit(1);
					}
					else 
					{				
						pPpt->Out("User answers YES. Setting num_branches = NPIPES.\n\n");
						num_branches = NPIPES;
					}
				}
				else
				{
					pPpt->Out("\nCannot continue. Please specify ");
					pPpt->Out(NPIPES - num_branches); pPpt->Out(" more branches in the input file. Exiting program.\n");
					exit(1);
				}
			}
			
			// Dimension intialisation parameters using nsections
			id_branch = new int [NPIPES];
			ref0_branch = new double [NPIPES];
			ref1_branch = new double [NPIPES];

			// By definition, the reference branches referenced to themselves implies
			ref0_branch[0] = 0;
			ref1_branch[0] = 0;
			ref1_branch[1] = 0;
		}

		// Model selection (in non-homentropic flow; only constant pressure model exists in homentropic flow)
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "CONSTP") == 0)
		{
			CONSTP = DoubleToBool(values[r]);
			if(!CONSTP && pPpt->HOMENTROPIC)
			{
				pPpt->Out("Warning: "); pPpt->Out(Identify()); pPpt->Out(" specified as a pressure loss junction (CONSTP = "); pPpt->Out(TrueOrFalse(CONSTP));
				pPpt->Out("). Can continue but must use constant pressure model in homentropic flow.\n");
				pPpt->Out("\nContinue? (y/n):");
				char cont; cin >> cont;	pPpt->Out("\n");
				if(cont=='n' || cont=='N'){pPpt->Out("User answers NO. Exiting program.\n\n"); exit(1);}
				else{pPpt->Out("User answers YES. Setting CONSTP = true.\n\n"); CONSTP = true;}
			}
		}
		if(strcmp(labels[r], "constp_until_time") == 0) constp_until_time = values[r];

			// If using constant pressure model, i.e., CONSTP == 1 == true
			// ----------------------------------------------------------------------------------------------------
			if(strcmp(labels[r], "tol_A_star") == 0) tol_A_star = values[r];

			// Else using pressure loss model, i.e., CONSTP == 0 == false
			// ----------------------------------------------------------------------------------------------------
	
			// Branch angles (used only in pressure loss models)
			// ---------------------------------------------------------------------------------------------------- 
			for(j=0; j<NPIPES; ++j)
			{
				if(strcmp(labels[r], ConstructString(pPpt, "id_branch", IntToString(j))) == 0 && NPIPES>j) id_branch[j] = int(values[r]);
				if(strcmp(labels[r], ConstructString(pPpt, "ref0_branch", IntToString(j))) == 0 && NPIPES>j && j>0) ref0_branch[j] = values[r];
				if(strcmp(labels[r], ConstructString(pPpt, "ref1_branch", IntToString(j))) == 0 && NPIPES>j && j>1) ref1_branch[j] = values[r];
			}

			// Algorithm control
			// ----------------------------------------------------------------------------------------------------
			if(strcmp(labels[r], "tol_mass_sum") == 0)
			{
				tol_mass_sum = values[r];
				if(tol_mass_sum < 1e-15)
				{
					pPpt->Out("Warning: ");	pPpt->Out(Identify()); 
					pPpt->Out(" specified with a continuity convergence tolerance of (tol_mass_sum = "); pPpt->Out(tol_mass_sum);
					pPpt->Out("). Can continue but must use less stringent tolerance.\n");
					pPpt->Out("\nContinue and set tol_mass_sum = 1e-15? (y/n):");
					char cont; cin >> cont; pPpt->Out("\n");
					if(cont=='n' || cont=='N'){pPpt->Out("User answers NO. Exiting program.\n\n"); exit(1);}
					else{pPpt->Out("User answers YES. Setting tol_mass_sum = 1e-15.\n\n"); tol_mass_sum = 1e-15;}
				}
			}
			if(strcmp(labels[r], "tol_del_star") == 0)
			{
				tol_del_star = values[r];
				if(tol_del_star < 1e-15)
				{
					pPpt->Out("Warning: "); pPpt->Out(Identify());
					pPpt->Out(" specified with a pressure loss convergence tolerance of (tol_del_star = "); pPpt->Out(tol_del_star);
					pPpt->Out("). Can continue but must use less stringent tolerance.\n");
					pPpt->Out("\nContinue and set tol_del_star = 1e-15? (y/n):");
					char cont; cin >> cont; pPpt->Out("\n");
					if(cont=='n' || cont=='N'){pPpt->Out("User answers NO. Exiting program.\n\n"); exit(1);}
					else{pPpt->Out("User answers YES. Setting tol_del_star = 1e-15.\n\n"); tol_del_star = 1e-15;}
				}
			}
			if(strcmp(labels[r], "loop_limit_main") == 0) loop_limit_main = values[r];

		// Deprecated
		// ----------------------------------------------------------------------------------------------------
/*		
		if(strcmp(labels[r], "type") == 0)
		{
			switch(int(values[r]))
			{
			case 0:
				type = HCPJID;
				break;
			case 1:
				type = NHCPJID;
				break;
			case 2:
				type = NHPLTJID;
				break;
			case 3:
				type = NHPLJID;
				break;
			case 4:
				type = NHPLPCID;
				break;
			default:
				type = NHCPJID; // Use non-homentropic constant pressure model as default
				if(EX) cout << "CJunction::ReadInput; Exhaust Junction [" << ID << "]: unknown junction model specifier, using NHCPJID\n";
				else cout << "CJunction::ReadInput; Intake Junction [" << ID << "]: unknown junction model specifier, using NHCPJID\n";
				break;
			}
		}
		if(strcmp(labels[r], "angle_deg") == 0)
		{
			angle_deg = values[r];
			angle_rad = angle_deg/180*PI;
		}
		if(strcmp(labels[r], "tol_main") == 0) tol_main = values[r];
		if(strcmp(labels[r], "tol_cont") == 0) tol_cont = values[r];
		if(strcmp(labels[r], "loop_limit_cont") == 0) loop_limit_cont = values[r];
		if(strcmp(labels[r], "loop_limit_switch") == 0) loop_limit_switch = values[r];
		//if(strcmp(labels[r], "PRINT_MOVIE_FILE") == 0) PRINT_MOVIE_FILE = bool(values[r]);
*/		
		// ====================================================================================================
		// End of file
		// ====================================================================================================
	}

	// Calculate and store angle between each branch and every other branch 
	// ----------------------------------------------------------------------------------------------------
	for(j=0; j<NPIPES; ++j){ref0_branch[j] *= PI/180; ref1_branch[j] *= PI/180;} // Convert all angles to radians
	
	// Assign temporary vector coordinates for each branch
	double** vector_branch;
	vector_branch = new double* [NPIPES];
	for(j=0; j<NPIPES; ++j) vector_branch[j] = new double [3]; // For each coordinate X, Y, Z
	int X, Y, Z; X = 0; Y = 1; Z = 2; // Labels
	double temp_ref0, temp_ref1;
	for(j=0; j<NPIPES; ++j)
	{
		temp_ref0 = ref0_branch[j];
		temp_ref1 = ref1_branch[j];
		
		// For -PI/2 <= temp_ref0 <= PI/2
		vector_branch[j][X] = cos(temp_ref0); 
		vector_branch[j][Y] = sin(temp_ref0)*cos(temp_ref1);
		vector_branch[j][Z] = sin(temp_ref0)*sin(temp_ref1);
	}

	// Work out angle from vectors using cosine rule
	ref_datum_branch = new double* [NPIPES];
	for(j=0; j<NPIPES; ++j) ref_datum_branch[j] = new double [NPIPES]; // NPIPES X NPIPES array

	double temp_cos_theta;
	for(j=0; j<NPIPES; ++j)
	{
		for(jj=0; jj<NPIPES; ++jj)
		{
			// Cosine of angle between branches [j] and [jj] is simply
			temp_cos_theta =				( vector_branch[j][X]*vector_branch[jj][X] 
											+ vector_branch[j][Y]*vector_branch[jj][Y] 
											+ vector_branch[j][Z]*vector_branch[jj][Z]
											)/
											(sqrt(pow(vector_branch[j][X],2)
												+ pow(vector_branch[j][Y],2)
												+ pow(vector_branch[j][Z],2)
												 )*
											 sqrt(pow(vector_branch[jj][X],2)
												+ pow(vector_branch[jj][Y],2)
												+ pow(vector_branch[jj][Z],2)
												 )
											 );
			if(fabs(temp_cos_theta - 1.) < 1e-12) temp_cos_theta = 1;
			ref_datum_branch[j][jj] = acos(temp_cos_theta);
		}
	}	

	for(j=0; j<NPIPES; ++j) delete [] vector_branch[j];
	delete [] vector_branch;
}

void CJunction::ListProperties(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ListProperties\n");}

	// ====================================================================================================
	// Parameter file for Exhaust Junction [0]
	// ====================================================================================================

	int j, jj;
	pPpt->Out(Underline(Identify(), "=", "\t", strDesc));
	pPpt->Out("\n");

	// Model selection (in non-homentropic flow; only constant pressure model exists in homentropic flow)
	// ----------------------------------------------------------------------------------------------------
	pPpt->Out(Underline("Model selection", "-", "\t"));
	if(pPpt->HOMENTROPIC)
	{
		CONSTP = true;
		pPpt->Out("\tHomentropic flow - must use constant pressure model, CONSTP\t=\t");
		pPpt->Out(TrueOrFalse(CONSTP)); pPpt->Out("\n");
	}
	else
	{
		if(CONSTP)
		{
			pPpt->Out("\tConstant pressure junction specified, CONSTP\t=\t");
			pPpt->Out(TrueOrFalse(CONSTP)); pPpt->Out("\n");
			pPpt->Out("\t- using non-homentropic constant pressure model\n");
			pPpt->Out("\tNHCPJ A_star convergence tolerance, tol_A_star\t=\t"); pPpt->Out(tol_A_star); pPpt->Out("\n");
		}
		else
		{
			pPpt->Out("\tPressure loss junction specified, CONSTP\t=\t");
			pPpt->Out(TrueOrFalse(CONSTP)); pPpt->Out("\n");
			pPpt->Out("\t- using Bassett-Winterbone pressure loss model\n");
		}
	}
	if(constp_until_time>0)
	{
		pPpt->Out("\tForce NHCPJ model until, constp_until_time\t=\t"); pPpt->Out(constp_until_time); pPpt->Out(" s\n");
		pPpt->Out("\tNHCPJ A_star convergence tolerance, tol_A_star\t=\t"); pPpt->Out(tol_A_star); pPpt->Out("\n");
	}
	pPpt->Out("\n");

	// Branch angles (used only in pressure loss models)
	// ---------------------------------------------------------------------------------------------------- 
	if(!CONSTP)
	{
		pPpt->Out(Underline("Branch angles (used only in pressure loss models)", "-", "\t"));
		pPpt->Out("\tNumber of branches, NPIPES\t\t=\t"); pPpt->Out(NPIPES); pPpt->Out("\n");
		pPpt->Out("\n");
		pPpt->Out("\tBranch no.\t\tPipe ID\t\t\tAngle ref. 1 (deg.)\tAngle ref. 2 (deg.)\n");
		pPpt->Out("\t----------\t\t-------\t\t\t-------------------\t-------------------\n");
		for(j=0; j<NPIPES; ++j)
		{
			pPpt->Out("\t"); pPpt->Out(j+1);
			if(j==0) pPpt->Out(" (1st ref. pipe)\t");
			else
			{
				if(j==1) pPpt->Out(" (2nd ref. pipe)\t");
				else pPpt->Out("\t\t\t"); 
			}
			if(EX) pPpt->Out("Exhaust "); else pPpt->Out("Intake "); 
			pPpt->Out("Pipe ["); pPpt->Out(id_branch[j]); pPpt->Out("]\t"); 
			pPpt->Out(ref0_branch[j]*180/PI); pPpt->Out("\t\t\t"); pPpt->Out(ref1_branch[j]*180/PI); pPpt->Out("\n");
		}
		pPpt->Out("\n");

		// Calculated angles
		// ----------------------------------------------------------------------------------------------------
		pPpt->Out(Underline("Calculated angle (degrees) between branches", "-", "\t"));
		//pPpt->Out("\n");
		pPpt->Out("\t"); pPpt->Out("Datum\\Branch"); for(jj=0; jj<NPIPES; ++jj){pPpt->Out("\t[");pPpt->Out(id_branch[jj]);/*pPpt->Out(jj);*/pPpt->Out("]");}
		pPpt->Out("\n");
		for(j=0; j<NPIPES; ++j)
		{
			pPpt->Out("\t["); pPpt->Out(id_branch[j]); /*pPpt->Out(j);*/ pPpt->Out("]\t");
	
			for(jj=0; jj<NPIPES; ++jj)
			{
				pPpt->Out("\t");
				pPpt->Out(ref_datum_branch[j][jj]*180/PI, 4);
			}
			pPpt->Out("\n");
		}
		pPpt->Out("\n");

		// Algorithm control
		// ----------------------------------------------------------------------------------------------------
		pPpt->Out(Underline("Algorithm control", "-", "\t"));
		pPpt->Out("\tContinuity convergence tolerance, tol_mass_sum\t=\t"); pPpt->Out(tol_mass_sum); pPpt->Out("\n");
		pPpt->Out("\tLoss convergence tolerance, tol_del_star\t=\t"); pPpt->Out(tol_del_star); pPpt->Out("\n");
		pPpt->Out("\tMax. pressure loss iterations,loop_limit_main\t=\t"); pPpt->Out(loop_limit_main); pPpt->Out("\n");
		pPpt->Out("\n");
	}

	// Deprecated
	// ----------------------------------------------------------------------------------------------------
	/*
	cout << "\tNumber of branches, NPIPES\t\t\t=\t" << NPIPES << endl;
	switch(type)
	{
	case HCPJID:
		cout << "\tUsing homentropic constant pressure junction model (HCPJID)\n";
		break;
	case NHCPJID:
		cout << "\tUsing non-homentropic constant pressure junction model (NHCPJID)\n";
		break;
	case NHPLTJID:
		cout << "\tUsing non-homentropic pressure loss tee-junction model (NHPLTJID)\n";
		break;
	case NHPLJID:
		cout << "\tUsing non-homentropic pressure loss junction model (NHPLJID)\n";
		break;
	case NHPLPCID:
		cout << "\tUsing non-homentropic pressure loss pulse converter model (NHPLPCID)\n";
		break;
	default:
		cout << "\tUnknown pressure loss junction model...exiting\n";
		exit(1);
		break;
	}
	cout << "\tangle_deg\t=\t" << angle_deg << Deg() << endl;
	cout << "\ttol_main\t=\t" << tol_main << " \n";
	cout << "\tloop_limit_main\t=\t" << loop_limit_main << " \n";
	cout << "\ttol_cont\t=\t" << tol_cont << " \n";
	cout << "\tloop_limit_cont\t=\t" << loop_limit_cont << " \n";
	cout << "\tloop_limit_switch\t=\t" << loop_limit_switch << " \n";
	if(PRINT_MOVIE_FILE) cout << "\tPrinting movie file\n";
	cout << endl;
	*/

	pPpt->Out("\n");
	
	// ====================================================================================================
	// End of file
	// ====================================================================================================
}
