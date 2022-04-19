// Boundary.cpp: implementation of the CBoundary class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Boundary.h"
#include "Tools.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CBoundary::CBoundary()
{

}

CBoundary::~CBoundary()
{

}

void CBoundary::InitialiseGen(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rPIPES, int** &rPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, int assyid, string calling_function_str, string parent_assy_res_dir)
// ====================================================================================================
// General boundary condition initialization routine
// ====================================================================================================
{
	int i;

	// Identification
	// ----------------------------------------------------------------------------------------------------
	ID = id;
	EX = ex;
	AssyID = assyid;
	NPIPES = npipes;
	if(pPpt->SHOW_calls){pPpt->Out("CBoundary.InitialiseGen ("); pPpt->Out(calling_function_str); pPpt->Out(")\n");}

	// Construct the correct assembly results directory
	// ----------------------------------------------------------------------------------------------------
	RES_DIR = StringToChar(parent_assy_res_dir);

	// Join the pipe pointer(s) for this boundary to the relevant pipe(s) on the correct assembly
	// ----------------------------------------------------------------------------------------------------
	pPipe = new CPipe* [NPIPES];
	for(i=0; i<NPIPES; ++i) pPipe[i] = pPipes[i];
	
	// Set the end correction(s) for the pipe(s) joining this boundary
	// ----------------------------------------------------------------------------------------------------
	end = new int [NPIPES];

	for(i=0; i<NPIPES; ++i) {
		if(rPIPES_ENDS[ID][i] == ODD) {
			end[i] = ODD;
			pPipe[i]->end_corr_odd_p = rENDCORR[ID];
		}
		else {
			if(rPIPES_ENDS[ID][i] == EVEN){
				end[i] = EVEN;
				pPipe[i]->end_corr_even_p = rENDCORR[ID];
			}
			else {
				cout << rPIPES_ENDS[ID][i] << " is not a valid pipe end specifier" << endl;
				cout << "End specifiers must be either " << ODD << " (ODD), or " << EVEN << " (EVEN)" << endl;
				exit(1);
			}
		}
	}
}

void CBoundary::Configure(CProperties* pPpt)
// ====================================================================================================
// General boundary condition configuration routine
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Configure\n");}

	// Join the boundary node pointer(s) for this boundary to the relevant pipe nodes
	// ----------------------------------------------------------------------------------------------------
	pBN = new CNode* [NPIPES];
	pCLIN = new double** [NPIPES];
	pCLOUT = new double** [NPIPES];
	CLIN_STAR = new double* [NPIPES];
	CLOUT_STAR = new double* [NPIPES];
	pend_flow = new int* [NPIPES];
	pPathLine = new CPathLine* [NPIPES];
	int i;
	for(i=0; i<NPIPES; ++i){

		CLIN_STAR[i] = new double[2];
		CLOUT_STAR[i] = new double[2];
		
		if(end[i] == ODD){
			pBN[i] = &(pPipe[i]->Node[0]);
			pCLIN[i] = &(pPipe[i]->Node[0].CL2);
			pCLOUT[i] = &(pPipe[i]->Node[0].CL1);
			CLIN_STAR[i][R] = pPipe[i]->Node[0].CL2[R] / pPipe[i]->Node[0].AA[R];
			CLOUT_STAR[i][R] = pPipe[i]->Node[0].CL1[R] / pPipe[i]->Node[0].AA[R];
			pend_flow[i] = &(pPipe[i]->odd_end_flow);
			pPathLine[i] = &(pPipe[i]->PathLine[0]);
		}
		else{
			if(end[i] == EVEN){
				pBN[i] = &(pPipe[i]->Node[pPipe[i]->N-1]);
				pCLIN[i] = &(pPipe[i]->Node[pPipe[i]->N-1].CL1);
				pCLOUT[i] = &(pPipe[i]->Node[pPipe[i]->N-1].CL2);
				CLIN_STAR[i][R] = pPipe[i]->Node[pPipe[i]->N - 1].CL1[R] / pPipe[i]->Node[pPipe[i]->N - 1].AA[R];
				CLOUT_STAR[i][R] = pPipe[i]->Node[pPipe[i]->N - 1].CL2[R] / pPipe[i]->Node[pPipe[i]->N - 1].AA[R];
				pend_flow[i] = &(pPipe[i]->even_end_flow);
				pPathLine[i] = &(pPipe[i]->PathLine[pPipe[i]->num_pathlines-1]);
			}
			else{
				cout << end[i] << " is not a valid pipe end specifier" << endl;
				cout << "End specifiers must be either " << ODD << " (ODD), or " << EVEN << " (EVEN)" << endl;
				exit(1);
			}
		}
		CLIN_STAR[i][R+1] = CLIN_STAR[i][R];
		CLOUT_STAR[i][R+1] = CLIN_STAR[i][R];
		pBN[i]->bc = NAME[i];
	}

	// Initialise member variables
	// ----------------------------------------------------------------------------------------------------
	pipe_flow = new int [NPIPES];
	pipe_flow_old = new int [NPIPES];
	for(i=0; i<NPIPES; ++i){
		pipe_flow[i] = *(pend_flow[i]);
		pipe_flow_old[i] = *(pend_flow[i]);
	}
}

void CBoundary::PrintConnections(CProperties* pPpt)
// ====================================================================================================
// Print boundary pipe connections
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".PrintConnections\n");}

	cout << Underline(Identify(), "-", "\t");
	for(int p=0; p<NPIPES; ++p) {
		pPpt->Out("\tconnection [");  pPpt->Out(p); pPpt->Out("] = "); pPpt->Out(pPipe[p]->Identify()); pPpt->Out("\n");
	}
	pPpt->Out("\n");
}

void CBoundary::SetupFiles(CProperties* pPpt, std::string res_str)
// ====================================================================================================
// General routine to open boundary condition results file
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".SetupFiles("); pPpt->Out(res_str); pPpt->Out(")\n");}
	OUTPUT_FILE = fopen(ConstructString(pPpt, RES_DIR, res_str, ID, pPpt->strFileExt), "w");
	if(PRINT_DEBUG_FILE) OUTPUT_FILE_DEBUG = fopen(ConstructString(pPpt, RES_DIR, res_str, ID, ".debug"), "w");

	std::string strMovExt;
	strMovExt = "_mov";
	strMovExt += pPpt->strFileExt; // Add the file extension
	if(PRINT_MOVIE_FILE) MOVIE_FILE = fopen(ConstructString(pPpt, RES_DIR, res_str, ID, strMovExt), "w");
}

void CBoundary::CloseFiles(CProperties* pPpt)
// ====================================================================================================
// General routine to close boundary condition results file
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".CloseFiles\n");}

	fclose(OUTPUT_FILE);
	if(PRINT_DEBUG_FILE) fclose(OUTPUT_FILE_DEBUG);
	if(PRINT_MOVIE_FILE) fclose(MOVIE_FILE);
}

char* CBoundary::Identify()
// ====================================================================================================
// Returns identification of the current object
// ====================================================================================================
{
	std::string sz;
	sz = "Assembly [";
	sz += IntToString(AssyID);
	sz += "], ";
	if(EX) sz += "Exhaust "; else sz += "Intake ";

	sz += GetBoundaryName(NAME[0]); sz += " ["; sz += IntToString(ID); sz += "]";

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str());
	return szz;
	delete [] szz;
}

void CBoundary::common_HI_code(CProperties* pPpt, double lambda_in, double &rlambda_out, bool &rCHOKED, double A0, double T)
// ====================================================================================================
// Inflow routine for open end, partially open end, valve, or port; homentropic flow
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".common_HI_code\n");}

	double A = (3-pPpt->gammaAir(T))/(pPpt->gammaAir(T)+1);
	double B = 4*(pPpt->gammaAir(T)-1)/(pPpt->gammaAir(T)+1);
	double C = ((pPpt->gammaAir(T)+1)/(3-pPpt->gammaAir(T)))*lambda_in; // Choked flow limit upon lambda_out
	
	rlambda_out = A*lambda_in + sqrt(B*pow(A0,2) - (1.0-pow(A,2))*pow(lambda_in,2));
	// This simplifies to rlambda_out=lambda_in for lambda_in=A0, i.e. closed end

	if(rlambda_out >= C){ // Test for choked flow
		rCHOKED = true;
		rlambda_out = C;	
	}
}

double CBoundary::common_HN_code(CProperties* pPpt, int timestep, double time, double lambda_in, double phi, double Pb, bool &rCHOKED, double T)
// ====================================================================================================
// Outflow routine for open end, partially open end, valve, or port; homentropic flow
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".common_HN_code\n");}

	double* A;
	A = new double [6];
	// A[0] is used at the end of the convergence loop for the ATT sum
	A[1] = (pPpt->gammaAir(T)+1.0)/(pPpt->gammaAir(T)-1.0);
	A[2] = 2.0/(pPpt->gammaAir(T)-1.0);
	A[3] = 2.0*A[2];
	A[4] = 2.0/(pPpt->gammaAir(T)+1.0);
	A[5] = 1.0/A[2];
	
	double K1 = (pPpt->gammaAir(T)-1.0)/(2.0*pPpt->gammaAir(T));
	
	double lambda_out;
	double AEST, DAEST, SUM;
	double UTEMP, UT;
	double ATT, DATT, Y;
	double CON, SKS;
	
	AEST = 0.5*(lambda_in+1.0);
	DAEST = 0.25*(lambda_in-1.0);

	do{
		SUM = pow(lambda_in - AEST, 2)*(pow(AEST, A[3]) - pow(phi, 2)) - (A[5]*pow(phi, 2))*(pow(AEST, 2) - 1.0); // Eqn (6.35)
		if(SUM<0.0)	AEST -= DAEST;
		else AEST += DAEST;
		DAEST /= 2;
	}
	while(DAEST>pPpt->NOZZ_TOL || fabs(SUM)>=pPpt->NOZZ_TOL);

	lambda_out = 2.0*AEST - lambda_in;
	UTEMP = (lambda_in - lambda_out)/(pPpt->gammaAir(T)-1.0);
	UT = pow(AEST, A[2])*UTEMP/phi;

	if(!(UT<=1.0))
	{
		rCHOKED = true;
		// Initial values of ATT and DATT
		ATT = 0.5 + 0.5*sqrt(A[4]);
		DATT = 0.25 - 0.25*sqrt(A[4]);
		do{
			Y = 1.0/ATT;
			SUM = pow(phi, 2) - (pow(Y, A[3]))*(A[1] - A[2]*pow(Y, 2)); // Eqn (6.38)
			if(SUM<0.0) ATT -= DATT;
			else ATT += DATT;
			DATT /= 2;
		}while(DATT>pPpt->NOZZ_TOL || fabs(SUM)>=pPpt->NOZZ_TOL);
		
		// DATT <= 0.0001, so just do:
		A[0] = pow(ATT, A[1]);
		CON = pPpt->K_entropy*A[5]*phi*A[0]; // Intermediate variable
		// K_entropy will be 1, apart from when mean entropy model 3 is used
		
		SKS = (1.0 - CON)/(1.0 + CON); // Equation 6.41
		lambda_out = SKS*lambda_in;
	}
	else rCHOKED = false;
	delete [] A;
	return lambda_out;
}

void CBoundary::common_NHI_code(CProperties* pPpt, double lambda_in_n, double lambda_out_n, double AA_n, double &rlambda_in_c, double &rlambda_out_c, double &rAA_c, double psi, double P0, double T0, bool &rCHOKED, bool &rSONIC, double T, int timestep, double time, bool INFLOW_CONST_P, double tolLoop)
// ====================================================================================================
// Inflow routine for open end, partially open end, valve, or port; non-homentropic flow
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".common_NHI_code\n");}

	double Ac, rc, AA_is, del_AA;
	double U, C, R, PpbyPc, AA_calc;
	double lambda_in_c, lambda_out_c, AA_est, lambda_in_c_prev, lambda_out_c_prev;
	bool CHOKED, SONIC;

	bool GREATER, GREATER_PREV;
	GREATER = false;
	GREATER_PREV = false;

	double k = pPpt->gammaAir(T);
		
	// Numbers in braces refer to steps in Benson's (1982, pg. 378) non-homentropic outflow to a pipe routine (Benson, 1982, pg. 378)
	// ----------------------------------------------------------------------------------------------------

	// (1) Enter with known characteristic values from boundary point
	// ----------------------------------------------------------------------------------------------------
	double A0 = sqrt(T0/(EX ? pPpt->TREFe : pPpt->TREFi));
	Ac = A0;	// Non-dimensional speed of sound in source (e.g., cylinder)
	
	// (2) Set uncorrected values
	// ----------------------------------------------------------------------------------------------------
	double PC = P0;
	AA_is = Ac*pow(pPpt->PREF/PC, (k-1)/(2*k));	// Entropy level for isentropic process
	lambda_in_c = lambda_in_n;					// Use uncorrected value as initial estimate
	lambda_out_c = lambda_out_n;				// Use previous result fot initial estimate
	
	// (3) Test rc
	// ----------------------------------------------------------------------------------------------------
	rc = PC/pPpt->PREF;
	if(rc < 1.0){AA_est = AA_is; del_AA = 0.5;}
	else{
		AA_est = AA_is; 
		del_AA = fabs(Ac - AA_is)/2;
		if(del_AA==0) del_AA = 1e-6; 
	}
	
	bool RESET;
	int loop = 0;
	do{
		RESET = false;
		++loop;

		// (4) Calculate new value of lambda_in from eqn (7.155)
		// ----------------------------------------------------------------------------------------------------
		lambda_in_c_prev = lambda_in_c;
		//lambda_in_c = 2*lambda_in_n*(AA_est/(AA_est + AA_n)) + lambda_out_c*((AA_est - AA_n)/(AA_est + AA_n));
		lambda_in_c = ((lambda_in_c + lambda_out_c)/2)*(1 - (AA_n/AA_est)) + lambda_in_n; // Eqn 7.155

		// (5) Test (Ac - lambda_in_c)
		// ----------------------------------------------------------------------------------------------------
		if(Ac - lambda_in_c < 0){
			RESET = true;
			AA_est = AA_is;
			del_AA /= 2;
			lambda_in_c_prev = lambda_in_c;
			//lambda_in_c = 2*lambda_in_n*(AA_est/(AA_est + AA_n)) + lambda_out_c*((AA_est - AA_n)/(AA_est + AA_n));
			lambda_in_c = ((lambda_in_c + lambda_out_c)/2)*(1 - (AA_n/AA_est)) + lambda_in_n; // Eqn 7.155
		}

		// (6) Calculate new value of lambda_out from eqn (7.154)
		// ----------------------------------------------------------------------------------------------------
		lambda_out_c_prev = lambda_out_c;
		lambda_out_c = ((3-k)/(k+1))*lambda_in_c + (2/(k+1))*sqrt((pow(k,2)-1)*pow(Ac,2) + 2*(1-k)*pow(lambda_in_c,2));

		// (7) Calculate U from eqn (7.149)
		// ----------------------------------------------------------------------------------------------------
		U = (lambda_out_c - lambda_in_c)/((k-1)*Ac);

		if(INFLOW_CONST_P){ // Inflow constant pressure (valve) model

			// (8) Calculate C from eqn (7.139)
			// ----------------------------------------------------------------------------------------------------
			C = (((k-1)/2)*pow(U,2))/pow((1 - ((k-1)/2)*pow(U,2)),2);

			// (9) Test 4C/(k^2 - 1) - psi^2
			// ----------------------------------------------------------------------------------------------------
			if((4*C)/(pow(k,2)-1) - pow(psi,2) >= 0){	// Sonic flow in the valve
//cout << "Sonic flow in the valve" << endl;
				// (10) If sonic flow, calculate pp/pc from eqn (7.142)
				// ----------------------------------------------------------------------------------------------------
				CHOKED = true;
				SONIC = true;
                if(U==0) PpbyPc = 1.0;
                else PpbyPc = psi*pow(2/(k+1), (k+1)/(2*(k-1)))*(1 - ((k-1)/2)*pow(U,2))*(1/U);
				if(PpbyPc==0) PpbyPc=1; // Can happen if U is very small but not zero
			}
			else{
				// (10) If subsonic flow calculate pp/pc fom eqn (7.138)
				// ----------------------------------------------------------------------------------------------------
				CHOKED = false;
				SONIC = false;
                if(U==0) PpbyPc = 1.0;
                else PpbyPc = pow((1/(2*C))*(psi*sqrt(pow(psi,2) + 4*C) - pow(psi,2)), k/(k-1));
				if(PpbyPc==0) PpbyPc=1; // Can happen if U is very small but not zero
			}
		}
		else{ // Inflow pressure drop (port) model
			// (8) Calculate R from eqn (7.157)
			// ----------------------------------------------------------------------------------------------------
			R = (pow(psi,2)*pow(1 - ((k-1)/2)*pow(U,2), (k-3)/(k-1)))/((k-1)*pow(U,2));

			// (9) Test (2/(k+1))*pow(1 - ((k-1)/2)*pow(U,2), 2/(k-1))*pow(U,2) - pow(psi,2)
			// ----------------------------------------------------------------------------------------------------
			if( (2/(k+1))*pow(1 - ((k-1)/2)*pow(U,2), 2/(k-1))*pow(U,2) - pow(psi,2) >= 0){	// Sonic flow in the ports
				// (10) If sonic flow, calculate pp/pc from eqn (7.159)
				// ----------------------------------------------------------------------------------------------------
				CHOKED = true;
				SONIC = true;
                if(U==0) PpbyPc = 1;
                else PpbyPc = psi*pow(2/(k+1), (k+1)/(2*(k-1)))*(1 - ((k-1)/2)*pow(U,2))*(1/U);
				if(PpbyPc==0) PpbyPc=1; // Can happen if U is very small but not zero
			}
			else{
				// (10) If subsonic flow calculate pp/pc fom eqn (7.156)
				// ----------------------------------------------------------------------------------------------------
				CHOKED = false;
				SONIC = false;
                if(U==0) PpbyPc = 1;
                else PpbyPc = pow(R*(pow(1 + (2/R)*(1 - ((k-1)/2)*pow(U,2)), 0.5) - 1), k/(k-1));
				if(PpbyPc==0) PpbyPc=1; // Can happen if U is very small but not zero
			}
		}

		// (11) Calculate AA from eqn (7.153)
		// ----------------------------------------------------------------------------------------------------
		AA_calc = ((lambda_in_c + lambda_out_c)/2)*pow(1/PpbyPc, (k-1)/(2*k))*pow(1/rc, (k-1)/(2*k));

		// (12) Test {(AA)calc - AAest}
		// ----------------------------------------------------------------------------------------------------
		if(loop>1) GREATER_PREV = GREATER;
		if(AA_calc - AA_est > 0){
			AA_est += del_AA;
			GREATER = true;
		}
		else{
			AA_est -= del_AA;
			GREATER = false;
		}
		if(GREATER_PREV != GREATER && loop > 1) del_AA /= 2;

		// Debug
		// ----------------------------------------------------------------------------------------------------
		/*
		if(time>0.0065){
		cout << "time = " << time << endl; 
		cout << "timestep = " << timestep << endl;
		cout << "U = " << U << endl;
		cout << "C = " << C << endl;
		cout << "PpbyPc = " << PpbyPc << endl;
		cout << "AA_calc = " << AA_calc << endl;
		cout << "AA_est = " << AA_est << endl;
		cout << "fabs(AA_calc - AA_est) = " << fabs(AA_calc - AA_est) << endl;
		cout << "lambda_in_c = " << lambda_in_c << endl;
		cout << "lambda_in_c_prev = " << lambda_in_c_prev << endl;
		cout << "fabs(lambda_in_c - lambda_in_c_prev) = " << fabs(lambda_in_c - lambda_in_c_prev) << endl;
		cout << "lambda_out_c = " << lambda_out_c << endl;
		cout << "lambda_out_c_prev = " << lambda_out_c_prev << endl;
		cout << "fabs(lambda_out_c - lambda_out_c_prev) = " << fabs(lambda_out_c - lambda_out_c_prev) << endl;
		cout << endl;
		}
		//*/
	}while(fabs(AA_calc - AA_est) >= tolLoop || fabs(lambda_in_c - lambda_in_c_prev) >= tolLoop || fabs(lambda_out_c - lambda_out_c_prev) >= tolLoop);
  rlambda_in_c = lambda_in_c;
	rlambda_out_c = lambda_out_c;
	rAA_c = AA_est;
	rCHOKED = CHOKED;
	rSONIC = SONIC;

	// Debug
	// ----------------------------------------------------------------------------------------------------
	/*
	if(time>0.0065){
	cout << "time = " << time << endl; 
	cout << "timestep = " << timestep << endl;
	cout << "Final loop = " << loop << endl;
	cout << "U = " << U << endl;
	cout << "C = " << C << endl;
	cout << "PpbyPc = " << PpbyPc << endl;
	cout << "AA_calc = " << AA_calc << endl;
	cout << "AA_est = " << AA_est << endl;
	cout << "fabs(AA_calc - AA_est) = " << fabs(AA_calc - AA_est) << endl;
	cout << "lambda_in_c = " << lambda_in_c << endl;
	cout << "lambda_in_c_prev = " << lambda_in_c_prev << endl;
	cout << "fabs(lambda_in_c - lambda_in_c_prev) = " << fabs(lambda_in_c - lambda_in_c_prev) << endl;
	cout << "lambda_out_c = " << lambda_out_c << endl;
	cout << "lambda_out_c_prev = " << lambda_out_c_prev << endl;
	cout << "fabs(lambda_out_c - lambda_out_c_prev) = " << fabs(lambda_out_c - lambda_out_c_prev) << endl;
	cout << "U = " << (lambda_in_c - lambda_out_c)/(k-1) << endl;
	cout << "u = " << ((lambda_in_c - lambda_out_c)/(k-1))*pPipe[ONE_SIDE]->AREF << endl;
	cout << endl;
	}
	//*/
}

double CBoundary::common_NHN_code(CProperties* pPpt, double lambda_in, double AA, double phi, double Pb, bool &rCHOKED, double T, double time)
// ====================================================================================================
// Outflow routine for open end, partially open end, valve, or port; non-homentropic flow
// ====================================================================================================
{	
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".common_NHN_code\n");}

	double k = pPpt->gammaAir(T);
	double rb = Pb/pPpt->PREF;
	double lambda_in_star = lambda_in/(AA*pow(rb, (k-1)/(2*k)));
	double lambda_out;
	double A_star = (lambda_in_star + 1)/2;		// A_star lies between 1.0 and lambda_in_star - start at halfway
	//double del_A_star = (lambda_in_star - 1)/4;	// Use a quarter of the difference the first time
	double del_A_star = fabs( (lambda_in_star - 1) / 4);	// Use a quarter of the difference the first time
	double SUM;
double A_star_prev;
	bool ALWAYS_HALVE = true; // Always halve step size even if there is no sign change in SUM (leads to less loops)
	bool GREATER = false;
	int loop = 0;
	do{
		A_star_prev = A_star;
		++loop;
		SUM = (pow(A_star, 4/(k-1)) - pow(phi,2))*pow((lambda_in_star - A_star),2) - ((k-1)/2)*pow(phi,2)*(pow(A_star,2) - 1); // Eqn (7.122)
		if(SUM<0.0){	
			A_star -= del_A_star;
			if(ALWAYS_HALVE || (GREATER && loop>1)) del_A_star /= 2;
			GREATER = false;
		}
		else{
			A_star += del_A_star;
			if(ALWAYS_HALVE || (!GREATER && loop>1)) del_A_star /= 2;
			GREATER = true;
		}
	}while(fabs(SUM)>pPpt->NOZZ_TOL || fabs(A_star - A_star_prev)>pPpt->NOZZ_TOL);

	double lambda_out_star = 2*A_star - lambda_in_star;				// Using eqn (7.118)
	double U_star_temp = (lambda_in_star - lambda_out_star)/(k-1);	// Using eqn (7.117)
	double Mt = U_star_temp*pow(A_star, 2/(k-1))/phi;				// Eqn (7.124)
	
	if(Mt>=1.0/*choked flow*/ || Mt<0/*unrealistic flow situation: inflow during outflow routine*/){
		rCHOKED = true;
		double A_star_cr = (1 + sqrt((k+1)/2))/2;		// From eqn (7.128) upper limit of A_star_cr is ((k+1)/2)^0.5, lower limit is 1, start at halfway
		double del_A_star_cr = (sqrt((k+1)/2) - 1)/4;	// Start with 1/4 of the difference
		double SUM_cr;
		GREATER = false;
		loop = 0;
		do{
			++loop;
			SUM_cr = pow(phi, 2) - ((k+1)/(k-1) - (2/(k-1))*pow(A_star_cr,2))*pow(A_star_cr, 4/(k-1)); // Eqn (7.128)
			if(SUM_cr>0.0){	
				A_star_cr -= del_A_star_cr;
				if(ALWAYS_HALVE || (!GREATER && loop>1)) del_A_star_cr /= 2;
				GREATER = true;
			}
			else{
				A_star_cr += del_A_star_cr;
				if(ALWAYS_HALVE || (GREATER && loop>1)) del_A_star_cr /= 2;
				GREATER = false;
			}
		}while(fabs(SUM_cr)>pPpt->NOZZ_TOL);
		lambda_out = (1 - ((k-1)/2)*phi*pow(1/A_star_cr, (k+1)/(k-1))) / (1 + ((k-1)/2)*phi*pow(1/A_star_cr, (k+1)/(k-1))) * lambda_in; // Using eqn (7.131)
	}
	else{
		rCHOKED = false;
		lambda_out = lambda_out_star*AA*pow(rb, (k-1)/(2*k));
	}
	return lambda_out;
}

void CBoundary::common_UPDATE_H(CProperties* pPpt, double* lambda_in_c, double* lambda_out_c, int* pipe_flow, bool* CHOKED)
// ====================================================================================================
// Update pipe variables with boundary condition results; homentropic flow
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".common_UPDATE_H\n");}

	for(int i=0; i<NPIPES; ++i) // Permits multiple sided boundary conditions
	{
		(*(pCLIN[i]))[R+1] = lambda_in_c[i];	// Not absolutely necessary as this won't change for homentropic methods
		(*(pCLOUT[i]))[R+1] = lambda_out_c[i];	// Always update lambda_out
		pipe_flow_old[i] = *(pend_flow[i]);		// Record previous flow directions
		*(pend_flow[i]) = pipe_flow[i];			// Update flow direction
		pBN[i]->CHOKED = CHOKED[i];				// Set choked condition of boundary node
	}
}

void CBoundary::common_UPDATE_NH(CProperties* pPpt, double* lambda_in_c, double* lambda_out_c, double* AA_c, int* pipe_flow, bool* CHOKED)
// ====================================================================================================
// Update pipe variables with boundary condition results; non-homentropic flow
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".common_UPDATE_NH\n");}

	for(int i=0; i<NPIPES; ++i) { // Permits multiple sided boundary conditions

		// Update for inflow and outflow
		// ----------------------------------------------------------------------------------------------------
		(*(pCLOUT[i]))[R+1] = lambda_out_c[i];	// Always update lambda_out
		pipe_flow_old[i] = *(pend_flow[i]);		// Record previous flow directions
		*(pend_flow[i]) = pipe_flow[i];			// Update flow direction
		
		if(*(pend_flow[i])==INFLOW)				// Only INFLOW creates pathlines
		{
			(*(pCLIN[i]))[R+1] = lambda_in_c[i];// Only need to update on INFLOW, because OUTFLOWs will maintain the uncorrected value
			if((pPipe[i]->METHOD==pPipe[i]->MMOC && !pPpt->HOMENTROPIC) || (pPipe[i]->METHOD==pPipe[i]->W_ALPHA_BETA && pPpt->COMBINED_WAB_MOC)) 
				(*(pPathLine[i])).AA = AA_c[i];		// Only INFLOW will have a new pathline at the end of the pipe ready to go
			
			pBN[i]->AA[R+1] = AA_c[i];			// This seems to be necessary to get the correct pipe pressure
			
			if(pipe_flow_old[i] != *(pend_flow[i]))
			// If old flow direction was not INFLOW then there was no new path line created
			// so adjust existing one by setting its XK value to the appropriate end
			{
				if((pPipe[i]->METHOD==pPipe[i]->MMOC && !pPpt->HOMENTROPIC) || (pPipe[i]->METHOD==pPipe[i]->W_ALPHA_BETA && pPpt->COMBINED_WAB_MOC))
				{
					//cout << "(*(pPathLine[i])).XK = " << (*(pPathLine[i])).XK << endl;
					//cout << "end[i] = " << end[i] << endl; 
					if(end[i]==ODD){(*(pPathLine[i])).XK = pPipe[i]->Node[0].X/*0*/; (*(pPathLine[i])).XK_old = (*(pPathLine[i])).XK;}
					else{(*(pPathLine[i])).XK = /*pPipe[i]->XPIPE*/pPipe[i]->Node[pPipe[i]->N-1].X; (*(pPathLine[i])).XK_old = (*(pPathLine[i])).XK;}
				}
			}
			else 
			// Here pipe_flow_old[i] == *(pend_flow[i])
			// Check that a new path line HAS been created correctly if INFLOW
			{
				if(pPpt->COMBINED_WAB_MOC)
				{
					if(*(pend_flow[i])==INFLOW && ((end[i]==ODD && (*(pPathLine[i])).XK!=pPipe[i]->Node[0].X/*0*/) || (end[i]==EVEN && (*(pPathLine[i])).XK!=/*pPipe->XPIPE*/pPipe[i]->Node[pPipe[i]->N-1].X)))
					{
						cout << "CBoundary::common_UPDATE_NH: should've created new path line but its XK is not right!!" << endl;
						if(end[i]==ODD)
							cout << "(*(pPathLine[i])).XK = " << (*(pPathLine[i])).XK << ", but it should be " << pPipe[i]->Node[0].X/*0*/ << endl;
						else cout << "(*(pPathLine[i])).XK = " << (*(pPathLine[i])).XK << ", but it should be " << /*pPipe->XPIPE*/pPipe[i]->Node[pPipe[i]->N-1].X << endl;
					}
				}
			}
		}
		pBN[i]->CHOKED = CHOKED[i];				// Set choked condition of boundary node

		// Finally provide starred lambda values at the boundary
		CLIN_STAR[i][R + 1] = lambda_in_c[i] / AA_c[i];
		CLOUT_STAR[i][R + 1] = lambda_out_c[i] / AA_c[i];
	}
}