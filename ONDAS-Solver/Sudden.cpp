// Sudden.cpp: implementation of the CSudden class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Sudden.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CSudden::CSudden()
{

}

CSudden::~CSudden()
{

}

void CSudden::Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rSPIPES, int** &rSPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, int assyid, string parent_assy_res_dir, string calling_object_str)
{
	if (pPpt->SHOW_calls) { pPpt->Out("CSudden.Initialise ("); pPpt->Out(calling_object_str); pPpt->Out(")\n"); }

	// Join the pipe pointer(s) to the relevant pipe(s)
	InitialiseGen(pPpt, pPipes, rPipe, rSPIPES, rSPIPES_ENDS, rENDCORR, id, ex, npipes, assyid, "CSudden", parent_assy_res_dir);

	// Boundary name
	NAME = new int [NPIPES];
	NAME[LEFT_SIDE] = SUDDEN_LEFT;	
	NAME[RIGHT_SIDE] = SUDDEN_RIGHT;

//	pipe_flow = new int [2];
//	pipe_flow_old = new int [2];

	LEFT_TO_RIGHT = true;
}

void CSudden::RunBoundary(CProperties* pPpt, double time, int timestep)
{
	//if(timestep>=504) cout << "top Timestep = " << timestep << endl;
	General(pPpt, time, timestep);
	//if(timestep>=504) cout << "bottom Timestep = " << timestep << endl;
}

void CSudden::General(CProperties* pPpt, double time, int timestep)
//--------------------------------------------------//
// Non-homentropic flow direction test				//
// -----------------------------------				//
// Determines flow direction and whether			//
// the flow experiences enlargement, contraction	//
// or none.											//
//--------------------------------------------------//
{
	int *pipe_flow; pipe_flow = new int [NPIPES];

	double lambda_in_c2, AA_c2, lambda_out1, lambda_out2;

	double* lambda_in_star_n;
	lambda_in_star_n = new double [2];

	lambda_in_star_n[LEFT_SIDE] = (*(pCLIN[LEFT_SIDE]))[R+1] / pBN[LEFT_SIDE]->AA[1];
	lambda_in_star_n[RIGHT_SIDE] = (*(pCLIN[RIGHT_SIDE]))[R+1] / pBN[RIGHT_SIDE]->AA[1];

//cout << "lambda_in_star_n[LEFT_SIDE] = " << lambda_in_star_n[LEFT_SIDE] << endl;
//cout << "lambda_in_star_n[RIGHT_SIDE] = " << lambda_in_star_n[RIGHT_SIDE] << endl;

	if(lambda_in_star_n[LEFT_SIDE] - lambda_in_star_n[RIGHT_SIDE] > 1e-12)
	{
		// Flow from left to right; left is upstream, right is downstream
//cout << "Flow is left to right" << endl;
//exit(1);
		UPSTREAM = LEFT_SIDE;
		DOWNSTREAM = RIGHT_SIDE;
		if(!LEFT_TO_RIGHT)
		{
//			cout << "Flow in area change is changing direction, to left-to-right\n";
//			cin >> pause;
		}
		LEFT_TO_RIGHT = true;
	}
	else
	{
		if(lambda_in_star_n[RIGHT_SIDE] - lambda_in_star_n[LEFT_SIDE] > 1e-12)
		{
			// Flow from right to left; left is downstream, right is upstream
//cout << "Flow is right to left\n";
			UPSTREAM = RIGHT_SIDE;
			DOWNSTREAM = LEFT_SIDE;
			if(LEFT_TO_RIGHT)
			{
//				cout << "Flow in area change is changing direction, to right-to-left\n";
//				cin >> pause;
			}
			LEFT_TO_RIGHT = false;
		}
		else
		{
//cout << "No flow in area change\n";
//			cin >> pause;
			// No flow - treat both as closed ends
			(*(pCLOUT[LEFT_SIDE]))[R+1] = (*(pCLIN[LEFT_SIDE]))[R+1];
			(*(pCLOUT[RIGHT_SIDE]))[R+1] = (*(pCLIN[RIGHT_SIDE]))[R+1];
			// No need to update lambda_in or AA
			return;
		}
	}
//if(timestep>=504) cout << "top Timestep = " << timestep << endl;
	if((pBN[UPSTREAM]->f_dash*pPpt->fref) - (pBN[DOWNSTREAM]->f_dash*pPpt->fref) > 1e-12)
	{
		// Upstream area is greater than downstream, hence sudden contraction
//if(ID==3)
//if(timestep>=504)
//cout << "Sudden contraction" << endl;
		TYPE = CONTRACTION;
		Contraction(pPpt, lambda_in_c2, AA_c2, lambda_out1, lambda_out2, pipe_flow);
	}
	else
	{
		if((pBN[DOWNSTREAM]->f_dash*pPpt->fref) - (pBN[UPSTREAM]->f_dash*pPpt->fref) > 1e-12)
		{
			// Downstream area is greater than upstream, hence sudden enlargement
//if(ID==3)
//if(timestep>=504)
//cout << "Sudden enlargement" << endl;
			TYPE = ENLARGEMENT;
			Enlargement(pPpt, lambda_in_c2, AA_c2, lambda_out1, lambda_out2, pipe_flow, time, timestep);
		}
		else
		{
			// Areas the same, but there is flow
			TYPE = EQUAL_AREAS;
//if(timestep>=504)
//cout << "Areas the same, but there is flow" << endl;
//cin >> pause;
			Enlargement(pPpt, lambda_in_c2, AA_c2, lambda_out1, lambda_out2, pipe_flow, time, timestep);
			//Contraction(pPpt, lambda_in_c2, AA_c2, lambda_out1, lambda_out2, pipe_flow);
		}
	}
//if(timestep>=504) cout << "bottom Timestep = " << timestep << endl;
	//////////////////////////////////////////
	double *lambda_in_c, *lambda_out_c, *AA_c;
	lambda_in_c = new double [NPIPES]; lambda_out_c = new double [NPIPES]; AA_c = new double [NPIPES];
			
	lambda_out_c[UPSTREAM] = lambda_out1;
	lambda_out_c[DOWNSTREAM] = lambda_out2;
	lambda_in_c[DOWNSTREAM] = lambda_in_c2; // lambda_in correction needed for DOWNSTREAM only
	AA_c[DOWNSTREAM] = AA_c2;

	bool* CHOKED; CHOKED = new bool [NPIPES];
	for(int p=0; p<NPIPES; ++p) CHOKED[p] = false;
	
	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out_c, AA_c, pipe_flow, CHOKED);
	//////////////////////////////////////////

	// Free dynamic memory
	// -------------------
	delete [] pipe_flow;
	delete [] lambda_in_star_n;
	delete [] lambda_in_c;
	delete [] lambda_out_c;
	delete [] AA_c;
	delete [] CHOKED;
}

void CSudden::Enlargement(CProperties* pPpt, double &rlambda_in_c2, double &rAA_c2, double &rlambda_out1, double &rlambda_out2, int* &rpipe_flow, double time, int timestep)
//--------------------------------------------------//
// Non-isentropic sudden enlargement				//
// ---------------------------------				//
// Developed from Benson's sudden enlargement model.//
//--------------------------------------------------//
{
	//double lambda_in_n[2];
	//double AA_n[2];
	//double lambda_in[2];
	//double lambda_out[2];
	//double A_star[2];
	//double U_star[2];

	// Allocate dynamic memory
	// -----------------------
	double* lambda_in_n;	lambda_in_n = new double [2];
	double*	AA_n;			AA_n = new double [2];
	double* lambda_in;		lambda_in = new double [2];
	double* lambda_out;		lambda_out = new double [2];
	double* A_star;			A_star = new double [2];
	double* U_star;			U_star = new double [2];

	double psi, RR, del_RR, S, S_sqrd, S_sign, N, T1, T2, T_ratio, temp_T2, del_T2, L, M, NN, sign, A1_over_A2, AA2_over_AA1;
	double lambda_in_star2, AA_1, lambda_in_d1, AA_c2, lambda_in_c2;
	bool T2_converged, T2_negative, GREATER;
	int counter;

	// Set flow directions
	// ===================
	rpipe_flow[UPSTREAM] = OUTFLOW;
	rpipe_flow[DOWNSTREAM] = INFLOW;

	// Calculate area ratio parameter for this flow direction (<1)
	// ===========================================================
	psi = (pBN[UPSTREAM]->f_dash*pPpt->fref)/(pBN[DOWNSTREAM]->f_dash*pPpt->fref);
	area_ratio = psi;

	// Enter inital characteristic values
	// ==================================
	lambda_in_n[UPSTREAM] = (*(pCLIN[UPSTREAM]))[1];
	lambda_in_n[DOWNSTREAM] = (*(pCLIN[DOWNSTREAM]))[1];
	AA_n[UPSTREAM] = (pBN[UPSTREAM])->AA[1];
	AA_n[DOWNSTREAM] = (pBN[DOWNSTREAM])->AA[1];

	lambda_in[UPSTREAM] = lambda_in_n[UPSTREAM];
	AA_1 = AA_n[UPSTREAM]; // Since no correction needed on this side

	// Initial values of RR and del_RR
	// ===============================
	RR = (1 + sqrt(2/(pPpt->gammaAir(pBN[UPSTREAM]->T)+1)))/2;		// Eq. 8.55
	del_RR = (1 - sqrt(2/(pPpt->gammaAir(pBN[UPSTREAM]->T)+1)))/4;	// Eq. 8.56

	// First part - loop to find limiting value of RR that gives N==0
	// ==============================================================
	GREATER = true;
	int counterRR = 0;
	do
	{
		++counterRR; if(counterRR>10000){pPpt->Out(Identify()); pPpt->Out(":Enlargement: counterRR = "); pPpt->Out(counterRR); pPpt->Out("\n");}
		S_sqrd = ((pPpt->gammaAir(pBN[UPSTREAM]->T)+1)/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1))*pow(RR,2) - 2/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1); // Eq. 8.52
		if(S_sqrd>0)
		{
			S = sqrt(S_sqrd); 
			S_sign = 1;
		}
		else
		{
			S = sqrt(-S_sqrd);
			S_sign = -1;
		}
		N = psi + ((RR*S_sign*S)/(1 + pPpt->gammaAir(pBN[UPSTREAM]->T)*S_sqrd*(RR - (S_sign*S)))); // Eq. 8.54
		if(N > 0)
		{
			if(!GREATER) del_RR = del_RR/2; // Halve step size only if error has changed sign
			RR = RR - del_RR;
			GREATER = true;
		}
		else
		{
			if(GREATER) del_RR = del_RR/2; // Halve step size only if error has changed sign
			RR = RR + del_RR;
			GREATER = false;
		}
	}while(fabs(N)>pPpt->tolSuddenEnlrgN);

	// Second part - loop to converge on T2
	// ====================================
	T2_converged = false;
	GREATER = false;
	del_T2 = S/4.5; // Better not to have del_T2 fit into T2 a whole number of times
	T2 = -0.0001;   // T2 must never be 0; leads to div by 0

	counter = 0;
	lambda_in_c2 = lambda_in_n[DOWNSTREAM];	// As initial estimate
	double T2_old = T2;
	do
	{
		++counter; if(counter>10000){pPpt->Out(Identify()); pPpt->Out(":Enlargement: T2 loop counter = "); pPpt->Out(counter); pPpt->Out("\n");}
		lambda_in[DOWNSTREAM] = lambda_in_c2;
		
		NN = (pPpt->gammaAir(pBN[UPSTREAM]->T)*psi - (pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/2);	// Eq. 8.59
/*
		// Equations as used in Benson's flow diagram
		if(NN<0) sign = 1; else sign = -1;
		L = ((1 + pPpt->gammaAir()*pow(T2,2))*psi)/(2*fabs(NN)*T2);				// Eq. 8.58
		M = (1 + ((pPpt->gammaAir()-1)/2)*pow(T2,2))/fabs(NN);						// Eq. 8.58	
		T1 = L + sqrt(pow(L,2) + sign*M);							// Eq. 8.57
		A1_over_A2 = sqrt(1 + ((pPpt->gammaAir()-1)/2)*(pow(T2,2) - pow(T1,2)));	// Eq. 8.60
		AA2_over_AA1 = pow( (pow(1/A1_over_A2, 2/(pPpt->gammaAir()-1)) * ((T2)/T1) * (1/psi)), (pPpt->gammaAir()-1)/(2*pPpt->gammaAir()) );	// Eq. 8.61
*/		

//		if(NN<0) sign = -1; else sign = 1;
		if(NN<0) sign = 1; else sign = -1;
		L = -( ((1 + pPpt->gammaAir(pBN[UPSTREAM]->T)*pow(T2,2))*psi)/(2*NN*T2) );					// Eq. 8.34
		M = (1 + ((pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/2)*pow(T2,2))/NN;							// Eq. 8.35

		if(pow(L,2) > M) // Check that we are not taking sqrt(negative)
		{
			T1 = L + sign*sqrt(pow(L,2) - M);							// Eq. 8.33
		}
		else
		{
			// Do not calculate new T1; go back and half step size for T2
			T1= T1;
			T2 = T2_old;
			del_T2 = del_T2/2;
		}
		A1_over_A2 = sqrt(1 + ((pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/2)*(pow(T2,2) - pow(T1,2)));	// Eq. 8.38

		// To prevent division by zero:
		if(fabs(T1)<1e-16) T_ratio = 1;
		else T_ratio = (fabs(T2)/T1);
		AA2_over_AA1 = pow( (pow(1/A1_over_A2, 2/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1)) * T_ratio * (1/psi)), (pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/(2*pPpt->gammaAir(pBN[UPSTREAM]->T)) );	// Eq. 8.39

		lambda_in_star2 = (lambda_in[DOWNSTREAM]/AA_1) * (1/AA2_over_AA1);		// Eq. 8.62
		A_star[DOWNSTREAM] = lambda_in_star2/(1 + ((pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/2)*T2);				// Eq. 8.63
		A_star[UPSTREAM] = (A1_over_A2)*(AA2_over_AA1)*A_star[DOWNSTREAM];		// Eq. 8.64
		U_star[UPSTREAM] = T1*(AA2_over_AA1)*A_star[DOWNSTREAM];				// Eq. 8.65

		lambda_in_d1 = AA_1*(A_star[UPSTREAM] + ((pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/2)*U_star[UPSTREAM]);	// Eq. 8.66

		T2_negative = false;
		T2_old = T2;
		if(lambda_in_d1>lambda_in[UPSTREAM])
		{
			if(!GREATER) del_T2 = del_T2/2;	// Halve step size only if error has changed sign
			do
			{
				temp_T2 = T2 + del_T2;
				if(temp_T2 >= 0) del_T2 = del_T2/2;	// T2 must always be -ve
				else T2_negative = true;
			}while(!T2_negative);
			T2 = temp_T2;
			GREATER = true;
		}
		else
		{
			if(GREATER) del_T2 = del_T2/2; // Halve step size only if error has changed sign
			do
			{
				temp_T2 = T2 - del_T2;
				if(temp_T2 >= 0) del_T2 = del_T2/2; // T2 must always be -ve
				else T2_negative = true;
			}while(!T2_negative);
			T2 = temp_T2;
			GREATER = false;
		}

		AA_c2 = (AA2_over_AA1)*AA_1;															// Eq. 8.68
		lambda_in_c2 = lambda_in_n[DOWNSTREAM] + A_star[DOWNSTREAM]*(AA_c2 - AA_n[DOWNSTREAM]);	// Eq. 8.67
		if(
			fabs(lambda_in[UPSTREAM] - lambda_in_d1)<pPpt->tolSuddenEnlrgLin
			|| fabs(del_T2)<1e-20
			) 
			T2_converged = true;
	}while(!T2_converged);

	lambda_out[UPSTREAM] = AA_1*(A_star[UPSTREAM] - ((pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/2)*U_star[UPSTREAM]); // Eq. 8.69
	U_star[DOWNSTREAM] = -1*T2*A_star[DOWNSTREAM];
	lambda_out[DOWNSTREAM] = AA_c2*(A_star[DOWNSTREAM] + ((pPpt->gammaAir(pBN[DOWNSTREAM]->T)-1)/2)*U_star[DOWNSTREAM]); // Eq. 8.70

	// Update lambda_in2 = lambda_in_c2, AA_2 = AA_c2, lambda_out2, lambda_out1, by reference
	// ======================================================================================
	rlambda_in_c2 = lambda_in_c2;
	rAA_c2 = AA_c2;
	rlambda_out1 = lambda_out[UPSTREAM];
	rlambda_out2 = lambda_out[DOWNSTREAM];

	// Update validation variables
	// ===========================
	V1 = U_star[UPSTREAM]/lambda_in_star2;
	del1 = A_star[UPSTREAM]/lambda_in_star2;
	del2 = A_star[DOWNSTREAM]/lambda_in_star2;

	// Free dynamic memory
	// -------------------
	delete [] lambda_in_n;
	delete [] AA_n;
	delete [] lambda_in;
	delete [] lambda_out;
	delete [] A_star;
	delete [] U_star;
	return;
}

void CSudden::Contraction(CProperties* pPpt, double &rlambda_in_c2, double &rAA_c2, double &rlambda_out1, double &rlambda_out2, int* &rpipe_flow)
//--------------------------------------------------//
// Isentropic sudden contraction					//
// -----------------------------					//
// Based on Benson's sudden contraction model,		//
// which assumes an isentropic change of state		//
// at the area change. Although there is no			//
// entropy change across the contraction there		//
// can be a temperature discontinuity, so a			//
// correction is still applied to the downstream	//
// lambda_in.										//
// The area ratio, phi, must be less than unity		//
// for the method to converge correctly.			//
//--------------------------------------------------//
{
	char pause;

	//double lambda_in_n[2];
	//double AA_n[2];
	//double lambda_in[2];
	//double lambda_out[2];
	//double A_star[2];
	//double U_star[2];
	//double A[2];
	//double U[2];

	// Allocate dynamic memory
	// -----------------------
	double* lambda_in_n;	lambda_in_n = new double [2];
	double*	AA_n;			AA_n = new double [2];
	double* lambda_in;		lambda_in = new double [2];
	double* lambda_out;		lambda_out = new double [2];
	double* A_star;			A_star = new double [2];
	double* U_star;			U_star = new double [2];
	double* A;				A = new double [2];
	double* U;				U = new double [2];

	double phi, RR, RR_old, N, N_old, grad, factor;
	double S, T2, del_S, T1;
	double lambda_in_star2, AA_1, lambda_in_d1, AA_c2, lambda_in_c2;
	double limit;
	
	bool changed_sign, GREATER, stop_RR, RR_converged, converged;
	
	limit = 1e-3;
	stop_RR = false;
	converged = false;

	rpipe_flow[UPSTREAM] = OUTFLOW;
	rpipe_flow[DOWNSTREAM] = INFLOW;

	phi = (pBN[DOWNSTREAM]->f_dash*pPpt->fref)/(pBN[UPSTREAM]->f_dash*pPpt->fref);
	this->area_ratio = phi;

	// Enter known characteristic values
	lambda_in_n[UPSTREAM] = (*(pCLIN[UPSTREAM]))[1];
	lambda_in_n[DOWNSTREAM] = (*(pCLIN[DOWNSTREAM]))[1];
	AA_n[UPSTREAM] = (pBN[UPSTREAM])->AA[1];
	AA_n[DOWNSTREAM] = (pBN[DOWNSTREAM])->AA[1];

	lambda_in[UPSTREAM] = lambda_in_n[UPSTREAM];
	lambda_in[DOWNSTREAM] = lambda_in_n[DOWNSTREAM];
	AA_1 = AA_n[UPSTREAM]; // Since no correction needed on this side

	// Initial values of RR and del_RR
	// ===============================

	// Iterate for RR == A_cr first
	RR_converged = false;
	changed_sign = false;
	factor = 1;
	
	RR = 1.5;
	N = pow(phi,2) - ( (pPpt->gammaAir(pBN[UPSTREAM]->T)+1)/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1) - (2/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1))*pow(RR, -2) )*pow(RR, -4/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1));
	RR_old = RR;
	RR = 1.25;
	do
	{
		N_old = N;
		N = pow(phi,2) - ( (pPpt->gammaAir(pBN[UPSTREAM]->T)+1)/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1) - (2/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1))*pow(1/RR, 2) )*pow(1/RR, 4/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1));	// Eq. 8.92
		grad = (N - N_old)/(RR - RR_old);
		RR_old = RR;
		
		// If sign of sum has changed, reduce the distance moved by half
		if( (N>0 && N_old<0) || (N<0 && N_old>0) ) changed_sign = true;
	
		if( fabs(N) < 1e-12 && fabs(RR - RR_old) < 1e-12 ) RR_converged = true;
		else
		{
			if(changed_sign)
			{ 
				factor *= 0.5;
				changed_sign = false;
			}
			RR = RR - factor*(N/grad);
			RR_converged = false;
		}
	}while(!RR_converged);

	// Second part
	// ===========
	double T1_sqrd, T2_sqrd;
	GREATER = false;
	S = 1/RR;
	del_S = fabs(0.5*(1/RR - 1));
	del_S = 0.5*(1/RR - 1);

	S = 1;
	del_S = 0.0001;

	
	bool T_sqrd_negative = true;

	lambda_in_c2 = lambda_in_n[DOWNSTREAM];	// As initial estimate
	int counter = 0;
	do
	{
		++ counter;
		lambda_in[DOWNSTREAM] = lambda_in_c2;
		T_sqrd_negative = true;
		do
		{
			T1_sqrd = ( (2/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1))*(1 - pow(S,2)) )/( 1 - (1/pow(phi,2))*pow(S,4/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1)) );
			T2_sqrd = (pow(S, 4/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1))/pow(phi,2))*( ((2/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1))*(1 - pow(S,2)))/(1 - (1/pow(phi,2))*pow(S,4/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1))) );

			//if(T1_sqrd<0 || T2_sqrd<0)
			if(T1_sqrd<=0 || T2_sqrd<=0)
			{
				//S = S - del_S; // Starting S value was the maximum, so reduce
				S = S + del_S;
				del_S = del_S/2;
			}
			else T_sqrd_negative = false;

		}while(T_sqrd_negative);

		T1 = sqrt(T1_sqrd);		// T1 must always be positive, take the positive root
		if(T1 >= phi*pow(RR, 2/(pPpt->gammaAir(pBN[UPSTREAM]->T)-1)))
		{
			cout << "Choked flow\n";
			cin >> pause;
		}
		T2 = -sqrt(T2_sqrd);	// T2 must always be negative, take the negative root

		AA_c2 = AA_1;

		lambda_in_star2 = (lambda_in[DOWNSTREAM]/AA_c2);

		A_star[DOWNSTREAM] = lambda_in_star2/(1 + ((pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/2)*T2);	// Eq. 8.89

		A[DOWNSTREAM] = A_star[DOWNSTREAM]*AA_c2;
		
		A[UPSTREAM] = S*A[DOWNSTREAM];
			A_star[UPSTREAM] = A[UPSTREAM]/AA_1;
		
		U[UPSTREAM] = T1*A[DOWNSTREAM];								// Eq. 8.90
			U_star[UPSTREAM] = U[UPSTREAM]/AA_1;	

		U[DOWNSTREAM] = T2*A[DOWNSTREAM];	
		
		lambda_in_d1 = A[UPSTREAM] + ((pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/2)*U[UPSTREAM];			// Eq. 8.91

		if(lambda_in_d1>lambda_in[UPSTREAM])
		{
			if(!GREATER) del_S = del_S/2;	// Has switched
			
			S = S - del_S;
			GREATER = true;
		}
		else
		{
			if(GREATER) del_S = del_S/2;		// Has switched

			S = S + del_S;
			GREATER = false;
		}

		AA_c2 = AA_1;
		lambda_in_c2 = lambda_in_n[DOWNSTREAM] + A_star[DOWNSTREAM]*(AA_c2 - AA_n[DOWNSTREAM]);	// Eq. 8.67
		if(fabs(lambda_in[UPSTREAM] - lambda_in_d1)<1e-6) converged = true;

	}while(!converged);

//	cout << "counter = " << counter << endl;

	//lambda_out[UPSTREAM] = AA_1*(A_star[UPSTREAM] - ((pPpt->gammaAir()-1)/2)*U_star[UPSTREAM]);				// Eq. 8.69
	lambda_out[UPSTREAM] = A[UPSTREAM] - ((pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/2)*U[UPSTREAM];
	U_star[DOWNSTREAM] = -1*T2*A_star[DOWNSTREAM];
	lambda_out[DOWNSTREAM] = AA_c2*(A_star[DOWNSTREAM] + ((pPpt->gammaAir(pBN[DOWNSTREAM]->T)-1)/2)*U_star[DOWNSTREAM]);			// Eq. 8.70

//	lambda_out[DOWNSTREAM] = A[DOWNSTREAM] + ((pPpt->gammaAir()-1)/2)*U[DOWNSTREAM];	



	// Update lambda_in2 = lambda_in_c2, AA_2 = AA_c2, lambda_out2, lambda_out1, by reference
	rlambda_in_c2 = lambda_in_c2;
	rAA_c2 = AA_c2;
	rlambda_out1 = lambda_out[UPSTREAM];
	rlambda_out2 = lambda_out[DOWNSTREAM];

	// Update validation variables
	V1 = U_star[UPSTREAM]/lambda_in_star2;
	del1 = A_star[UPSTREAM]/lambda_in_star2;
	del2 = A_star[DOWNSTREAM]/lambda_in_star2;
	
	// Free dynamic memory
	// -------------------
	delete [] lambda_in_n;
	delete [] AA_n;
	delete [] lambda_in;
	delete [] lambda_out;
	delete [] A_star;
	delete [] U_star;
	delete [] A;
	delete [] U;

	return;
}