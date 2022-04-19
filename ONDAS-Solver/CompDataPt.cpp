// CompDataPt.cpp: implementation of the CCompDataPt class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "CompDataPt.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCompDataPt::CCompDataPt()
{

}

CCompDataPt::~CCompDataPt()
{

}

// Copy constructor
CCompDataPt::CCompDataPt(const CCompDataPt& inCompDataPt)
{
	sp = inCompDataPt.sp;
	mfp = inCompDataPt.mfp;
	pr = inCompDataPt.pr;
	eta = inCompDataPt.eta;

	T2_over_T1 = inCompDataPt.T2_over_T1;
	M2 = inCompDataPt.M2;
	lambda_in_star = inCompDataPt.lambda_in_star;
	lambda_out_star = inCompDataPt.lambda_out_star;
	AA2_over_AA1 = inCompDataPt.AA2_over_AA1;

	lambda_in_star_max_comp = inCompDataPt.lambda_in_star_max_comp;
//	M2_max = inCompDataPt.M2_max;
	KS_comp = inCompDataPt.KS_comp;
	lambda_in_star_min_comp = inCompDataPt.lambda_in_star_min_comp;
}

CCompDataPt& CCompDataPt::operator=(const CCompDataPt& inCompDataPt)
{
	if(this != &inCompDataPt)
	{
		sp = inCompDataPt.sp;
		mfp = inCompDataPt.mfp;
		pr = inCompDataPt.pr;
		eta = inCompDataPt.eta;

		T2_over_T1 = inCompDataPt.T2_over_T1;
		M2 = inCompDataPt.M2;
		lambda_in_star = inCompDataPt.lambda_in_star;
		lambda_out_star = inCompDataPt.lambda_out_star;
		AA2_over_AA1 = inCompDataPt.AA2_over_AA1;

		lambda_in_star_max_comp = inCompDataPt.lambda_in_star_max_comp;
//		M2_max = inCompDataPt.M2_max;
		KS_comp = inCompDataPt.KS_comp;
		lambda_in_star_min_comp = inCompDataPt.lambda_in_star_min_comp;
	}
	return *this;
}

CCompDataPt CCompDataPt::operator+ (const CCompDataPt& rhsCompDataPt)
{
	// this + rhsCompDataPt
		
	CCompDataPt temp;
	temp = *this; // Set temp equal to this, then add rhs values

	temp.sp += rhsCompDataPt.sp;
	temp.mfp += rhsCompDataPt.mfp;
	temp.pr += rhsCompDataPt.pr;
	temp.eta += rhsCompDataPt.eta;

	temp.T2_over_T1 += rhsCompDataPt.T2_over_T1;
	temp.M2 += rhsCompDataPt.M2;
	temp.lambda_in_star += rhsCompDataPt.lambda_in_star;
	temp.lambda_out_star += rhsCompDataPt.lambda_out_star;
	temp.AA2_over_AA1 += rhsCompDataPt.AA2_over_AA1;

	temp.lambda_in_star_max_comp += rhsCompDataPt.lambda_in_star_max_comp;
//	temp.M2_max += rhsCompDataPt.M2_max;
	temp.KS_comp += rhsCompDataPt.KS_comp;
	temp.lambda_in_star_min_comp += rhsCompDataPt.lambda_in_star_min_comp;

	return temp;
}

CCompDataPt CCompDataPt::operator- (const CCompDataPt& rhsCompDataPt)
{
	// this - rhsCompDataPt
		
	CCompDataPt temp;
	temp = *this; // Set temp equal to this, then subtract rhs values

	temp.sp -= rhsCompDataPt.sp;
	temp.mfp -= rhsCompDataPt.mfp;
	temp.pr -= rhsCompDataPt.pr;
	temp.eta -= rhsCompDataPt.eta;

	temp.T2_over_T1 -= rhsCompDataPt.T2_over_T1;
	temp.M2 -= rhsCompDataPt.M2;
	temp.lambda_in_star -= rhsCompDataPt.lambda_in_star;
	temp.lambda_out_star -= rhsCompDataPt.lambda_out_star;
	temp.AA2_over_AA1 -= rhsCompDataPt.AA2_over_AA1;

	temp.lambda_in_star_max_comp -= rhsCompDataPt.lambda_in_star_max_comp;
//	temp.M2_max -= rhsCompDataPt.M2_max;
	temp.KS_comp -= rhsCompDataPt.KS_comp;
	temp.lambda_in_star_min_comp -= rhsCompDataPt.lambda_in_star_min_comp;
	
	return temp;
}

CCompDataPt CCompDataPt::operator/ (const CCompDataPt& denCompDataPt)
{
	// this / denCompDataPt
		
	CCompDataPt temp;
	temp = *this; // Set temp equal to this, then divide by denominator values

	temp.sp /= denCompDataPt.sp;
	temp.mfp /= denCompDataPt.mfp;
	temp.pr /= denCompDataPt.pr;
	temp.eta /= denCompDataPt.eta;

	temp.T2_over_T1 /= denCompDataPt.T2_over_T1;
	temp.M2 /= denCompDataPt.M2;
	temp.lambda_in_star /= denCompDataPt.lambda_in_star;
	temp.lambda_out_star /= denCompDataPt.lambda_out_star;
	temp.AA2_over_AA1 /= denCompDataPt.AA2_over_AA1;

	temp.lambda_in_star_max_comp /= denCompDataPt.lambda_in_star_max_comp;
//	temp.M2_max /= denCompDataPt.M2_max;
	temp.KS_comp /= denCompDataPt.KS_comp;
	temp.lambda_in_star_min_comp /= denCompDataPt.lambda_in_star_min_comp;

	return temp;
}

CCompDataPt CCompDataPt::operator* (const CCompDataPt& rhsCompDataPt)
{
	// this * rhsCompDataPt
		
	CCompDataPt temp;
	temp = *this; // Set temp equal to this, then multiply by rhs values

	temp.sp *= rhsCompDataPt.sp;
	temp.mfp *= rhsCompDataPt.mfp;
	temp.pr *= rhsCompDataPt.pr;
	temp.eta *= rhsCompDataPt.eta;

	temp.T2_over_T1 *= rhsCompDataPt.T2_over_T1;
	temp.M2 *= rhsCompDataPt.M2;
	temp.lambda_in_star *= rhsCompDataPt.lambda_in_star;
	temp.lambda_out_star *= rhsCompDataPt.lambda_out_star;
	temp.AA2_over_AA1 *= rhsCompDataPt.AA2_over_AA1;

	temp.lambda_in_star_max_comp *= rhsCompDataPt.lambda_in_star_max_comp;
//	temp.M2_max *= rhsCompDataPt.M2_max;
	temp.KS_comp *= rhsCompDataPt.KS_comp;
	temp.lambda_in_star_min_comp *= rhsCompDataPt.lambda_in_star_min_comp;

	return temp;
}

CCompDataPt CCompDataPt::operator* (const double value)
{
	// this * some double
		
	CCompDataPt temp;
	temp = *this; // Set temp equal to this, then multiply by rhs values

	temp.sp *= value;
	temp.mfp *= value;
	temp.pr *= value;
	temp.eta *= value;

	temp.T2_over_T1 *= value;
	temp.M2 *= value;
	temp.lambda_in_star *= value;
	temp.lambda_out_star *= value;
	temp.AA2_over_AA1 *= value;

	temp.lambda_in_star_max_comp *= value;
//	temp.M2_max *= value;
	temp.KS_comp *= value;
	temp.lambda_in_star_min_comp *= value;

	return temp;
}

void CCompDataPt::ProcessCompDataPtVar(CProperties* pPpt, double* f, double* C, double T)
// ----------------------------------------------------------------------//
// Processes input PR and MFP data into Riemann variables.		//																	//
// ----------------------------------------------------------------------//
{
	double c, M1, fM1, del_M1;
	double p01_over_p2, p01_over_p1, T01_over_T1, T2_over_T01;

	int UPSTREAM = 0;
	int DOWNSTREAM = 1;

	p01_over_p2 = this->pr;

	if(this->mfp<=0.0) // Zero or reverse flow
	{
		this->G1 = 0;
		this->T2_over_T1 = 1.0;
		if(fabs(this->mfp) < 1e-6)
		{
			this->p1_over_p2 = p01_over_p2;
			this->lambda_in_star_ratio = pow(p1_over_p2, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T)))
											* ( (1 + C[UPSTREAM]*G1)/(1 - C[DOWNSTREAM]*G1*p1_over_p2*sqrt(T2_over_T1)) );
		}
	}
	else if (this->mfp>0.0)
	{
		c = pow( ( (this->mfp/f[UPSTREAM])*sqrt(pPpt->R_air/pPpt->gammaAir(T)) ), (2*(pPpt->gammaAir(T)-1)/(pPpt->gammaAir(T)+1)) );
		M1 = 0.001;

		do
		{
			fM1 = pow( c + ((pPpt->gammaAir(T)-1)/2)*c*pow(M1,2), (pPpt->gammaAir(T)+1)/(2*(pPpt->gammaAir(T)-1)) );
			del_M1 = /*fabs*/(M1 - fM1);
			M1 = fM1;
		}
		while(del_M1>0.0001);

		p01_over_p1 = pow( 1 + ((pPpt->gammaAir(T)-1)/2)*pow(M1,2), pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1) );
		T01_over_T1 = 1 + ((pPpt->gammaAir(T)-1)/2)*pow(M1,2);
		T2_over_T01 = 1 - eta*(1 - pow(1/p01_over_p2, (pPpt->gammaAir(T)-1)/pPpt->gammaAir(T)));

		this->G1 = mfp*p01_over_p1*sqrt(1/T01_over_T1);
		this->p1_over_p2 = (1/p01_over_p1)*p01_over_p2;
		this->T2_over_T1 = T2_over_T01*T01_over_T1;

		this->lambda_in_star_ratio = pow(p1_over_p2, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T)))
											* ( (1 + C[UPSTREAM]*G1)/(1 - C[DOWNSTREAM]*G1*p1_over_p2*sqrt(T2_over_T1)) );
		// Store lambda_in_star_ratio, G1, p1_over_p2, T2_over_T1 as member variables
	}
}