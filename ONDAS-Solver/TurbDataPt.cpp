// TurbDataPt.cpp: implementation of the CTurbDataPt class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "TurbDataPt.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTurbDataPt::CTurbDataPt()
{

}

CTurbDataPt::~CTurbDataPt()
{

}

// Copy constructor
CTurbDataPt::CTurbDataPt(const CTurbDataPt& inTurbDataPt)
{
	lambda_in_star = inTurbDataPt.lambda_in_star;
	lambda_out_star = inTurbDataPt.lambda_out_star;
	sp = inTurbDataPt.sp;
	mfp = inTurbDataPt.mfp;
	pr = inTurbDataPt.pr;

	eta = inTurbDataPt.eta;
	G1 = inTurbDataPt.G1;
	p1_over_p2 = inTurbDataPt.p1_over_p2;
	T2_over_T1 = inTurbDataPt.T2_over_T1;	
	lambda_in_star_ratio = inTurbDataPt.lambda_in_star_ratio;
}

CTurbDataPt& CTurbDataPt::operator=(const CTurbDataPt& inTurbDataPt)
{
	if(this != &inTurbDataPt)
	{
		lambda_in_star = inTurbDataPt.lambda_in_star;
		lambda_out_star = inTurbDataPt.lambda_out_star;
		sp = inTurbDataPt.sp;
		mfp = inTurbDataPt.mfp;
		pr = inTurbDataPt.pr;

		eta = inTurbDataPt.eta;
		G1 = inTurbDataPt.G1;
		p1_over_p2 = inTurbDataPt.p1_over_p2;
		T2_over_T1 = inTurbDataPt.T2_over_T1;	
		lambda_in_star_ratio = inTurbDataPt.lambda_in_star_ratio;
	}
	return *this;
}

CTurbDataPt CTurbDataPt::operator+ (const CTurbDataPt& rhsTurbDataPt)
{
	// this + rhsTurbDataPt
		
	CTurbDataPt temp;
	temp = *this; // Set temp equal to this, then add rhs values

	temp.lambda_in_star += rhsTurbDataPt.lambda_in_star;
	temp.lambda_out_star += rhsTurbDataPt.lambda_out_star;
	temp.sp += rhsTurbDataPt.sp;
	temp.mfp += rhsTurbDataPt.mfp;
	temp.pr += rhsTurbDataPt.pr;

	temp.eta += rhsTurbDataPt.eta;
	temp.G1 += rhsTurbDataPt.G1;
	temp.p1_over_p2 += rhsTurbDataPt.p1_over_p2;
	temp.T2_over_T1 += rhsTurbDataPt.T2_over_T1;
	temp.lambda_in_star_ratio += rhsTurbDataPt.lambda_in_star_ratio;

	return temp;
}

CTurbDataPt CTurbDataPt::operator- (const CTurbDataPt& rhsTurbDataPt)
{
	// this - rhsTurbDataPt
		
	CTurbDataPt temp;
	temp = *this; // Set temp equal to this, then subtract rhs values

	temp.lambda_in_star -= rhsTurbDataPt.lambda_in_star;
	temp.lambda_out_star -= rhsTurbDataPt.lambda_out_star;
	temp.sp -= rhsTurbDataPt.sp;
	temp.mfp -= rhsTurbDataPt.mfp;
	temp.pr -= rhsTurbDataPt.pr;

	temp.eta -= rhsTurbDataPt.eta;
	temp.G1 -= rhsTurbDataPt.G1;
	temp.p1_over_p2 -= rhsTurbDataPt.p1_over_p2;
	temp.T2_over_T1 -= rhsTurbDataPt.T2_over_T1;
	temp.lambda_in_star_ratio -= rhsTurbDataPt.lambda_in_star_ratio;

	return temp;
}

CTurbDataPt CTurbDataPt::operator/ (const CTurbDataPt& denTurbDataPt)
{
	// this / denTurbDataPt
		
	CTurbDataPt temp;
	temp = *this; // Set temp equal to this, then divide by denominator values

	temp.lambda_in_star /= denTurbDataPt.lambda_in_star;
	temp.lambda_out_star /= denTurbDataPt.lambda_out_star;
	temp.sp /= denTurbDataPt.sp;
	temp.mfp /= denTurbDataPt.mfp;
	temp.pr /= denTurbDataPt.pr;

	temp.eta /= denTurbDataPt.eta;
	temp.G1 /= denTurbDataPt.G1;
	temp.p1_over_p2 /= denTurbDataPt.p1_over_p2;
	temp.T2_over_T1 /= denTurbDataPt.T2_over_T1;
	temp.lambda_in_star_ratio /= denTurbDataPt.lambda_in_star_ratio;
	
	return temp;
}

CTurbDataPt CTurbDataPt::operator* (const CTurbDataPt& rhsTurbDataPt)
{
	// this * rhsTurbDataPt
		
	CTurbDataPt temp;
	temp = *this; // Set temp equal to this, then multiply by rhs values

	temp.lambda_in_star *= rhsTurbDataPt.lambda_in_star;
	temp.lambda_out_star *= rhsTurbDataPt.lambda_out_star;
	temp.sp *= rhsTurbDataPt.sp;
	temp.mfp *= rhsTurbDataPt.mfp;
	temp.pr *= rhsTurbDataPt.pr;

	temp.eta *= rhsTurbDataPt.eta;
	temp.G1 *= rhsTurbDataPt.G1;
	temp.p1_over_p2 *= rhsTurbDataPt.p1_over_p2;
	temp.T2_over_T1 *= rhsTurbDataPt.T2_over_T1;
	temp.lambda_in_star_ratio *= rhsTurbDataPt.lambda_in_star_ratio;
	
	return temp;
}

CTurbDataPt CTurbDataPt::operator* (const double value)
{
	// this * some double
		
	CTurbDataPt temp;
	temp = *this; // Set temp equal to this, then multiply by rhs values

	temp.lambda_in_star *= value;
	temp.lambda_out_star *= value;
	temp.sp *= value;
	temp.mfp *= value;
	temp.pr *= value;

	temp.eta *= value;
	temp.G1 *= value;
	temp.p1_over_p2 *= value;
	temp.T2_over_T1 *= value;
	temp.lambda_in_star_ratio *= value;

	return temp;
}


void CTurbDataPt::ProcessTurbDataPtVar(CProperties* pPpt, double* f, double* C, double T)
// ----------------------------------------------------------------------//
// Processes input PR and MFP data into Riemann variables.		//
// pr == p01/p2																		//
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
